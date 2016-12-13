#!/usr/bin/env python

import os, sys
import tableio
from pyfits import getheader, getval
import numpy as numpy
from math import log10
import time
import Numeric
import pyfits
import string


class combcat:
    ''' Combine, swarp and get catalogs '''

    def __init__(self,
                 assocfile,
                 datapath='',
                 outpath='',
                 pixscale=0.2666,
                 dryrun=None,
                 verb='yes'):

        self.assocfile = assocfile
        self.datapath = datapath
        #self.tilename  = os.path.basename(assocfile).split('.')[0]
        self.tilename = os.path.basename(os.path.splitext(assocfile)[0])
        self.outpath = outpath
        self.dryrun = dryrun
        self.verb = verb
        self.DetImage = None
        self.centered = None
        self.got_zeropt = False

        # Check for environ vars
        if not os.getenv('BCSPIPE'):
            os.environ['BCSPIPE'] = os.path.join(os.environ['HOME'], 'BCSPIPE')
        self.BCSPIPE = os.getenv('BCSPIPE')

        # Set the dust_directory
        #os.environ['DUST_DIR'] = os.path.join(self.BCSPIPE,"LIB")

        # Set the pixel scale
        self.pixscale = pixscale
        self.read_assoc()

        # Initialize and set to zero the dust exctiction corrections
        self.XCorr = {}
        self.XCorrError = {}
        for filter in self.filters:
            self.XCorr[filter] = 0.0
            self.XCorrError[filter] = 0.0

        return

    def read_assoc(self):
        ''' Read the association files for a given tile'''

        self.filters = []
        self.filelist = []
        self.infiles = []

        self.filter = {}
        self.exptime = {}
        self.airmass = {}

        # Make image list per filter
        self.files = {}
        # Exptimes per filter
        self.exptimes = {}

        print "# Will read %s" % self.assocfile

        # Read in the assoc file
        for line in open(self.assocfile).readlines():

            if line[0] == "#":
                continue

            vals = line.split()
            fname = os.path.basename(vals[0])
            infile = vals[0]
            self.filter[fname] = vals[1]
            self.exptime[fname] = float(vals[2])
            self.airmass[fname] = float(vals[3])

            print fname, self.filter[fname]

            # A list of the files
            self.filelist.append(fname)
            self.infiles.append(infile)
            if vals[1] not in self.filters:
                self.filters.append(vals[1])

        self.filters.sort()
        for filter in self.filters:
            self.files[filter] = []
            self.exptimes[filter] = []

        # Read in the assoc file
        for line in open(self.assocfile).readlines():

            if line[0] == "#":
                continue
            vals = line.split()
            fname = os.path.basename(vals[0])
            infile = vals[0]
            filtername = vals[1]

            # Loop over filters
            for filter in self.filters:
                if filter == filtername:
                    self.files[filter].append(fname)
                    self.exptimes[filter].append(self.exptime[fname])

        return

    def read_assoc_old(self):
        ''' Read the association files for a given tile'''

        self.filters = []
        self.filelist = []
        self.infiles = []

        self.filter = {}
        self.exptime = {}
        self.airmass = {}

        print "# Will read %s" % self.assocfile

        # Read in the assoc file
        for line in open(self.assocfile).readlines():

            if line[0] == "#":
                continue

            vals = line.split()
            fname = os.path.basename(vals[0])
            infile = vals[0]
            self.filter[fname] = vals[1]
            self.exptime[fname] = float(vals[2])
            self.airmass[fname] = float(vals[3])

            # A list of the files
            self.filelist.append(fname)
            self.infiles.append(infile)
            if vals[1] not in self.filters:
                self.filters.append(vals[1])

        self.filters.sort()

        # Make image list per filter
        self.files = {}

        # Exptimes per filter
        self.exptimes = {}

        for filter in self.filters:

            self.files[filter] = []
            self.exptimes[filter] = []

            for file in self.filelist:

                if filter == self.filter[file]:
                    self.files[filter].append(file)

                if filter == self.filter[file]:
                    self.exptimes[filter].append(self.exptime[file])

            print "# %s %s" % (filter, self.files[filter])

        return

    def copyfiles(self, copy="yes"):
        ''' Copy files from the remote location to the current dir'''

        # Fits make sure that the output dir exists
        self.outdir = os.path.join(self.outpath, self.tilename)
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        # Go to the directory
        os.chdir(self.outdir)

        print >> sys.stderr, "We are in:", self.outdir

        for file in self.infiles:
            cmd1 = "rsync -av -e ssh --progress %s ." % os.path.join(
                self.datapath, file)
            print cmd1
            if copy:
                os.system(cmd1)

        return

    def center_dither(self, conf="SWarp-center.conf"):
        ''' Center the dither pattern for SWARP '''

        check_exe("swarp")

        # The configuration file
        center_conf = os.path.join(self.BCSPIPE, 'LIB/pars', conf)

        # First we need to get the center for all the files, swarp them all
        opts = {}
        opts["IMAGEOUT_NAME"] = os.path.join(
            self.outdir, "SWarp-%s-center.fits" % (self.tilename))
        opts["PIXEL_SCALE"] = self.pixscale
        opts["RESAMPLING_TYPE"] = "LANCZOS3"
        opts["CENTER_TYPE"] = "ALL"
        opts["NTHREADS"] = "0"

        #cmd  = "swarp %s*[0-9][0-9][0-9].fits -c %s " % (self.tilename,center_conf)
        cmd = "swarp %s -c %s " % (string.join(self.filelist), center_conf)
        for param, value in opts.items():
            cmd = cmd + "-%s %s " % (param, value)

        print cmd
        if not self.dryrun:
            print >> sys.stderr, "Centering mosaic"
            os.system(cmd)

        # Read in the header
        header = getheader(opts["IMAGEOUT_NAME"])
        nx = header["NAXIS1"]
        ny = header["NAXIS2"]
        x_center = header["CRVAL1"]
        y_center = header["CRVAL2"]

        x_center = dec2sex(x_center / 15)
        y_center = dec2sex(y_center)

        print "\tImage Size:  %s x %s" % (nx, ny)
        print "\tCentered on: %s   %s" % (x_center, y_center)

        self.nx = nx
        self.ny = ny
        self.xo = x_center
        self.yo = y_center

        # Store wether center was done....
        self.centered = True

        return

    def get_FLXSCALE(self, magbase=30):

        print "# Computing FLXSCALE for magbase=%s" % magbase

        self.magbase = magbase
        self.flxscale = {}
        for filter in self.filters:

            self.flxscale[filter] = []

            for fname in self.files[filter]:
                header = getheader(fname)
                zp = header['MAGZERO'] + 2.5 * log10(header['EXPTIME'])
                flxscale = 10.0**(0.4 * (magbase - zp))
                self.flxscale[filter].append(flxscale)

        return

    def swarp_files(self,
                    conf="SWarp-common.conf",
                    dryrun=None,
                    combtype="MEDIAN",
                    reSWarp=None,
                    filters=None):

        if not filters:
            filters = self.filters

        # Get the dither centroid
        if not self.centered:
            self.center_dither()

        self.combtype = combtype

        # Compare and get the effective/combined zeropt
        # self.get_zeropt() -- # We do it at SEx time instead

        self.get_FLXSCALE(magbase=30)

        # Keys to keep
        #keywords = "OBJECT,EXPTIME,AIRMASS,TIMESYS,DATE-OBS,TIME-OBS,OBSTYPE,OBSERVAT,TELESCOP,HA,ZD,DETECTOR,DARKTIME"
        keywords = "OBJECT,OBSTYPE"
        keys = keywords.split(",")

        common_conf = os.path.join(self.BCSPIPE, 'LIB/pars', conf)
        pars = {}
        pars["IMAGE_SIZE"] = "%s,%s" % (self.nx, self.ny)
        pars["CENTER_TYPE"] = "MANUAL"
        pars["CENTER"] = "%s,%s" % (self.xo, self.yo)
        pars["PIXEL_SCALE"] = self.pixscale
        pars["PIXELSCALE_TYPE"] = "MANUAL"
        #pars["FSCALE_KEYWORD"]  = "FLXSCALE"
        pars["COPY_KEYWORDS"] = keywords
        pars["COMBINE_TYPE"] = combtype
        pars["NTHREADS"] = "0"

        # The options
        opts = ""
        for param, value in pars.items():
            opts = opts + "-%s %s " % (param, value)

        self.combima = {}
        self.weightima = {}
        self.comb_exptime = {}
        cmd = ""
        for filter in filters:

            pars[""] = "@%s" % (self.flxscale[filter])

            outimage = "%s%s.fits" % (self.tilename, filter)
            # We only want the i-band weight to save space
            if filter == 'i':
                outweight = "%s%s_weight.fits" % (self.tilename, filter)
            else:
                outweight = "coadd.weight.fits"

            # Store the names
            self.combima[filter] = "%s%s" % (self.tilename, filter)
            self.weightima[filter] = "%s%s_weight.fits" % (self.tilename,
                                                           filter)

            filelist = string.join(self.files[filter])
            cmd = "swarp %s -c %s -IMAGEOUT_NAME %s -WEIGHTOUT_NAME %s " % (
                filelist, common_conf, outimage, outweight)
            cmd = cmd + " -FSCALE_DEFAULT %s " % (
                ",".join(map(str, self.flxscale[filter])))
            cmd = cmd + opts

            if not dryrun:
                print "# Will run:\n\t%s" % cmd
                os.system(cmd)
            else:
                print cmd

        return

    # Put correction
    def header_FSCALE(self, filters=None, dm=0.05):

        self.do_level = {}
        self.filters_level = []

        if not filters:
            filters = self.filters

        for filter in filters:

            if len(self.files[filter]) < 2:
                print "Only one frame, no FLUX SCALE for %s %s filter" % (
                    self.tilename, filter)
                continue

            # First Check if we should correct
            self.do_level[filter] = False
            for file in self.files[filter]:

                if abs(2.5 * log10(self.FS_MAG[file])) >= 0.05:
                    #if self.MCorr[file] > dm:
                    self.do_level[filter] = True
                    self.filters_level.append(filter)

                    # Skip when values are < dm
            if not self.do_level[filter]:
                print >> sys.stderr, "# No need to correct frames %s for %s" % (
                    self.tilename, filter)
                continue

            # Other wise call put_FSCALE
            for file in self.files[filter]:
                self.put_FSCALE(file)

        return

    def swarp_weight(self,
                     conf="SWarp-common.conf",
                     dryrun=None,
                     combtype="SUM",
                     filters=None):

        if dryrun:
            print "Skipping weight generation for %s" % self.tilename
            return

        # Get the dither centroid
        if not self.centered:
            self.center_dither()

        common_conf = os.path.join(self.BCSPIPE, 'LIB/pars', conf)
        pars = {}
        pars["IMAGE_SIZE"] = "%s,%s" % (self.nx, self.ny)
        pars["CENTER_TYPE"] = "MANUAL"
        pars["CENTER"] = "%s,%s" % (self.xo, self.yo)
        pars["PIXEL_SCALE"] = self.pixscale
        pars["PIXELSCALE_TYPE"] = "MANUAL"
        pars["FSCALE_KEYWORD"] = "NONE"
        pars["FSCALASTRO_TYPE"] = "NONE"
        pars["COMBINE_TYPE"] = combtype
        pars['SUBTRACT_BACK'] = "N"

        # The options
        opts = ""
        for param, value in pars.items():
            opts = opts + "-%s %s " % (param, value)
        cmd = ""

        # Only on selected filters
        if not filters:
            filters = self.filters

        for filter in filters:

            outimage = "%s%s_nweight.fits" % (self.tilename, filter)
            outweight = "weight.fits"

            outfiles = []
            for file in self.files[filter]:
                outfits = "rep_" + os.path.basename(file)
                replace_vals_image(file, outfits, repval=1)
                outfiles.append(outfits)

            filelist = string.join(outfiles)
            cmd = "swarp %s -c %s -IMAGEOUT_NAME %s -WEIGHTOUT_NAME %s " % (
                filelist, common_conf, outimage, outweight)
            cmd = cmd + opts

            if not dryrun:
                print cmd
                os.system(cmd)
                chtype_fits(outimage, type='UInt8', verb=self.verb)
            else:
                print cmd

            # Clean up files
            os.system("rm %s" % outweight)
            for file in outfiles:
                print "Removing %s" % file
                os.system("rm %s" % file)

        return

    def swarp_BPM(self,
                  conf="SWarp-common.conf",
                  dryrun=None,
                  combtype="SUM",
                  filters=None):

        # Get the dither centroid
        if not self.centered:
            self.center_dither()

        common_conf = os.path.join(self.BCSPIPE, 'LIB/pars', conf)
        pars = {}
        pars["IMAGE_SIZE"] = "%s,%s" % (self.nx, self.ny)
        pars["CENTER_TYPE"] = "MANUAL"
        pars["CENTER"] = "%s,%s" % (self.xo, self.yo)
        pars["PIXEL_SCALE"] = self.pixscale
        pars["PIXELSCALE_TYPE"] = "MANUAL"
        pars["FSCALE_KEYWORD"] = "FLXSCALE"
        pars["FSCALASTRO_TYPE"] = "NONE"
        pars["COMBINE_TYPE"] = combtype
        pars['SUBTRACT_BACK'] = "N"

        # The options
        opts = ""
        for param, value in pars.items():
            opts = opts + "-%s %s " % (param, value)
        cmd = ""

        # Only on selected filters
        if not filters:
            filters = self.filters

        for filter in filters:

            outimage = "%s%s_bpmask.fits" % (self.tilename, filter)
            outweight = "weight.fits"

            bpmfiles = []
            for file in self.files[filter]:
                bpm_mef = pl2fits(file, N=16, verb=self.verb)
                bpmfiles.append(bpm_mef)

            filelist = string.join(bpmfiles)
            cmd = "swarp %s -c %s -IMAGEOUT_NAME %s -WEIGHTOUT_NAME %s " % (
                filelist, common_conf, outimage, outweight)
            cmd = cmd + opts
            os.system(cmd)
            chtype_replace_nonzero(outimage, repval=1, verb=self.verb)

            # Clean up files
            os.system("rm %s" % outweight)
            for file in bpmfiles:
                print "Removing %s" % file
                os.system("rm %s" % file)

        return

        # Generate the mask from the weight
    def generate_masks(self, filters, dryrun):

        self.mask = {}
        for filter in filters:

            self.mask[filter] = "%s%s_mask.fits" % (self.tilename, filter)
            weight = self.weightima[filter]
            mask = self.mask[filter]

            if dryrun:
                print >> sys.stderr, "# Skipping Mask Generation..."
                continue

            print >> sys.stderr, "# Generating mask image from %s --> %s" % (
                weight, mask)
            mask_from_weight(weight, mask, value=0)

        return

    # Make the detection Image using the mask
    def makeDetectionIma(self, filter='i'):

        mask = self.mask[filter]
        image = self.combima[filter] + ".fits"
        DetImage = "%s_detection.fits" % self.tilename

        print >> sys.stderr, "# Generating Detection image for filter:%s" % filter
        print >> sys.stderr, "# \t\t%s x %s --> %s" % (image, mask, DetImage)

        # Read in the mask and the data
        m_data, m_hdr = pyfits.getdata(mask, header=True)
        i_data, i_hdr = pyfits.getdata(image, header=True)

        # Multiply to get the detection
        d_data = m_data * i_data

        # Make a fits file out of it
        newfits = pyfits.HDUList()
        hdu = pyfits.PrimaryHDU()
        hdu.data = d_data
        hdu.header = i_hdr

        # Remove old version of file before
        if os.path.isfile(DetImage):
            os.remove(DetImage)

        newfits.append(hdu)
        newfits.writeto(DetImage)
        newfits.close

        # Make it visible outside
        self.DetImage = DetImage

        return

    # Run SExtractor
    def SEx(self, det_filter='i'):

        check_exe("sex")

        # Compare and get the effective/combined zeropt
        #self.compare_zeropts()
        #self.get_zeropt() # Maybe  we can run this at swarp time ?

        # list of output catalog names
        self.combcat = {}
        self.SExinpar = os.path.join(os.environ['BCSPIPE'],
                                     'LIB/pars/bcs_Catalog.inpar')

        # The detection image that we'll use
        if self.DetImage:
            det_ima = self.DetImage
            print >> sys.stderr, "# Will use default Detection Image:%s " % self.DetImage
        else:
            det_ima = self.combima[det_filter] + ".fits"
            print >> sys.stderr, "# Will use %s band Detection Image:%s " % (
                det_filter, det_ima)

        self.getbpz = 1  # This var is not really used...

        for filter in self.filters:

            self.combcat[filter] = self.combima[filter] + ".cat"
            input = self.combima[filter] + ".fits"
            output = self.combcat[filter]

            opts = ''
            if filter == 'i':
                backima = self.combima[filter] + "_BACK.fits"
                opts = " -CHECKIMAGE_TYPE BACKGROUND_RMS -CHECKIMAGE_NAME %s " % backima

            #opts = opts + ' -GAIN %s ' % header['GAIN']
            opts = opts + ' -WEIGHT_TYPE MAP_WEIGHT,BACKGROUND '  # weight map for detection, and BACKGROUND for meassurement
            opts = opts + ' -WEIGHT_IMAGE %s%s_weight.fits ' % (self.tilename,
                                                                det_filter)

            # Do the SEx
            cmd = "sex %s,%s -CATALOG_NAME %s -MAG_ZEROPOINT %s -c %s %s 1>&2" % (
                det_ima, input, output, self.magbase, self.SExinpar, opts)

            print cmd
            if not self.dryrun:
                os.system(cmd)

        return

    # Build th color catalog to use when computing the photo-z
    # Adapted from JHU APSIS pipeline
    def BuildColorCat(self):

        # Change accordingly
        zp_error = 0.05

        # The default output names
        self.colorCat = self.tilename + ".color"
        self.columnsFile = self.tilename + ".columns"

        print >> sys.stderr, 'Processing catalogs... for: ', self.tilename

        flux = {}
        fluxerr = {}

        m = {}
        em = {}

        # Get the detection catalog required columns
        outColumns = ['NUMBER', 'X_IMAGE', 'Y_IMAGE']
        detCatalog = self.combcat['i']
        detcols = SEx_head(detCatalog, verb=None)
        detectionList = []
        for key in outColumns:
            detectionList.append(detcols[key])
        detectionColumns = tuple(
            detectionList)  # the get_data function requires a tuple
        detection_variables = tableio.get_data(detCatalog, detectionColumns)

        # Read in the MAG_ISO and MAG_ISOERR from each catalog
        for filter in self.filters:

            # Get the columns
            sexcols = SEx_head(self.combcat[filter], verb=None)

            ## Info for flux columns
            fluxList = []
            fluxList.append(sexcols['FLUX_ISO'])
            fluxList.append(sexcols['FLUXERR_ISO'])
            fluxColumns = tuple(
                fluxList)  # the get_data function interface requires a tuple

            # Get the array using tableio
            flux[filter], fluxerr[filter] = tableio.get_data(
                self.combcat[filter], fluxColumns)
            m[filter] = flux[filter] * 0.0
            em[filter] = flux[filter] * 0.0

            # Fix the NAN values
            flux[filter] = deNAN(flux[filter])

            # Those objects with flux equal or less than 0 are assigned a magnitude of 99
            # and a limiting magnitude equal to their SExtractor photometric error. This
            # is interpreted by BPZ as a nondetection with zero flux and 1-sigma error
            # equal to the limiting magnitude

            nondetected = Numeric.less_equal(
                flux[filter], 0.0) * Numeric.greater(fluxerr[filter], 0.0)

            # Those objects with error flux and flux equal to 0 are assigned a magnitude of -99
            # and a flux of 0, which is interpreted by SExtractor as a non-observed object

            nonobserved = Numeric.less_equal(fluxerr[filter], 0.0)

            # When flux error > 100*(flux), mark as nonobserved (Benitez, 24-Oct-03).

            # Fix for fc11 -- y[:] has change meaning
            #nonobserved = Numeric.where(fluxerr[filter] > 100*(abs(flux[filter])),1.0,nonobserved[:])
            nonobserved = Numeric.where(fluxerr[filter] > 100 *
                                        (abs(flux[filter])), 1.0, nonobserved)

            detected = Numeric.logical_not(nonobserved + nondetected)

            # Get the zero point for the final magnitudes
            zpoint = self.magbase

            print filter, zpoint

            flux[filter] = Numeric.clip(flux[filter], 1e-100, 1e100)
            m[filter] = Numeric.where(
                detected, -2.5 * Numeric.log10(abs(flux[filter])) + zpoint -
                self.XCorr[filter], m[filter])
            m[filter] = Numeric.where(nondetected, 99.0, m[filter])
            m[filter] = Numeric.where(nonobserved, -99.0, m[filter])

            # clip values from being too small or large, i.e. 0 or inf.
            fluxerr[filter] = Numeric.clip(fluxerr[filter], 1e-100, 1e100)
            em[filter] = Numeric.where(
                detected,
                2.5 * Numeric.log10(1.0 + abs(fluxerr[filter] / flux[filter]))
                + self.XCorrError[filter], em[filter])
            em[filter] = Numeric.where(
                nondetected,
                2.5 * Numeric.log10(abs(fluxerr[filter])) - zpoint, em[filter])
            em[filter] = Numeric.where(nonobserved, 0.0, em[filter])

            #outColumns.append(filter +'_SDSS_MAG_ISO')
            #outColumns.append(filter +'_SDSS_MAGERR_ISO')
            outColumns.append(filter + '_MOSAICII_MAG_ISO')
            outColumns.append(filter + '_MOSAICII_MAGERR_ISO')

        # Prepare the header
        header = \
               '## ' + time.ctime() + '\n'+\
               '## BPZ Catalog file for Observation: ' + self.tilename + '\n'+\
               '## (This file was generated automatically by the BCS Rutgers pipeline)\n##\n'
        for i in range(len(outColumns)):
            header = header + '# ' + str(i + 1) + '\t' + outColumns[i] + '\n'

            # Prepare the data
        vars = list(detection_variables)
        for filter in self.filters:
            vars.append(m[filter])
            vars.append(em[filter])

        variables = tuple(vars)
        format = '%i\t %10.2f %10.2f' + '%10.4f  ' * (len(variables) - 3)
        print >> sys.stderr, 'Writing data to multicolor catalog...'
        tableio.put_data(self.colorCat,
                         variables,
                         header=header,
                         format=format,
                         append='no')
        print >> sys.stderr, 'Multicolor catalog complete.'

        # And now write .columns file
        cfile = open(self.columnsFile, 'w')
        cfile.write('## ' + time.ctime() + '\n')
        cfile.write('## ' + 'BPZ' + ' .columns file for Observation: ' +
                    self.tilename + '\n')
        cfile.write(
            '## (This file was generated automatically by the BCS Rutgers pipeline)\n##\n')

        i = len(detection_variables)
        for filter in self.filters:

            if filter == 'i':
                n_mo = str(i + 1)
            colmag = i + 1
            colmagerr = i + 2
            cfile.write('%s_MOSAICII\t %s,%s\t AB\t %.2f\t 0.0\n' %
                        (filter, i + 1, i + 2, zp_error))
            i = i + 2

        cfile.write('M_0\t%s\n' % n_mo)
        cfile.close()
        return

    # Run Benitez BPZ
    def runBPZ(self):
        """Runs BPZ on the multicolor catalog file using the .columns """

        print >> sys.stderr, 'Starting photometric redshift determination...'
        bpz = os.path.join(os.environ['P_BPZPATH'], 'bpz.py ')
        bpzcat = self.tilename + ".bpz"
        bpzprobs = self.tilename + ".probs"

        #cmd = 'python ' + bpz + self.colorCat + ' -ZMAX 6.0 -VERBOSE no -INTERP 2 -DZ 0.005 -SPECTRA CWWSB_Benitez2003.list'
        #cmd = 'python ' + bpz + self.colorCat + ' -ZMAX 6.0 -VERBOSE no -INTERP 0 -DZ 0.01 -SPECTRA CWWSB_Benitez2003.list'
        #cmd = 'python ' + bpz + self.colorCat + ' -ZMAX 1.8 -VERBOSE no -INTERP 2 -DZ 0.01 -SPECTRA CWWSB_Benitez2003.list -PRIOR hdfn'
        cmd = 'python ' + bpz + self.colorCat + ' -ZMAX 1.8 -VERBOSE no -INTERP 2 -DZ 0.01 -SPECTRA CWWSB_Benitez2003.list -PRIOR full -PROBS_LITE ' + bpzprobs

        if not self.dryrun:
            print cmd
            print >> sys.stderr, "Running full prior"
            os.system(cmd)
            print >> sys.stderr, "Photo-z ready"
        else:
            print cmd
        return

    def mean_airmass(self, filter):

        x = []
        for file in self.files[filter]:
            #print file,self.airmass[file]
            x.append(self.airmass[file])
        x = numpy.asarray(x)
        return x.mean()


# Create a mask file from the weights
def mask_from_weight(infile, outfile, value=0):
    import pyfits

    (data, hdr) = rfits(infile)

    newdata = numpy.where(data > 0, 1, 0)

    newfits = pyfits.HDUList()
    hdu = pyfits.PrimaryHDU()
    #hdu.data   = newdata.astype("Int16") # Just as integer
    hdu.data = newdata  #.astype("UInt8") # Just as ushort integer
    hdu.header = hdr

    # Remove old version of file before
    if os.path.isfile(outfile):
        os.remove(outfile)

    newfits.append(hdu)
    newfits.writeto(outfile)
    newfits.close
    return


def replace_vals_image(infits, outfits, repval=1):
    # Replace all values in MOSAIC multi-extension fits file by
    # repvalue

    import pyfits

    hdulist = pyfits.open(infits)
    Nima = len(hdulist)
    print >> sys.stderr, "Replacing with n=%s on %s images, %s --> %s" % (
        repval, Nima - 1, infits, outfits)
    for i in range(Nima)[1:]:
        hdulist[i].data = (hdulist[i].data * 0 + 1).astype('UInt8')
    hdulist.verify('fix')
    hdulist.writeto(outfits)
    hdulist.close()
    return


def add_airmass(filename, airmass):

    import pyfits
    c = pyfits.Card("AIRMASS", airmass, "Computed Stacked Mean Airmass")
    f = pyfits.open(filename, mode="update")
    #f[0].header.ascardlist().append(c)
    f[0].header.update(c.key, c.value, c.comment)
    f.close()


def put_exptime(file, exptime):

    import pyfits

    # Open the file, update mode
    f = pyfits.open(file, mode="update")  # open a FITS file
    hdr = f[0].header  # the primary HDU header

    print >> sys.stderr, "# Updating %s with EXPTIME=%s after SWarp" % (
        file, exptime)
    c = pyfits.Card("EXPTIME", exptime, "After SWarp equivalent exptime (sec)")
    c.verify()
    hdr.update('EXPTIME', c.value, c.comment)  #,after=after)

    # Close the file
    f.verify('fix')
    f.close()
    return


def put_scaled_head(file):

    import pyfits

    # Open the file, update mode
    f = pyfits.open(file, mode="update")  # open a FITS file
    hdr = f[0].header  # the primary HDU header

    #print >>sys.stderr,"# Updating %s with EXPTIME=%s after SWarp" % (file,exptime)
    c = pyfits.Card("FSCALED", True,
                    "Has the flux been scaled in individual frames")
    c.verify()
    hdr.update('FSCALED', c.value, c.comment)  #,after=after)

    # Close the file
    f.verify('fix')
    f.close()
    return


def put_zeropt(file, zeropt, photo='yes'):
    import pyfits

    # Open thr child info, update mode
    f = pyfits.open(file, mode="update")  # open a FITS file
    hdr = f[0].header  # the primary HDU header

    print >> sys.stderr, "# Updating %s with ZP=%s as in SExtractor, photo:%s" % (
        file, zeropt, photo)

    after = "DATE"

    c = {}
    if photo == 'yes':
        c['ZEROPT'] = pyfits.Card("ZEROPT", zeropt,
                                  "Computed ZEROPOINT AB mags/sec")
        c['PHOTOM'] = pyfits.Card("PHOTOM", 1,
                                  "Photometric quality 1=photo/0=not")
    else:
        c['ZEROPT'] = pyfits.Card("ZEROPT", zeropt,
                                  "Non-photo ZEROPOINT AB mags/sec")
        c['PHOTOM'] = pyfits.Card("PHOTOM", 0,
                                  "Photometric quality 1=photo/0=not")

    for key in c.keys():
        c[key].verify()
        hdr.update(key, c[key].value, c[key].comment, after=after)

    # Close the file
    #f.verify('silentfix')
    f.verify('fix')
    #f.verify('ignore')
    f.close()
    return


def chtype_replace_nonzero(infits, repval=1, out=None, verb=None):

    # Replace all values in !=0 in fits file to 1
    # certain value file by repvalue and changes the type to UInt8

    import pyfits
    import numpy

    if verb:
        print >> sys.stderr, "\tChange type UInt8 on %s --- Replacing pixels not eq 0  --> %s" % (
            infits, repval)

    hdulist = pyfits.open(infits, mode='update')
    Nima = len(hdulist)
    if Nima == 1:
        ima = hdulist[0].data
        ima = numpy.where(ima != 0, repval, ima).astype('UInt8')
        hdulist[0].data = ima
    else:
        for i in range(Nima)[1:]:
            ima = hdulist[i].data
            ima = numpy.where(ima != 0, repval, ima).astype('UInt8')
            hdulist[i].data = ima
    hdulist.verify('silentfix')
    hdulist.flush()
    hdulist.close()
    return


def chtype_fits(infits, type='UInt8', out=None, verb=None):

    # Change type of fits file
    import pyfits
    import numpy

    if verb:
        print >> sys.stderr, "\tChange type on %s to ---> %s" % (infits, type)

    hdulist = pyfits.open(infits, mode='update')
    Nima = len(hdulist)
    if Nima == 1:
        ima = (hdulist[0].data).astype(type)
        hdulist[0].data = ima
    else:
        for i in range(Nima)[1:]:
            ima = (hdulist[i].data).astype(type)
            hdulist[i].data = ima
    hdulist.verify('silentfix')
    hdulist.flush()
    hdulist.close()
    return


# Clean up the FLXCORR key in the fits header
def clean_FLXCORR(file, N=16):

    from pyraf import iraf
    from pyfits import getheader

    print "Analizing %s ..." % file
    header = getheader(file, 1)

    try:
        FLXCORR = header['FLXCORR']
        print "Cleaning up %s" % file
        for i in range(N):
            n = i + 1
            image = "%s[%s]" % (file, n)
            iraf.hedit(image,
                       'FLXCORR',
                       delete='yes',
                       verify='no',
                       update='yes',
                       show='yes')
            # Remove key from [0] extension
            #print image
    except:
        print "NO FLXCORR key found"

    return


# Checks if program is in PATH
def check_exe(exe, verb="yes"):

    path = os.environ['PATH'].split(':')
    for p in path:

        f = os.path.join(p, exe)
        if os.path.isfile(f):
            if verb:
                print >> sys.stderr, "# Found %s in %s" % (exe, f)
            return f
    print >> sys.stderr, "# ERROR: Couldn't find %s" % exe
    return


def elapsed_time(t1, text=''):
    import time
    t2 = time.time()
    hh = int((t2 - t1) / 3600.)
    mm = int((t2 - t1) / 60 - hh * 60)
    ss = (t2 - t1) - 60 * mm - 3600 * hh
    print "Elapsed time: %dh %dm %2.2fs %s" % (hh, mm, ss, text)
    print >> sys.stderr, "Elapsed time: %dh %dm %2.2fs %s" % (hh, mm, ss, text)
    #print >>sys.stderr,"Elapsed time: %dm %2.2fs" % ( int( (t2-t1)/60.), (t2-t1) - 60*int((t2-t1)/60.))
    return


def elapsed_time_str(t1):
    import time
    t2 = time.time()
    hh = int((t2 - t1) / 3600.)
    mm = int((t2 - t1) / 60 - hh * 60)
    ss = (t2 - t1) - 60 * mm - 3600 * hh
    return "Elapsed time: %dh %dm %2.4fs" % (hh, mm, ss)


# Parses the header of a SExtractor catalog and stores the column
# number (in python' order starting from zero) for every keyword columns
# It returs a dictionary with the column number for evry keyword.
def SEx_head(catalog, verb='yes'):

    if verb:
        print >> sys.stderr, "\r Parsing SEx head for:", catalog

    # Dictionary with column numbers
    SExcols = {}

    # Read the SExtractor catalog
    for line in open(catalog).readlines():

        if line[0] <> '#':
            break

        if line[:2] == "##":
            continue

        try:
            line = line.strip()
            vals = line.split()
            col = vals[1]
            key = vals[2]
            SExcols[key] = int(col) - 1
            if verb:
                print >> sys.stderr, "# %-20s %s" % (key, SExcols[key] + 1)
        except:
            continue

    return SExcols

########################################
# Read and write fits simple functions
#######################################


# To read in simple fits file
def rfits(filename):
    import pyfits
    import numpy

    ff = pyfits.open(filename, "readonly")
    data = ff[0].data
    hdr = ff[0].header
    ff.close()
    return data, hdr


# To write in a simple file
def writefits(array, filename):
    import numpy
    import pyfits, os
    # Writes fits files of a given array
    newfits = pyfits.HDUList()
    hdu = pyfits.PrimaryHDU()
    hdu.data = array
    newfits.append(hdu)

    # Remove old version of file before
    if os.path.isfile(filename):
        os.remove(filename)

    newfits.writeto(filename)
    newfits.close
    return


wfits = writefits


def cmdline():

    from optparse import OptionParser

    # Read in the command line options
    USAGE = " usage:\t %prog <assoc> <in_datapath> <out_datapath> [options] \n i.e.: %prog ~/bcs-inventory/BCS0519-5448.assoc ~/PROC ~/BCS/PROC"
    parser = OptionParser(usage=USAGE)

    parser.add_option("--noCopy",
                      action="store_true",
                      dest="noCopy",
                      default=0,
                      help="No Copy files")

    parser.add_option("--noSWarp",
                      action="store_true",
                      dest="noSWarp",
                      default=0,
                      help="No SWarp")

    parser.add_option("--noBPZ",
                      action="store_true",
                      dest="noBPZ",
                      default=0,
                      help="No BPZ")

    parser.add_option("--noSEx",
                      action="store_true",
                      dest="noSEx",
                      default=0,
                      help="No SEx")

    parser.add_option("--Dust",
                      action="store_true",
                      dest="Dust",
                      default=0,
                      help="Dust Extinction Correction")

    parser.add_option("--WeightOnly",
                      action="store_true",
                      dest="WeightOnly",
                      default=0,
                      help="Create custom weights from images and exits")

    parser.add_option("--Weight",
                      action="store_true",
                      dest="Weight",
                      default=0,
                      help="Create custom weights from images")

    parser.add_option("--noMask",
                      action="store_true",
                      dest="noMask",
                      default=0,
                      help="Do not generate the mask from weight")

    parser.add_option("--useMask",
                      action="store_true",
                      dest="useMask",
                      default=0,
                      help="Use Trimming mask")

    parser.add_option("--dryrun",
                      action="store_true",
                      dest="dryrun",
                      default=0,
                      help="Dry Run (only SExtractor remains)")

    parser.add_option("--combtype",
                      dest="combtype",
                      default='MEDIAN',
                      help="SWarp COMBINE TYPE")

    (options, args) = parser.parse_args()

    if len(args) < 3:
        parser.error(USAGE + "\nMust supply at least one argument required")

    # Dry run turns off everything... almost
    if options.dryrun:
        options.noCopy = 1
        options.noSWarp = 1
        options.noBPZ = 1
        options.noMask = 1
        options.Dust = 0

    # Same for WeightOnly
    if options.WeightOnly:
        options.noCopy = 1
        options.noSWarp = 1
        options.noBPZ = 1
        options.noMask = 1
        options.Dust = 0

    return options, args


# Transform form Dec to Sexagesimal hh:mm:ss format
def dec2sex(hms):

    hh = int(hms)
    mm = int((hms - hh) * 60)
    ss = ((hms - hh) * 60 - mm) * 60

    if abs(mm) < 10:
        mm = "0%d" % abs(mm)
    else:
        mm = "%2d" % abs(mm)

    if abs(ss) < 10:
        ss = "0%.2f" % abs(ss)
    else:
        ss = "%5.2f" % abs(ss)

        #return "%2d:%d:%.2f" % (hh,abs(mm),abs(ss))
    return "%s:%s:%s" % (hh, mm, ss)

# Taken from APSIS fUtil
import Numeric


def deNAN(a, value=0.0):
    nans = Numeric.logical_not(Numeric.less(a, 0.0) + Numeric.greater_equal(
        a, 0.0))
    return Numeric.where(nans, value, a)


def main():

    # The start time
    tstart = time.time()

    # Get the command line options
    opt, arg = cmdline()

    assocfile = arg[0]
    inpath = arg[1]
    outpath = arg[2]

    # Read in the command line options
    #USAGE = " usage:\t"+os.path.split(sys.argv[0])[1]+" <assoc> <in_datapath> <out_datapath> [dryrun]"
    # Example
    # ./combcat.py ~/bcs-inventory/BCS0519-5448.assoc ~/PROC ~/BCS/PROC

    # Init the class
    c = combcat(assocfile,
                datapath=inpath,
                outpath=outpath,
                verb='yes',
                dryrun=opt.dryrun)

    # Copy the files into dir
    if opt.noCopy:
        c.copyfiles(copy=None)
    else:
        c.copyfiles(copy='yes')

    # SWarp
    c.swarp_files(dryrun=opt.noSWarp,
                  conf="SWarp-common.conf",
                  combtype=opt.combtype)  #,filters=opt.frames_filters)

    try:
        c.makeDetectionIma(filter='i')
    except:
        print "# Could not make Detection image with i-band"

    # Make the detection image
    if opt.useMask:
        # Make the mask from the weights
        c.generate_masks(filters=('i', ), dryrun=opt.noMask)

    # SExtractor
    if not opt.noSEx:
        c.SEx()

        # Dust Extinction Correction
        if opt.Dust:
            c.DustCorrection()
            #sys.exit()

            # Cats if allowed
        if c.getbpz and not opt.noBPZ:
            c.BuildColorCat()
            c.runBPZ()

    elapsed_time(tstart, c.tilename)
    return


main()
