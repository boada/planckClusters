#!/usr/bin/env python3

import os
import sys
import tableio
import numpy as np
from math import log10
import time
import string
from astropy.io import fits
from astropy.io.fits import getheader

class combcat:
    ''' Combine, swarp and get catalogs '''

    def __init__(self,
                 assocfile,
                 datapath='',
                 outpath='',
                 pixscale=0.2666,
                 dryrun=None,
                 noSWarp=False,
                 verb='yes'):

        self.assocfile = assocfile
        self.datapath = datapath
        self.tilename = os.path.basename(os.path.splitext(assocfile)[0])
        self.outpath = outpath
        self.dryrun = dryrun
        self.noSWarp = noSWarp
        self.verb = verb
        self.DetImage = None
        self.centered = None
        self.got_zeropt = False

        # Check for environ vars
        if not os.getenv('BCSPIPE'):
            os.environ['BCSPIPE'] = os.path.join(os.environ['HOME'],
                                                 'BCSPIPE')
        self.BCSPIPE = os.getenv('BCSPIPE')

        # Set the dust_directory
        #os.environ['DUST_DIR'] = os.path.join(self.BCSPIPE,"LIB")

        # Set the pixel scale
        self.pixscale = pixscale

        # Read the association file
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

        self.exptime = {}
        self.airmass = {}

        # Make image list per filter
        self.files = {}
        self.files_weight = {}
        # Exptimes per filter
        self.exptimes = {}

        # The the DQmaks files
        self.dqmask = {}

        print("# Will read %s" % self.assocfile)

        with open(self.assocfile, 'r') as assocfile:
            lines = assocfile.readlines()
            for line in lines:
                if line[0] == "#":
                    continue

                vals = line.split()
                fname = os.path.basename(vals[0])

                # hack to deal with compressed files in the association
                if '.fz' in fname:
                    fname = fname.rstrip('.fz')
                    if not os.path.isfile(fname) and not self.noSWarp:
                        os.system('funpack -v {}.fz'.format(fname))

                filtername = vals[1]
                self.exptime[fname] = float(vals[2])
                self.airmass[fname] = float(vals[3])

                # Figure out the dqmask
                if 'tu' in fname:
                    nid = int(os.path.splitext(fname)[0][2:]) + 1
                    ext = os.path.splitext(fname)[1]
                    pre = fname[0:2]
                    dqmask = "%s%s%s" % (pre, nid, ext)
                    if not os.path.isfile(dqmask) and not self.noSWarp:
                        os.system('funpack -v {}.fz'.format(dqmask))
                elif 'k4' in fname:
                    dqmask = fname.replace('opi', 'opd')
                    if not os.path.isfile(dqmask) and not self.noSWarp:
                        os.system('funpack -v {}.fz'.format(dqmask))

                # make a list of the file names
                self.filelist.append(fname)
                self.dqmask[fname] = dqmask
                # this is just a list of all files used
                self.infiles.append(fname)
                self.infiles.append(dqmask)

                # update the list of filters
                if filtername not in self.filters:
                    self.filters.append(vals[1])
                    self.files[filtername] = []
                    self.files_weight[filtername] = []
                    self.exptimes[filtername] = []

                # append the filename to the right filter
                self.files[filtername].append(fname)
                self.files_weight[filtername].append("%s.weight.fits'[0]'" %
                            os.path.splitext(fname)[0])
                self.exptimes[filtername].append(self.exptime[fname])

        return

    def copyfiles(self, copy=False):
        ''' Copy files from the remote location to the current dir'''

        if not copy:
            self.outdir = os.path.join(self.outpath, self.tilename)
            if not os.path.exists(self.outdir):
                os.mkdir(self.outdir)
            return

        # Fits make sure that the output dir exists
        self.outdir = os.path.join(self.outpath, self.tilename)
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        # Go to the directory
        os.chdir(self.outdir)

        print("We are in:", self.outdir, file=sys.stderr)

        for file in self.infiles:
            cmd1 = "rsync -avhP -e ssh %s ." % os.path.join(
                self.datapath, file)
            print(cmd1)
            if copy:
                os.system(cmd1)

        return

    def update_header_projection(self, proj='TAN'):
        for fname in self.infiles:
            print('# Updating {} '.format(fname), end='...')
            with fits.open(fname, mode='update') as hdulist:
                for i, hdu in enumerate(hdulist):
                    if not i:
                        continue
                    elif proj in hdu.header['CTYPE1']:
                        continue
                    else:
                        hdu.header['CTYPE1'] = 'RA---{}'.format(proj)
                        hdu.header['CTYPE2'] = 'DEC--{}'.format(proj)
            print(' to {}'.format(proj))
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

        cmd = "swarp %s -c %s " % (' '.join(self.filelist), center_conf)
        for param, value in list(opts.items()):
            cmd += "-%s %s " % (param, value)

        print(cmd)
        if not self.dryrun:
            print(os.getcwd())
            print("Centering mosaic", file=sys.stderr)
            os.system(cmd)

        # Read in the header
        print(opts['IMAGEOUT_NAME'])
        header = getheader(opts["IMAGEOUT_NAME"])
        nx = header["NAXIS1"]
        ny = header["NAXIS2"]
        x_center = header["CRVAL1"]
        y_center = header["CRVAL2"]

        x_center = dec2sex(x_center / 15)
        y_center = dec2sex(y_center)

        print("\tImage Size:  %s x %s" % (nx, ny))
        print("\tCentered on: %s   %s" % (x_center, y_center))

        self.nx = nx
        self.ny = ny
        self.xo = x_center
        self.yo = y_center

        # Store wether center was done....
        self.centered = True

        return

    def get_FLXSCALE(self, magbase=30):
        ''' I don't really know what this function is supposed to do. '''

        print("# Computing FLXSCALE for magbase=%s" % magbase)
        self.magbase = magbase
        self.flxscale = {}
        for filter in self.filters:

            self.flxscale[filter] = []

            for fname in self.files[filter]:
                header = getheader(fname)
                #zp = header['MAGZERO'] + 2.5*log10(header['EXPTIME'])
                zp = header['MAGZERO']
                flxscale = 10.0**(0.4 * (magbase - zp))
                self.flxscale[filter].append(flxscale)

        return

    def swarp_files(self,
                    conf="SWarp-common.conf",
                    dryrun=None,
                    combtype="MEDIAN",
                    reSWarp=None,
                    filters=None):
        ''' SWARPS all the files given in the association file together to make
        some moscaics.

        '''

        if not filters:
            filters = self.filters

        # update the projections
        #self.update_header_projection()

        # Get the dither centroid
        if not self.centered:
            self.center_dither()

        self.make_swarp_input_weights(clobber=False)
        self.combtype = combtype
        self.get_FLXSCALE(magbase=26)

        # Keys to keep
        keywords = ['OBJECT', 'EXPTIME', 'AIRMASS', 'TIMESYS', 'DATE-OBS',
                    'TIME-OBS', 'OBSTYPE', 'OBSERVAT', 'TELESCOP', 'HA', 'ZD',
                    'DETECTOR', 'DARKTIME', 'RA', 'DEC', 'MJD-OBS', 'INSTRUME',
                    'FILTER']

        pars = {}
        pars["IMAGE_SIZE"] = "%s,%s" % (self.nx, self.ny)
        pars["CENTER_TYPE"] = "MANUAL"
        pars["CENTER"] = "%s,%s" % (self.xo, self.yo)
        pars["PIXEL_SCALE"] = self.pixscale
        pars["PIXELSCALE_TYPE"] = "MANUAL"
        #pars["FSCALE_KEYWORD"]  = "FLXSCALE"
        pars["COPY_KEYWORDS"] = ','.join(keywords)
        pars["COMBINE_TYPE"] = combtype
        pars["NTHREADS"] = "0"
        pars["WEIGHT_TYPE"] = "MAP_WEIGHT"

        # The options
        opts = ""
        for param, value in list(pars.items()):
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

            filelist = ' '.join(self.files[filter])

            cmd = "swarp %s -IMAGEOUT_NAME %s -WEIGHTOUT_NAME %s " %\
                (filelist, outimage, outweight)
            cmd += " -WEIGHT_IMAGE %s " % (
                ",".join(self.files_weight[filter]))
            #cmd += " -WEIGHT_SUFFIX .weight.fits'[0]'"
            cmd = cmd + " -FSCALE_DEFAULT %s " % (
                ",".join(map(str, self.flxscale[filter])))
            cmd = cmd + opts

            if not dryrun:
                print("# Will run:\n\t%s" % cmd)
                os.system(cmd)
            else:
                print(cmd)

        return

    def swarp_newfirm(self):
        ''' This runs swarp on the stacked newfirm images to make sure they are
        in the same projection as the mosaic images. Right now, it needs to
        have the stacked data from the VO. In the future,  I might create my
        own stacks similar to the method above, but probably not. I think the
        stacked images from the VO are good enough to do what we want. '''

        from glob import glob

        # get the center file created above
        center = "{}/SWarp-{}-center.fits".format(self.outdir, self.tilename)

        # link this file into the newfirm directory and change the name
        # check first
        newfirm_dir = '../../newfirm/stacked/'
        if not os.path.isdir(newfirm_dir):
            return

        # make sure the images aren't compressed
        imgs = glob('{}*.fz'.format(newfirm_dir))
        print(imgs)
        if len(imgs) < 1:
            # check for uncompressed images
                imgs = glob('{}*.fits'.format(newfirm_dir))
                if len(imgs) < 1:
                    print('No NEWFIRM images found!')
                    return
        else:
            for img in imgs:
                with fits.open(img) as kimg:
                    try:
                        prod_type = kimg[0].header['prodtype']
                    except KeyError:
                        try:
                            prod_type = kimg[1].header['prodtype']
                        except KeyError:
                            print('Something is wrong with the images!')
                            return
                if 'image' in prod_type:
                    # use funpack to decompress items
                    check_exe('funpack')
                    os.system('funpack -v {}'.format(img))
                else:
                    continue

            try:
                relpath = os.path.relpath('./', newfirm_dir)
                print(relpath)
                os.symlink('{}/{}'.format(relpath, center),
                        '{}{}.head'.format(newfirm_dir, self.tilename))
            except FileExistsError:
                pass

            # build the swarp command
            check_exe('swarp')
            imgs = glob('{}*.fits'.format(newfirm_dir))
            for img in imgs:
                cmd = 'swarp {} '.format(img)
                cmd += '-IMAGEOUT_NAME {}{}k.fits '.format(newfirm_dir,
                                                        self.tilename)
                cmd += '-SUBTRACT_BACK N -WRITE_XML N'

                print(cmd)
                os.system(cmd)

                # clean up decompressed files
                os.remove('{}'.format(img))

        return

    # Put correction
    def header_FSCALE(self, filters=None, dm=0.05):

        self.do_level = {}
        self.filters_level = []

        if not filters:
            filters = self.filters

        for filter in filters:
            if len(self.files[filter]) < 2:
                print("Only one frame, no FLUX SCALE for %s %s filter" % (
                    self.tilename, filter))
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
                print("# No need to correct frames %s for %s" % (
                    self.tilename, filter), file=sys.stderr)
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
            print("Skipping weight generation for %s" % self.tilename)
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
        for param, value in list(pars.items()):
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
                print(cmd)
                os.system(cmd)
                chtype_fits(outimage, type='UInt8', verb=self.verb)
            else:
                print(cmd)

            # Clean up files
            os.system("rm %s" % outweight)
            for file in outfiles:
                print("Removing %s" % file)
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
                print("# Skipping Mask Generation...", file=sys.stderr)
                continue

            print("# Generating mask image from %s --> %s" % (
                weight, mask), file=sys.stderr)
            mask_from_weight(weight, mask, value=0)

        return

    # Generate the mask from the weight
    def make_swarp_input_weights(self, clobber=True):
        ''' Generate the masks from the weight images provided by the VO '''

        for fname, dqmask in list(self.dqmask.items()):
            outname = "%s.weight.fits" % os.path.splitext(fname)[0]
            print("# Generating weight image from %s --> %s" % (
                dqmask, outname), file=sys.stderr)
            weight_from_dqfile(dqmask, outname, clobber=clobber)
        return

    # Make the detection Image using the mask
    def makeDetectionIma(self, filter='i'):

        mask = self.mask[filter]
        image = self.combima[filter] + ".fits"
        DetImage = "%s_detection.fits" % self.tilename

        print("# Generating Detection image for filter:%s" % filter,
              file=sys.stderr)
        print("# \t\t%s x %s --> %s" % (image, mask, DetImage),
              file=sys.stderr)

        # Read in the mask and the data
        m_data, m_hdr = fits.getdata(mask, header=True)
        i_data, i_hdr = fits.getdata(image, header=True)

        # Multiply to get the detection
        d_data = m_data * i_data

        # Make a fits file out of it
        newfits = fits.HDUList()
        hdu = fits.PrimaryHDU()
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

    def get_zeropt(self):
        ''' This is going to call the photometrypipline script that lives in
        the projects directory. It should both astrometrically correct the
        mosaics and calculate the overall zeropoint of the mosaic. Currently,
        this has only been tested on the MOSAIC camera, and should NOT work
        with NEWFIRM, yet. It's on my TODO list.

        '''

        check_exe('pp_run')

        for filter in self.filters:
            mosaic = '{}.fits'.format(self.combima[filter])

            cmd = 'pp_run {}'.format(mosaic)

            print(cmd)
            if not self.dryrun:
                pass
                os.system(cmd)

            if os.path.isfile('photometry_control_star.dat'):
                os.rename('photometry_control_star.dat',
                          'photometry_control_star_{}.dat'.format(filter))
            else:
                print('Photometric calibration failed')

            # correct the header information to make sure floats are floats and
            # not strings
            with fits.open('{}.fits'.format(self.combima[filter]),
                                            mode='update') as f:
                header = f[0].header
                for key, val in list(header.items()):
                    if 'CD1_' in key or 'CD2_' in key or \
                        'CRVAL' in key or 'CRPIX' in key or \
                            'EQUINOX' in key:
                        f[0].header[key] = float(val)

        return

    # Run SExtractor
    def SEx(self, det_filter='i'):
        ''' Runs SEXTRACTOR on the mosaicked images created by swarp. It should
        use the zero point created by the photometrypipline.

        '''

        check_exe("sex")

        # list of output catalog names
        self.combcat = {}
        self.SExinpar = os.path.join(os.environ['BCSPIPE'],
                                     'LIB/pars/bcs_Catalog.inpar')

        # The detection image that we'll use
        if self.DetImage:
            det_ima = self.DetImage
            print("# Will use default Detection Image:%s " % self.DetImage,
                  file=sys.stderr)
        else:
            det_ima = self.combima[det_filter] + ".fits"
            print("# Will use %s band Detection Image:%s " % (
                det_filter, det_ima), file=sys.stderr)

        self.getbpz = 1  # This var is not really used...

        # make the zeropoint files with PP.
        self.get_zeropt()

        for filter in self.filters:

            self.combcat[filter] = self.combima[filter] + ".cat"
            input = self.combima[filter] + ".fits"
            output = self.combcat[filter]
            photocal = 'photometry_control_star_{}.dat'.format(filter)
            try:
                _tmp = np.genfromtxt(photocal, names=True, dtype=None)
                zeropt = _tmp['ZP']
            except OSError:
                print('WARNING!: Photometric calibration not set!')
                zeropt = self.magbase

            print(zeropt)

            opts = ''
            if filter == 'i':
                backima = self.combima[filter] + "_BACK.fits"
                opts = " -CHECKIMAGE_TYPE BACKGROUND_RMS"
                opts += " -CHECKIMAGE_NAME {} ".format(backima)

            # weight map for detection, and BACKGROUND for meassurement
            opts += ' -WEIGHT_TYPE MAP_WEIGHT,BACKGROUND '
            opts += ' -WEIGHT_IMAGE %s%s_weight.fits ' % (self.tilename,
                                                          det_filter)

            # Do the SEx
            cmd = "sex {},{} -CATALOG_NAME {}".format(det_ima, input, output)
            cmd += " -MAG_ZEROPOINT {} -c {} {}1>&2".format(zeropt,
                                                    self.SExinpar, opts)

            print(cmd)
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

        print('Processing catalogs... for: ', self.tilename, file=sys.stderr)

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
        detection_variables = tableio.get_data(detCatalog,
                                               detectionColumns)

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

            # Those objects with flux equal or less than 0 are assigned a
            # magnitude of 99 and a limiting magnitude equal to their
            # SExtractor photometric error. This is interpreted by BPZ as a
            # nondetection with zero flux and 1-sigma error equal to the
            # limiting magnitude

            nondetected = np.less_equal(
                flux[filter], 0.0) * np.greater(fluxerr[filter], 0.0)

            # Those objects with error flux and flux equal to 0 are assigned a
            # magnitude of -99
            # and a flux of 0, which is interpreted by SExtractor as a
            # non-observed object

            nonobserved = np.less_equal(fluxerr[filter], 0.0)

            # When flux error > 100*(flux), mark as nonobserved (Benitez,
            # 24-Oct-03).

            # Fix for fc11 -- y[:] has change meaning
            #nonobserved = np.where(fluxerr[filter] >
            #100*(abs(flux[filter])),1.0,nonobserved[:])
            nonobserved = np.where(fluxerr[filter] > 100 *
                                        (abs(flux[filter])), 1.0,
                                   nonobserved)

            detected = np.logical_not(nonobserved + nondetected)

            # Get the zero point for the final magnitudes
            zpoint = self.magbase

            print(filter, zpoint)

            flux[filter] = np.clip(flux[filter], 1e-100, 1e100)
            m[filter] = np.where(detected,
                                      -2.5 * np.log10(abs(flux[filter])) +
                                      zpoint - self.XCorr[filter],
                                 m[filter])
            m[filter] = np.where(nondetected, 99.0, m[filter])
            m[filter] = np.where(nonobserved, -99.0, m[filter])

            # clip values from being too small or large, i.e. 0 or inf.
            fluxerr[filter] = np.clip(fluxerr[filter], 1e-100, 1e100)
            em[filter] = np.where(
                detected,
                2.5 * np.log10(1.0 + abs(fluxerr[filter] / flux[filter])) +
                self.XCorrError[filter], em[filter])
            em[filter] = np.where(
                nondetected,
                2.5 * np.log10(abs(fluxerr[filter])) - zpoint, em[filter])
            em[filter] = np.where(nonobserved, 0.0, em[filter])

            #outColumns.append(filter +'_SDSS_MAG_ISO')
            #outColumns.append(filter +'_SDSS_MAGERR_ISO')
            outColumns.append(filter + '_MOSAICII_MAG_ISO')
            outColumns.append(filter + '_MOSAICII_MAGERR_ISO')

        # Prepare the header
        header = \
               '## ' + time.ctime() + '\n' + \
               '## BPZ Catalog file for Observation: ' + self.tilename + \
                '\n' + \
                '## (This file was generated automatically by' + \
                'the BCS Rutgers pipeline)\n##\n'
        for i in range(len(outColumns)):
            header = header + '# ' + str(i + 1) + '\t' + outColumns[i] + '\n'

            # Prepare the data
        vars = list(detection_variables)
        for filter in self.filters:
            vars.append(m[filter])
            vars.append(em[filter])

        variables = tuple(vars)
        format = '%i\t %10.2f %10.2f' + '%10.4f  ' * (len(variables) - 3)
        print('Writing data to multicolor catalog...', file=sys.stderr)
        tableio.put_data(self.colorCat,
                         variables,
                         header=header,
                         format=format,
                         append='no')
        print('Multicolor catalog complete.', file=sys.stderr)

        # And now write .columns file
        cfile = open(self.columnsFile, 'w')
        cfile.write('## ' + time.ctime() + '\n')
        cfile.write('## ' + 'BPZ' + ' .columns file for Observation: ' +
                    self.tilename + '\n')
        cfile.write(
            '## (This file was generated by the BCS Rutgers pipeline)\n##\n')

        i = len(detection_variables)
        for filter in self.filters:

            if filter == 'i':
                n_mo = str(i + 1)
            cfile.write('%s_MOSAICII\t %s,%s\t AB\t %.2f\t 0.0\n' %
                        (filter, i + 1, i + 2, zp_error))
            i = i + 2

        cfile.write('M_0\t%s\n' % n_mo)
        cfile.close()
        return

    # Run Benitez BPZ
    def runBPZ(self):
        """Runs BPZ on the multicolor catalog file using the .columns """

        print('Starting photometric redshift determination...',
              file=sys.stderr)
        bpz = os.path.join(os.environ['P_BPZPATH'], 'bpz.py ')
        #bpzcat = self.tilename + ".bpz"
        bpzprobs = self.tilename + ".probs"

        #cmd = 'python ' + bpz + self.colorCat + \
        #' -ZMAX 6.0 -VERBOSE no -INTERP 2 -DZ 0.005 -SPECTRA
        #CWWSB_Benitez2003.list'
        #cmd = 'python ' + bpz + self.colorCat + \
        #' -ZMAX 6.0 -VERBOSE no -INTERP 0 -DZ 0.01 -SPECTRA
        #CWWSB_Benitez2003.list'
        #cmd = 'python ' + bpz + self.colorCat + \
        #' -ZMAX 1.8 -VERBOSE no -INTERP 2 -DZ 0.01 -SPECTRA
        #CWWSB_Benitez2003.list -PRIOR hdfn'
        cmd = 'python ' + bpz + self.colorCat + \
            ' -ZMAX 1.8 -VERBOSE no -INTERP 2 -DZ 0.01 -SPECTRA ' + \
            'CWWSB_Benitez2003.list -PRIOR full -PROBS_LITE ' + bpzprobs

        if not self.dryrun:
            print(cmd)
            print("Running full prior", file=sys.stderr)
            os.system(cmd)
            print("Photo-z ready", file=sys.stderr)
        else:
            print(cmd)
        return

    def mean_airmass(self, filter):

        x = []
        for file in self.files[filter]:
            #print file,self.airmass[file]
            x.append(self.airmass[file])
        x = np.asarray(x)
        return x.mean()

    def cleanup_files(self):
        ''' Removes all of the intermediate files created during the
        processing. Basically that is every file in the association file.

        '''

        print('cleaning up files', end='.')
        for f in self.infiles:
                print('', end='.')
                os.remove(f)
        for filtername in self.filters:
            for f in self.files_weight[filtername]:
                print('', end='.')
                os.remove(f.rstrip("'[0]'"))
        print('')
        return

    def make_RGB(self, kband=False):

        try:
            check_exe('stiff')
        except FileNotFoundError:
            return

        # input files
        if kband:
            try:
                red = '../../newfirm/stacked/{}{}.fits'.format(self.tilename, 'k')
            except FileNotFoundError:
                print('k-band file not found, restoring defaults')
            green = './{}{}.fits'.format(self.tilename, 'i')
            blue = './{}{}.fits'.format(self.tilename, 'r')
        else:
            red = './{}{}.fits'.format(self.tilename, 'i')
            green = './{}{}.fits'.format(self.tilename, 'r')
            blue = './{}{}.fits'.format(self.tilename, 'g')

        # output file
        output = '{}.tiff'.format(self.tilename)

        # options
        opts = ['-MIN_LEVEL', '0.001',
                '-MAX_LEVEL', '0.999',
                '-MAX_TYPE', 'QUANTILE',
                '-DESCRIPTION', "'{} RGB'".format(self.tilename),
                '-WRITE_XML', 'N',
                '-COPYRIGHT', "'Steven Boada'"]

        # build the command -- a space is always last
        cmd = 'stiff {} {} {} '.format(red, green, blue)
        cmd += '-OUTFILE_NAME {} '.format(output)
        # append the options
        for opt in opts:
            cmd += '{} '.format(opt)

        print(cmd)
        os.system(cmd)
        return

# Create a mask file from the weights
def mask_from_weight(infile, outfile, value=0):

    (data, hdr) = rfits(infile)

    newdata = np.where(data > 0, 1, 0)

    newfits = fits.HDUList()
    hdu = fits.PrimaryHDU()
    #hdu.data   = newdata.astype("Int16") # Just as integer
    hdu.data = newdata.astype("UInt8")  # Just as ushort integer
    hdu.header = hdr

    # Remove old version of file before
    if os.path.isfile(outfile):
        os.remove(outfile)

    newfits.append(hdu)
    newfits.writeto(outfile)
    newfits.close
    return


# Create a mask file from the weights
def weight_from_dqfile(infile, outfile, clobber=False):
    # Remove old version of file before
    if os.path.isfile(outfile):
        if clobber:
            os.remove(outfile)
        else:
            print("# Skipping creation, %s image exists" % outfile)
            return

    inHDU = fits.open(infile, "readonly")
    data = inHDU[1].data
    newdata = np.where(data > 0, 0, 1)
    inHDU[1].data = newdata.astype("UInt8")  # Just as ushort integer
    newhdu = fits.HDUList()
    newhdu.append(inHDU[1])
    newhdu.writeto(outfile, clobber=clobber)
    newhdu.close()
    inHDU.close()
    return

def replace_vals_image(infits, outfits, repval=1):
    # Replace all values in MOSAIC multi-extension fits file by
    # repvalue

    hdulist = fits.open(infits)
    Nima = len(hdulist)
    print("Replacing with n=%s on %s images, %s --> %s" % (
        repval, Nima - 1, infits, outfits), file=sys.stderr)
    for i in range(Nima)[1:]:
        hdulist[i].data = (hdulist[i].data * 0 + 1).astype('UInt8')
    hdulist.verify('fix')
    hdulist.writeto(outfits)
    hdulist.close()
    return


def add_airmass(filename, airmass):

    c = fits.Card("AIRMASS", airmass, "Computed Stacked Mean Airmass")
    f = fits.open(filename, mode="update")
    #f[0].header.ascardlist().append(c)
    f[0].header.update(c.key, c.value, c.comment)
    f.close()


def put_exptime(file, exptime):

    # Open the file, update mode
    f = fits.open(file, mode="update")  # open a FITS file
    hdr = f[0].header  # the primary HDU header

    print("# Updating %s with EXPTIME=%s after SWarp" % (
        file, exptime), file=sys.stderr)
    c = fits.Card("EXPTIME", exptime, "After SWarp equivalent exptime (sec)")
    c.verify()
    hdr.update('EXPTIME', c.value, c.comment)  # ,after=after)

    # Close the file
    f.verify('fix')
    f.close()
    return


def put_scaled_head(file):

    # Open the file, update mode
    f = fits.open(file, mode="update")  # open a FITS file
    hdr = f[0].header  # the primary HDU header

    #print >>sys.stderr,"# Updating %s with EXPTIME=%s after SWarp" %
    #(file,exptime)
    c = fits.Card("FSCALED", True,
                    "Has the flux been scaled in individual frames")
    c.verify()
    hdr.update('FSCALED', c.value, c.comment)  # ,after=after)

    # Close the file
    f.verify('fix')
    f.close()
    return


def put_zeropt(file, zeropt, photo='yes'):

    # Open thr child info, update mode
    f = fits.open(file, mode="update")  # open a FITS file
    hdr = f[0].header  # the primary HDU header

    print("# Updating %s with ZP=%s as in SExtractor, photo:%s" % (
        file, zeropt, photo), file=sys.stderr)

    after = "DATE"

    c = {}
    if photo == 'yes':
        c['ZEROPT'] = fits.Card("ZEROPT", zeropt,
                                  "Computed ZEROPOINT AB mags/sec")
        c['PHOTOM'] = fits.Card("PHOTOM", 1,
                                  "Photometric quality 1=photo/0=not")
    else:
        c['ZEROPT'] = fits.Card("ZEROPT", zeropt,
                                  "Non-photo ZEROPOINT AB mags/sec")
        c['PHOTOM'] = fits.Card("PHOTOM", 0,
                                  "Photometric quality 1=photo/0=not")

    for key in list(c.keys()):
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

    if verb:
        print("\tChange type UInt8 on %s -- Replacing pixels not eq 0 --> %s"
              % (infits, repval), file=sys.stderr)

    hdulist = fits.open(infits, mode='update')
    Nima = len(hdulist)
    if Nima == 1:
        ima = hdulist[0].data
        ima = np.where(ima != 0, repval, ima).astype('UInt8')
        hdulist[0].data = ima
    else:
        for i in range(Nima)[1:]:
            ima = hdulist[i].data
            ima = np.where(ima != 0, repval, ima).astype('UInt8')
            hdulist[i].data = ima
    hdulist.verify('silentfix')
    hdulist.flush()
    hdulist.close()
    return


def chtype_fits(infits, type='UInt8', out=None, verb=None):

    if verb:
        print("\tChange type on %s to ---> %s" % (infits, type),
              file=sys.stderr)

    hdulist = fits.open(infits, mode='update')
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

    print("Analizing %s ..." % file)
    #header = getheader(file, 1)

    try:
        #FLXCORR = header['FLXCORR']
        print("Cleaning up %s" % file)
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
        print("NO FLXCORR key found")

    return

def check_exe(exe, verb="yes"):
    ''' Checks to make sure we have the appropriate system command available.
    If we don't it raises an exception.

    '''

    path = os.environ['PATH'].split(':')
    for p in path:
        f = os.path.join(p, exe)
        if os.path.isfile(f):
            if verb:
                print("# Found %s in %s" % (exe, f), file=sys.stderr)
            return True
    # it wasn't found
    print("# ERROR: Couldn't find %s" % exe, file=sys.stderr)
    raise FileNotFoundError(exe)
    return False

def elapsed_time(t1, text=''):
    import time
    t2 = time.time()
    hh = int((t2 - t1) / 3600.)
    mm = int((t2 - t1) / 60 - hh * 60)
    ss = (t2 - t1) - 60 * mm - 3600 * hh
    print("Elapsed time: %dh %dm %2.2fs %s" % (hh, mm, ss, text))
    print("Elapsed time: %dh %dm %2.2fs %s" % (hh, mm, ss, text),
          file=sys.stderr)
    #print >>sys.stderr,"Elapsed time: %dm %2.2fs" % ( int( (t2-t1)/60.),
    #(t2-t1) - 60*int((t2-t1)/60.))
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
        print("\r Parsing SEx head for:", catalog, file=sys.stderr)

    # Dictionary with column numbers
    SExcols = {}

    # Read the SExtractor catalog
    for line in open(catalog).readlines():

        if line[0] != '#':
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
                print("# %-20s %s" % (key, SExcols[key] + 1), file=sys.stderr)
        except:
            continue

    return SExcols

########################################
# Read and write fits simple functions
#######################################


# To read in simple fits file
def rfits(filename, hdu=0):

    ff = fits.open(filename, "readonly")
    data = ff[hdu].data
    hdr = ff[hdu].header
    ff.close()
    return data, hdr


# To write in a simple file
def writefits(array, filename):
    import os
    # Writes fits files of a given array
    newfits = fits.HDUList()
    hdu = fits.PrimaryHDU()
    hdu.data = array
    newfits.append(hdu)

    # Remove old version of file before
    if os.path.isfile(filename):
        os.remove(filename)

    newfits.writeto(filename)
    newfits.close
    return


wfits = writefits

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
def deNAN(a, value=0.0):
    nans = np.logical_not(np.less(a, 0.0) + np.greater_equal(
        a, 0.0))
    return np.where(nans, value, a)

def cmdline():

    from optparse import OptionParser

    # Read in the command line options
    USAGE = '''usage:\t %prog <assoc> <in_datapath> <out_datapath> [options]
    i.e.: %prog ~/bcs-inventory/BCS0519-5448.assoc ~/PROC ~/BCS/PROC'''

    parser = OptionParser(usage=USAGE)

    parser.add_option("--Copy",
                      action="store_true",
                      dest="Copy",
                      default=0,
                      help="Copy files")

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

    parser.add_option("--noCleanUP",
                      action='store_true',
                      dest='noCleanUP',
                      default=0,
                      help='Whether or not to remove uncompressed files')

    parser.add_option("--noRGB",
                      action='store_true',
                      dest='noRGB',
                      default=0,
                      help='Whether or not to create RGB images from mosaics')

    parser.add_option("--noNEWFIRM",
                      action='store_true',
                      dest='noNEWFIRM',
                      default=0,
                      help='Whether or not to swarp the newfirm mosaic data.')

    (options, args) = parser.parse_args()

    if len(args) < 3:
        parser.error(USAGE + "\nMust supply at least one argument required")

    # Dry run turns off everything... almost
    if options.dryrun:
        options.Copy = 0
        options.noSWarp = 1
        options.noBPZ = 1
        options.noMask = 1
        options.Dust = 0

    # Same for WeightOnly
    if options.WeightOnly:
        options.Copy = 0
        options.noSWarp = 1
        options.noBPZ = 1
        options.noMask = 1
        options.Dust = 0

    return options, args


def main():

    # The start time
    tstart = time.time()

    # Get the command line options
    opt, arg = cmdline()

    assocfile = arg[0]
    inpath = arg[1]
    outpath = arg[2]

    # Read in the command line options
    #USAGE = " usage:\t"+os.path.split(sys.argv[0])[1]+" <assoc> <in_datapath>
    #<out_datapath> [dryrun]"
    # Example
    # ./combcat.py ~/bcs-inventory/BCS0519-5448.assoc ~/PROC ~/BCS/PROC

    # Init the class
    c = combcat(assocfile,
                datapath=inpath,
                outpath=outpath,
                verb='yes',
                dryrun=opt.dryrun,
                noSWarp=opt.noSWarp)

    # Copy the files into dir
    if opt.Copy:
        c.copyfiles(copy=True)
    else:
        c.copyfiles(copy=False)

    # SWarp
    if not opt.noSWarp:
        c.swarp_files(dryrun=opt.noSWarp,
                  conf="SWarp-common.conf",
                  combtype=opt.combtype)

    # swarp NEWFIRM images (if they are there)
    if not opt.noNEWFIRM:
        c.swarp_newfirm()

    # Make the detection image
    if opt.useMask:
        # Make the mask from the weights
        c.generate_masks(filters=('i', ), dryrun=opt.noMask)
        c.makeDetectionIma(filter='i')

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

    # make RGB images (pngs)
    if not opt.noRGB:
        print('make rgb')
        c.make_RGB(kband=True)

    # cleanup
    if opt.noCleanUP:
        pass
    else:
        print("CLEANUP!")
        c.cleanup_files()

    elapsed_time(tstart, c.tilename)
    return


main()
