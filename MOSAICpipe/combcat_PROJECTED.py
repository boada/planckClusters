#!/usr/bin/env python3

import os
import sys
import tableio
import numpy as np
from math import log10
import time
import string
import subprocess
import shlex
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

        self.pipeline = '/home/boada/Projects/planckClusters/MOSAICpipe/'

        # Check for environ vars
        if not os.getenv('PIPE'):
            os.environ['PIPE'] = os.path.join(self.pipeline)

        if not os.getenv('BPZPATH'):
            os.environ['BPZPATH'] = os.path.join(self.pipeline + 'bpz-1.99.3')

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
            subprocs = []
            for line in lines:
                if line[0] == "#":
                    continue

                vals = line.split()
                fname = os.path.basename(vals[0])

                # hack to deal with compressed files in the association
                if '.fz' in fname:
                    fname = fname.rstrip('.fz')
                    if not os.path.isfile(fname) and not self.noSWarp:
                        subprocs.append(subprocess.Popen(
                            shlex.split('funpack -v {}.fz'.format(fname))))

                # Figure out the dqmask
                if 'tu' in fname:
                    i = 1
                    while True:
                        nid = int(os.path.splitext(fname)[0][2:]) + i
                        ext = os.path.splitext(fname)[1]
                        pre = fname[0:2]
                        dqmask = "%s%s%s" % (pre, nid, ext)
                        if os.path.isfile('{}.fz'.format(dqmask)):
                            break
                        else:
                            print('Looking...')
                            i += 1
                    if not os.path.isfile(dqmask) and not self.noSWarp:
                        subprocs.append(subprocess.Popen(
                            shlex.split('funpack -v {}.fz'.format(dqmask))))
                elif 'k4' in fname:
                    dqmask = fname.replace('opi', 'opd')
                    if not os.path.isfile(dqmask) and not self.noSWarp:
                        subprocs.append(subprocess.Popen(
                            shlex.split('funpack -v {}.fz'.format(dqmask))))

            [p.wait() for p in subprocs]
            [i.kill() for i in subprocs]

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
                    #if not os.path.isfile(fname) and not self.noSWarp:
                    #os.system('funpack -v {}.fz'.format(fname))

                filtername = vals[1]
                self.exptime[fname] = float(vals[2])
                self.airmass[fname] = float(vals[3])

                # Figure out the dqmask
                if 'tu' in fname:
                    i = 1
                    while True:
                        nid = int(os.path.splitext(fname)[0][2:]) + i
                        ext = os.path.splitext(fname)[1]
                        pre = fname[0:2]
                        dqmask = "%s%s%s" % (pre, nid, ext)
                        if os.path.isfile('{}.fz'.format(dqmask)):
                            break
                        else:
                            print('Looking...')
                            i += 1
                    if not os.path.isfile(dqmask) and not self.noSWarp:
                        pass
                        #os.system('funpack -v {}.fz'.format(dqmask))
                        #subprocs.append('funpack -v {}.fz'.format(dqmask))
                elif 'k4' in fname:
                    dqmask = fname.replace('opi', 'opd')
                    if not os.path.isfile(dqmask) and not self.noSWarp:
                        pass
                        #os.system('funpack -v {}.fz'.format(dqmask))
                        #subprocs.append('funpack -v {}.fz'.format(dqmask))

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
                if not self.noSWarp:
                    with fits.open(dqmask, mode='readonly') as inHDU:
                        if len(inHDU) > 2:
                            self.files_weight[filtername].append(
                                "%s.weight.fits" % os.path.splitext(fname)[0])
                        else:
                            self.files_weight[filtername].append(
                                "%s.weight.fits'[0]'" %
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

    def center_dither(self, conf="SWarp-center.conf", filter='i', dryrun=False):
        ''' Center the dither pattern for SWARP '''

        check_exe("swarp")

        # The configuration file
        center_conf = os.path.join(self.pipeline, 'confs', conf)

        # First we need to get the center for all the files, swarp them all
        opts = {}

        #opts["PIXEL_SCALE"] = self.pixscale
        opts["PIXELSCALE_TYPE"] = "MAX"
        opts["RESAMPLING_TYPE"] = "LANCZOS3"
        opts["CENTER_TYPE"] = "ALL"
        opts["NTHREADS"] = "0"

        if not filter:
            cmd = "swarp %s -c %s " % (' '.join(self.filelist), center_conf)
        else:
            try:
                cmd = "swarp %s -c %s " % (' '.join(self.files[filter]),
                                            center_conf)
            except KeyError:
                filter = 'r'
                try:
                    cmd = "swarp %s -c %s " % (' '.join(self.files[filter]),
                                            center_conf)
                except KeyError:
                    filter = 'K'
                    cmd = "swarp %s -c %s " % (' '.join(self.files[filter]),
                                            center_conf)

        if not filter:
            opts["IMAGEOUT_NAME"] = os.path.join(
                self.outdir, "SWarp-%s-center.fits" % (self.tilename))
        else:
            opts["IMAGEOUT_NAME"] = os.path.join(
                self.outdir, "SWarp-%s-center_%s.fits" % (self.tilename, filter))

        for param, value in list(opts.items()):
            cmd += "-%s %s " % (param, value)

        print(cmd)
        if not self.dryrun or dryrun:
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

        pixscale = abs(header['CD1_1'] * 3600)

        print("\tImage Size:  %s x %s" % (nx, ny))
        print("\tCentered on: %s   %s" % (x_center, y_center))
        print("\tPixelscale: %.4f" % (pixscale))

        self.nx = nx
        self.ny = ny
        self.xo = x_center
        self.yo = y_center
        self.pixscale = float('{0:.4f}'.format(pixscale))

        # Store wether center was done....
        if not filter:
            self.centered = True
        else:
            self.centered = True

        return

    def get_FLXSCALE(self, magbase=30, filename=None):
        ''' I don't really know what this function is supposed to do. '''

        print("# Computing FLXSCALE for magbase=%s" % magbase)
        self.magbase = magbase
        self.flxscale = {}
        for filter in self.filters:

            self.flxscale[filter] = []

            if not filename:
                for fname in self.files[filter]:
                    header = getheader(fname)
                    #zp = header['MAGZERO'] + 2.5*log10(header['EXPTIME'])
                    zp = header['MAGZERO']
                    flxscale = 10.0**(0.4 * (magbase - zp))
                    self.flxscale[filter].append(flxscale)
            if filename:
                header = getheader(filename)
                zp = header['MAGZERO']
                flxscale = 10.0**(0.4 * (magbase - zp))

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

        _conf = os.path.join(self.pipeline, 'confs', conf)

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
                    'FILTER', 'MAGZERO', 'MAGZSIG']

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

            cmd = "swarp %s -c %s -IMAGEOUT_NAME %s -WEIGHTOUT_NAME %s " %\
                (filelist, _conf, outimage, outweight)
            cmd += " -WEIGHT_IMAGE %s " % (
                ",".join(self.files_weight[filter]))
            cmd += " -FSCALE_DEFAULT %s " % (
                ",".join(map(str, self.flxscale[filter])))
            cmd += opts

            if not dryrun:
                print("# Will run:\n\t%s" % cmd)
                os.system(cmd)
                AM = self.calc_airmass(filter)
                put_airmass(outimage, AM)
                ET = self.calc_exptime(filter)
                put_exptime(outimage, ET)

                # make sure the header keywords have been propigated. This is
                # important for the astro and flux calibration
                put_headerKeywords(self.files[filter][0], outimage, keywords,
                                self.xo, self.yo)

            else:
                print(cmd)

        return

    def get_filenames(self):
        ''' This script will get called if we choose not to swarp the images.
        This makes sure that all of the variables are defined.

        '''
        # this is to convert the hms/dms coordinates back into decimal degrees.
        # I don't really want to use this, but I am not rewriting things.
        from astLib import astCoords

        self.combima = {}
        self.combcat = {}
        self.weightima = {}

        self.center_dither(dryrun=True)
        self.xo = astCoords.hms2decimal(self.xo, ':')
        self.yo = astCoords.dms2decimal(self.yo, ':')

        for filter in self.filters:

            # Store the names
            self.combima[filter] = "%s%s" % (self.tilename, filter)
            self.weightima[filter] = "%s%s_weight.fits" % (self.tilename,
                                                           filter)
            # make the combined catalog filenames
            if os.path.isfile('{}{}.cat'.format(self.tilename, filter)):
                self.combcat[filter] = '{}{}.cat'.format(self.tilename, filter)
            elif os.path.isfile('{}{}_cal.cat'.format(self.tilename, filter)):
                self.combcat[filter] = '{}{}_cal.cat'.format(self.tilename,
                                                             filter)
        # The default output names
        self.colorCat = self.tilename + ".color"
        self.columnsFile = self.tilename + ".columns"

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

        common_conf = os.path.join(self.pipeline, 'confs', conf)
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

    def get_astrometry(self):
        ''' This calls the pp script to astrometrically correct the images. I'm
        breaking it appart from the rest of the script so I can control when
        things happen. It is mostly for testing.

        '''

        # we have to prepare the images first.
        check_exe('pp_prepare')
        check_exe('pp_register')

        # little patch to keep it from crashing
        try:
            os.mkdir('.diagnostics')
        except:
            pass

        # make the diagnostics file too
        open('diagnostics.html', 'a').close()

        try:
            if ':' in self.xo or ':' in self.yo:
                from astLib import astCoords
                self.xo = astCoords.hms2decimal(self.xo, ':')
                self.yo = astCoords.dms2decimal(self.yo, ':')
        except TypeError:
            pass

        subprocs = []
        for filter in self.filters:
            mosaic = '{}.fits'.format(self.combima[filter])
            if self.pixscale == 0.25:
                cmd = 'pp_prepare -ra {} -dec {} -telescope {} {}'.format(
                    self.xo, self.yo, 'KPNO4MOS3', mosaic)
            elif self.pixscale == 0.2666:
                cmd = 'pp_prepare -ra {} -dec {} -telescope {} {}'.format(
                    self.xo, self.yo, 'KPNO4MOS1', mosaic)
            else:
                cmd = 'pp_prepare -ra {} -dec {} {}'.format(
                    self.xo, self.yo, mosaic)

            if not self.dryrun:
                print(cmd)
                subprocs.append(subprocess.Popen(shlex.split(cmd)))

        [p.wait(timeout=600) for p in subprocs]
        [i.kill() for i in subprocs]

        for filter in self.filters:
            mosaic = '{}.fits'.format(self.combima[filter])

            cmd = 'pp_register -snr 10 -minarea 12 {}'.format(mosaic)

            if not self.dryrun:
                subprocs.append(subprocess.Popen(shlex.split(cmd)))

        [p.wait(timeout=600) for p in subprocs]
        [i.kill() for i in subprocs]

        return

    def get_zeropt(self):
        from astropy.io import fits
        ''' This is going to call the photometrypipline script that lives in
        the projects directory. It should both astrometrically correct the
        mosaics and calculate the overall zeropoint of the mosaic. Currently,
        this has only been tested on the MOSAIC camera, and should NOT work
        with NEWFIRM, yet. It's on my TODO list.

        '''

        check_exe('pp_photometry')
        check_exe('pp_calibrate')
        check_exe('pp_distill')

        try:
            os.mkdir('.diagnostics')
        except:
            pass

        # make the diagnostics file too
        with open('diagnostics.html', 'a') as f:
            pass

        subprocs = []

        for filter in self.filters:
            if os.path.isfile('photometry_control_star_{}.dat'.format(filter)):
                print('# remove old photometry')
                os.remove('photometry_control_star_{}.dat'.format(filter))

            mosaic = '{}.fits'.format(self.combima[filter])
            # check to make sure it has the right header keywords
            with fits.open(mosaic, mode='update') as mos:
                keywords = ['TEL_KEYW', 'TELINSTR', 'MIDTIMJD', 'SECPIXY',
                            'SECPIXX', 'PHOT_K']

                if self.pixscale == 0.25:
                    values = ['KPNO4MOS3', 'KPNO4m/MOSAIC', 0.0, 0.25, 0.25,
                                  0.5]
                elif self.pixscale == 0.2666:
                    values = ['KPNO4MOS1', 'KPNO4m/MOSAIC', 0.0, 0.2666, 0.2666,
                                  0.5]
                else:
                    values = ['KPNO4NEWF', 'KPNO4m/NEWFIRM', 0.0, 0.4, 0.4,
                                  0.5]

                for kw, val in zip(keywords, values):
                    # only update if they don't exist
                    try:
                        mos[0].header[kw]
                    except KeyError:
                        mos[0].header[kw] = val

                #####################################################
                ##### HACK TO MAKE THE ZEROPOINT ERRORS SMALLER #####
                #####################################################
                if filter != 'K':
                    mos[0].header['GAIN'] = 1.

            if filter != 'K':
                cmd = 'pp_photometry -snr 10 -aprad 5.7 {}'.format(mosaic)
            elif self.pixscale == 0.25:
                cmd = 'pp_photometry -snr 10 -aprad 16 {}'.format(mosaic)
            elif self.pixscale == 0.2666:
                cmd = 'pp_photometry -snr 10 -aprad 15 {}'.format(mosaic)
            elif self.pixscale == 0.4:
                cmd = 'pp_photometry -snr 10 -aprad 10 {}'.format(mosaic)

            if not self.dryrun:
                subprocs.append(subprocess.Popen(shlex.split(cmd)))

        [p.wait() for p in subprocs]

        for filter in self.filters:
            mosaic = '{}.fits'.format(self.combima[filter])
            cmd = 'pp_calibrate {}'.format(mosaic)

            if not self.dryrun:
                os.system(cmd)
                #subprocs.append(subprocess.Popen(shlex.split(cmd)))

        #[p.wait() for p in subprocs]
        #[i.kill() for i in subprocs]

        for filter in self.filters:
            mosaic = '{}.fits'.format(self.combima[filter])
            # the positions file fails -- it's just used to make sure it doesn't
            # query JPL because our stuff isn't an asteroid.
            cmd = 'pp_distill {}'.format(mosaic)
            os.system(cmd) # gotta call it to make it work

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
                    if 'PV' in key:
                        f[0].header[key] = str(val)

        return

    # run sextractor
    def SEx(self, det_filter='i'):
        ''' Runs SEXTRACTOR on the mosaicked images created by swarp. It should
        use the zero point created by the photometrypipline.

        '''

        check_exe("sex")

        # list of output catalog names
        self.combcat = {}

        # The configuration file
        self.SExinpar = os.path.join(self.pipeline, 'confs',
                                     'bcs_Catalog.inpar')

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

        for filter in self.filters:

            input = self.combima[filter] + ".fits"
            photocal = 'photometry_control_star_{}.dat'.format(filter)
            try:
                _tmp = np.genfromtxt(photocal, names=True, dtype=None)
                zeropt = _tmp['ZP']
                if zeropt == 0.0:
                    print('WARNING!: Photometric calibration not set!')
                    try:
                        zeropt = self.magbase
                        self.combcat[filter] = self.combima[filter] + ".cat"
                    except AttributeError:
                        zeropt = 26
                        self.combcat[filter] = self.combima[filter] + ".cat"
                else:
                    self.combcat[filter] = self.combima[filter] + "_cal.cat"
            except OSError:
                print('WARNING!: Photometric calibration not set!')
                try:
                    zeropt = self.magbase
                    self.combcat[filter] = self.combima[filter] + ".cat"
                except AttributeError:
                    zeropt = 26
                    self.combcat[filter] = self.combima[filter] + ".cat"

            output = self.combcat[filter]

            print(zeropt)

            opts = ''
            if filter == 'i':
                backima = self.combima[filter] + "_SEG.fits"
                #opts = " -CHECKIMAGE_TYPE BACKGROUND_RMS"
                opts = " -CHECKIMAGE_TYPE SEGMENTATION"
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
    def BuildColorCat(self, newfirm=False):

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
        # the get_data function requires a tuple
        detectionColumns = tuple(detectionList)
        detection_variables = tableio.get_data(detCatalog,
                                               detectionColumns)

        # Read in the MAG_ISO and MAG_ISOERR from each catalog
        for filter in self.filters:
            if not newfirm and filter == 'K':
                continue
            # get the zeropoint Info
            tmp = np.genfromtxt('photometry_control_star_{}.dat'.format(
                                filter), names=True, dtype=None)
            zpoint = tmp['ZP']
            zp_error = tmp['ZP_sig']

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

            #nondetected = np.less_equal(flux[filter], 0.0) * \
            #              np.greater(fluxerr[filter], 0.0)

            # update: There are a lot of really small positive values. I am
            # going to modify this to look for values really close to zero.

            nondetected = (flux[filter] < 0.0) | (abs(flux[filter]) < 1E-3)

            # Those objects with error flux and flux equal to 0 are assigned a
            # magnitude of -99
            # and a flux of 0, which is interpreted by SExtractor as a
            # non-observed object

            nonobserved = np.less_equal(fluxerr[filter], 0.0)

            # When flux error > 100*(flux), mark as nonobserved (Benitez,
            # 24-Oct-03).

            nonobserved = np.where(fluxerr[filter] > 100 *
                                        (abs(flux[filter])), True,
                                   nonobserved)

            detected = np.logical_not(nonobserved + nondetected)

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
            if filter == 'K':
                outColumns.append(filter + '_KittPeak_MAG_ISO')
                outColumns.append(filter + '_KittPeak_MAGERR_ISO')
            else:
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
            if not newfirm and filter == 'K':
                continue
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
        with open(self.columnsFile, 'w') as cfile:
            cfile.write('## ' + time.ctime() + '\n')
            cfile.write('## BPZ .columns file for Observation: {}\n'.format(
                self.tilename))
            cfile.write('## (This file was generated by the '
                    'BCS Rutgers pipeline)\n##\n')

            i = len(detection_variables)
            for filter in self.filters:
                if not newfirm and filter == 'K':
                    continue
                # Get the zeropoint information
                tmp = np.genfromtxt('photometry_control_star_{}.dat'.format(
                                    filter), names=True, dtype=None)
                zpoint = tmp['ZP']
                zp_error = tmp['ZP_sig']

                if filter == 'i':
                    n_mo = str(i + 1)

                if filter == 'K':
                    cfile.write('%s_KittPeak\t %s,%s\t AB\t %.2f\t 0.0\n' %
                            (filter, i + 1, i + 2, zp_error))
                else:
                    cfile.write('%s_MOSAICII\t %s,%s\t AB\t %.2f\t 0.0\n' %
                            (filter, i + 1, i + 2, zp_error))
                i += 2

            cfile.write('M_0\t%s\n' % n_mo)
        return

    def cleanCatalogs(self):
        ''' This function takes all of the sextractor and newly created color
        catalogs and cleans them to make sure there are no horrible detections.
        Horrible detections include things where the detection image doesn't
        overlap with the other image. This causes us to get crazy magnitudes
        that don't make any sense. So throw those away. I am going to clean
        both the color catalog, and the Sextractor catalogs. The color catalogs
        will get -99 (for now), but it might be better to just remove them all
        together.

        '''

        from astropy.io import ascii

        color = ascii.read('{}.color'.format(self.tilename))
        for filter in self.filters:
            sex_cat = ascii.read(self.combcat[filter])
            # find all the spots where the FWHM == 0.
            fwhm_idx = np.where(sex_cat['FWHM_IMAGE'] == 0.0)[0]
            color['{}_MOSAICII_MAG_ISO'.format(filter)][fwhm_idx] = -99.
            color['{}_MOSAICII_MAGERR_ISO'.format(filter)][fwhm_idx] = 0.

    # Run Benitez BPZ
    def runBPZ(self, Specz=True):
        """Runs BPZ on the multicolor catalog file using the .columns """

        # first we update with Specz's if we want to
        match_SEx(self.tilename, self.filters)
        add_Speczs(self.tilename)
        #add_extra_photometry(self.tilename)

        print('Starting photometric redshift determination...',
              file=sys.stderr)
        bpz = os.path.join(os.environ['BPZPATH'], 'bpz.py ')
        #bpzcat = self.tilename + ".bpz"
        bpzprobs = self.tilename + ".probs"

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

    def calc_airmass(self, filter, combtype='median'):
        ''' Calculates the mean airmass for a set of observations. '''

        x = []
        for file in self.files[filter]:
            #print file,self.airmass[file]
            x.append(self.airmass[file])
        x = np.asarray(x)
        if combtype == 'mean':
            return x.mean()
        elif combtype == 'median':
            return np.median(x)
        else:
            print('combtype not understood. Check your call.')
            raise TypeError

    def calc_exptime(self, filter, combtype='median'):
        ''' Calculates the mean airmass for a set of observations. '''

        x = []
        for file in self.files[filter]:
            #print file,self.airmass[file]
            x.append(self.exptime[file])
        x = np.asarray(x)
        if combtype == 'mean':
            return x.mean()
        elif combtype == 'median':
            return np.median(x)
        else:
            print('combtype not understood. Check your call.')
            raise TypeError

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

    def make_RGB(self, kband=False, conf='stiff-common.conf'):

        try:
            check_exe('stiff')
        except FileNotFoundError:
            return

        _conf = os.path.join(self.pipeline, 'confs', conf)

        # input files
        if kband:
            red = '../../newfirm/stacked/{}{}.fits'.format(self.tilename, 'K')
            if os.path.isfile(red):
                green = './{}{}.fits'.format(self.tilename, 'i')
                blue = './{}{}.fits'.format(self.tilename, 'r')
                bands = 'Kir'
            else:
                print('k-band file not found, restoring defaults')
                red = './{}{}.fits'.format(self.tilename, 'i')
                green = './{}{}.fits'.format(self.tilename, 'r')
                blue = './{}{}.fits'.format(self.tilename, 'g')
                bands = 'irg'
        else:
            red = './{}{}.fits'.format(self.tilename, 'i')
            green = './{}{}.fits'.format(self.tilename, 'r')
            blue = './{}{}.fits'.format(self.tilename, 'g')
            bands = 'irg'

        # output file
        output = '{}{}.tiff'.format(self.tilename, bands)

        # options
        opts = ['-MIN_LEVEL', '0.01',
                '-MAX_LEVEL', '0.993',
                '-MIN_TYPE', 'GREYLEVEL',
                '-MAX_TYPE', 'QUANTILE',
                '-COLOUR_SAT', '2.0',
                '-DESCRIPTION', "'{} RGB'".format(self.tilename),
                '-WRITE_XML', 'N',
                '-COPYRIGHT', "'Steven Boada'"]

        # build the command -- a space is always last
        cmd = 'stiff {} {} {} -c {} '.format(red, green, blue, _conf)
        cmd += '-OUTFILE_NAME {} '.format(output)
        # append the options
        for opt in opts:
            cmd += '{} '.format(opt)

        print(cmd)
        if not self.dryrun:
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

    with fits.open(infile, mode='readonly') as inHDU:
        if len(inHDU) > 2:
            newhdu = fits.HDUList()
            newhdu.append(inHDU[0])
            for hdu in inHDU[1:]:
                data = hdu.data
                newdata = np.where(data > 0, 0, 1)
                hdu.data = newdata.astype("UInt8")  # Just as ushort integer
                newhdu.append(hdu)
        else:
            data = inHDU[-1].data
            newdata = np.where(data > 0, 0, 1)
            inHDU[-1].data = newdata.astype("UInt8")  # Just as ushort integer
            newhdu = fits.HDUList()
            newhdu.append(inHDU[-1])

        newhdu.writeto(outfile, clobber=clobber)
        newhdu.close()
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

def put_headerKeywords(origFile, newFile, keywords, ra, dec):
    ''' This takes all of the keywords from the original fits image and makes
    sure they are present in the new fits image. This is important for the
    mosaic3 data. For that data it doesn't seem to propigate the header
    keywords correctly. This function makes sure those keywords are present.

    '''

    print(origFile, newFile)
    with fits.open(origFile) as oF:
        with fits.open(newFile, mode='update') as nF:
            oHDR = oF[0].header # original header
            nHDR = nF[0].header # new header

            print('# Updating {} with new keywords'.format(newFile))
            for kw in keywords:
                try:
                    nHDR[kw]
                except KeyError: # only update if it doesn't exist
                    try:
                        if kw == 'RA':
                            nHDR[kw] = ra
                        elif kw == 'DEC':
                            nHDR[kw] = dec
                        else:
                            nHDR[kw] = oHDR[kw]
                    except KeyError: # sometimes the original keys aren't there
                        print('#{} Not Found'.format(kw))
    return

def put_airmass(filename, airmass):

    # Open the file, update mode
    with fits.open(filename, mode="update") as f: # open a FITS file
        hdr = f[0].header  # the primary HDU header

        print("# Updating %s with AIRMASS=%s after SWarp" % (
            filename, airmass), file=sys.stderr)
        hdr['AIRMASS'] = (airmass, 'After SWarp equivalent airmass')

        # check the file
        f.verify('fix')

    return

def put_exptime(filename, exptime):

    # Open the file, update mode
    with fits.open(filename, mode="update") as f: # open a FITS file
        hdr = f[0].header  # the primary HDU header

        print("# Updating %s with EXPTIME=%s after SWarp" % (
            filename, exptime), file=sys.stderr)
        hdr['EXPTIME'] = (exptime, 'After SWarp equivalent exptime (sec)')

        # check the file
        f.verify('fix')
    return

def put_scaled_head(file):
    ''' Currently not used. '''

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
    ''' Currently not used. '''

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
    ''' Don't really know what this is supposed to do. '''

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
    ''' Currently not used. '''

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

def match_SEx(tilename, filters):
    from astropy.io import ascii
    from astropy import wcs
    from astropy.table import Column
    from astropy.coordinates import SkyCoord
    from astropy import units as u

    # get the wcs info from the i-band because it is the dectection image
    if os.path.isfile('{}i.fits'.format(tilename)):
        w = wcs.WCS('{}i.fits'.format(tilename))
    else:
        print('# No i-band detection image aborting...')
        return

    # read in the sdss catalog. We are doing this first because we only need to
    # do it once
    try:
        sdss_cat = ascii.read('/home/boada/Projects/'
                         'planckClusters/data/extern/SDSS/{}/'
                         '{}_SDSS_catalog.csv'.format(tilename, tilename))
        if len(sdss_cat) < 2:
            print('# SDSS TOO SHORT!')
            sdss = False
        else:
            sdss = True
    except FileNotFoundError:
        print('# SDSS CATALOG NOT FOUND!')
        sdss = False

    try:
        ps1_cat = ascii.read('/home/boada/Projects/'
                         'planckClusters/data/extern/PS1/{}/'
                         '{}_PS1_catalog.csv'.format(tilename, tilename))
        if len(ps1_cat) < 2:
            print('# PS1 TOO SHORT!')
        else:
            ps1 = True
    except FileNotFoundError:
        print('# PS1 CATALOG NOT FOUND!')
        ps1 = False

    try:
        twoMASS_cat = ascii.read('/home/boada/Projects/'
                         'planckClusters/data/extern/2MASS/{}/'
                         '{}_2MASS_catalog.csv'.format(tilename, tilename))
        if len(twoMASS_cat) < 2:
            print('# 2MASS TOO SHORT!')
            twoMASS = False
        else:
            twoMASS = True
    except FileNotFoundError:
        print('# 2MASS CATALOG NOT FOUND!')
        twoMASS = False

    # need these coordinates for the matching
    if sdss:
        s_coord = SkyCoord(ra=sdss_cat['ra'] * u.degree, dec=sdss_cat['dec'] *
                       u.degree)
    if ps1:
        p_coord = SkyCoord(ra=ps1_cat['ramean'] * u.degree,
                        dec=ps1_cat['decmean'] * u.degree)
    if twoMASS:
        tm_coord = SkyCoord(ra=twoMASS_cat['ra'] * u.degree,
                    dec=twoMASS_cat['dec'] * u.degree)

    for filter in filters:
        # only do the Kband when we get there
        if filter == 'K':
            if os.path.isfile('{}{}_cal.cat'.format(tilename, filter)):
                kband = True
            elif os.path.isfile('{}{}.cat'.format(tilename, filter)):
                kband = True
            else:
                kband = False
                continue
        else:
            kband = False

        try:
            cat = ascii.read('{}{}_cal.cat'.format(tilename, filter))
            cal = True
        except FileNotFoundError:
            cat = ascii.read('{}{}.cat'.format(tilename, filter))
            cal = False

        # calulate the RA/DEC of all of the detections
        ra, dec = w.all_pix2world(cat['X_IMAGE'], cat['Y_IMAGE'], 0)
        ra = Column(ra, name='RA')
        dec = Column(dec, name='DEC')
        try:
            cat.add_column(ra)
        except ValueError:
            pass
        try:
            cat.add_column(dec)
        except ValueError:
            pass

        # need these coordinates for the matching
        c_coord = SkyCoord(ra=cat['RA'] * u.degree, dec=cat['DEC'] * u.degree)

        if sdss and not kband:
            # match the two catalogs -- idxc for cat, idxs for sdss
            idxc, idxs, d2d, d3d = s_coord.search_around_sky(c_coord, 1 * u.arcsec)

            # make some data to catch
            d = np.ones(len(cat)) * 99.0 # 99 is the non-detection value in SEx...

            # add some extra info from SDSS -- Specz, Photoz
            cols = []
            cols.append(Column(d, name='sdss_{}'.format(filter)))
            cols.append(Column(d, name='sdss_{}_err'.format(filter)))
            cols.append(Column(d, name='sdss_petro_{}'.format(filter)))
            cols.append(Column(d, name='sdss_petro_{}_err'.format(filter)))
            cols.append(Column(d, name='sdss_objid', dtype=np.long))
            cols.append(Column(d, name='sdss_specz'))
            cols.append(Column(d, name='sdss_specz_err'))
            cols.append(Column(d, name='sdss_photoz'))
            cols.append(Column(d, name='sdss_photoz_err'))
            cols.append(Column(d, name='sdss_type'))
            for col in cols:
                try:
                    cat.add_column(col)
                except ValueError:
                    pass

            # merge the matches
            cat['sdss_objid'][idxc] = sdss_cat['objid'][idxs]
            cat['sdss_{}'.format(filter)][idxc] = sdss_cat[
                                                'fiberMag_{}'.format(filter)][idxs]
            cat['sdss_{}_err'.format(filter)][idxc] = sdss_cat[
                                            'fiberMagErr_{}'.format(filter)][idxs]
            cat['sdss_petro_{}'.format(filter)][idxc] = \
                                    sdss_cat['petroRad_{}'.format(filter)][idxs]
            cat['sdss_petro_{}_err'.format(filter)][idxc] = \
                                    sdss_cat['petroRadErr_{}'.format(filter)][idxs]
            cat['sdss_specz'][idxc] = sdss_cat['specz'][idxs]
            cat['sdss_specz_err'][idxc] = sdss_cat['specz_err'][idxs]
            cat['sdss_photoz'][idxc] = sdss_cat['photoz'][idxs]
            cat['sdss_photoz_err'][idxc] = sdss_cat['photoz_err'][idxs]
            cat['sdss_type'][idxc] = sdss_cat['type'][idxs]

        #### NOW THE PANSTARRS DATA ####
        if ps1 and not kband:
            idxc, idxp, d2d, d3d = p_coord.search_around_sky(c_coord, 1 * u.arcsec)
            # make some data to catch
            d = np.ones(len(cat)) * 99.0 # 99 is the non-detection value in SEx...
            col = Column(d, name='ps1_{}'.format(filter))
            col_err = Column(d, name='ps1_{}_err'.format(filter))
            try:
                cat.add_column(col)
                cat.add_column(col_err)
            except ValueError:
                pass
            # merge the matches
            cat['ps1_{}'.format(filter)][idxc] = ps1_cat[
                                            '{}meanpsfmag'.format(filter)][idxp]
            cat['ps1_{}_err'.format(filter)][idxc] = ps1_cat[
                                            '{}meanpsfmagerr'.format(filter)][idxp]

        #### NOW THE 2MASS DATA -- IF KBAND ####
        if kband and twoMASS:
            idxc, idxp, d2d, d3d = tm_coord.search_around_sky(c_coord, 1 * u.arcsec)
            # make some data to catch
            d = np.ones(len(cat)) * 99.0 # 99 is the non-detection value in SEx...
            col = Column(d, name='2MASS_{}'.format(filter))
            col_err = Column(d, name='2MASS_{}_err'.format(filter))
            try:
                cat.add_column(col)
                cat.add_column(col_err)
            except ValueError:
                pass
            # merge the matches
            cat['2MASS_{}'.format(filter)][idxc] = twoMASS_cat[
                                            '{}mag'.format(filter)][idxp]
            cat['2MASS_{}_err'.format(filter)][idxc] = twoMASS_cat[
                                            'e_{}mag'.format(filter)][idxp]

        # write out the results
        cat.write('tmp.color', format='ascii.commented_header', overwrite=True)

        # mv the old color catalog
        if cal:
            os.rename('{}{}_cal.cat'.format(tilename, filter),
                      '{}{}_cal.cat.orig'.format(tilename, filter))
        else:
            os.rename('{}{}.cat'.format(tilename, filter),
                      '{}{}.cat.orig'.format(tilename, filter))

        with open('tmp.color') as f:
            line = f.readline()
            line = line.split()
            if cal:
                file = '{}{}_cal.cat'.format(tilename, filter)
            else:
                file = '{}{}.cat'.format(tilename, filter)
            with open(file, 'wt') as f2:
                f2.write('## ' + time.ctime() + '\n')
                f2.write('## Matched Catalog file for Observation: '
                            '{}\n'.format(tilename))
                f2.write('## (This file was generated by the '
                        'BCS Rutgers pipeline)\n##\n')
                for i, l in enumerate(line[1:]):
                    f2.write('{} {} {}\n'.format(line[0], i + 1, l))
                line = f.readline()
                while line:
                    f2.write(line)
                    line = f.readline()

        os.remove('tmp.color')

def add_Speczs(tilename):
    ''' This function adds the spec-z's to the color catalog. The previous
    function adds all of the extra information to the sextractor catalogs. So
    you can run either one, or both. Because they don't do the same thing even
    though it looks like it.

    '''

    from astropy.io import ascii
    from astropy import wcs
    from astropy.table import Column
    from astropy.coordinates import SkyCoord
    from astropy import units as u

    w = wcs.WCS('{}i.fits'.format(tilename))
    cat = ascii.read('{}.color'.format(tilename))
    ra, dec = w.all_pix2world(cat['X_IMAGE'], cat['Y_IMAGE'], 0)
    ra = Column(ra, name='RA')
    dec = Column(dec, name='DEC')
    cat.add_column(ra)
    cat.add_column(dec)

    # match the catalogs
    try:
        sdss_cat = ascii.read('/home/boada/Projects/'
                         'planckClusters/data/extern/SDSS/{}/'
                         '{}_SDSS_catalog.csv'.format(tilename, tilename))
        if len(sdss_cat) < 2:
            return
    except FileNotFoundError:
        print('SDSS CATALOG NOT FOUND!')
        return
    c_coord = SkyCoord(ra=cat['RA'] * u.degree, dec=cat['DEC'] * u.degree)
    s_coord = SkyCoord(ra=sdss_cat['ra'] * u.degree, dec=sdss_cat['dec'] *
                       u.degree)
    idxc, idxs, d2d, d3d = s_coord.search_around_sky(c_coord, 1 * u.arcsec)

    # set the fill value
    sdss_cat['specz'].fill_value = -99.0

    # make a new column to catch the results
    zspec = np.ones(len(cat)) * -99.0
    zspec[idxc] = sdss_cat['specz'].filled()[idxs]
    zspec = Column(zspec, name='Z_S')
    cat.add_column(zspec)
    cat.write('tmp.color', format='ascii.commented_header', overwrite=True)

    # mv the old color catalog
    os.rename('{}.color'.format(tilename), '{}.color.orig'.format(tilename))

    with open('tmp.color') as f:
        line = f.readline()
        line = line.split()
        with open('{}.color'.format(tilename), 'wt') as f2:
            f2.write('## ' + time.ctime() + '\n')
            f2.write('## BPZ .columns file for Observation: '
                        '{}\n'.format(tilename))
            f2.write('## (This file was generated by the '
                      'BCS Rutgers pipeline)\n##\n')
            for i, l in enumerate(line[1:]):
                f2.write('{} {} {}\n'.format(line[0], i + 1, l))
            line = f.readline()
            while line:
                f2.write(line)
                line = f.readline()

    # fix the columns file
    with open('{}.columns'.format(tilename), 'at') as f:
        f.write('Z_S {}'.format(i + 1))

    os.remove('tmp.color')

    return

def add_extra_photometry(tilename, filters=['u']):
    ''' This function adds the spec-z's to the color catalog. The previous
    function adds all of the extra information to the sextractor catalogs. So
    you can run either one, or both. Because they don't do the same thing even
    though it looks like it.

    '''

    from astropy.io import ascii
    from astropy.table import Column
    from astropy.coordinates import SkyCoord
    from astropy import units as u

    # read in the sdss catalog. We are doing this first because we only need to
    # do it once
    try:
        cat = ascii.read('/home/boada/Projects/'
                         'planckClusters/data/extern/SDSS/{}/'
                         '{}_SDSS_catalog.csv'.format(tilename, tilename))
        sdss = True
        ps1 = False
        if len(cat) < 2:
            print('# SDSS CATALOG TOO SHORT! -- TRYING PS1')
            sdss = False
    except FileNotFoundError:
        sdss = False
        print('# SDSS CATALOG NOT FOUND! -- TRYING PS1')
    if not sdss:
        try:
            cat = ascii.read('/home/boada/Projects/'
                             'planckClusters/data/extern/PS1/{}/'
                             '{}_PS1_catalog.csv'.format(tilename, tilename))
            ps1 = True
            if len(cat) < 2:
                print('# PS1 CATALOG TOO SHORT!')
                return
        except FileNotFoundError:
            print('# PS1 CATALOG NOT FOUND!')
            return
    # need these coordinates for the matching
    if sdss:
        s_coord = SkyCoord(ra=cat['ra'] * u.degree, dec=cat['dec'] * u.degree)
    elif ps1:
        s_coord = SkyCoord(ra=cat['ramean'] * u.degree, dec=cat['decmean'] *
                           u.degree)
    else:
        return

    color = ascii.read('{}.color'.format(tilename))
    c_coord = SkyCoord(ra=color['RA'] * u.degree, dec=color['DEC'] * u.degree)

    # match the catalogs
    idxc, idxs, d2d, d3d = s_coord.search_around_sky(c_coord, 1 * u.arcsec)

    for filter in filters:
        if sdss:
            # set the fill value
            cat[filter].fill_value = -99.0
            cat['{}_err'.format(filter)].fill_value = 0.

            # make a new column to catch the results
            newPhot = np.ones(len(color)) * -99.0
            newPhot_err = np.zeros(len(color))
            newPhot[idxc] = cat[filter].filled()[idxs]
            newPhot_err[idxc] = cat['{}_err'.format(filter)].filled()[idxs]

            newPhot = Column(newPhot, name='sdss_{}'.format(filter))
            newPhot_err = Column(newPhot_err, name='sdss_{}_err'.format(filter))
            color.add_column(newPhot)
            color.add_column(newPhot_err)
        elif ps1:
            cat[filter].fill_value = -99.0
            cat['{}_err'.format(filter)].fill_value = 0.

            # make a new column to catch the results
            newPhot = np.ones(len(color)) * -99.0
            newPhot_err = np.zeros(len(color))
            newPhot[idxc] = cat['{}meanpsfmag'.format(filter)].filled()[idxs]
            newPhot_err[idxc] = cat['{}meanpsfmagerr'.format(
                                                    filter)].filled()[idxs]

            newPhot = Column(newPhot, name='ps1_{}'.format(filter))
            newPhot_err = Column(newPhot, name='ps1_{}_err'.format(filter))
            color.add_column(newPhot)
            color.add_column(newPhot_err)
        else:
            return

    color.write('tmp.color', format='ascii.commented_header', overwrite=True)

    with open('tmp.color') as f:
        line = f.readline()
        line = line.split()
        with open('{}.color'.format(tilename), 'wt') as f2:
            f2.write('## ' + time.ctime() + '\n')
            f2.write('## BPZ .columns file for Observation: '
                        '{}\n'.format(tilename))
            f2.write('## (This file was generated by the '
                      'BCS Rutgers pipeline)\n##\n')
            for i, l in enumerate(line[1:]):
                f2.write('{} {} {}\n'.format(line[0], i + 1, l))
            line = f.readline()
            while line:
                f2.write(line)
                line = f.readline()

    # fix the columns file
    with open('{}.columns'.format(tilename)) as f:
        line = f.readline()
        with open('tmp.columns', 'w') as f2:
            while line:
                if 'M_0' in line:
                    f2.write('SLOAN-SDSS.{}\t {},{}'.format(filter, i, i + 1))
                    f2.write('\t AB 0.05\t 0.0\n')
                f2.write(line)
                line = f.readline()

    # mv the old color catalog
    os.rename('tmp.columns'.format(tilename), '{}.columns'.format(tilename))
    os.remove('tmp.color')

    return

def write_ds9(d):
    with open('sdss_specz2.reg', 'wt') as f:
        f.write('# Region file formart: DS9 version 4.1\n')
        f.write('# Filename: sdss_specz.reg\n')
        mask = np.where(d['Z_S'] > -1)[0]
        for ra, dec, z in zip(d['RA'][mask], d['DEC'][mask], d['Z_S'][mask]):
            f.write('fk5;circle({},{},2")'.format(ra, dec))
            f.write('# width=2 text={' + str(z) + '} \n')

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

    parser.add_option("--noAstro",
                      action='store_true',
                      dest='noAstro',
                      default=0,
                      help='Whether or not to astro calibrate the mosiacs.')

    parser.add_option("--noPhoto",
                      action='store_true',
                      dest='noPhoto',
                      default=0,
                      help='Whether or not to photo calibrate the mosaics.')

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

    if opt.noSWarp:
        c.get_filenames()

    if not opt.noAstro:
        c.get_astrometry()

    if not opt.noPhoto:
        c.get_zeropt()

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

    if not opt.noBPZ:
        c.BuildColorCat()
        c.runBPZ()

    # make RGB images (pngs)
    if not opt.noRGB:
        print('make rgb')
        c.make_RGB(kband=False)

    # cleanup
    if opt.noCleanUP or opt.noSWarp:
        pass
    else:
        print("CLEANUP!")
        c.cleanup_files()

    elapsed_time(tstart, c.tilename)
    return


main()
