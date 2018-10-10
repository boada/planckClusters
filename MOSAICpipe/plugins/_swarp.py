import os
import sys
from astropy.io import fits
from astropy.io.fits import getheader

# get the utils from the parent directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import (check_exe, dec2sex, weight_from_dqfile, put_airmass,
                    put_exptime, put_headerKeywords)


def center_dither(self, conf="SWarp-center.conf", filter='i', dryrun=False):
    ''' Center the dither pattern for SWARP '''

    print()
    check_exe("swarp")

    # The configuration file
    center_conf = os.path.join(self.pipeline, 'confs', conf)
    self.outdir = os.path.join(self.outpath, self.tilename)
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

    if self.dryrun or not dryrun:
        print("# Centering mosaic", file=sys.stderr)
        print(cmd)
        os.system(cmd)
    else:
        print('# Mosaic already centered', file=sys.stderr)
        if not os.path.isfile(opts['IMAGEOUT_NAME']):
            opts['IMAGEOUT_NAME'] = '{}{}.fits'.format(self.tilename, filter)

    # Read in the header
    print(opts['IMAGEOUT_NAME'])
    header = getheader(opts["IMAGEOUT_NAME"])
    nx = header["NAXIS1"]
    ny = header["NAXIS2"]
    x_center = header["CRVAL1"]
    y_center = header["CRVAL2"]

    x_center = dec2sex(float(x_center) / 15)
    y_center = dec2sex(float(y_center))

    pixscale = abs(float(header['CD1_1']) * 3600)

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
        #if filter == 'i':
        outweight = "%s%s_weight.fits" % (self.tilename, filter)
        #else:
        #    outweight = "coadd.weight.fits"

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

def swarp_extras(self,
                conf="SWarp-common.conf",
                dryrun=None,
                combtype="MEDIAN",
                reSWarp=None,
                filters=None,
                newfirm=False):
    # Make the detection image an colors

    # Keys to keep
    keywords = ['OBJECT', 'EXPTIME', 'AIRMASS', 'TIMESYS', 'DATE-OBS',
                'TIME-OBS', 'OBSTYPE', 'OBSERVAT', 'TELESCOP', 'HA', 'ZD',
                'DETECTOR', 'DARKTIME', 'RA', 'DEC', 'MJD-OBS', 'INSTRUME',
                'FILTER', 'MAGZERO', 'MAGZSIG']

    if not filters:
        filters = self.filters

    _conf = os.path.join(self.pipeline, 'confs', conf)

    pars = {}
    pars["IMAGE_SIZE"] = "%s,%s" % (self.nx, self.ny)
    pars["CENTER_TYPE"] = "MANUAL"
    pars["CENTER"] = "%s,%s" % (self.xo, self.yo)
    pars["PIXEL_SCALE"] = self.pixscale
    pars["PIXELSCALE_TYPE"] = "MANUAL"
    pars["COPY_KEYWORDS"] = ','.join(keywords)
    pars["NTHREADS"] = "0"
    pars["WEIGHT_TYPE"] = "MAP_WEIGHT"
    #pars["COMBINE_TYPE"]    = combtype

    # Color Channels
    filters = {}

    try:
        self.filters.index('K')

        if not newfirm:
            raise(ValueError)

        filters["Red"] = ['K']
        filters["Green"] = ['z', 'i']
        filters["Blue"] = ['g', 'r']
        filters["Detec"] = ['i', 'K']
        ir = True
    except ValueError:
        filters["Red"] = ['i']
        filters["Green"] = ['r']
        filters["Blue"] = ['g']
        filters["Detec"] = ['i', 'z']
        ir = False

    for color in list(filters.keys()):

        outimage = "%s%s.fits" % (self.tilename, color)
        outweight = "%s%s_weight.fits" % (self.tilename, color)
        # Store the names
        self.combima[color] = "%s%s" % (self.tilename, color)
        # Different combination for detection image"
        if color == "Detec":
            pars["COMBINE_TYPE"] = "CHI_OLD"
            self.DetImage = outimage
        else:
            pars["COMBINE_TYPE"] = combtype

        # The options for SWarp
        opts = ""
        for param, value in list(pars.items()):
            opts += "-%s %s " % (param, value)

        print(color, pars["COMBINE_TYPE"])

        # Make the list of images to combine
        try:
            del filelist, weightlist
        except UnboundLocalError:
            pass
        for filter in filters[color]:
            if not os.path.isfile("%s%s.fits" % (self.tilename, filter)):
                continue
            try:
                filelist += ' ' + ' '.join(["%s%s.fits" % (
                                                    self.tilename, filter)])
                weightlist += ',' + ','.join(["%s%s_weight.fits" % (
                                                    self.tilename, filter)])
            except UnboundLocalError:
                filelist = ' '.join(["%s%s.fits" % (
                                                    self.tilename, filter)])
                weightlist = ','.join(["%s%s_weight.fits" % (
                                                    self.tilename, filter)])

        print(color, filelist)
        cmd = "swarp %s -c %s -IMAGEOUT_NAME %s -WEIGHTOUT_NAME %s " %\
                        (filelist, _conf, outimage, outweight)
        cmd += " -WEIGHT_IMAGE %s " % (weightlist)
        cmd += opts

        if not dryrun:
            if color == 'Red':
                cmd = 'cp {} {}'.format(filelist, outimage)
            elif not ir and color != 'Detec':
                cmd = 'cp {} {}'.format(filelist, outimage)
            print(cmd)
            os.system(cmd)
        else:
            print(cmd)

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
            try:
                os.system("rm %s" % file)
            except FileNotFoundError:
                pass

    return

# Generate the mask from the weight
def make_swarp_input_weights(self, clobber=True):
    ''' Generate the masks from the weight images provided by the VO '''

    for fname, dqmask in list(self.dqmask.items()):
        outname = "%s.weight.fits" % os.path.splitext(fname)[0]
        print("# Generating weight image from %s --> %s" % (
            dqmask, outname), file=sys.stderr)
        weight_from_dqfile(dqmask, outname, clobber=clobber)

        # Here we are going to make a fix for the bad k-band data
        # until the VO can actually provide a usuable fix for things.
        # THis should be commented out when not working specifically with
        # the bad k-band data. I'll make that easy with a if True: statement

        if False:
            with fits.open(fname, mode='readonly') as orig:
                if orig[0].header['instrume'] == 'newfirm':

                    # make a science image mask
                    sci_data = orig[0].data
                    vgap = sci_data[:, 2000:2250]
                    mask = (-0.12 < vgap) & (vgap < 0.12)
                    vgap[mask] = 0 # the vertical chip gap
                    hgap = sci_data[2000:2250, :]
                    mask = (-0.12 < hgap) & (hgap < 0.12)
                    hgap[mask] = 0 # the horizontal chip gap
                    sci_data[2000:2250, :] = hgap
                    sci_data[:, 2000:2250] = vgap

                    print('Updating the weight map.. {}'.format(outname))
                    with fits.open(outname, mode='readonly') as hdu:
                        data = hdu[0].data
                        data = np.where(sci_data == 0, 0, data)
                        nhdu = fits.PrimaryHDU(data, header=hdu[0].header)
                    nhdu.writeto(outname, clobber=True)

    return
