import os
import sys
from astropy.io.fits import getheader

# get the utils from the parent directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import check_exe

# run sextractor
def SEx(self, det_filter='Detec', deblend=False):
    ''' Runs SEXTRACTOR on the mosaicked images created by swarp. It should
    use the zero point created by the photometrypipline.

    '''

    print()
    check_exe("sex")

    # list of output catalog names
    self.combcat = {}

    # The configuration file

    if not deblend:
        self.SExinpar = os.path.join(self.pipeline, 'confs',
                                 'bcs_Catalog.inpar')
    else:
        self.SExinpar = os.path.join(self.pipeline, 'confs',
                                 'bcs_Catalog_deblend.inpar')

    # The detection image that we'll use
    if self.DetImage:
        det_ima = self.DetImage
        print("# Will use default Detection Image:%s " % self.DetImage,
              file=sys.stderr)
    else:
        det_ima = '{}{}.fits'.format(self.tilename, det_filter)
        #det_ima = self.combima[det_filter] + ".fits"
        print("# Will use %s band Detection Image:%s " % (
            det_filter, det_ima), file=sys.stderr)

    self.getbpz = 1  # This var is not really used...

    for filter in self.filters:

        input = self.combima[filter] + ".fits"
        hdr = getheader(input)
        try:
            zeropt = hdr['MAGZERO']
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
        opts += ' -WEIGHT_TYPE MAP_WEIGHT,MAP_WEIGHT '
        opts += ' -WEIGHT_IMAGE %s%s_weight.fits,' % (self.tilename,
                                                      det_filter)
        opts += '%s%s_weight.fits ' % (self.tilename, filter)

        # Do the SEx
        cmd = "sex {},{} -CATALOG_NAME {}".format(det_ima, input, output)
        cmd += " -MAG_ZEROPOINT {} -c {} {} 1>&2".format(zeropt,
                                                self.SExinpar, opts)

        print(cmd)
        if not self.dryrun:
            os.system(cmd)

    return
