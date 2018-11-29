import os
import sys
from astropy.io.fits import getheader
import numpy as np

# get the utils from the parent directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import check_exe


def calc_airmass(self, filter, combtype='median'):
    ''' Calculates the mean airmass for a set of observations. '''

    x = [self.airmass[file] for file in self.files[filter]]
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

    x = [self.exptime[file] for file in self.files[filter]]
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
            try:
                os.remove(f)
            except FileNotFoundError:
                continue
            print('', end='.')
    for filtername in self.filters:
        for f in self.files_weight[filtername]:
            try:
                os.remove(f.rstrip("'[0]'"))
            except FileNotFoundError:
                continue
            print('', end='.')
    print('')
    return

def make_RGB(self, newfirm=False, conf='stiff-common.conf'):

    try:
        print()
        check_exe('stiff')
    except FileNotFoundError:
        return

    _conf = os.path.join(self.pipeline, 'confs', conf)

    if newfirm:
        red = './{}{}.fits'.format(self.tilename, 'Red')
        green = './{}{}.fits'.format(self.tilename, 'Green')
        blue = './{}{}.fits'.format(self.tilename, 'Blue')
        bands = 'rgb'
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
