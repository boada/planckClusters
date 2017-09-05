from astropy.io import fits
from glob import glob
import sys
import os

if os.path.isdir(sys.argv[1]):
    files = glob('{}/*.fits.fz'.format(sys.argv[1]))

    for f in files:
        with fits.open(f, mode='update') as hdus:
            for hdu in hdus[1:]:
                del hdu.header['PV*']
                if 'opd' in f:
                    hdu.header['PRODTYPE'] = 'dqmask'


