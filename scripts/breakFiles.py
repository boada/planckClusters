from __future__ import print_function
from astropy.io import fits
from glob import glob

files = glob('k4m*.fits.fz')

for f in files:
    hdulist = fits.open(f)
    for i, hdu in enumerate(hdulist):
        if i == 0:
            continue
        print(f)
        f = f.rstrip('.fits.fz')
        hdu.writeto(f+'_ccd{}.fits.fz'.format(str(i)))
