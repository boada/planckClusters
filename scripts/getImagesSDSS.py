import numpy as np
from astropy import coordinates as coords
from astroquery.sdss import SDSS
import os
from itertools import cycle

# get file data
data = np.genfromtxt('../../PSZ2_unconfirmed_catalog - Master.csv',
           delimiter=',', names=True, dtype=None)

# check to see if we are continuing
try:
    with open('imagesGot.txt', 'r') as f:
        done = int(f.readline())
except FileNotFoundError:
    done = -np.inf

for i, (ra, dec, name) in enumerate(zip(data['RA'], data['Dec'],
                                        data['Name'])):
    # make string comparisons work
    name = name.decode()

    # here's the little bit for persistance (see below)
    if i <= done:
        continue
    # only do the things in the footprint
    if not data['SDSS_Footprint'][i].decode() == 'TRUE':
        continue
    print(name)
    if not os.path.isdir(name):
        os.mkdir(name)

    pos = coords.SkyCoord(ra, dec, frame='icrs', unit='deg')
    print('fetching images....')
    imgs = SDSS.get_images(pos, radius="20'", band='ugriz', data_release=12)

    print('writing %s' % name, end='')
    counter = 0
    for HDU, band in zip(imgs, cycle('ugriz')):
        print('.', end='')
        if os.path.isfile(f'./{name}/{name}_sdss_{band}_{counter}.fits'):
            HDU.writeto(f'./{name}/{name}_sdss_{band}_{counter + 1}.fits',
                        clobber=True)
            if band == 'z':
                counter += 1
        else:
            HDU.writeto(f'./{name}/{name}_sdss_{band}_{counter}.fits',
                        clobber=True)

    # add a little thing to keep track how many we have done. Poor man's
    # persistance
    with open('imagesGot.txt', 'w') as f:
        f.writelines(f'{i}')

    print('')


