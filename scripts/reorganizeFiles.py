from astropy.io import fits
from glob import glob
import os

dataDir = '../data'
files = glob('{}/2*/*/*.fits.fz'.format(dataDir))

for f in files:
    fName = f.split('/')[-1]
    dirName = f.split('/')[-2]

    print(f)
    with fits.open(f) as oimg:
        try:
            obj = oimg[1].header['object']
            tel = oimg[1].header['instrume'].split('_')[0]
        except IndexError:
            print('{} only has one extension'.format(f))
            continue
        except KeyError:
            print("{} isn't a image file?".format(f))
            continue
        # mk directories recursively
        if not os.path.isdir('./{}/{}/{}'.format(dataDir, obj, tel)):
            os.makedirs('./{}/{}/{}'.format(dataDir, obj, tel))
        try:
            os.symlink('../../{}/{}'.format(dirName, fName),
                   './{}/{}/{}/{}'.format(dataDir, obj, tel, fName))
        except FileExistsError:
            pass
