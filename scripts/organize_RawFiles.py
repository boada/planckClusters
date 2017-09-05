from astropy.io import fits
from glob import glob
import os

''' This creates a directory that has all of the raw files (links to) sorted
into directories based on their type. This doesn't really do anything useful
except make it easy to find a specific file based on it's object name.

'''


dataDir = '.'
files = glob('*.fits.fz')

for f in files:
    print(f)
    with fits.open(f) as oimg:
        wanted = False
        try:
            obj = oimg[0].header['object']
            if 'flat' in obj.lower():
                pass
            elif 'test' in obj.lower():
                continue
            elif 'focus' in obj.lower():
                continue
            elif 'dark' in obj.lower():
                pass
            else:
                pass
        except IndexError:
            print('{} only has one extension'.format(f))
            continue
        except KeyError:
            print("{} isn't a image file?".format(f))
            continue
        # mk directories recursively
        if not os.path.isdir('./sorted/{}'.format(obj)):
            os.makedirs('./sorted/{}'.format(obj))
        try:
            os.symlink('../../{}'.format(f),
                   './sorted/{}/{}'.format(obj, f))
        except FileExistsError:
            pass
