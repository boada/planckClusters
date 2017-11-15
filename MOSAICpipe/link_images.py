import os
from glob import glob

files = glob('./PS*/**/PS*.tiff', recursive=True)

for f in files:
    target = f[:-8] + '.tiff'
    try:
        os.symlink(f.split('/')[-1], target)
    except FileExistsError:
        pass

