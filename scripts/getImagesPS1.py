from time import sleep
import subprocess
import shlex
import os
import sys
from astLib.astCoords import decimal2hms, decimal2dms
sys.path.append(f'{os.environ["HOME"]}/Projects/planckClusters/catalogs')
from load_catalogs import load_PSZcatalog

def check_redo(outpath, name):
    ''' check to see if we really need to do anything '''

    if not os.path.isdir(f'{outpath}/{name}/'):
        return 1

    files = [f for f in os.walk(f'{outpath}/{name}/')]

    if not len(files[0][2]) == 5:
        return 1

    for f in files[0][2]:
        if 'PSZ' not in f:
            return 1

    return 0

def work(outpath, ra, dec, name):

    ra_sex = decimal2hms(ra, ':')
    dec_sex = decimal2dms(dec, ':')

    cmd = 'panstamps --width=40 --filters grizy '
    cmd += f'--downloadFolder={outpath}/{name} stack {ra_sex} {dec_sex}'
    print(cmd)
    p = subprocess.Popen(shlex.split(cmd))

    p.wait()

def rename(outpath, name):

    if not os.path.isdir(f'{outpath}/{name}/'):
        return

    files = [f for f in os.walk(f'{outpath}/{name}/')]

    if not len(files[0][2]):
        print(f'{name} has no data!')
        return

    for f in files[0][2]:
        if 'PSZ' in f:
            # don't need to rename
            continue
        filt = f[6]

        print(filt)

        # make sure it's in our list
        if filt not in ['g', 'r', 'i', 'z', 'y']:
            print('problem')
        else:
            os.rename(f'{outpath}/{name}/{f}',
                      f'{outpath}/{name}/{name}_PS1stack_{filt}.fits')


data = load_PSZcatalog()

outpath = './PS1'
datapath = '../data/proc2'

for i, row in data.iterrows():
    n = row.NAME.replace(' ', '_')
    #print(n)
    if os.path.isdir(f'{datapath}/{n}'):
        name_us = n
    else:
        try:
            n_psz1 = row.NAME_PSZ1.replace(' ', '_')
        except AttributeError:
            continue
        if os.path.isdir(f'{datapath}/{n_psz1}'):
            name_us = n_psz1
        else:
            continue

    ra = float(row.RA)
    dec = float(row.DEC)

#    if dec < -31:
#        continue

#    if not check_redo(outpath, f'{name_us}'):
#        continue

    if not os.path.isfile(f'{datapath}/{name_us}/{name_us}K.fits'):
        continue

    if os.path.isdir(f'{outpath}/{n}'):
        continue

    work(outpath, ra, dec, f'{n}')
    rename(outpath, f'{n}')
    sleep(1)
