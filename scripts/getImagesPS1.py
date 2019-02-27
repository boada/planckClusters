from time import sleep
import subprocess
import shlex
import os
from astropy.table import Table
from numpy import append as npappend
from astLib.astCoords import decimal2hms, decimal2dms

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

    cmd = 'panstamps --width=10 --filters grizy '
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

    filters = ['g', 'r', 'i', 'z', 'y']

    for filt, f in zip(filters, files[0][2]):
        if 'PSZ' in f:
            # don't need to rename
            continue
        elif filt in f:
            os.rename(f'{outpath}/{name}/{f}',
                      f'{outpath}/{name}/{name}_PS1stack_{filt}.fits')
        else:
            print(f'Missing the {filt} filter!')

def load_PSZcatalog():

    datapath = './../../planckClusters/catalogs/'

    ps1 = Table.read(f'{datapath}/PSZ1v2.1.fits')
    ps2 = Table.read(f'{datapath}/PSZ2v1.fits')

    # convert to pandas
    df1 = ps1.to_pandas()
    df2 = ps2.to_pandas()

    # clean up strings -- not required
    df1 = df1.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)
    df2 = df2.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)

    # merge the catalogs together
    df_m = df1.merge(
        df2,
        how='outer',
        left_on='INDEX',
        right_on='PSZ',
        suffixes=('_PSZ1', '_PSZ2'))

    # get the columns that we want
    cols = df_m.columns[[0, 1, 4, 5, 8, 29, 33, 34, 37, 38, 40, 51]]
    df_final = df_m[cols]

    # remerge to find bits that were missing
    df_final_bigger = df_final.merge(
        df2, how='left', left_on='INDEX_PSZ1', right_on='PSZ')
    # fill in nans
    for col in ['NAME', 'RA', 'DEC', 'SNR', 'REDSHIFT', 'INDEX']:
        df_final_bigger[col + '_PSZ2'] = df_final_bigger[col + '_PSZ2'].fillna(
            df_final_bigger[col])
    # fill in nans
    for col in ['NAME', 'RA', 'DEC', 'SNR', 'REDSHIFT', 'INDEX']:
        df_final_bigger[col + '_PSZ2'] = df_final_bigger[col + '_PSZ2'].fillna(
            df_final_bigger[col])
    for col in ['NAME', 'RA', 'DEC']:
        df_final_bigger[col] = df_final_bigger[col + '_PSZ2'].fillna(
            df_final_bigger[col + '_PSZ1'])

    df_final_bigger = df_final_bigger[npappend(
        df_final_bigger.columns[:12].values, ['NAME', 'RA', 'DEC'])]

    return df_final_bigger


data = load_PSZcatalog()

outpath = './../data/extern/PS1'

for i, (ra, dec, name) in enumerate(
        zip(data['RA'], data['DEC'], data['NAME'])):

    print(name)
    name = name.replace(' ', '_')

    if dec < -31:
        continue

    if not check_redo(outpath, f'{name}'):
        continue

    work(outpath, ra, dec, f'{name}')
    rename(outpath, f'{name}')
    sleep(1)

