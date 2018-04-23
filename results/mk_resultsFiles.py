from astropy.io import ascii
from astropy.table import vstack
from numpy import int64
from glob import glob
from sys import argv


def mk_results(round=1, user='boada'):
    ''' This aggregates the *.info files into single files so we don't have to
    read the individual files in every single time. This file should only need
    to be run once per user. If there has been a significant amount of changes,
    then you should probably run it again.

    >>> python3 mk_resultsFiles.py 1 boada

    '''

    files = glob('{}/{}/**/*.info'.format(round, user), recursive=True)

    t = []
    for f in files:
        d = ascii.read(f)
        if isinstance(d['ID_BCG'][0], int64):
            d['ID_BCG'] = d['ID_BCG'].astype('U23')
            d['ID_BCG'] = '{}_0'.format(f.split('/')[0])
            d['RA'] = d['RA'].astype('U10')
            d['DEC'] = d['DEC'].astype('U10')
        t.append(d)

    table = vstack(t)

    table.write('{}/{}/{}_results.csv'.format(round, user, user),
                format='ascii.csv')


if __name__ == "__main__":
    try:
        round = argv[1]
        user = argv[2]
    except IndexError:
        print('specify the `round` and `user`.')
        print('>>> python3 mk_resultsFiles.py 1 boada')
