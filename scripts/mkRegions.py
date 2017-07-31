#def mk_regions(coords):
import os
from numpy import genfromtxt
import sys


def mk_regions(dir):
    if not os.path.isdir(dir):
        print('argument is not a directory')
        return

    for coords in os.listdir(dir):
        if coords.endswith('.csv'):
            print(coords)
            d = genfromtxt('{}/{}'.format(dir, coords),
                           names=True,
                           dtype=None,
                           delimiter=',')
            with open(coords.rstrip('csv') + 'reg', 'wt') as f1:
                f1.writelines('# Region file formart: DS9 version 4.1\n')
                f1.writelines('# Filename: ' + coords.rstrip('txt') + 'reg\n')
                f1.writelines('global color=cyan\n')

                for ra, dec in zip(d['ramean'], d['decmean']):
                    f1.writelines('fk5;circle({}, {}, 5")'.format(ra, dec))
                    f1.writelines('# width=2')
                    #f1.writelines('tag={' + coords.rstrip('') + '}\n')

            f1.close()


def mk_regions_array(output, ra, dec, info, tag='regions'):
    with open(output, 'wt') as f:
        f.writelines('# Region file formart: DS9 version 4.1\n')
        f.writelines('# Filename: {}'.format(output))
        f.writelines('global color=cyan\n')

        for r, d, i in zip(ra, dec, info):
            f.writelines('fk5;circle({}, {}, 5")'.format(r, d))
            f.writelines('# width=2 ')
            f.writelines('text={{{:.3f}}} '.format(i))
            f.writelines('tag={{{}}}\n'.format(tag))


if __name__ == "__main__":
    mk_regions(sys.argv[1])
