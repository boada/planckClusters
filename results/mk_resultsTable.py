from astropy.io import ascii
from astropy.table import vstack
from numpy import sort

def main():
    ''' This creates a simple latex table with the results using the columns
    specified below. A few things will need to be done by hand after this is
    created. The corrected errors to redshifts, and corrected Ngals will need
    to be manually added by looking at the results images. The cluster finder
    doesn't save those columns by default.

    '''

    high_conf = ['PSZ2_G145.25+50.84',
                'PSZ2_G120.76+44.14',
                'PSZ2_G305.76+44.79',
                'PSZ2_G029.66-47.63',
                'PSZ2_G173.76+22.92',
                'PSZ1_G224.82+13.62',
                'PSZ2_G048.47+34.86',
                'PSZ2_G106.11+24.11',
                'PSZ1_G084.62-15.86',
                'PSZ2_G125.55+32.72',
                'PSZ2_G043.44-41.27',
                'PSZ2_G096.43-20.89']

    results_dir = '/home/boada/Projects/planckClusters/results/boada'

    hc = sort(high_conf)

    results = [ascii.read('{}/{}/{}/{}.info'.format(results_dir, c, c, c)) for c in
                        hc]

    table = vstack(results)
    # get the columns we want
    t = table[['RA', 'DEC', 'zBCG', 'z_cl', 'Ngal']]
    # add the cluster names
    t['Cluster'] = hc

    # reorder things
    t = t[['Cluster', 'RA', 'DEC', 'zBCG', 'z_cl', 'Ngal']]


if __name__ == "__main__":
    main()
