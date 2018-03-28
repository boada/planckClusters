from astropy.io import ascii
from astropy.table import Table, vstack
from numpy import sort
from get_results import loadClusters

def main():
    ''' This creates a simple latex table with the results using the columns
    specified below. A few things will need to be done by hand after this is
    created. The corrected errors to redshifts, and corrected Ngals will need
    to be manually added by looking at the results images. The cluster finder
    doesn't save those columns by default.

    '''

    # the confirmed = True gets the 12 confirmed clusters
    results = loadClusters(confirmed=True)

    # load the master spreadsheet
    t_ss = Table.read('../catalogs/PSZ2_unconfirmed_catalog - current.csv')
    df_ss = t_ss.to_pandas()

    observed = df_ss.loc[~df_ss['MOSAIC Imaging'].isnull()]

    confirmed = observed.merge(results, left_on='Name', right_on='Cluster',
                               how='left')


if __name__ == "__main__":
    main()
