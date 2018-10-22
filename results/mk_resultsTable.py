from astropy.table import Table
from get_results import loadClusters, loadMembers
from astLib import astCoords
from numpy import nan

def main(round=3):
    ''' This creates a simple latex table with the results using the columns
    specified below. A few things will need to be done by hand after this is
    created. The corrected errors to redshifts, and corrected Ngals will need
    to be manually added by looking at the results images. The cluster finder
    doesn't save those columns by default.

    '''

    # the confirmed = True gets the 12 confirmed clusters
    results = loadClusters(round=round, confirmed=True)

    # load the master spreadsheet
    t_ss = Table.read('../catalogs/PSZ2_unconfirmed_catalog - current.csv')
    df_ss = t_ss.to_pandas()

    observed = df_ss.loc[~df_ss['MOSAIC Imaging'].isnull()]

    confirmed = observed.merge(results, left_on='Name', right_on='Cluster',
                               how='left')

    # load the PSZ1 catalog
    t_1 = Table.read('../catalogs/PSZ1v2.1.fits')
    df1 = t_1.to_pandas()

    # now we add all of the extra info
    # Berrena paper
    Bpaper = Table.read('../papers/1803.05764/Barrena_tbl3.csv')
    df_paper = Bpaper.to_pandas()

    complete = confirmed.merge(df_paper, left_on='Name', right_on='Planck Name',
                               how='left')

    # tack on the PSZ1 catalog
    complete = complete.merge(df1, left_on='PSZ1 Indx', right_on='INDEX',
                              how='left')

    # combine some columns to get extra info
    complete['z_extern'] = complete['PSZ1 Redshift'].fillna(complete['z_cl'])
    complete['S/N'] = complete['SNR_PSZ2'].fillna(complete['SNR_PSZ1'])

    # get the columns we want
    cols = ['Name', 'PSZ1 Indx', 'PSZ2 Indx', 'S/N', 'RA_SEX', 'DEC_SEX',
            'dist_BCG', 'z_cl_boada', 'z_clerr_boada', 'Ngal_c_boada', 'Cmag_i',
            'z_extern', 'REDSHIFT_SOURCE', 'BCG_boada']

    c = complete.loc[:, cols]

    c.loc[:, ('RA_SEX', 'DEC_SEX')] = nan

    c.loc[(~c['z_extern'].isnull()) & (c['REDSHIFT_SOURCE'] == -1.0),
          'REDSHIFT_SOURCE'] = 99

    # let's the the RA/DEC of our BCGs
    m = c.loc[c['BCG_boada'].notnull()]

    for i, row in m.iterrows():
        mems = loadMembers('boada', row['Name'], round=round)
        try:
            ra = mems.loc[mems['ID'] == row['BCG_boada'], 'RA'].values[0]
        except IndexError:
            continue
        dec = mems.loc[mems['ID'] == row['BCG_boada'], 'DEC'].values[0]
        # convert to sexigesimal
        ra_sex = astCoords.decimal2hms(ra, ':')
        dec_sex = astCoords.decimal2dms(dec, ':')

        # compute the distance to the PSZ position
        psz_ra = complete.iloc[i]['RA_x']
        psz_dec = complete.iloc[i]['DEC_x']

        m.loc[i, 'dist_BCG'] = astCoords.calcAngSepDeg(psz_ra,
                                                       psz_dec, ra, dec) * 60

        # write it back into the main frame
        m.loc[i, 'RA_SEX'] = ra_sex
        m.loc[i, 'DEC_SEX'] = dec_sex

    return m

    #m.to_latex('results_table.tex', index=False, float_format='%0.2f')


if __name__ == "__main__":
    main()
