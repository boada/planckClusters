from astropy.table import Table
from numpy import append as npappend
import os
from pandas import to_numeric

def load_PSZcatalog(unconf=False, full=False, extras=False, **kwargs):
    ''' Load the PSZ catalog data into a pandas dataframe. This is useful for
    getting the catalog data into other scripts in an easy way.

    By default, the script loads all *unique* entries in the combined PSZ1 and
    PSZ2 catalogs. The objects are updated to the PSZ2 values if they appear in
    both catalogs. This should be good enough for most applications where we
    want to include both confirmed and unconfirmed objects.

    Key options:

    unconf = True -- Gives *only* the unconfirmed objects in the PSZ catalogs
    full = True -- Gives the full catalogs instead of just the names and basic
                infomation
    extras = True -- Loads extra information from either (or both) the Barrena
    et al catalog, and denotes where we have mosaic/newfirm imaging.

    **kwargs -- whether to load Barrena AND/OR our catalog. Options are
    `barrena = True` and `us = True`

    returns a pandas dataframe.

    '''

    datapath = f'{os.environ["HOME"]}/Projects/planckClusters/catalogs'

    ps1 = Table.read(f'{datapath}/PSZ1v2.1.fits')
    ps2 = Table.read(f'{datapath}/PSZ2v1.fits')

    # convert to pandas
    df1 = ps1.to_pandas()
    df2 = ps2.to_pandas()

    if unconf:
        # only get unconfirmed sources
        df1 = df1.loc[df1['VALIDATION'] <= 3]
        df2 = df2.loc[df2['VALIDATION'] == -1]

    # clean up strings -- not required
    df1 = df1.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)
    df2 = df2.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)

    # merge the catalogs together
    df_m = df1.merge(df2, how='outer', left_on='INDEX', right_on='PSZ',
                     suffixes=('_PSZ1', '_PSZ2'))

    # get the columns that we want
    if full:
        df_final = df_m
    else:
        cols = df_m.columns[[0, 1, 4, 5, 8, 29, 33, 34, 37, 38, 40, 51]]
        df_final = df_m[cols]

    # remerge to find bits that were missing
    df_final_bigger = df_final.merge(df2, how='left', left_on='INDEX_PSZ1',
                                 right_on='PSZ')
    # fill in nans
    for col in ['NAME', 'RA', 'DEC', 'SNR', 'REDSHIFT', 'INDEX']:
        df_final_bigger[col + '_PSZ2'] = \
            df_final_bigger[col + '_PSZ2'].fillna(df_final_bigger[col])
    # fill in nans
    for col in ['NAME', 'RA', 'DEC', 'SNR', 'REDSHIFT', 'INDEX']:
        df_final_bigger[col + '_PSZ2'] = \
                df_final_bigger[col + '_PSZ2'].fillna(df_final_bigger[col])

    for col in ['NAME', 'RA', 'DEC', 'SNR']:
        df_final_bigger[col] = \
                df_final_bigger[col + '_PSZ2'].fillna(df_final_bigger[
                    col + '_PSZ1'])

    df_final_bigger = df_final_bigger[npappend(
                df_final_bigger.columns[:13].values, ['NAME', 'RA', 'DEC', 'SNR'])]

    if extras:
        df_final_bigger = load_extras(df_final_bigger, **kwargs)

    return df_final_bigger

def load_extras(df, barrena=False, us=False):
    if us:
        df['mosaic'] = False
        df['newfirm'] = False
        for i, n in enumerate(df['NAME']):
            n = n.replace(' ', '_')
            if os.path.isfile(f'../data/proc2/{n}/{n}i.fits'):
                df.iloc[i]['mosaic'] = True
            elif os.path.isfile(f'../data/proc2/'
                                f'{df.iloc[i]["NAME_PSZ1"]}/'
                                f'{df.iloc[i]["NAME_PSZ1"]}i.fits'):
                df.iloc[i]['mosaic'] = True

            if os.path.isfile(f'../data/proc2/{n}/{n}K.fits'):
                df.iloc[i]['newfirm'] = True
            elif os.path.isfile(f'../data/proc2/'
                                f'{df.iloc[i]["NAME_PSZ1"]}/'
                                f'{df.iloc[i]["NAME_PSZ1"]}K.fits'):
                df.iloc[i]['newfirm'] = True

    if barrena:
        datapath = f'{os.environ["HOME"]}/Projects/planckClusters/catalogs'
        # load Barrena table -- Barrena_tbl3.csv
        t = Table.read(f'{datapath}/Barrena_tbl3.csv')
        df_b = t.to_pandas()

        # add a column for and deal with multiBCGs
        df_b['multiBCG'] = False
        df_b.loc[df_b.ID.str.contains('-A', regex=False), 'multiBCG'] = True

        df_b.loc[df_b.ID.str.contains(
            '-A', regex=False), 'ID'] = df_b.loc[df_b.ID.str.contains(
                '-A', regex=False), 'ID'].str.replace('-A', '')

        df_b.drop(df_b.loc[df_b.ID.str.contains('-B', regex=False)].index,
                inplace=True)
        df_b.drop(df_b.loc[df_b.ID.str.contains('-C', regex=False)].index,
                inplace=True)

        # merge the dataframes together
        df_b['ID'] = to_numeric(df_b['ID'])
        df_extra = df.merge(df_b, how='left', left_on='INDEX_PSZ1',
                            right_on='ID', suffixes=('', '_Barrena'))

        # clean things up a bit
        df_b.drop(['Planck Name', 'SZ S/N', 'Notes', 'ID'], axis=1,
                  inplace=True)
    return df_extra
