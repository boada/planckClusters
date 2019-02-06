import os

def load_PSZcatalog(unconf=False):
    from astropy.table import Table
    from numpy import append as npappend

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
    for col in ['NAME', 'RA', 'DEC']:
        df_final_bigger[col] = \
            df_final_bigger[col + '_PSZ2'].fillna(df_final_bigger[col + '_PSZ1'])

    df_final_bigger = \
        df_final_bigger[npappend(df_final_bigger.columns[:12].values, ['NAME',
                                                                   'RA',
                                                                   'DEC'])]

    return df_final_bigger


base_string = ('https://ned.ipac.caltech.edu/cgi-bin/objsearch?in_csys='
               'Equatorial&in_equinox=J2000.0&lon={:.6f}d&lat={:.6f}d'
               '&radius=5&hconst=73&omegam=0.27&omegav=0.73&corr_z=1'
               '&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&in_'
               'objtypes1=Galaxies&in_objtypes1=GPairs&in_objtypes1=GTriples&'
               'in_objtypes1=GGroups&in_objtypes1=GClusters&search_type='
               'Near+Position+Search&nmp_op=ANY&out_csys=Equatorial&out_'
               'equinox=J2000.0&obj_sort=Distance+to+search+center&of=pre_text'
               '&zv_breaker=30000.0&list_limit=5&img_stamp=YES')

data = load_PSZcatalog(True)
data['NED'] = ''


for index, row in data.iterrows():
    data['NED'][index] = base_string.format(row['RA'], row['DEC'])

data.to_csv('updated.csv')
