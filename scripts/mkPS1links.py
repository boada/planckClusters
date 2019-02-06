import os
from astropy import coordinates

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


base_string = ('http://ps1images.stsci.edu/cgi-bin/ps1cutouts?'
               'pos={}%2C{}&'
               'filter=color&filter=g&filter=r&filter=i&filetypes=stack&'
               'auxiliary=data&size=1200&output_size=0&verbose=0&'
               'autoscale=99.500000&catlist=')

data = load_PSZcatalog(True)
data['PanSTARRS'] = ''

for index, row in data.iterrows():
    if row.DEC < -31:
        continue
    else:
        data['PanSTARRS'][index] = base_string.format(row.RA, row.DEC)
data.to_csv('updated_withPS1.csv', index=False)








