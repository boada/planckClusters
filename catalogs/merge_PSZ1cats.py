from astropy.table import Table
from numpy import append as npappend

ps1 = Table.read('./PSZ1v2.1.fits')
ps2 = Table.read('./PSZ2v1.fits')

# convert to pandas
df1 = ps1.to_pandas()
df2 = ps2.to_pandas()

# clean up strings -- not required
df1 = df1.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)
df2 = df2.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)

# only get unconfirmed sources
df1_uc = df1.loc[df1['VALIDATION'] <= 3]
df2_uc = df2.loc[df2['VALIDATION'] == -1]


#df_m = df2_uc.merge(df1_uc, how='outer', left_on='PSZ', right_on='INDEX',
#                    suffixes=('_PSZ2', '_PSZ1'))

# merge
df_m = df1_uc.merge(df2_uc, how='outer', left_on='INDEX', right_on='PSZ',
                    suffixes=('_PSZ1', '_PSZ2'))

# get the columns we want
# PSZ1 first
cols = df_m.columns[[0, 1, 4, 5, 8, 29, 33, 34, 37, 38, 40, 51]]

# PSZ2 first
#cols = df_m.columns[[0, 1, 4, 5, 7, 30, 31, 34, 35, 59]]

df_final = df_m[cols]

# remerge to find bits that were missing
df_final_bigger = df_final.merge(df2, how='left', left_on='INDEX_PSZ1',
                                 right_on='PSZ')

# fill in nans
for col in ['NAME', 'RA', 'DEC', 'SNR', 'REDSHIFT', 'INDEX']:
    df_final_bigger[col+'_PSZ2'] =\
        df_final_bigger[col+'_PSZ2'].fillna(df_final_bigger[col])

# finally make single columns column

for col in ['NAME', 'RA', 'DEC']:
    df_final_bigger[col] =\
        df_final_bigger[col+'_PSZ2'].fillna(df_final_bigger[col+'_PSZ1'])

df_final_bigger = df_final_bigger[npappend(df_final_bigger.columns[:12].values,
                                           ['NAME', 'RA', 'DEC'])]


