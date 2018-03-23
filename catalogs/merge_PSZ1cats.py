from astropy.table import Table

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

# merge
df_m = df1_uc.merge(df2_uc, how='outer', left_on='INDEX', right_on='PSZ',
                    suffixes=('_PSZ1', '_PSZ2'))

# get the columns we want
cols = df_m.columns[[0, 1, 4, 5, 29, 33, 34, 37, 38, 40]]
df_final = df_m[cols]


