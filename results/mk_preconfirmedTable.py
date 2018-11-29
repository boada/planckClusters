from astropy.table import Table
from get_results import loadClusters
from astLib import astCoords
from numpy import nan
import pandas as pd

''' This creates a simple latex table with the results using the columns
specified below. A few things will need to be done by hand after this is
created. The corrected errors to redshifts, and corrected Ngals will need
to be manually added by looking at the results images. The cluster finder
doesn't save those columns by default.

'''

# the confirmed = True gets the 12 confirmed clusters
results = loadClusters(round=2)

# load the master spreadsheet
t_ss = Table.read('../catalogs/PSZ2_unconfirmed_catalog - current.csv')
df_ss = t_ss.to_pandas()

observed = df_ss.loc[df_ss['MOSAIC Imaging'].notnull()]

confirmed = observed.merge(results, left_on='Name', right_on='Cluster',
                            how='left')

# load the PSZ1 catalog
t_1 = Table.read('../catalogs/PSZ1v2.1.fits')
df1 = t_1.to_pandas()

# now we add all of the extra info
# Berrena paper
Bpaper = Table.read('../papers/1803.05764/Barrena_tbl3.csv')
df_paper = Bpaper.to_pandas()

# do this to not match on values of NAN
df_paper = df_paper.loc[df_paper['Planck Name'].notnull()]

complete = confirmed.merge(df_paper, left_on='PSZ1 Name', right_on='Planck Name',
                            how='left')

# tack on the PSZ1 catalog
complete = complete.merge(df1, left_on='PSZ1 Indx', right_on='INDEX',
                            how='left')

# combine some columns to get extra info
complete['z_extern'] = complete['PSZ1 Redshift'].fillna(complete['z_cl'])
complete['z_extern'] = complete['z_extern'].fillna(complete['zphot'])
complete['S/N'] = complete['SNR_PSZ2'].fillna(complete['SNR_PSZ1'])

t_pir = Table.read('../catalogs/PIRXXXVI.fits')
df_pir = t_pir.to_pandas()

df_pir['ID'] += 1
complete = complete.merge(df_pir, left_on='PSZ1 Indx', right_on='ID',
                        suffixes=('_cat', '_pir'), how='left')

m = complete.loc[complete['z_extern'].notnull()]

# PSZ info
table = pd.DataFrame(index=m.index)
table['Name'] = m['Name']
table[['RA PSZ', 'DEC PSZ']] = m[['RA_SEX', 'DEC_SEX']]

# our info
table['RA'] = nan
table['DEC'] = nan
table['Dist'] = nan
table['z'] = m['zBCG_boada']
table['Mag Lim'] = m['Cmag_i']
#table['N'] = nan

# extern info
table['RA EX'] = m['RA_y']
table['DEC EX'] = m['DEC_y']
table['Dist EX'] = nan
table['z EX'] = m['z_extern']
#table['N EX'] = m['R']
table['Flag'] = m['Flag']
table['Ref'] = m['REDSHIFT_SOURCE']

fixed1 = m.loc[pd.notnull(m['RAJ2000']),
                'RAJ2000'].apply(lambda x: astCoords.decimal2hms(x, ':'))
fixed2 = m.loc[pd.notnull(m['DEJ2000']),
                'DEJ2000'].apply(lambda x: astCoords.decimal2dms(x, ':'))

table.loc[fixed1.index, 'RA EX'] = fixed1
table.loc[fixed2.index, 'DEC EX'] = fixed2

for i, row in table.iterrows():
    table.loc[i, 'Dist EX'] = astCoords.calcAngSepDeg(row['RA PSZ'],
                                                      row['DEC PSZ'],
                                                      row['RA EX'],
                                                      row['DEC EX']) * 60

print('done')

