import pandas as pd
import numpy as np

r''' Reads in the latex table from this paper. You have to clean up the
latex before you try to read it. Basically delete all of the extra stuff
that isn't the actual data... \hline, \begin{..}, etc. Delete all of the
\cr and footnotes. Then things will read and be really nice.

'''

df = pd.read_csv('table2_psz1_2.tex',
                sep='&',
                header=None,
                engine='python', comment='%')

for i in range(12):
    df[i] = df[i].map(lambda x: x.replace('$', ''))
    df[i] = df[i].map(lambda x: x.replace('*', ''))
    df[i] = df[i].map(lambda x: x.strip())

df.columns = ['ID', 'Planck Name', 'SZ S/N', 'RA', 'DEC', 'Dist', 'z_cl',
                'Nspec', 'zphot', 'R', 'Flag', 'Notes']

df['Planck Name'] = df['Planck Name'].map(lambda x: x.replace('Z1 G', 'Z1_G'))
df['RA'] = df['RA'].map(lambda x: x.replace(' ', ':'))
df['DEC'] = df['DEC'].map(lambda x: x.replace(' ', ':'))

for c in df.columns:
    df.loc[df.loc[:, c].str.strip() == "-", c] = np.nan

df.insert(7, 'z_bcg', np.nan)
df.insert(9, 'zphot_err', np.nan)
df.insert(10, 'R_err', np.nan)

for i, j in df['z_cl'].iteritems():
    try:
        parts = j.split('~;~')
    except AttributeError:
        continue
    df.loc[i, 'z_cl'] = parts[0]
    df.loc[i, 'z_bcg'] = parts[-1]

for i, j in df['zphot'].iteritems():
    try:
        parts = j.split('\pm')
    except AttributeError:
        continue
    df.loc[i, 'zphot'] = parts[0]
    df.loc[i, 'zphot_err'] = parts[-1]

for i, j in df['R'].iteritems():
    try:
        parts = j.split('\pm')
    except AttributeError:
        continue
    df.loc[i, 'R'] = parts[0]
    df.loc[i, 'R_err'] = parts[-1]

if True:
    df.to_csv('Barrena_tbl3.csv', index=False)
