import pandas as pd

# read in the catalogs
psz1 = pd.read_csv('./PSZ2_unconfirmed_catalog - PSZ1 catalog.csv')
psz2 = pd.read_csv('./PSZ2_unconfirmed_catalog - PSZ2 catalog.csv')
swift = pd.read_csv('./PSZ2_unconfirmed_catalog - SWIFT program.csv')

# merge psz2 into psz1 based on the psz1 name
psz12 = pd.merge(psz2, psz1, left_on='PSZ1 Name', right_on='Name', how='outer')

# drop the old PSZ1 Name
del psz12['PSZ1 Name']

# merge all the main names
psz12['Name_x'].fillna(psz12['Name_y'], inplace=True)

# merge RA and DEC
psz12['RA_x'].fillna(psz12['RA_y'], inplace=True)
psz12['DEC'].fillna(psz12['Dec'], inplace=True)

# merge rass rates
psz12['RASS rate_x'].fillna(psz12['RASS rate_y'], inplace=True)
psz12['RASS S/N_x'].fillna(psz12['RASS S/N_y'], inplace=True)
psz12['ROSAT_x'].fillna(psz12['ROSAT_y'], inplace=True)

# merge more columns
psz12['PSZ1 Indx'].fillna(psz12['PSZ1 Index'], inplace=True)
psz12['Comments'].fillna(psz12['Notes'], inplace=True)
psz12['New S/N'].fillna(psz12['SNR_x'], inplace=True)
psz12['Old S/N'].fillna(psz12['SNR_y'], inplace=True)

# imaging columns
psz12['MOSAIC Imaging'].fillna(psz12['Optical run1'], inplace=True)
psz12['Optical run2_x'].fillna(psz12['Optical run2_y'], inplace=True)
psz12['NEWFIRM Imaging'].fillna(psz12['IR run1'], inplace=True)

# drop extra columns
cols = ['Name_y', 'SNR_x', 'SNR_y', 'RA_y', 'Dec', 'RASS rate_y', 'RASS S/N_y',
        'ROSAT_y', 'Notes', 'PSZ1 Index']
for col in cols:
    del psz12[col]

# drop blank rows
psz12 = psz12.drop(psz12.index[psz12['Name_x'].isnull()])

psz12.rename(columns={'New S/N': 'SNR_PSZ2',
                      'Old S/N': 'SNR_PSZ1',
                      'RASS rate_x': 'RASS rate',
                      'RASS S/N_x': 'RASS S/N',
                      'ROSAT_X': 'ROSAT',
                      'Optical run2_x': 'Optical run2',
                      'Name_x': 'Name',
                      'RA_x': 'RA'}, inplace=True)

psz12.to_csv('test.csv')
# done with the two main catalogs -- now we add the swift spreadsheet.
# I don't think I will combine this catalog

# merge psz2 into psz1 based on the psz1 name
psz12s = pd.merge(psz12, swift, left_on='Name', right_on='Target1', how='outer')

# combine a few of the columns
psz12s['RA_x'].fillna(psz12s['RA_y'], inplace=True)
psz12s['DEC'].fillna(psz12s['Dec'], inplace=True)
psz12s['RA_SEX_x'].fillna(psz12s['RA_SEX_y'], inplace=True)
psz12s['DEC_SEX_x'].fillna(psz12s['DEC_SEX_y'], inplace=True)
psz12s['Redshift_x'].fillna(psz12s['Redshift_y'], inplace=True)

# drop extra columns
cols = ['RA_y', 'Dec', 'RA_SEX_y', 'DEC_SEX_y', 'Redshift_y', 'PSZ2 ID#']
for col in cols:
    del psz12s[col]

#psz12s.to_csv('test.csv')

