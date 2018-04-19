from numpy import genfromtxt
import calendar

# read the catalog
d = genfromtxt('./PSZ2_unconfirmed_catalog - current.csv', delimiter=',',
               names=True, dtype=None)
# make a new catalog with selected columns
d2 = d[['Name', 'RA_SEX', 'DEC_SEX', 'SNR_PSZ1', 'SNR_PSZ2', 'MOSAIC_Imaging']]
# mask out the data we want
mask = d2['MOSAIC_Imaging'] != b''

# get the dates
mi = [date.decode() for date in d2['MOSAIC_Imaging'][mask]]

# Change the xx/xx/xx format to Month, Year
mi = ['{}, {}'.format(calendar.month_abbr[int(i.split('/')[0])],
                      i.split('/')[-1]) if i != '' else r'\nd' for i in mi]

d2 = d2[mask]

for i, _ in enumerate(d2):
    print((r'{} & {} & {} & {:.2f} & {:.2f} & {}'
          r'\\'.format(d2['Name'][i].decode().replace('_', ' '),
                                                        d2['RA_SEX'][i].decode(),
                                                        d2['DEC_SEX'][i].decode(),
                                                        d2['SNR_PSZ1'][i],
                                                        d2['SNR_PSZ2'][i],
                                                        mi[i])))
