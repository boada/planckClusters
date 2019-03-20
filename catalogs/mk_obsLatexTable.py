from load_catalogs import load_PSZcatalog
import calendar
import os
from astropy.io.fits import getheader
from astropy.coordinates import SkyCoord
import astropy.units as u

data = load_PSZcatalog()

# sort the names
data = data.sort_values('NAME')

datapath = '../data/proc2'

for i, row in data.iterrows():
    n = row.NAME.replace(' ', '_')
    #print(n)
    if os.path.isdir(f'{datapath}/{n}'):
        name_us = n
    else:
        try:
            n_psz1 = row.NAME_PSZ1.replace(' ', '_')
        except AttributeError:
            continue
        if os.path.isdir(f'{datapath}/{n_psz1}'):
            name_us = n_psz1
        else:
            continue

    try:
        hdr = getheader(f'{datapath}/{name_us}/{name_us}K.fits')
    except FileNotFoundError:
        continue

    date = hdr['DATE-OBS']
    date = date[:9]

    # Change the xx-xx-xx format to Month, Year

    date = (f"{calendar.month_abbr[int(date.split('-')[1])]},"
            f" {date.split('-')[0]}")

    coord = SkyCoord(ra=row.RA, dec=row.DEC, unit=u.deg)
    ra = coord.ra.to_string(sep=':', unit=u.hour, precision=2)
    dec = coord.dec.to_string(sep=':', precision=2)

    print((f'{row.NAME} & '
           f'{ra} & '
           f'{dec} & '
           f'{row.INDEX_PSZ2} ({row.INDEX_PSZ1}) & '
           f'{row.SNR_PSZ2:.3} ({row.SNR_PSZ1:.3}) & '
           f'{date} '
            r'\\'))
