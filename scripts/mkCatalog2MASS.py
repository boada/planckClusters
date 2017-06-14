from astroquery.vizier import Vizier
import astropy.coordinates as coord
import astropy.units as u
from astropy.table import Column
from astropy.io import ascii
import numpy
from time import sleep

def work(vquery, ra, dec, target):

    field = coord.SkyCoord(ra=ra,
                        dec=dec,
                        unit=(u.deg, u.deg),
                        frame='icrs')

    data = vquery.query_region(field,
                            radius=("%fd" % 0.5),
                            catalog="II/246/out")[0]

    ### rename column names using PP conventions
    data.remove_column('_RAJ2000')
    data.remove_column('_DEJ2000')
    data.rename_column('_2MASS', '2MASS')
    data.rename_column('RAJ2000', 'ra')
    data.rename_column('DEJ2000', 'dec')

    ### determine RA and Dec positional uncertainties and
    #   add respective columns
    data['errPA'][data['errPA'] == 0] = 1  # workaround

    arc_xopt = numpy.arctan(-data['errMin'] / data['errMaj'] *
                            numpy.tan(data['errPA'].to(u.rad)))
    ra_err = abs(data['errMaj'] * numpy.cos(arc_xopt) *
                numpy.cos(data['errPA'].to(u.rad)) - data['errMin'] *
                numpy.sin(arc_xopt) * numpy.sin(data['errPA'].to(u.rad)))

    data.add_column(Column(data=ra_err, name='e_ra.deg', unit=u.arcsec),
                    index=2)

    arc_yopt = numpy.arctan(data['errMin'] / data['errMaj'] *
                            numpy.cos(data['errPA'].to(u.rad)) /
                            numpy.sin(data['errPA'].to(u.rad)))
    dec_err = abs(data['errMaj'] * numpy.cos(arc_yopt) *
                numpy.sin(data['errPA'].to(u.rad)) + data['errMin'] *
                numpy.sin(arc_yopt) * numpy.cos(data['errPA'].to(u.rad)))
    data.add_column(Column(data=dec_err, name='e_dec.deg',
                        unit=u.arcsec), index=3)

    # remove error ellipse columns
    data.remove_column('errMaj')
    data.remove_column('errMin')
    data.remove_column('errPA')

    return data


max_mag = 21
max_sources = 1e5

vquery = Vizier(columns=['2MASS', 'RAJ2000', 'DEJ2000', 'errMaj', 'errMin',
                         'errPA', 'Jmag', 'e_Jmag', 'Hmag', 'e_Hmag', 'Kmag',
                         'e_Kmag', 'Qflg'],
                column_filters={"Jmag": ("<%f" % max_mag)},
                row_limit=max_sources)

data = numpy.genfromtxt('../catalogs/PSZ2_unconfirmed_catalog - Master.csv',
                     delimiter=',', names=True, dtype=None)

for i, (ra, dec, name) in enumerate(zip(data['RA'], data['DEC'],
                                        data['Name'])):
    print(data['Name'][i])

    result = work(vquery, ra, dec, '{}'.format(name))
    ascii.write(result, format='csv', output='./2MASS/%s_2MASS_catalog.csv' %
                (name.decode()))
    sleep(1)

