import numpy
import pandas as pd
from scipy import interpolate
import os

def calc_completeness(field, c_lvl=0.8, filter='i', Niter=4, Ngal=100, m1=20,
                      m2=26, dm=0.2):
    ''' Calculates the magnitude for a given completeness, c_lvl, for a
    specific field and filter.

    The other parameters are from the types of simulations done. So those
    should be taken from there.

    '''

    # input values
    mag = numpy.arange(m1, m2, dm)
    frac = numpy.zeros_like(mag)
    Ngal = Ngal * Niter
    for i, m in enumerate(mag):
        cmd = "cat %s/%s_m%.2f_%s_%s.dat | wc" % (path, field,
                                                    mag[i],
                                                    filter, '*')

        # Do simple cat + wc and redirect and split stdout
        Nobs = os.popen(cmd).readline().split()[0]
        frac[i] = float(Nobs) / Ngal
    # plot the 80% completeness lines
    func = interpolate.interp1d(frac, mag)
    try:
        return c_lvl, func(c_lvl)
    except ValueError:
        # return the highest fraction if not c_lvl
        return numpy.max(frac), func(frac[numpy.argmax(frac)])


if __name__ == "__main__":
    path = '/home/boada/Projects/planckClusters/data/sims/Catalogs_Gal'
    spreadsheet = pd.read_csv(
                            '../catalogs/PSZ2_unconfirmed_catalog - proc2.csv')

    new_sheet = pd.DataFrame(spreadsheet['Name'])

    filters = ['g', 'r', 'i', 'z', 'K']

    for filter in filters:
        # initialize the new columns
        new_sheet['Cmag_{}'.format(filter)] = 0.0
        new_sheet['Clvl_{}'.format(filter)] = 0.0
        for i, n in new_sheet['Name'].iteritems():
            if os.path.exists('{}/{}{}.mch'.format(path, n, filter)):
                lvl, mag = calc_completeness(n, filter=filter)
                new_sheet.loc[i, 'Cmag_{}'.format(filter)] = mag
                new_sheet.loc[i, 'Clvl_{}'.format(filter)] = lvl


