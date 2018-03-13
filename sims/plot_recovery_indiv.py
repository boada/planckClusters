#!/usr/bin/env python3

#from astropy.io import ascii
#from astropy.table import vstack
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from glob import glob
from matplotlib import rc
import os

rc('axes', labelsize=20)
rc('axes', titlesize=20)
rc('xtick', labelsize=18)
rc('ytick', labelsize=18)

# input values
nfields = 1
Ngal_o = 100
m1 = 20.0
m2 = 26.0
dm = 0.2
Niter = 4
filter = 'i'
path = '/home/boada/Projects/planckClusters/data/sims/Catalogs_Gal'
files = glob('{}/PSZ*{}.mch'.format(path, filter))
fields = [f.split('/')[-1][:-5] for f in files]

c = ['#348abd', '#7a68a6', '#467821', '#cf4451', '#a60628']

for f in fields:
    mag = numpy.arange(m1, m2, dm)
    for ii, filter in enumerate(['g', 'r', 'i', 'z', 'K']):
        frac = numpy.zeros_like(mag)
        d_mag = numpy.zeros_like(mag)
        Ngal = Ngal_o * Niter
        for i, m in enumerate(mag):
            cmd = "cat %s/%s_m%.2f_%s_%s.dat | wc" % (path, f,
                                                        mag[i],
                                                        filter, '*')
#            fnames = [ascii.read('%s/%s_m%.2f_%s_%s.dat' % (path, f,  m, filter,
#                                                str(n).zfill(3))) for n in
#                    range(1, Niter + 1)]

            # stacked data table
#            d = vstack(fnames)

            # Do simple cat + wc and redirect and split stdout
            Nobs = os.popen(cmd).readline().split()[0]

#           d_mag[i] = numpy.std(d['col6']) # the recovery mag column
            frac[i] = float(Nobs) / Ngal

        plt.errorbar(mag,
                frac,
                xerr=d_mag * 0.0,
                #linestyle='dashed',
                color=c[ii],
                marker='o',
                ms=5,
                mec='k',
                mfc='0.5',
                label=filter)

        xmin = m1
        xmax = m2

        # make a legend
        plt.legend()

        # 80%,50% lines
        xo = numpy.array([m1, m2])
        y1 = numpy.array([1.0, 1.0])
        y2 = numpy.array([0.8, 0.8])
        y3 = numpy.array([0.5, 0.5])

        plt.plot(xo, y1, 'k:')
        plt.plot(xo, y2, 'k:')
        plt.plot(xo, y3, 'k:')
        plt.text(xmax - 0.1, 0.5 + 0.02, "50%", ha='right')
        plt.text(xmax - 0.1, 0.8 + 0.02, "80%", ha='right')

        plt.xlim(xmin, xmax)
        plt.ylim(0.001, 1.01)

        plt.ylabel('Recovery Fraction')
        #pylab.xlabel('$r$-band apparent magnitude')
        plt.xlabel('magnitude')

        plt.tight_layout()
        #pylab.savefig('recovery_i.pdf')
        plt.savefig('plots_gal/{}_recovery.png'.format(f), bbox='tight')
    plt.close()
