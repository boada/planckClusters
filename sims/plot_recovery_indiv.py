#!/usr/bin/env python3

from scipy import interpolate
from astropy.io import ascii
from astropy.table import vstack
import numpy
from glob import glob
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

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

    fig, axPlot = plt.subplots(figsize=(5.5, 5.5))

    # now we add the top "histogram"
    divider = make_axes_locatable(axPlot)
    axHistx = divider.append_axes("top", 1.2, pad=0.1, sharex=axPlot)

    # set the top xlabel to invisible
    plt.setp(axHistx.get_xticklabels(), visible=False)

    for ii, filter in enumerate(['g', 'r', 'i', 'z', 'K']):
        frac = numpy.zeros_like(mag)
        d_mag = numpy.zeros_like(mag)
        Ngal = Ngal_o * Niter
        for i, m in enumerate(mag):
            cmd = "cat %s/%s_m%.2f_%s_%s.dat | wc" % (path, f,
                                                        mag[i],
                                                        filter, '*')

            # Do simple cat + wc and redirect and split stdout
            Nobs = os.popen(cmd).readline().split()[0]

            frac[i] = float(Nobs) / Ngal

            # now we compute the delta mag
            fnames = []
            for n in range(1, Niter + 1):
                try:
                    t = ascii.read('%s/%s_m%.2f_%s_%s.dat' % (path, f,
                                                                m, filter,
                                                            str(n).zfill(3)))
                    fnames.append(t)
                except (ValueError, FileNotFoundError):
                    continue

            # stacked data table
            try:
                d = vstack(fnames)
                #d_mag[i] = numpy.std(d['col6']) # the recovery mag column
                d_mag[i] = m - numpy.mean(d['col6']) # the recovery mag column
            except TypeError:
                d_mag[i] = 0.0

        # plot the crap
        axPlot.plot(mag, frac, color=c[ii], marker='o', ms=5, label=filter)
        axHistx.scatter(mag, d_mag, color=c[ii], marker='o', s=5, label=filter)

        # plot the 80% completeness lines
        func = interpolate.interp1d(frac, mag)
        axPlot.axvline(func(0.8), lw=1, color=c[ii], zorder=0)

    xmin = m1
    xmax = m2

    # make a legend
    axPlot.legend()

    # 80%,50% lines
    axPlot.axhline(1, ls=':', c='k')
    axPlot.axhline(0.8, ls=':', c='k')
    axPlot.axhline(0.5, ls=':', c='k')

    axPlot.text(xmax - 0.1, 0.5 + 0.02, "50%", ha='right')
    axPlot.text(xmax - 0.1, 0.8 + 0.02, "80%", ha='right')

    axPlot.set_xlim(xmin, xmax)
    axPlot.set_ylim(0.001, 1.01)
    axHistx.set_ylim(-0.7, 2)

    axPlot.set_ylabel('Recovery Fraction')
    #pylab.xlabel('$r$-band apparent magnitude')
    axPlot.set_xlabel('magnitude')
    axHistx.set_ylabel('m - <$m_{rec}$>', fontsize=18)

    plt.tight_layout()
    #pylab.savefig('recovery_i.pdf')
    plt.savefig('plots_gal/{}_recovery.png'.format(f), bbox='tight')
    plt.close()
    #plt.show()
