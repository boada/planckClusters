#!/usr/bin/env python3

import numpy
import os
import pylab
from glob import glob
from matplotlib import rc
rc('axes', labelsize=20)
rc('axes', titlesize=20)
rc('xtick', labelsize=18)
rc('ytick', labelsize=18)

# NTT recovered values
nfields = 1
Ngal = 100
m1 = 20.0
m2 = 25.0
dm = 0.5
rh = 3  # kpc
Niter = 4
filter = 'g'
NTTpath = '/home/boada/Projects/planckClusters/data/sims/Catalogs'
files = glob('{}/PSZ*{}.mch'.format(NTTpath, filter))
fields = [f.split('/')[-1][:-5] for f in files]
print(fields)

for f in fields:
    mag_NTT = numpy.arange(m1, m2, dm)
    frac_NTT = mag_NTT * 0.0
    Ngal_NTT = Ngal * nfields * Niter

    for i in range(len(mag_NTT)):
        # Do simple cat + wc and redirect and split stdout
        cmd = "cat %s/%s_m%.2f_%s_%s.dat | wc" % (NTTpath, f,
                                                      mag_NTT[i],
                                                      filter, '*')
        Nobs = os.popen(cmd).readline().split()[0]
        frac_NTT[i] = float(Nobs) / Ngal_NTT
        print(f, mag_NTT[i], frac_NTT[i])

    pylab.plot(mag_NTT,
            frac_NTT,
            linestyle='dashed',
            color='0.5',
            marker='o',
            ms=5,
            mec='k',
            mfc='0.5')

xmin = m1
xmax = m2
# 80%,50% lines
xo = numpy.array([m1, m2])
y1 = numpy.array([1.0, 1.0])
y2 = numpy.array([0.8, 0.8])
y3 = numpy.array([0.5, 0.5])

pylab.plot(xo, y1, 'k:')
pylab.plot(xo, y2, 'k:')
pylab.plot(xo, y3, 'k:')
pylab.text(xmax - 0.1, 0.5 + 0.02, "50%", ha='right')
pylab.text(xmax - 0.1, 0.8 + 0.02, "80%", ha='right')

pylab.xlim(xmin, xmax)
pylab.ylim(0.001, 1.01)

pylab.ylabel('Recovery Fraction')
#pylab.xlabel('$r$-band apparent magnitude')
pylab.xlabel('i magnitude')

#pylab.savefig('recovery_i.pdf')
#pylab.savefig('recovery_i.eps')
pylab.show()
