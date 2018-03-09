#!/usr/bin/env python

import numpy
import os, sys
import pylab

# NTT recovered values
nfields = 10
Ngal = 20
m1 = 22.0
m2 = 26.0
dm = 0.1
rh = 0.5  # kpc
Niter = 15
NTTpath = os.path.join(os.environ['HOME'], "NTT-data/Sim/Catalogs_r")
mag_NTT = numpy.arange(m1, m2, dm)
frac_NTT = mag_NTT * 0.0
Ngal_NTT = Ngal * nfields * Niter

recov = open('recovery_fraction_NTT_r.dat', 'w')
for i in range(len(mag_NTT)):
    # Do simple cat + wc and redirect and split stdout
    cmd = "cat %s/iter_m%.2f_rh%.1f_%s.dat | wc" % (NTTpath, mag_NTT[i], rh,
                                                    '*')
    Nobs = os.popen(cmd).readline().split()[0]
    frac_NTT[i] = float(Nobs) / Ngal_NTT
    recov.write("%8.3f %8.4f\n" % (mag_NTT[i], frac_NTT[i]))
recov.close()

# SOAR values
nfields = 10
Ngal = 25
m1 = 22.0
m2 = 28.0
dm = 0.1
rh = 0.5  # kpc
Niter = 10
SOARpath = os.path.join(os.environ['HOME'], "SOAR-data/Sim/Catalogs_r")
mag_SOAR = numpy.arange(m1, m2, dm)
frac_SOAR = mag_SOAR * 0.0
Ngal_SOAR = Ngal * nfields * Niter

recov = open('recovery_fraction_SOAR_r.dat', 'w')
for i in range(len(mag_SOAR)):
    # Do simple cat + wc and redirect and split stdout
    cmd = "cat %s/iter_m%.2f_rh%.1f_%s.dat | wc" % (SOARpath, mag_SOAR[i], rh,
                                                    '*')
    Nobs = os.popen(cmd).readline().split()[0]
    frac_SOAR[i] = float(Nobs) / Ngal_SOAR
    recov.write("%8.3f %8.4f\n" % (mag_SOAR[i], frac_SOAR[i]))
recov.close()

#print frac_NTT.max()
#print frac_SOAR.max()

frac_SOAR = frac_SOAR * 1.03
pylab.figure(1)
pylab.plot(mag_NTT,
           frac_NTT,
           linestyle='dashed',
           color='0.5',
           marker='o',
           ms=5,
           mec='k',
           mfc='0.5')
pylab.plot(mag_SOAR,
           frac_SOAR,
           linestyle='solid',
           color='0.0',
           marker='o',
           ms=5,
           mec='k',
           mfc='0.0')
pylab.legend(('NTT/EFOSC', 'SOI/SOAR'),
             'lower left',
             fancybox=True,
             shadow=True,
             borderpad=0.8)

xmin = 22
xmax = 26
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
pylab.xlabel('r magnitude')

pylab.savefig('recovery_r.pdf')
pylab.savefig('recovery_r.eps')
pylab.show()
