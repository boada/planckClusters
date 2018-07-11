import matplotlib.pyplot as plt
import os
import numpy
from glob import glob
from astLib import astWCS
from scipy import interpolate


def calc_completeness():
    ''' This is taken from the plot_mag_z_relation.py file in the sims folder.
    I've just copy and pasted it here for ease of use.

    '''

    Ngal_o = 100
    m1 = 20.0
    m2 = 26.0
    dm = 0.2
    Niter = 4
    filter = 'i'
    path = '../data/sims/Catalogs_Gal_small/'
    mag = numpy.arange(m1, m2, dm)

    completeness = []
    for f in fields:
        frac = numpy.zeros_like(mag)
        Ngal = Ngal_o * Niter
        for i, m in enumerate(mag):
            cmd = "cat %s/%s_m%.2f_%s_%s.dat | wc" % (path, f,
                                                        mag[i],
                                                        filter, '*')

            # Do simple cat + wc and redirect and split stdout
            Nobs = os.popen(cmd).readline().split()[0]

            frac[i] = float(Nobs) / Ngal
        # figure out the 80% completeness
        func = interpolate.interp1d(frac, mag)
        try:
            completeness.append(func(0.8))
        except ValueError:
            completeness.append(mag[numpy.argmax(frac)])

    return completeness


# find all of the fields we have hunted
imgs = glob('../cluster_search/round2/PSZ*/**/*A.png', recursive=True)
fields = [i.split('/')[-2] for i in imgs]

# make some stuff before we get going
fig, ax = plt.subplots(1)
data_dir = '../data/proc2'

# where the magic happens
numPerArcmin = []
for f in fields:
    print(f)
    # get the size of the catalog
    cat = '{}/{}/{}i_cal.cat'.format(data_dir, f, f)
    cmd = 'cat {} | wc -l'.format(cat)
    Nobs = os.popen(cmd).readline().split()[0]
    Nobs = int(Nobs) - 47 # subtract off the sextractor header

    # get the size of the image
    wcs = astWCS.WCS('{}/{}/{}Detec.fits'.format(data_dir, f, f))
    size = wcs.getFullSizeSkyDeg()

    numPerArcmin.append(Nobs / (size[0] * size[1] * 60**2))

# confirmed clusters
high_conf = ['PSZ1_G206.45+13.89',
            'PSZ1_G224.82+13.62',
            'PSZ2_G029.66-47.63',
            'PSZ2_G043.44-41.27',
            'PSZ2_G096.43-20.89',
            'PSZ2_G120.76+44.14',
            'PSZ2_G125.55+32.72',
            'PSZ2_G137.24+53.93',
            'PSZ2_G305.76+44.79',
            'PSZ2_G107.83-45.45',
            'PSZ2_G098.38+77.22',
            'PSZ1_G084.62-15.86',
            'PSZ2_G106.11+24.11',
            'PSZ2_G173.76+22.92',
            'PSZ2_G191.82-26.64']

# get the density for the confirmed fields.
conf = [numPerArcmin[fields.index(hc)] for hc in high_conf]

# all
n, bins, patches = ax.hist(numPerArcmin, bins=20)
# confirmed
n2, bins, patches = ax.hist(conf, bins=bins)

ax.set_xlabel('Object Density (arcmin$^{-2}$)')
ax.set_ylabel('N$_{fields}$')

