import matplotlib.pyplot as plt
import os
import numpy
from glob import glob
from astLib import astWCS
from scipy import interpolate
from scipy.stats import linregress

def calc_completeness_dndm():
    ''' Calculates the completeness using a histogram. '''

    from astropy.io import ascii

    bins = numpy.arange(15, 30, 0.5)
    centers = (bins[:-1] + bins[1:]) / 2

    data_dir = '../data/proc2_small/'

    completeness_dndm = []
    for f in fields:
        cat = '{}/{}/{}i_cal.cat'.format(data_dir, f, f)
        cat = ascii.read(cat)
        cat = cat.to_pandas()
        cat = cat.loc[cat.MAG_AUTO < 40]
        cat = cat.loc[cat.CLASS_STAR < 0.8]
        n, bins_ = numpy.histogram(cat['MAG_AUTO'], bins=bins)

        # figure out the completeness
        completeness_dndm.append(centers[numpy.argmax(n) + 1])

        mag = centers[numpy.argmax(n) + 1]

        print(f, mag)

        # make a bunch of figures
        fig, ax = plt.subplots(1)
        n, bins, patches = ax.hist(cat['MAG_AUTO'], bins=bins)
        ax.axvline(mag)
        ax.text(mag + 0.1, ax.get_ylim()[1] * 0.8,
                '{:.2f}'.format(mag))
        ax.set_xlabel('i Magnitude')
        ax.set_ylabel('N')
        ax.set_title(f)

        plt.tight_layout()

        plt.savefig('./completeness_plots/{}_dndm.png'.format(f), bbox='tight')

        plt.close()

    return completeness_dndm


def calc_completeness_hist():
    ''' Calculates the completeness using a histogram. '''

    from astropy.io import ascii

    bins = numpy.arange(15, 30, 0.5)
    centers = (bins[:-1] + bins[1:]) / 2

    data_dir = '../data/proc2_small/'

    completeness_hist = []
    for f in fields:
        cat = '{}/{}/{}i_cal.cat'.format(data_dir, f, f)
        cat = ascii.read(cat)
        cat = cat.to_pandas()
        cat = cat.loc[cat.MAG_AUTO < 40]
        cat = cat.loc[cat.CLASS_STAR < 0.8]

        # make a bunch of figures
        fig, ax = plt.subplots(1)
        n, bins_, patches = ax.hist(cat['MAG_AUTO'], bins=bins)

        # make it a log plot
        logn = numpy.log10(n)

        # find the peak
        peak = numpy.argmax(logn)

        # make a model from mag 18.5 - 21.5
        model = linregress(centers[peak - 5:peak], logn[peak - 5:peak])

        # convert the linear model in lin-log space to log in linear space
        #  and figure out where 80% completeness is
        # see https://en.wikipedia.org/wiki/Semi-log_plot
        y = n / (10**model.intercept * 10**(centers * model.slope))
        x = centers

        # plot(y, x) to see how the ratio curve goes.
        func = interpolate.interp1d(x, y)

        # the interpolate wasn't doing very well...
        # when just asked what is 80%
        mags = numpy.arange(centers[0], centers[-1], 0.1)
        magdiff = 0.8 - func(mags)

        # find the last bin where the difference is negative
        # this is the bin, with the highest magnitude, where we go from having
        # more observed objects to more objects in the model.
        mag_idx = numpy.where(magdiff < 0)[0][-1]

        completeness_hist.append(mags[mag_idx])

        print(f, mags[mag_idx])

        ax.axvline(mags[mag_idx])
        ax.text(mags[mag_idx] + 0.1, ax.get_ylim()[1] * 0.8,
                '{:.2f}'.format(mags[mag_idx]))

        # add the power law
        ax.set_ylim(ax.get_ylim()) # keep the limits normal
        ax.plot(centers,
                     10**model.intercept * 10**(centers * model.slope),
                     marker='o')

        ax.set_xlabel('i Magnitude')
        ax.set_ylabel('N')
        ax.set_title(f)

        plt.tight_layout()

        plt.savefig('./completeness_plots/{}_model.png'.format(f), bbox='tight')

        plt.close()

    return completeness_hist


def calc_completeness():
    ''' This is taken from the plot_mag_z_relation.py file in the sims folder.
    I've just copy and pasted it here for ease of use.

    '''

    Ngal_o = 100
    m1 = 20.0
    m2 = 25.0
    dm = 0.2
    Niter = 10
    filter = 'i'
    path = '../data/sims/Catalogs_Gal_small/'
    #path = '../data/sims/Catalogs_Gal_big/'
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
data_dir = '../data/proc2'
f = plt.figure(figsize=(7, 7 * (numpy.sqrt(5.) - 1.0) / 2.0))
ax = plt.subplot2grid((1, 3), (0, 0), colspan=2)
axs = plt.subplot2grid((1, 3), (0, 2))

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
n, bins, patches = axs.hist(numPerArcmin, bins=20, orientation='horizontal',
                           color='#348abd')
# confirmed
n2, bins, patches = axs.hist(conf, bins=bins, orientation='horizontal',
                            color='#e24a33')

ax.set_ylabel('Object Density (arcmin$^{-2}$)')
ax.set_xlabel('Limiting i Magnitude')
axs.set_xlabel('N$_{fields}$')

