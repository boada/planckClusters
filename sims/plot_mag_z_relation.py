#!/usr/bin/env python

from scipy import interpolate
from glob import glob
import os
import cosmology
import numpy
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import linregress
Polygon = matplotlib.patches.Polygon

def mag_z(axes):

    LBCG = 4.0
    z = numpy.arange(0.01, 1.5, 0.01)
    mstar = mi_star_evol(z)
    mstar_sub = mstar - 2.5 * numpy.log10(0.4)
    BCG = mstar - 2.5 * numpy.log10(LBCG)

    axes.plot(z, mstar_sub, 'k-', linewidth=0.5, label='$0.4L_\star$ galaxy')
    axes.plot(z, mstar, 'k-', linewidth=1.5, label='$L_\star$ galaxy')
    axes.plot(z, BCG, 'k-', linewidth=2.5, label='$%dL_\star$ (BCG)' % LBCG)
    axes.set_xlabel('Redshift')

    axes.legend(loc='lower right', fancybox=True, shadow=True)

    axes.set_xlim(0.05, 1.5)
    axes.set_ylim(20, 26)

    return

def calc_completeness_model(fields):
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
        n, bins_ = numpy.histogram(cat['MAG_AUTO'], bins=bins)

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

        print(f, '{:.3f}'.format(mags[mag_idx]))

    return completeness_hist

def mag_lim_hist(axes):
    Ngal_o = 100
    m1 = 20.0
    m2 = 25.0
    dm = 0.2
    Niter = 10
    filter = 'i'
    path = '../data/sims/Catalogs_Gal_small/'

    # find all of the fields we have hunted
    imgs = glob('../cluster_search/round2/PSZ*/**/*A.png', recursive=True)
    fields = [i.split('/')[-2] for i in imgs]

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

    axes.hist(completeness, bins=mag, color='#348abd', orientation='horizontal')
    #axes.hist(completeness_low, bins=mag, color='#348abd')
    #axes.set_ylabel('$N_{fields}$')
    axes.set_xlabel('$N_{fields}$')

    return axes

def mag_lim_hist_model(axes):
    m1 = 20.0
    m2 = 25.0
    dm = 0.2
    mag = numpy.arange(m1, m2, dm)

    # find all of the fields we have hunted
    imgs = glob('../cluster_search/round2/PSZ*/**/*A.png', recursive=True)
    fields = [i.split('/')[-2] for i in imgs]

    completeness = calc_completeness_model(fields)

    axes.hist(completeness, bins=mag, color='#348abd',
              orientation='horizontal', histtype='stepfilled')
    # flip the axis
    axes.invert_xaxis()
    axes.set_ylim(20, 26)
    axes.set_xlabel('$N_{fields}$')
    axes.set_ylabel('Limiting i Magnitude')

    return axes, fields, completeness

# observed mi_star as a function of redshift
def mi_star_evol(z, h=0.7, cosmo=(0.3, 0.7, 0.7)):

    # Blanton's number i.e. M* - 1.5 mags
    BCSPIPE = '/home/boada/Projects/planckClusters/MOSAICpipe'
    evolfile = "1_0gyr_hr_m62_salp.color"
    evolfile = os.path.join(BCSPIPE, "LIB/evol", evolfile)

    k, ev, c = KEfit(evolfile)
    dlum = cosmology.dl(z, cosmology=cosmo)
    # Blanton M*
    Mi_star = -21.22 - 5 * numpy.log10(h)  # + self.evf['i'](z)[0]
    dlum = cosmology.dl(z, cosmology=cosmo)
    DM = 25.0 + 5.0 * numpy.log10(dlum)
    mx = Mi_star + DM + k['i'](z) + ev['i'](z) - ev['i'](0.1)
    return mx

    ##################################################################
    # Read both kcorrection k(z) and evolution ev(z) from BC03 model
    ##################################################################
def KEfit(modelfile):

    import scipy
    import scipy.interpolate
    import tableio

    print("# Getting K(z) and Ev(z) corrections from file:  %s\n" % modelfile)

    e = {}
    k = {}
    c = {}

    (z, c_gr, c_ri, c_iz, k_g, k_r, k_i, k_z, e_g, e_r, e_i,
     e_z) = tableio.get_data(
         modelfile, cols=(0, 3, 4, 5, 10, 11, 12, 13, 14, 15, 16, 17))

    # K-only correction at each age SED,
    k['g'] = scipy.interpolate.interp1d(z, k_g)
    k['r'] = scipy.interpolate.interp1d(z, k_r)
    k['i'] = scipy.interpolate.interp1d(z, k_i)
    k['z'] = scipy.interpolate.interp1d(z, k_z)

    # Evolution term alone
    e['g'] = scipy.interpolate.interp1d(z, e_g)
    e['r'] = scipy.interpolate.interp1d(z, e_r)
    e['i'] = scipy.interpolate.interp1d(z, e_i)
    e['z'] = scipy.interpolate.interp1d(z, e_z)

    # Color redshift
    c['gr'] = scipy.interpolate.interp1d(z, c_gr)
    c['ri'] = scipy.interpolate.interp1d(z, c_ri)
    c['iz'] = scipy.interpolate.interp1d(z, c_iz)

    return k, e, c

def add_z_cl(ax, fields, completeness):
    # fix the path to get the results
    import sys
    sys.path.append('../results/')
    from get_results import loadClusters

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
    depth = [completeness[fields.index(hc)] for hc in numpy.sort(high_conf)]

    # the confirmed = True gets the 15 confirmed clusters
    results = loadClusters(round=3, confirmed=True)

    # sort the results
    results.sort_values('Cluster', inplace=True)

    ax.scatter(results['z_cl_boada'], depth, s=150, marker='*', color='#e24a33')

    return ax

if __name__ == "__main__":

    f = plt.figure(figsize=(7, 7 * (numpy.sqrt(5.) - 1.0) / 2.0))
    ax = plt.subplot2grid((1, 4), (0, 0), colspan=2)
    axs = plt.subplot2grid((1, 4), (0, 2), colspan=2)

    ax, fields, completeness = mag_lim_hist_model(ax)
    mag_z(axs)

    add_z_cl(axs, fields, completeness)

    plt.tight_layout()
    plt.show()
