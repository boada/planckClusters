#!/usr/bin/env python

from scipy import interpolate
from glob import glob
import os
import cosmology
import numpy
import matplotlib.pyplot as plt
import matplotlib
Polygon = matplotlib.patches.Polygon

def mag_z(axes):

    LBCG = 4.0
    z = numpy.arange(0.01, 1.5, 0.01)
    mstar = mi_star_evol(z)
    mstar_sub = mstar - 2.5 * numpy.log10(0.4)
    BCG = mstar - 2.5 * numpy.log10(LBCG)

    axes.plot(z, mstar, 'k-', linewidth=3, label='$L_\star$ galaxy')
    axes.plot(z, mstar_sub, 'k--', linewidth=1, label='$0.4L_\star$ galaxy')
    axes.plot(z, BCG, 'k-', linewidth=1, label='$%dL_\star$ (BCG)' % LBCG)
    axes.set_xlabel('Redshift')
    axes.set_ylabel('i magnitude')

    axes.legend(loc='lower right', fancybox=True, shadow=True)

    axes.set_xlim(0.05, 1.5)
    axes.set_ylim(16.5, 26)

    return

def mag_lim_hist(axes):
    Ngal_o = 100
    m1 = 20.0
    m2 = 26.0
    dm = 0.2
    Niter = 4
    filter = 'i'
    path = '../data/sims/Catalogs_Gal_small/'
    files = glob('{}/PSZ*{}.mch'.format(path, filter))
    fields = [f.split('/')[-1][:-5] for f in files]

    mag = numpy.arange(m1, m2, dm)

    completeness = []
    completeness_low = []
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
            completeness_low.append(mag[numpy.argmax(frac)])

    axes.hist(completeness, bins=mag)
    axes.hist(completeness_low, bins=mag)
    axes.set_ylabel('$N_{fields}$')
    axes.set_xlabel("$i'$ 80% Limit")

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


if __name__ == "__main__":

    fig, axes = plt.subplots(ncols=2,
                             squeeze=True,
                             figsize=(7, 7 * (numpy.sqrt(5.) - 1.0) / 2.0))
    mag_lim_hist(axes[0])
    mag_z(axes[1])

    plt.tight_layout()
    plt.show()
