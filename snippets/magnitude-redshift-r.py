#!/usr/bin/env python

import os, sys
import cosmology
import cosmopy
import numpy
import math
import pylab
import matplotlib
#from matplotlib import rc
Polygon = matplotlib.patches.Polygon


def main():

    #rc('axes',labelsize=14)
    #rc('axes',titlesize=14)
    #rc('xtick',labelsize=12)
    #rc('ytick',labelsize=12)

    z = numpy.arange(0.01, 1.5, 0.01)
    mstar = mr_star(z)
    mstar_sub = mstar - 2.5 * math.log10(0.4)
    BCG = mstar - 1.5

    # 5-sigma values
    mlim_SOAR = 24.39
    mlim_NTT = 23.60

    # 80% recovery
    mlim_SOAR = 24.60
    mlim_NTT = 23.60

    x = numpy.array([0, 1.5])
    y = numpy.array([mlim_SOAR, mlim_SOAR])

    ax = pylab.figure(1, figsize=(8, 6))
    pylab.plot(z, mstar, 'k-', linewidth=3)
    pylab.plot(z, mstar_sub, 'k--', linewidth=1)
    pylab.plot(z, BCG, 'k-', linewidth=1)
    ax = pylab.axes()
    pylab.xlabel('Redshift')
    pylab.ylabel('r magnitude')

    # Add the Blakeless points now
    z_850 = numpy.array([21.1, 21.2, 21.4])
    r_pts = z_850 + 2.0  # r_SDSS-z_850
    z_pts = r_pts * 0.0 + 1.24
    pylab.plot(z_pts, r_pts, 'o', mfc='0.7',
               mec='k', ms=6)  #mec='k',mfc='0.7',ms=5)

    # Add the Mullis points now
    z_VLT = numpy.array([20.3, 21.0, 21.2])
    r_pts = z_VLT + 2.2  # r-z
    z_pts = r_pts * 0.0 + 1.40
    pylab.plot(z_pts, r_pts, 's', mec='k', mfc='0.9', ms=5)

    # Add the Standford clusters
    z_850 = numpy.array([22.80, 22.84, 22.78])
    r_pts = z_850 + 2.4
    z_pts = r_pts * 0.0 + 1.46

    pylab.plot(z_pts, r_pts, 'x', mec='k', mfc='0.9', ms=5)

    #pylab.legend( ('$L^*$ galaxy', '$0.4L^*$ galaxy', '$4L^*$ (BCG)','Mag lim $r=24.5$'),'lower right' )
    pylab.legend(('$L^*$ galaxy', '$0.4L^*$ galaxy', '$4L^*$ (BCG)'),
                 'lower right',
                 fancybox=True,
                 shadow=True)

    pylab.text(1.24, 22.5, "CL-1252.9-2927", size=8, ha='center')
    pylab.text(1.36, 21.7, "XMMU-2235.3-2557", size=8, ha='center')

    # Make detection region
    x = numpy.array([0, 1.5])
    y = numpy.array([mlim_NTT, mlim_NTT])
    pylab.plot(x, y, 'k:')
    X1 = 0
    X2 = 1.5
    Y1 = mlim_NTT
    Y2 = mlim_NTT + 2.0
    verts = [(X1, Y1), (X1, Y2), (X2, Y2), (X2, Y1)]
    poly = Polygon(verts, facecolor='0.95', edgecolor='0.85')
    ax.add_patch(poly)
    pylab.text(0.05, mlim_NTT + 0.1, "NTT/EFOSC")

    x = numpy.array([0, 1.5])
    y = numpy.array([mlim_SOAR, mlim_SOAR])
    pylab.plot(x, y, 'k:')
    Y1 = mlim_SOAR
    Y2 = mlim_SOAR + 2.0
    verts = [(X1, Y1), (X1, Y2), (X2, Y2), (X2, Y1)]
    poly = Polygon(verts, facecolor='0.70', edgecolor='0.85')
    ax.add_patch(poly)
    pylab.text(0.05, mlim_SOAR + 0.1, "SOAR/SOI")

    pylab.xlim(0, 1.5)
    pylab.ylim(14.5, 26)
    pylab.savefig('mr_z.pdf')
    pylab.show()

    return


# observed mi_star as a function of redshift
def mi_star(z, h=0.7, cosmo=(0.3, 0.7, 0.7)):
    #dlum = self.c.dlum(z)[0]
    dlum = cosmology.dl(z, cosmology=cosmo)

    # Red galaxies fit from Brown et al
    Mb_star = -19.43 - 1.01 * z
    Mi_star = cosmology.reobs('El_Benitez2003',
                              m=Mb_star,
                              oldfilter="B_Johnson",
                              newfilter="i_MOSAICII")

    # Alternative -- Paolillo et al. (2001) LF for clusters
    #Mr_star = -21.53 + self.evf['r'](z)[0]
    #Mi_star = cosmology.reobs('El_Benitez2003',m=Mr_star, oldfilter="R_Cousins", newfilter="i_MOSAICII")

    # Blanton M*
    #Mi_star = -21.22 - 5*math.log10(h) + self.evf['i'](z)[0]

    return Mi_star + 5.0 * math.log10(dlum) + 25


# observed mi_star as a function of redshift
def mr_star(z, h=0.7, cosmo=(0.3, 0.7, 0.7)):
    #dlum = self.c.dlum(z)[0]
    dlum = cosmology.dl(z, cosmology=cosmo)

    # Red galaxies fit from Brown et al
    Mb_star = -19.43 - 1.01 * z
    Mr_star = cosmology.reobs('El_Benitez2003',
                              m=Mb_star,
                              oldfilter="B_Johnson",
                              newfilter="r_SDSS")

    # Alternative -- Paolillo et al. (2001) LF for clusters
    #Mr_star = -21.53 + self.evf['r'](z)[0]
    #Mi_star = cosmology.reobs('El_Benitez2003',m=Mr_star, oldfilter="R_Cousins", newfilter="i_MOSAICII")

    # Blanton M*
    Mi_star = -21.22 - 5 * math.log10(h)  #+ self.evf['i'](z)[0]
    Mr_star = cosmology.reobs('El_Benitez2003',
                              m=Mi_star,
                              oldfilter="i_SDSS",
                              newfilter="r_SDSS")
    return Mr_star + 5.0 * numpy.log10(dlum) + 25


def mr_star_new(z):

    # Blanton's number i.e. M* - 1.5 mags
    BCSPIPE = os.getenv('BCSPIPE')
    evolfile = "0_1gyr_hr_m62_salp.color"
    evolfile = os.path.join(BCSPIPE, "LIB/evol", evolfile)

    k, ev = KEfit(evolfile)
    h = 0.7
    flat = (0.3, 0.7, 0.7)
    cset = cosmopy.set(flat)
    #z    = numpy.arange(x1*0.95,x2*1.01,0.01)
    Mev = -20.44 + 5 * math.log10(h)  # + ev['r'](z)
    dlum = cset.dlum(z)
    DM = 25.0 + 5.0 * numpy.log10(dlum)
    mx = Mev + DM + k['r'](z) + ev['r'](z)
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

    (z, k_g, k_r, k_i, k_z, e_g, e_r, e_i,
     e_z) = tableio.get_data(modelfile,
                             cols=(0, 10, 11, 12, 13, 14, 15, 16, 17))

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

    return k, e


main()
