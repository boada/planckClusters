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
from matplotlib import rc
#rc('text', usetex=True)
#rc('font',**{'family':'sans-serif','sans-serif':['VeraSe']})
rc('axes', labelsize=20)
rc('axes', titlesize=20)
rc('xtick', labelsize=18)
rc('ytick', labelsize=18)


def main():

    LBCG = 4.0
    z = numpy.arange(0.01, 1.5, 0.01)
    #mstar     = mi_star(z)
    mstar = mi_star_evol(z)
    mstar_sub = mstar - 2.5 * math.log10(0.4)
    BCG = mstar - 2.5 * math.log10(LBCG)

    # 80% recovery - iband
    mlim_SOAR = 24.30
    mlim_NTT = 23.50

    x = numpy.array([0, 1.5])
    y = numpy.array([mlim_SOAR, mlim_SOAR])

    ax = pylab.figure(1, figsize=(8, 6))
    pylab.plot(z, mstar, 'k-', linewidth=3)
    pylab.plot(z, mstar_sub, 'k--', linewidth=1)
    pylab.plot(z, BCG, 'k-', linewidth=1)
    ax = pylab.axes()
    pylab.xlabel('Redshift')
    pylab.ylabel('i magnitude')

    # For the evolution of colors
    BCSPIPE = os.getenv('BCSPIPE')
    evolfile = "1_0gyr_hr_m62_salp.color"
    evolfile = os.path.join(BCSPIPE, "LIB/evol", evolfile)
    k, ev, colorZ = KEfit(evolfile)
    SED = 'El_Benitez2003'
    #SED = 'El_cww'

    # Add the Blakeless points now
    z_850 = numpy.array([21.1, 21.2, 21.4])
    i_pts = z_850 + cosmology.color_z(SED, 'i_SDSS', 'z_WFC', 1.24)
    #i_pts = z_850 + colorZ['iz'](1.24)
    z_pts = i_pts * 0.0 + 1.24
    pylab.plot(z_pts, i_pts, 'o', mfc='0.7',
               mec='k', ms=6)  #mec='k',mfc='0.7',ms=5)

    # Add the Mullis points now
    z_VLT = numpy.array([20.3, 21.0, 21.2])
    i_pts = z_VLT + cosmology.color_z(SED, 'i_SDSS', 'z_SDSS', 1.40)
    #i_pts = z_VLT + colorZ['iz'](1.40)
    z_pts = i_pts * 0.0 + 1.40
    pylab.plot(z_pts, i_pts, 's', mec='k', mfc='0.9', ms=5)

    # Add the Standford clusters
    z_850 = numpy.array([22.80, 22.84, 22.78])
    i_pts = z_850 + cosmology.color_z(SED, 'i_SDSS', 'z_WFC', 1.46)
    #i_pts = z_850 + colorZ['iz'](1.46)
    z_pts = i_pts * 0.0 + 1.46
    pylab.plot(z_pts, i_pts, 'x', mec='k', mfc='0.9', ms=5)

    pylab.legend(('$L^*$ galaxy', '$0.4L^*$ galaxy', '$%dL^*$ (BCG)' % LBCG),
                 'lower right',
                 fancybox=True,
                 shadow=True)

    pylab.text(1.24, 21.7, "CL-1252.9-2927", size=8, ha='center')
    pylab.text(1.36, 20.7, "XMMU-2235.3-2557", size=8, ha='center')
    pylab.text(1.36, 23.95, "XMM-2215.9-1738", size=8, ha='center')

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
    pylab.text(0.1, mlim_NTT + 0.1, "NTT/EFOSC")

    x = numpy.array([0, 1.5])
    y = numpy.array([mlim_SOAR, mlim_SOAR])
    pylab.plot(x, y, 'k:')
    Y1 = mlim_SOAR
    Y2 = mlim_SOAR + 2.0
    verts = [(X1, Y1), (X1, Y2), (X2, Y2), (X2, Y1)]
    poly = Polygon(verts, facecolor='0.70', edgecolor='0.85')
    ax.add_patch(poly)
    pylab.text(0.1, mlim_SOAR + 0.1, "SOAR/SOI")

    pylab.xlim(0.05, 1.5)
    pylab.ylim(16.5, 26)
    pylab.savefig('mi_z.pdf')
    pylab.savefig('mi_z.eps')
    pylab.show()

    return


# observed mi_star as a function of redshift
def mr_star(z, h=0.7, cosmo=(0.3, 0.7, 0.7)):

    dlum = cosmology.dl(z, cosmology=cosmo)

    # Red galaxies fit from Brown et al
    #Mb_star = -19.43 - 1.01*z
    #Mr_star = cosmology.reobs('El_Benitez2003',m=Mb_star, oldfilter="B_Johnson", newfilter="r_SDSS")

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


# observed mi_star as a function of redshift
def mi_star(z, h=0.7, cosmo=(0.3, 0.7, 0.7), SED='El_Benitez2003'):
    dlum = cosmology.dl(z, cosmology=cosmo)

    # Red galaxies fit from Brown et al
    #Mb_star = -19.43 - 1.01*z
    #Mr_star = cosmology.reobs('El_Benitez2003',m=Mb_star, oldfilter="B_Johnson", newfilter="r_SDSS")

    # Alternative -- Paolillo et al. (2001) LF for clusters
    #Mr_star = -21.53 + self.evf['r'](z)[0]
    #Mi_star = cosmology.reobs('El_Benitez2003',m=Mr_star, oldfilter="R_Cousins", newfilter="i_MOSAICII")

    # Blanton M*
    #Mi_star = -21.22 - 5*math.log10(h) #+ self.evf['i'](z)[0]
    Mi_star = -21.22 - 5 * math.log10(h)  # cosmology.kcor(z,SED,'i_SDSS')
    #Mr_star = cosmology.reobs('El_Benitez2003',m=Mi_star, oldfilter="i_SDSS", newfilter="r_SDSS")
    return Mi_star + 5.0 * numpy.log10(dlum) + 25


def mi_star_evol(z, h=0.7, cosmo=(0.3, 0.7, 0.7)):

    # Blanton's number i.e. M* - 1.5 mags
    BCSPIPE = os.getenv('BCSPIPE')
    evolfile = "1_0gyr_hr_m62_salp.color"
    evolfile = os.path.join(BCSPIPE, "LIB/evol", evolfile)

    k, ev, c = KEfit(evolfile)
    dlum = cosmology.dl(z, cosmology=cosmo)
    # Blanton M*
    Mi_star = -21.22 - 5 * math.log10(h)  #+ self.evf['i'](z)[0]
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

    print "# Getting K(z) and Ev(z) corrections from file:  %s\n" % modelfile

    e = {}
    k = {}
    c = {}

    (z, c_gr, c_ri, c_iz, k_g, k_r, k_i, k_z, e_g, e_r, e_i,
     e_z) = tableio.get_data(modelfile,
                             cols=(0, 3, 4, 5, 10, 11, 12, 13, 14, 15, 16, 17))

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


main()
