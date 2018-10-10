#!/usr/bin/env python
#

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import zip
from builtins import map
from builtins import range
from builtins import object
from past.utils import old_div
import os
import sys
from . import extras
import numpy.random_array as nrandom
import re, time
import tableio
import cosmopy
import math
import cosmology  # from bpz
from . import astrometry
import fitsio
import bpz_mix
import numpy

# matplotlib stuff
import pylab
from matplotlib import rc
import MLab
import random
import aux

import scipy
import scipy.interpolate
import scipy.special
import scipy.misc.pilutil as pilutil
erf = scipy.special.erf  # normal distribution error function

sout = sys.stderr
#logfile = "logfile_%s.log" % time.strftime("%d%b%Y_%H:%M", time.localtime())
#sout = open(logfile,"w")
#print >>sys.stderr,"# Will write to %s" % logfile
land = numpy.logical_and
lge = numpy.greater_equal
lle = numpy.less_equal
lor = numpy.logical_or


class BCGfinder(object):
    """ A class to find clusters in the BCS Survey"""

    def __init__(self, catsfile,probsfile,
                 maglim=24.0,
                 zlim =1.2,
                 dz=0.05,
                 cosmo=(0.3,0.7,0.7),
                 zuse="ZB", # Use ZB (Bayesian) or ML (Max Like)
                 outpath='plots',
                 path = "/home/felipe/COMB-07",
                 evolfile = "0_1gyr_hr_m62_salp.color"):

        # Check for environ vars
        if not os.getenv('BCSPIPE'):
            os.environ['BCSPIPE'] = '/home/felipe/BCSPIPE'
        self.BCSPIPE = os.getenv('BCSPIPE')

        self.catsfile = catsfile
        self.probsfile = probsfile
        self.maglim = maglim
        self.zlim = zlim
        self.cosmo = cosmo
        self.evolfile = os.path.join(self.BCSPIPE, "LIB/evol", evolfile)
        self.dz = dz
        self.zuse = zuse
        self.outpath = outpath
        self.path = path

        # Set the cosmology now
        self.cset = cosmopy.set(self.cosmo)
        self.Om = cosmo[0]
        self.OL = cosmo[1]
        self.h = cosmo[2]
        self.Ho = self.h * 100.0

        self.read_cat()  # Read catalogs avoding, faint, high-z and 99 objects
        self.read_probs()  # Read probs function of objects from catalogs
        self.get_absmags()  # We compute Abs Mags for each object

        # Not need, done by self.get_absmags()
        #self.get_evol()    # Compute evol(z) for each object
        # Set the BCG masked to False, so we select BCGs on the firts run
        self.BCG_masked = False
        self.BCG_probs = False
        self.pixscale = 0.266

        # Check if the destination folder exists
        if os.path.exists(self.outpath):
            sout.write("# Will put files to: %s\n" % self.outpath)
        else:
            sout.write("# Will create new folder: %s" % self.outpath)
            os.mkdir(self.outpath)

        return

    ##################################################################
    # Define the sub-sample for BCGs candidates around a position
    ##################################################################
    def get_BCG_candidates_radec(self,
                                 ID,
                                 ra,
                                 dec,
                                 zo,
                                 Mr_limit=-22.71,
                                 p_lim=1e-4,
                                 i_lim=24,
                                 plot='yes',
                                 D_factor=1.0,
                                 dz=0.1):

        t0 = time.time()
        RA = astrometry.dec2deg(old_div(ra, 15.))
        DEC = astrometry.dec2deg(dec)

        # Make some variables part of the class
        self.ra0 = ra
        self.dec0 = dec
        self.RA = RA
        self.DEC = DEC
        self.cID = ID  # Candidate's ID
        self.z_c = zo

        # The Abs mag limit @ z=0.1 in the i-band
        Mi_limit = cosmology.reobs('El_Benitez2003',
                                   m=Mr_limit,
                                   oldfilter="r_MOSAICII",
                                   newfilter="i_MOSAICII")

        #dz  = 0.1
        z1 = zo - dz
        z2 = zo + dz
        star_lim = 0.50
        Dmin = D_factor * 250.0  # in kpc

        # Ccompute the largerst angular distance closer to the lower redshift limit of the interval
        if zo <= 0.1:  # to avoid inf
            zcenter = zo
        else:
            zcenter = zo - dz
        dmin = old_div(
            astrometry.kpc2arc(zcenter, Dmin, self.cosmo),
            3600.)  # in degrees.
        print("# Will Use Dmin: %.2f @ zo: %s" % (dmin * 60.0, zcenter))

        # Evaluate the genertic mask for BCG only onece
        if not self.BCG_masked:

            # We also might want to take into account the error in the
            # photo-z (around ~0.05) that can change the value of the
            # distance modulus DM = 25 + 5*log10(dl) significantly,
            # about 1Mag a more at z=0.1 and below. The factor is :
            # 5log10(dL(z+dz)/dL(z)), where z is the objects's photo-z and
            # dz the error. So the limit is 5log10(dL(z+dz)/dL(z)) magnitude
            # fainter.
            # dL_up = self.cset.dlum(self.z_ph+self.dz)
            # dL_lo = self.cset.dlum(self.z_ph)
            # self.DM_factor = 5*numpy.log10(dL_up/dL_lo)
            #
            # We get the limit at the z_ph of each candidate, corrected by z=0.1
            Mr_BCG_limit = Mr_limit + self.ev_r - self.evf['r'](
                0.1)  #+ self.DM_factor
            Mi_BCG_limit = Mi_limit + self.ev_i - self.evf['i'](
                0.1)  #+ self.DM_factor
            # Evaluate the BCG Probability function, we get the limit for each object
            self.p = p_BCG(self.Mr, Mr_BCG_limit)

            sout.write(
                "# Selecting BCG candidates from %s objects... will do this only once\n"
                % (len(self.ra)))
            mask_p = numpy.where(self.p >= p_lim, 1, 0)
            mask_g = numpy.where(self.g < i_lim + 5, 1, 0)
            mask_r = numpy.where(self.r < i_lim + 2, 1, 0)
            mask_i = numpy.where(self.i < i_lim, 1, 0)
            mask_z = numpy.where(self.z < i_lim + 1, 1, 0)

            # Avoid freakishly bright objects, 2.5 mags brighter than the M_BCG_limit
            mask_br = numpy.where(self.Mr > Mr_BCG_limit - 2.5, 1, 0)
            mask_bi = numpy.where(self.Mi > Mi_BCG_limit - 2.5, 1, 0)

            # Put a more strict cut in class_star for bcg candidates
            sout.write("# Avoiding CLASS_STAR > %s in BGCs\n" % star_lim)
            mask_star = numpy.where(self.class_star <= star_lim, 1, 0)

            # Construct the final mask now
            self.mask_BCG = mask_g * mask_r * mask_i * mask_z * mask_br * mask_bi * mask_p
            self.BCG_masked = True

            # Model color only once
            self.zx = numpy.arange(0.01, self.zlim, 0.01)
            self.gr_model = cosmology.color_z(sed='El_Benitez2003',
                                              filter_new='g_MOSAICII',
                                              filter_old='r_MOSAICII',
                                              z=self.zx,
                                              calibration='AB')
            self.ri_model = cosmology.color_z(sed='El_Benitez2003',
                                              filter_new='r_MOSAICII',
                                              filter_old='i_MOSAICII',
                                              z=self.zx,
                                              calibration='AB')
            self.iz_model = cosmology.color_z(sed='El_Benitez2003',
                                              filter_new='i_MOSAICII',
                                              filter_old='z_MOSAICII',
                                              z=self.zx,
                                              calibration='AB')

            sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t0))

        # Now we select based on position and redshift.
        # These are the only two masks that depend on the position and
        # redshift and will be compurted on each call.
        mask_zph = numpy.where(land(self.z_ph >= z1, self.z_ph <= z2), 1, 0)

        # Mask in positions
        #distance = numpy.sqrt( (ra - self.ra)**2 + (dec - self.dec)**2)
        distance = astrometry.circle_distance(ra,
                                              dec,
                                              self.ra,
                                              self.dec,
                                              units='deg')
        mask_pos = numpy.where(distance <= dmin, 1, 0)

        # Select the candidates now
        idx = numpy.where(mask_zph * mask_pos * self.mask_BCG == 1)

        # And pass up to to class
        self.idx_BCG = idx
        self.id_BCG = self.id[idx]
        self.ra_BCG = self.ra[idx]
        self.dec_BCG = self.dec[idx]
        self.p_BCG = self.p[idx]
        self.z_BCG = self.z_ph[idx]
        self.t_BCG = self.type[idx]
        self.N_BCG = len(idx[0])
        self.Mi_BCG = self.Mi[idx]
        self.Mr_BCG = self.Mr[idx]
        self.DM_BCG = self.DM[idx]  # distance modulus
        self.dang_BCG = self.dang[idx]  # distance modulus

        self.zml_BCG = self.z_ml[idx]
        self.tml_BCG = self.t_ml[idx]
        self.zb_BCG = self.z_b[idx]
        self.tb_BCG = self.t_b[idx]
        self.class_BCG = self.class_star[idx]
        self.a_BCG = self.a_image[idx]
        self.b_BCG = self.b_image[idx]
        self.theta_BCG = self.theta[idx]

        # r,i-band stuff
        self.r_BCG = self.r[idx]
        self.i_BCG = self.i[idx]

        # Get the 1-sigma intervals
        self.z1_BCG = self.z1[idx]
        self.z2_BCG = self.z2[idx]

        # Get the tile name closet to that postion
        iclose = numpy.argmin(distance)
        IDclose = self.id[iclose]
        self.tile = IDclose.split('_')[0]

        # The distance to the candidate's position for each BCG, in arcmin
        self.d_BCG = distance[idx] * 60.0

        # The radius used in the search in arcsec
        self.dmin = dmin * 60.0

        # The r-band Luminosity of the BCGs
        self.LBCG = self.Lr[idx]

        sout.write(
            "# Found %s BCG candidates with i<%s at zo:%s, (z1,z2) = %s-%s\n" %
            (self.N_BCG, i_lim, zo, z1, z2))
        # Optional, plot to see that we are getting the right stuff
        if plot:
            self.BCG_plots(ID, zo, Mr_limit)
        return

    ##################################################################
    # Define the sub-sample for BCGs candidates around a position
    ##################################################################
    def get_BCG_ID(self, ID, Mr_limit=-22.71, p_lim=1e-4, SCSname=None):

        t0 = time.time()

        # The Abs mag limit @ z=0.1 in the i-band
        Mi_limit = cosmology.reobs('El_Benitez2003',
                                   m=Mr_limit,
                                   oldfilter="r_MOSAICII",
                                   newfilter="i_MOSAICII")

        # Evaluate the genertic mask for BCG only onece
        if not self.BCG_probs:

            # We get the limit at the z_ph of each candidate, corrected by z=0.1
            Mr_BCG_limit = Mr_limit + self.ev_r - self.evf['r'](
                0.1)  #+ self.DM_factor
            Mi_BCG_limit = Mi_limit + self.ev_i - self.evf['i'](
                0.1)  #+ self.DM_factor
            # Evaluate the BCG Probability function, we get the limit for each object
            self.p = p_BCG(self.Mr, Mr_BCG_limit)
            self.BCG_probs = True
            # Model color only once
            self.zx = numpy.arange(0.01, self.zlim, 0.01)
            self.gr_model = cosmology.color_z(sed='El_Benitez2003',
                                              filter_new='g_MOSAICII',
                                              filter_old='r_MOSAICII',
                                              z=self.zx,
                                              calibration='AB')
            self.ri_model = cosmology.color_z(sed='El_Benitez2003',
                                              filter_new='r_MOSAICII',
                                              filter_old='i_MOSAICII',
                                              z=self.zx,
                                              calibration='AB')
            self.iz_model = cosmology.color_z(sed='El_Benitez2003',
                                              filter_new='i_MOSAICII',
                                              filter_old='z_MOSAICII',
                                              z=self.zx,
                                              calibration='AB')
            sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t0))

        # And pass up to to class
        idx = numpy.where(self.id == ID)
        # The index number
        iBCG = idx  #[0]
        self.idx_BCG = idx  #[0]
        self.id_BCG = self.id[iBCG]
        self.ra_BCG = self.ra[iBCG]
        self.dec_BCG = self.dec[iBCG]
        self.p_BCG = self.p[iBCG]
        self.z_BCG = self.z_ph[iBCG]
        self.t_BCG = self.type[iBCG]
        self.N_BCG = len(idx[0])
        self.Mi_BCG = self.Mi[iBCG]
        self.Mr_BCG = self.Mr[iBCG]
        self.DM_BCG = self.DM[iBCG]  # distance modulus
        self.dang_BCG = self.dang[iBCG]  # distance modulus

        self.zml_BCG = self.z_ml[iBCG]
        self.tml_BCG = self.t_ml[iBCG]
        self.zb_BCG = self.z_b[iBCG]
        self.tb_BCG = self.t_b[iBCG]
        self.class_BCG = self.class_star[iBCG]
        self.a_BCG = self.a_image[iBCG]
        self.b_BCG = self.b_image[iBCG]
        self.theta_BCG = self.theta[iBCG]

        # r,i-band stuff
        self.r_BCG = self.r[iBCG]
        self.i_BCG = self.i[iBCG]

        # Get the 1-sigma intervals
        self.z1_BCG = self.z1[iBCG]
        self.z2_BCG = self.z2[iBCG]

        # Get the tile name closet to that postion
        self.tile = ID.split('_')[0]

        print(old_div(self.ra_BCG, 15))
        print(old_div(self.dec_BCG, 15))

        self.RA = astrometry.dec2deg(old_div(self.ra_BCG, 15.))[0]
        self.DEC = astrometry.dec2deg(self.dec_BCG)[0]

        # The SCS rootname for the files
        if SCSname:
            self.SCSname = SCSname
        else:
            RA = astrometry.dec2deg(old_div(self.ra_BCG[0], 15), sep="")
            DEC = astrometry.dec2deg(self.dec_BCG[0], sep="")
            # self.SCSname = "SCSO_J%s%s" % (RA[0:4],DEC[0:5])
            self.SCSname = "SCSO_J%s%s" % (RA[0:6],
                                           DEC[0:7])  # to avoid confusion

        # The r-band Luminosity of the BCGs
        self.LBCG = self.Lr[idx]

        # The distance to the candidate's position for each BCG, in arcmin
        self.d_BCG = 0.0
        sout.write("# Found BCG candidates for %s @ i:%s \n" % (ID, iBCG[0]))

    #######################
    # BCG plot diagnostics
    #######################
    def BCG_plots(self, ID, zo, Mr_limit):

        params = {
            'axes.labelsize': 12,
            'text.fontsize': 10,
            'legend.fontsize': 10,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            # The figure subplot parameters.  All dimensions are fraction of the
            # figure width or height
            'figure.subplot.left':
            0.10,  # the left side of the subplots of the figure
            'figure.subplot.right':
            0.95,  # the right side of the subplots of the figure
            'figure.subplot.bottom':
            0.05,  # the bottom of the subplots of the figure
            'figure.subplot.top':
            0.95,  # the top of the subplots of the figure
            'figure.subplot.wspace':
            0.2,  # the amount of width reserved for blank space between subplots
            'figure.subplot.hspace':
            0.3,  # the amount of height reserved for white space between subplots
            #'font.size': 8,
            #'backend': 'ps',
            #'text.usetex': True,
            'figure.figsize': (11.5, 9)
        }

        pylab.rcParams.update(params)
        zx = self.zx
        idx = self.idx_BCG
        #######################
        # p_BCG(z) plot at zo
        #######################
        Mr_limit_zo = Mr_limit + self.evf['r'](zo) - self.evf['r'](0.1)
        M1 = Mr_limit_zo - 3
        M2 = Mr_limit_zo + 3
        M_x = numpy.arange(M1, M2, 0.05)
        P_x = p_BCG(M_x, Mr_limit_zo)
        pylab.figure(1)
        pylab.subplot(321)
        pylab.plot(M_x, P_x, 'k-')
        pylab.plot(self.Mr_BCG, self.p_BCG, 'ro')
        xx = numpy.asarray([Mr_limit_zo, Mr_limit_zo])
        yy = numpy.asarray([-3, 3])
        pylab.plot(xx, yy, '--')
        dx = [-2, -1, 1, +2]  # make dotted lines at N +/- mags
        for d in dx:
            pylab.plot(xx + d, yy, 'k:')
        pylab.xlabel("r-band Abs Mag")
        pylab.ylabel("BCG Probability p(M)")
        pylab.ylim(-0.02, 1.02)
        pylab.xlim(M1, M2)
        ################
        # M vs m plot
        ################
        dlum = self.cset.dlum(zx)
        dm = 25.0 + 5.0 * numpy.log10(dlum)
        Mo = Mr_limit + self.evf['r'](zx) - self.evf['r'](0.1)
        mo = Mo + dm + self.kcorr['r'](zx)
        M1 = Mo + 1.0
        M2 = Mo - 1.0
        pylab.subplot(323)
        pylab.plot(mo, Mo, 'k')
        pylab.plot(mo, M1, 'k:')
        pylab.plot(mo, M2, 'k:')
        pylab.plot(self.r_BCG, self.Mr_BCG, 'ro')
        pylab.xlabel("r-band magnitude")
        pylab.ylabel("Abs Magnitude")
        ################
        # m vs z  plot
        ################
        pylab.subplot(325)
        pylab.plot(zx, mo, 'k--')
        pylab.plot(self.z_BCG, self.r_BCG, 'ro')
        pylab.ylabel("r-band magnitude")
        pylab.xlabel("Redshift")
        pylab.xlim(0.05, 0.9)
        pylab.ylim(13.5, 24.1)
        #############################
        # Color - redshift plots
        #############################
        gr = self.g_bpz[idx] - self.r_bpz[idx]
        ri = self.r_bpz[idx] - self.i_bpz[idx]
        iz = self.i_bpz[idx] - self.z_bpz[idx]
        # (g-r)
        pylab.subplot(322)
        pylab.plot(self.z_BCG, gr, 'ro')
        pylab.plot(zx, self.gr_model, 'k--')
        pylab.plot(zx, self.gr_model - 0.3, 'k:')
        pylab.plot(zx, self.gr_model + 0.3, 'k:')
        pylab.ylabel("g-r")
        pylab.xlim(0.05, 0.9)
        # (r-i)
        pylab.subplot(324)
        pylab.plot(self.z_BCG, ri, 'ro')
        pylab.plot(zx, self.ri_model, 'k--')
        pylab.plot(zx, self.ri_model - 0.3, 'k:')
        pylab.plot(zx, self.ri_model + 0.3, 'k:')
        pylab.ylabel("r-i")
        pylab.xlim(0.05, 0.9)
        # (i-z)
        pylab.subplot(326)
        pylab.plot(self.z_BCG, iz, 'ro')
        pylab.plot(zx, self.iz_model, 'k--')
        pylab.plot(zx, self.iz_model + 0.3, 'k:')
        pylab.plot(zx, self.iz_model - 0.3, 'k:')
        pylab.xlabel("Redshift")
        pylab.ylabel("i-z")
        pylab.xlim(0.05, 0.9)
        # Make the ps file
        pylab.savefig(os.path.join(self.outpath, "%s_plots.eps" % ID))
        pylab.close()
        return

    ##################################################################
    # Select neighbors around ra,dec and at the right brightness @ zo
    ##################################################################
    def get_BCG_neighbors(self, i, Mi_lim=-20.25, dz=0.05):

        t0 = time.time()
        sout.write("# Selecting Neighbors ")
        # Get the relevant info for ith BCG
        zo = self.z_BCG[i]
        ra0 = self.ra_BCG[i]
        dec0 = self.dec_BCG[i]
        Mi_BCG = self.Mi_BCG[i]
        DM = self.DM_BCG[i]

        # 1 - Select in position around ra0,dec0
        # Define 1h^-1 Mpc radius in degress @ zo
        R1Mpc = 1000 * 1.0 / self.h  # in kpc
        rmin = old_div(
            astrometry.kpc2arc(zo, R1Mpc, self.cosmo), 3600.)  # in degrees.
        #dist  = numpy.sqrt( (ra0 - self.ra)**2 + (dec0 - self.dec)**2)
        dist = astrometry.circle_distance(ra0,
                                          dec0,
                                          self.ra,
                                          self.dec,
                                          units='deg')

        mask_pos = numpy.where(dist < rmin, 1, 0)

        # 2 - Select in redshift
        #dz = self.dz*(1 + zo)
        z1 = zo - dz
        z2 = zo + dz
        mask_z = numpy.where(land(self.z_ph >= z1, self.z_ph <= z2), 1, 0)

        # 3 - Select in brightness
        Mi_lim_zo = Mi_lim + self.evf['i'](zo) - self.evf['i'](0.1)
        mask_L1 = numpy.where(self.Mi <= Mi_lim_zo, 1, 0)  # Faint  cut > 0.4L*
        mask_L2 = numpy.where(self.Mi >= Mi_BCG, 1, 0)  # Bright cut < L_BCG

        # The final selection mask, position x redshift x Luminosity
        mask_sel = mask_pos * mask_L1 * mask_L2 * mask_z

        idx = numpy.where(mask_sel == 1)
        sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t0))
        return idx

    #########################################
    # Read in the big catalog of photometry
    #########################################
    def read_cat(self):

        cols = (1,
                2,
                23,
                27,
                26,
                28,
                29,
                30,
                3,
                4,  #5,
                6,
                7,  #8,
                9,
                10,  #11,
                12,
                13,  #14,
                15,
                16,
                17,
                18,
                19,
                20,
                21,
                22,
                31,
                32,
                33,
                34)

        t1 = time.time()
        sout.write("# Reading cols:%s\n# Reading cats from: %s... \n" %
                   (cols, self.catsfile))

        (ra,
         dec,
         z_b,
         odds,
         t_b,
         z_ml,
         t_ml,
         chi,
         g,
         g_err,  #g_sn,
         r,
         r_err,  #r_sn,
         i,
         i_err,  #i_sn,
         z,
         z_err,  #z_sn,
         g_bpz,
         g_berr,
         r_bpz,
         r_berr,
         i_bpz,
         i_berr,
         z_bpz,
         z_berr,  #) = tableio.get_data(self.catsfile,cols=cols)
         class_star,
         a_image,
         b_image,
         theta) = tableio.get_data(self.catsfile, cols=cols)

        (id) = tableio.get_str(self.catsfile, cols=(0, ))

        ############################################
        # Choose the photo-z to use, ml or bayesian
        ############################################
        sout.write("# Will use %s redshifts\n" % self.zuse)
        if self.zuse == "ML":
            z_ph = z_ml
            t = t_ml
        elif self.zuse == "ZB":
            z_ph = z_b
            t = t_b

        i_lim = self.maglim
        odds_lim = 0.80
        star_lim = 0.80

        # Clean up according to BPZ
        sout.write("# Avoiding magnitudes -99 and 99 in BPZ \n")
        g_mask = numpy.where(lor(g_bpz == 99, g_bpz == -99), 0, 1)
        r_mask = numpy.where(lor(r_bpz == 99, r_bpz == -99), 0, 1)
        i_mask = numpy.where(lor(i_bpz == 99, i_bpz == -99), 0, 1)
        z_mask = numpy.where(lor(z_bpz == 99, z_bpz == -99), 0, 1)
        bpz_mask = g_mask * r_mask * i_mask * z_mask

        # Clean up to avoid 99 values and very faint i_mag values
        sout.write("# Avoiding magnitudes 99 in MAG_AUTO \n")
        g_mask = numpy.where(g >= 99, 0, 1)
        r_mask = numpy.where(r >= 99, 0, 1)
        i_mask = numpy.where(i >= i_lim, 0, 1)
        z_mask = numpy.where(z >= 99, 0, 1)
        sout.write("# Avoiding magnitudes i > %s in MAG_AUTO \n" % i_lim)

        # Clean by class_star
        sout.write("# Avoiding CLASS_STAR > %s \n" % star_lim)
        mask_star = numpy.where(class_star > star_lim, 0, 1)

        # Clean up by odds
        sout.write("# Avoiding ODDS < %s in BPZ \n" % odds_lim)
        odds_mask = numpy.where(odds > odds_lim, 1, 0)

        # Avoid z> zlim objects too.
        sout.write("# Avoiding objects with z > %s " % self.zlim)
        zp_mask = numpy.where(z_ph > self.zlim, 0, 1)

        # The final 'good' mask
        mask_good = g_mask * r_mask * i_mask * z_mask * zp_mask * odds_mask * mask_star
        idx = numpy.where(mask_good == 1)

        # Make ids a Char String in numarray
        self.id = nstr.array(id)[idx]

        # Only keep the 'good' one, avoid -99 and 99 values in BPZ mags
        self.ra = ra[idx]
        self.dec = dec[idx]
        self.z_b = z_b[idx]
        self.odds = odds[idx]

        self.z_ml = z_ml[idx]
        self.t_ml = t_ml[idx]
        self.t_b = t_b[idx]
        self.t_ml = t_ml[idx]

        ############################################
        # Choose the photo-z to use, ml or bayesian
        ############################################
        if self.zuse == "ML":
            self.z_ph = self.z_ml
            self.type = self.t_ml
        elif self.zuse == "ZB":
            self.z_ph = self.z_b
            self.type = self.t_b

        self.g = g[idx]
        self.r = r[idx]
        self.i = i[idx]
        self.z = z[idx]
        self.g_err = g_err[idx]
        self.r_err = r_err[idx]
        self.i_err = i_err[idx]
        self.z_err = z_err[idx]

        self.g_bpz = g_bpz[idx]
        self.r_bpz = r_bpz[idx]
        self.i_bpz = i_bpz[idx]
        self.z_bpz = z_bpz[idx]
        self.g_berr = g_berr[idx]
        self.r_berr = r_berr[idx]
        self.i_berr = i_berr[idx]
        self.z_berr = z_berr[idx]

        self.class_star = class_star[idx]
        self.a_image = a_image[idx]
        self.b_image = b_image[idx]
        self.theta = theta[idx]

        # Color of selected galaxies
        self.gr = self.g_bpz - self.r_bpz
        self.ri = self.r_bpz - self.i_bpz
        self.iz = self.i_bpz - self.z_bpz

        # Min and and max values in RA/DEC
        self.ramin = self.ra.min()
        self.ramax = self.ra.max()
        self.decmin = self.dec.min()
        self.decmax = self.dec.max()

        self.idx_cat = idx

        sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t1))
        return

    ####################################
    # Read in the the probabilty file
    ####################################
    def read_probs(self):

        # The reg expresion to compile
        regexp_point = re.compile(r"arange\("
                                  r"(?P<z1>[0-9]+.[0-9]+),"
                                  r"(?P<z2>[0-9]+.[0-9]+),"
                                  r"(?P<dz>[0-9]+.[0-9]+)\)")
        t0 = time.time()
        sout.write("# Reading probs from :%s... " % self.probsfile)

        # probability arrays
        probs = []
        for line in open(self.probsfile).readlines():

            fields = line.split()
            if fields[0][0] == "#":
                point = regexp_point.search(line)
                # Extract the information if a point was selected
                if point:
                    z1 = float(point.group('z1'))
                    z2 = float(point.group('z2'))
                    dz = float(point.group('dz'))
                    zx = numpy.arange(z1, z2, dz)
                continue
            ID = fields[0]
            probs.append(numpy.asarray(list(map(float, fields[1:]))))

        # Transform the list into an N array
        p_z = numpy.asarray(probs)

        # select same galaxies as in catalogs we just read
        self.p_z = p_z[self.idx_cat][:]
        self.zx = zx
        sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t0))

        t1 = time.time()
        # Get the 1-sigma z1, z2 limits for each galaxy
        # Cumulatibe P(<z) function for each selected galaxy
        self.Psum = numpy.cumsum(self.p_z, axis=1)
        sout.write("# Getting +/- 1sigma (z1,z2) limits for each galaxy ")
        self.z1 = self.ra * 0.0
        self.z2 = self.ra * 0.0

        # One by one in the list
        for i in range(len(self.ra)):
            i1 = numpy.where(self.Psum[i, :] >= 0.159)[0][0]
            i2 = numpy.where(self.Psum[i, :] > 0.842)[0][0]
            self.z1[i] = self.zx[i1]
            self.z2[i] = self.zx[i2]

        sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t1))
        return

    ################################################
    # Get the absolute magnitudes for each object
    ################################################
    def get_absmags(self):

        # Distance modulus, dlum and dangular
        self.dlum = self.cset.dlum(self.z_ph)
        self.dang = self.cset.dang(self.z_ph)
        self.DM = 25.0 + 5.0 * numpy.log10(self.dlum)

        t0 = time.time()
        # Get the absolute magnitudes, *** not including evolution ***, only Kcorr
        # We use a BPZ E's template for Kcorr
        #sout.write("# Computing absolute magnitudes interpolating Kcorr ")
        #k = Kcorr_fit(sed='El_Benitez2003')

        # Alternatibely we can get both the kcorr and the evol from
        # the *.color file from BC03 *.ised file
        sout.write(
            "# Computing absolute magnitudes interpolating konly from BC03 model \n")
        k, ev = KEfit(self.evolfile)

        self.Mg = self.g - self.DM - k['g'](self.z_ph)
        self.Mr = self.r - self.DM - k['r'](self.z_ph)
        self.Mi = self.i - self.DM - k['i'](self.z_ph)
        self.Mz = self.z - self.DM - k['z'](self.z_ph)

        sout.write("# Computing evolution ev(z) for each galaxy ")
        self.ev_g = ev['g'](self.z_ph)
        self.ev_r = ev['r'](self.z_ph)
        self.ev_i = ev['i'](self.z_ph)
        self.ev_z = ev['z'](self.z_ph)

        # Also get the luminosities in Msun
        # taken from http://www.ucolick.org/~cnaw/sun.html
        self.Msun = {}
        self.Msun['g'] = 5.11
        self.Msun['r'] = 4.65
        self.Msun['i'] = 4.54
        self.Msun['z'] = 4.52

        # Mags k-corrected to z=0.25 as done in Reyes el al 2009
        Mg = self.g - self.DM - k['g'](self.z_ph) + k['g'](0.25)
        Mr = self.r - self.DM - k['r'](self.z_ph) + k['r'](0.25)
        Mi = self.i - self.DM - k['i'](self.z_ph) + k['i'](0.25)
        Mz = self.z - self.DM - k['z'](self.z_ph) + k['z'](0.25)
        self.Lg = 10.0**(-0.4 * (Mg - self.Msun['g']))
        self.Lr = 10.0**(-0.4 * (Mr - self.Msun['r']))
        self.Li = 10.0**(-0.4 * (Mi - self.Msun['i']))
        self.Lz = 10.0**(-0.4 * (Mz - self.Msun['z']))
        self.Lg_err = self.Lg * self.g_err / 1.0857
        self.Lr_err = self.Lr * self.r_err / 1.0857
        self.Li_err = self.Li * self.i_err / 1.0857
        self.Lz_err = self.Lz * self.z_err / 1.0857

        # Pass it up to the class
        self.kcorr = k
        self.evf = ev
        sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t0))
        return

# #########################################################################
#     Not needed, all done now at self.get_absmag()
#
#     ###########################################################
#     # Get the evolutionay correction for each object's redshift
#     ###########################################################
#     def get_evol(self):
#         t0  = time.time()
#         sout.write("# Computing evolution")
#         ev = evolfit(self.evolfile)
#         self.evf = ev # pass the function to the class
#         self.ev_g = numpy.asarray(ev['g'](self.z_ph))
#         self.ev_r = numpy.asarray(ev['r'](self.z_ph))
#         self.ev_i = numpy.asarray(ev['i'](self.z_ph))
#         self.ev_z = numpy.asarray(ev['z'](self.z_ph))
#         sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t0))
############################################################################

###########################################################################
# intergrates the p_z Bayesian probability between zlow and zhigh interval
###########################################################################

    def p_z_int(self, z1, z2, idx=None):
        #sout.write("# Integrating p(z) between %s -- %s ..." % (z1,z2))
        # the limits of integration of the probabilty
        i1 = numpy.where(self.zx >= z1)[0][0]
        i2 = numpy.where(self.zx >= z2)[0][0]
        #dz = abs(self.zx[1]-self.zx[0])
        if idx:
            p_int = numpy.sum(self.p_z[idx][:, i1:i2], axis=1)
        else:
            p_int = numpy.sum(self.p_z[:, i1:i2], axis=1)
        #sout.write(" Done\n")
        return p_int

    ################################################
    # Set the coeffs of P(r) and normalize it unity
    ################################################
    def set_Pr(self, zo, dang, ro=0.180, n=1.e5):  # ro in [Mpc]

        # Set the cosmology
        #dang = self.cset.dang(zo)[0]
        #dang = 1.0

        # Set rc and rmax
        #scale    = 180.*60/math.pi
        scale = old_div(180.0, math.pi)  # Mpc to degrees
        self.rc = ro * scale / dang
        #self.rmax = 10.0*self.rc
        self.rmax = 1.0 * scale / dang

        # Normalize the function to unity
        #sout.write("# Normalizing P(r,z) at z=%s... " % zo)
        r1 = 0.
        r2 = self.rmax * 1.1
        dr = old_div((r2 - r1), n)  # arcmins in 10^3 steps
        rx = numpy.arange(r1, r2, dr)
        self.Pn = 1.0  # reset normalization before normalizing
        self.Pn = 2.0 * math.pi * (self.P_r(rx) * rx * dr).sum()

        #print self.Pn
        #print self.P_sum()
        #sout.write(" Done\n")
        return

    #####################################
    # Profile weight, r is in arcmins
    #####################################
    def P_r(self, r):

        rc = self.rc
        rmax = self.rmax

        #r    = extras.asnumpy.r)
        #r    = asarray(r) # fix to work with numpy.and numpy
        r = numpy.asarray(r)
        Pr = r * 0.0

        idx = numpy.where(r <= rmax)  # P(r) = 0  for r > rmax
        ri = r[idx]
        Pr[idx] = old_div(1.0, numpy.sqrt(1 + (old_div(ri, rc))**2)) - old_div(
            1.0, numpy.sqrt(1 + (old_div(rmax, rc))**2))
        return old_div(Pr, self.Pn)

    ##############################################
    # Normalize phi and set the Lum weight, mstar
    ##############################################
    def set_L(self, zo, alpha=-0.5, n=1.e5):

        self.alpha = alpha
        self.mstar = mi_star(zo, self.cosmo)

        # intergration limits for normalization
        m2 = self.maglim
        m1 = 10.0
        dm = old_div(abs(m2 - m1), n)

        # Normalize the schecter luminosity function -- GET phi*
        sout.write("# Normalizing L(m,z) at z=%s... " % zo)
        m = numpy.arange(m1, m2, dm)
        phi = PHI(m, self.mstar, self.alpha)
        self.phinorm = (phi * dm).sum()
        sout.write(" Done\n")
        return

    ####################
    # Luminosity weight
    ####################
    def L(self, m):

        # Background number of galaxies from number counts Yasuda et al (2001)
        b = 0.537 * 10**(-0.4 * (m - 20.))
        # Postman Luminosity weight
        phi = PHI(m, self.mstar, self.alpha)
        L = phi / b / self.phinorm
        return L

    # Just to check the normalization
    def P_sum(self):

        r1 = 0.
        r2 = self.rmax * 1.1
        dr = old_div((r2 - r1), 1000.)  # arcmins
        rx = numpy.arange(r1, r2, dr)
        return 2 * math.pi * (self.P_r(rx) * rx * dr).sum()

    ############################################################################
    # make the jpeg for each candidate and label each accordingly, A, B, C, etc.
    ############################################################################
    def make_jpeg(self, k=None):

        t0 = time.time()

        # Order is important for color image
        filters = ('i', 'r', 'g')

        # Use the candidates (ra,dec) as position, default
        if k == None:
            xo = self.RA
            yo = self.DEC
            #zo = self.z_c
            zo = self.zcl
        # Use the ranked (A) BCG as the center
        else:
            xo = astrometry.dec2deg(old_div(self.ra_BCG[k], 15.))
            yo = astrometry.dec2deg(self.dec_BCG[k])
            zo = self.z_BCG[k]
            zo = self.zcl

        # Find out tile name for ranking BCG ID name
        tile = self.tile
        # Size of the jpeg in pixels depending on the redshift
        size = old_div(1000, self.h)  # 1h^1Mpc [in kpc]
        dx = old_div(astrometry.kpc2arc(zo, size, self.cosmo), self.pixscale)
        # Avoid crashing getfits, it fails for size ~ 1900 pix and above
        if dx > 1900:
            dx = 1900
            # Call imhead instead
        dx = int(dx)
        dy = dx
        sout.write("# %.3f Mpc @ z=%s is %s x %s pixels\n" %
                   (old_div(size, 1000.0), zo, dx, dy))
        # Cut the fits files using getfits
        files = ""
        for filter in filters:
            fitsfile = os.path.join(self.path, tile,
                                    tile + filter + "_ext.fits")
            fitsout = os.path.join(self.outpath,
                                   "tmp_%s%s.fits" % (tile, filter))
            cmd = "getfits %s %s %s %s %s -o %s > /dev/null 2>&1 " % (
                fitsfile, xo, yo, dx, dy, fitsout)
            os.system(cmd)
            files = files + fitsout + " "

        # Keep names in the class
        self.jpeg_name = os.path.join(self.outpath, "%s.jpg" % self.cID)
        self.tiff_name = os.path.join(self.outpath, "%s.tif" % self.cID)
        self.fitsfile = fitsout
        self.tmp_fits = files
        self.dx = dx
        self.dy = dy
        # Make the color tiff
        conf = os.path.join(os.environ['BCSPIPE'], 'LIB/stiff.conf')
        opts = {}
        opts["OUTFILE_NAME"] = self.tiff_name
        opts["BINNING"] = 1
        opts['GAMMA'] = 2.0 + (zo - 0.2)  # For better contrast at higher z
        opts['MAX_LEVEL'] = 0.99 - old_div((zo - 0.2), 50.)
        opts['VERBOSE_TYPE'] = "QUIET"
        cmd = "stiff "
        cmd = cmd + " %s  -c %s " % (files, conf)
        for param, value in list(opts.items()):
            cmd = cmd + "-%s %s " % (param, value)
        os.system(cmd)
        # Make the color jpg
        cmd = "convert %s %s" % (self.tiff_name, self.jpeg_name)
        os.system(cmd)
        return

    ############################################################################
    # make the jpeg for each candidate and label each accordingly, A, B, C, etc.
    ############################################################################
    def make_jpeg_ID(self, k=0):

        t0 = time.time()

        # Order is important for color image
        filters = ('i', 'r', 'g')

        xo = astrometry.dec2deg(old_div(self.ra_BCG[k], 15.))
        yo = astrometry.dec2deg(self.dec_BCG[k])
        #zo = self.z_BCG[k]
        zo = self.zcl

        # Find out tile name for ranking BCG ID name
        tile = self.tile
        # Size of the jpeg in pixels depending on the redshift
        size = old_div(1000, self.h)  # 1h^1Mpc [in kpc]
        dx = 2 * astrometry.kpc2arc(zo, size, self.cosmo) / self.pixscale
        dx = int(dx)
        dy = dx

        sout.write("# 2 x %.3f Mpc @ z=%s is %s x %s pixels\n" %
                   (old_div(size, 1000.0), zo, dx, dy))
        # Cut the fits files using getfits
        files = ""
        for filter in filters:
            fitsfile = os.path.join(self.path, tile,
                                    tile + filter + "_ext.fits")
            fitsout = os.path.join(self.outpath,
                                   "tmp_%s%s.fits" % (tile, filter))
            files = files + fitsout + " "
            xpix, ypix = astrometry.rd2xy(self.ra_BCG[k], self.dec_BCG[k],
                                          fitsfile)
            cutfits(fitsfile, xpix, ypix, dx, dy, fitsout)

            # Keep names in the class
        self.jpeg_name = os.path.join(self.outpath, "%s.jpg" % self.SCSname)
        self.tiff_name = os.path.join(self.outpath, "%s.tif" % self.SCSname)
        self.fitsfile = fitsout
        self.tmp_fits = files
        self.dx = dx
        self.dy = dy
        # Make the color tiff
        conf = os.path.join(os.environ['BCSPIPE'], 'LIB/stiff.conf')
        opts = {}
        opts["OUTFILE_NAME"] = self.tiff_name
        opts["BINNING"] = 1
        opts['GAMMA'] = 2.0 + (zo - 0.2)  # For better contrast at higher z
        opts['MAX_LEVEL'] = 0.99 - old_div((zo - 0.2), 50.)
        opts['VERBOSE_TYPE'] = "QUIET"
        cmd = "stiff "
        cmd = cmd + " %s  -c %s " % (files, conf)
        for param, value in list(opts.items()):
            cmd = cmd + "-%s %s " % (param, value)
        os.system(cmd)
        # Make the color jpg
        cmd = "convert %s %s" % (self.tiff_name, self.jpeg_name)
        os.system(cmd)
        return

    ##############################
    # Draw the BGCs in the jpeg
    ##############################
    def ellipse_BCGs(self, index=None):

        sout.write("# Drawing ellipses + info for %s\n" % self.cID)
        self.jpg_array = pilutil.imread(self.jpeg_name)
        self.jpg_region = self.jpg_array  #[y1:y2, x1:x2, :]
        (ny, nx, nz) = self.jpg_array.shape

        x1 = 0.05
        y1 = 0.05
        dsx = old_div((0.95 - x1), 3.)
        dsy = old_div((0.99 - y1), 4.)
        dx = 3.0 * dsx
        dy = 3.0 * dsy
        y1 = dsy
        fig = pylab.figure(1, figsize=(9, 10))
        ax1 = pylab.axes([x1, y1, dx, dy])
        pylab.imshow(self.jpg_region)
        # Change ax to arcmin
        self.ax_to_arcmin(ds=1.0)
        pylab.xlabel("x[arcmin]")
        pylab.ylabel("y[arcmin]")

        # Add a green cross at the center of the candidate's position
        pylab.plot([old_div(nx, 2.0)], [old_div(ny, 2.0)],
                   'gx',
                   markersize=15,
                   markeredgewidth=1.0)

        # Radius of search in from arcmin --> pixels
        radius = self.dmin * 60.0 / self.pixscale

        # And a circle in the search area
        C = PCircle((old_div(nx, 2.0), old_div(ny, 2.0)),
                    radius,
                    resolution=80,
                    fill=0,
                    edgecolor="white",
                    linestyle='dashed',
                    linewidth=0.3)
        ax1.add_patch(C)

        # Label with radius search size and singal of detection
        ax1.text(
            old_div(nx, 20.),
            old_div(ny, 10.),
            "Dm=%.2f'\nSn=%.2f" % (self.dmin, self.Sn),
            color='white',
            #family='monospace',
            horizontalalignment='left',
            fontsize=11)

        # If no BCGs, plot simple and exit
        if index == None:
            RA = self.RA
            DEC = self.DEC
            #zo  = self.z_c
            zo = self.zc
            pylab.title("%s ---  %s, %s -- zc:%.1f (NO BCGs)" %
                        (self.cID, RA, DEC, zo),
                        fontsize=10)
            pylab.savefig(os.path.join(self.outpath, "%s_finder.png" %
                                       self.cID))
            pylab.close()
            # Clean up some files
            os.system("rm %s %s" % (self.jpeg_name, self.tiff_name))
            # Remove temporary fits files
            os.system("rm %s" % self.tmp_fits)
            return

        RA = astrometry.dec2deg(old_div(self.ra_BCG[index[0]], 15.0))
        DEC = astrometry.dec2deg(self.dec_BCG[index[0]])
        zo = self.z_BCG[index[0]]
        pylab.title("%s ---  %s, %s -- z:%.3f (A)" % (self.cID, RA, DEC, zo),
                    fontsize=10)

        # construct the ellipse for the current display
        label = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L',
                 'M', 'N', 'O', 'Q', 'R', 'S', 'T']
        j = 0
        for i in index[0:7]:
            ra = self.ra_BCG[i]
            dec = self.dec_BCG[i]
            a = self.a_BCG[i]
            b = self.b_BCG[i]
            theta = self.theta_BCG[i]  #*math.pi/180.0
            (xo, yo) = astrometry.rd2xy(ra, dec, self.fitsfile)
            # Change the referece pixel to reflect jpg standards where the
            # origin is at (0,ny), is the upper left corner
            yo = ny - yo
            E = PEllipse((xo, yo), (a, b),
                         resolution=80,
                         angle=theta,
                         fill=0,
                         edgecolor="red",
                         linewidth=0.5)
            ax1.add_patch(E)

            #angle = random.randint(0,360)
            angle = theta
            xsh = 3 * a * math.cos(angle * math.pi / 180.)
            ysh = 3 * b * math.sin(angle * math.pi / 180.)

            # Draw labels, A,B,C, etc.
            ax1.annotate(label[j],
                         xy=(xo, yo),
                         xycoords='data',
                         xytext=(xsh, ysh),
                         textcoords='offset points',
                         arrowprops=dict(arrowstyle="->",
                                         connectionstyle="arc3,rad=.2",
                                         edgecolor='white',
                                         shrinkA=0.0,
                                         linewidth=0.5),
                         bbox=dict(boxstyle="round",
                                   fc="1.0",
                                   edgecolor=(1., .5, .5),
                                   fill=1,
                                   alpha=0.4))
            #bbox=dict(boxstyle="round", fc=(1.0, 0.7, 0.7), ec=(1., .5, .5),alpha=0.5),
            #bbox=dict(boxstyle="round", fc=(1.0, 0.0, 0.0), ec=(1., .5, .5),alpha=0.5),
            #arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2",edgecolor='red'),color='white')
            j = j + 1

        # The sub-plot with the info text
        ya = 0
        yb = 10
        x1 = 0.05
        y1 = 0.01
        dx = dx
        dy = dsy - 0.05
        ax2 = pylab.axes([x1, y1, dx, dy])
        # Header with info
        x = 0
        y = yb - 2

        header = " %-2s %6s %6s %5s %5s %6s %7s %7s %8s %8s %8s %6s" % (
            '', 'z_B', 'z_ML', 'Ngal', 'N200', 'R200', 'M_r', 'm_i', 'p(BCG)',
            'P(color)', 'P(total)', 'D\"')
        ax2.text(x,
                 y,
                 header,
                 family='monospace',
                 horizontalalignment='left',
                 fontsize=11)
        j = 0
        format = " %-2s %6.3f %6.3f %5d %5d %6.2f %7.2f %7.2f %8.3f %8.1f %8.1f %6.2f"
        for i in index[0:7]:
            y = yb - 3 - j
            vars = (label[j], self.zb_BCG[i], self.zml_BCG[i], self.Ngal[i],
                    self.N200[i], self.R200[i], self.Mr_BCG[i], self.i_BCG[i],
                    self.p_BCG[i], self.P_color[i], self.P_total[i],
                    self.d_BCG[i])
            texto = format % vars
            ax2.text(x,
                     y,
                     texto,
                     family='monospace',
                     horizontalalignment='left',
                     fontsize=11)

            j = j + 1

        pylab.axis('off')
        pylab.ylim(ya, yb)
        pylab.savefig(os.path.join(self.outpath, "%s_finder.png" % self.cID))
        pylab.close()

        # Clean up some files
        os.system("rm %s %s" % (self.jpeg_name, self.tiff_name))
        # Remove temporary fits files
        os.system("rm %s" % self.tmp_fits)

        return

    ##########################################
    # Write the BCGs candidates information
    ##########################################
    def write_BCGs(self, k):

        n = len(k)
        label = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L',
                 'M', 'N', 'O', 'Q', 'R', 'S', 'T']
        # Make ids a Char String in numarray
        lab = nstr.array(label)
        header = "# %2s %20s %6s %11s %11s %5s %5s %6s %7s %7s %8s %8s %7s %7s %6s" % (
            '', 'ID', 'z', 'RA', 'DEC', 'Ngal', 'N200', 'R200\'', 'p(BCG)',
            'P_z', 'P(color)', 'P(total)', 'M_r', 'm_i', 'D\"')
        format = "%4s %20s %6.3f %11.6f %11.6f %5d %5d %6.2f %7.3f %7.1f %8.1f %8.1f %7.2f %7.2f %6.2f"
        vars = (lab[0:n], self.id_BCG[k], self.z_BCG[k], self.ra_BCG[k],
                self.dec_BCG[k], self.Ngal[k], self.N200[k], self.R200[k],
                self.p_BCG[k], self.P_z[k], self.P_color[k], self.P_total[k],
                self.Mr_BCG[k], self.i_BCG[k], self.d_BCG[k])

        P_file = os.path.join(self.outpath, "%s_Ptotal.dat" % self.cID)
        tableio.put_data(P_file, vars, format=format, header=header)
        sout.write("# BCGs info in %s\n" % P_file)
        return

    #####################################
    # Draw the BGCs and cluster members
    #####################################
    def ellipse_members(self, k=0):

        sout.write("# Drawing ellipses + info for %s\n" % self.SCSname)
        self.jpg_array = pilutil.imread(self.jpeg_name)
        self.jpg_region = self.jpg_array  #[y1:y2, x1:x2, :]
        (ny, nx, nz) = self.jpg_array.shape

        x1 = 0.05
        y1 = 0.05
        dsx = old_div((0.95 - x1), 3.)
        dsy = old_div((0.99 - y1), 4.)
        dx = 3.0 * dsx
        dy = 3.0 * dsy
        y1 = dsy
        fig = pylab.figure(1, figsize=(9, 10))
        ax1 = pylab.axes([x1, y1, dx, dy])
        pylab.imshow(self.jpg_region)
        # Change ax to arcmin
        self.ax_to_arcmin(ds=1.0)
        pylab.xlabel("x[arcmin]")
        pylab.ylabel("y[arcmin]")

        RA = astrometry.dec2deg(old_div(self.ra_BCG[k], 15.0))
        DEC = astrometry.dec2deg(self.dec_BCG[k])
        zo = self.z_BCG[k]
        pylab.title("%s -  %s, %s - z:%.3f Ngal:%d" %
                    (self.SCSname, RA, DEC, zo, self.N1Mpc),
                    fontsize=10)

        # construct the ellipses for each members
        for i in self.iR1Mpc[0]:
            ra = self.ra[i]
            dec = self.dec[i]
            a = self.a_image[i]
            b = self.b_image[i]
            theta = self.theta[i]  #*math.pi/180.0
            (xo, yo) = astrometry.rd2xy(ra, dec, self.fitsfile)
            # Change the referece pixel to reflect jpg standards where the
            # origin is at (0,ny), is the upper left corner
            yo = ny - yo
            if i == self.idx_BCG[0]:
                ec = 'blue'
            else:
                ec = 'green'
            E = PEllipse((xo, yo), (a, b),
                         resolution=80,
                         angle=theta,
                         fill=0,
                         edgecolor=ec,
                         linewidth=0.5)
            ax1.add_patch(E)

        # And a circle of 1Mpc/h radiys
        radius = self.r1Mpc * 3600.0 / self.pixscale
        C = PCircle((old_div(nx, 2.0), old_div(ny, 2.0)),
                    radius,
                    resolution=80,
                    fill=0,
                    edgecolor="white",
                    linestyle='dashed',
                    linewidth=0.5)
        ax1.add_patch(C)

        # The sub-plot with the info text
        ya = 0
        yb = 10
        x1 = 0.05
        y1 = 0.01
        dx = dx
        dy = dsy - 0.05
        ax2 = pylab.axes([x1, y1, dx, dy])
        # Header with info
        x = 0
        y = yb - 2

        header = "%6s %6s %5s %5s %6s %6s %7s %7s %8s %8s %8s" % (
            'z_B', 'z_ML', 'Ngal', 'N200', 'R200', 'r200', 'M_r', 'm_i',
            'p(BCG)', 'P(color)', 'P(total)')
        ax2.text(x,
                 y,
                 header,
                 family='monospace',
                 horizontalalignment='left',
                 fontsize=11)
        format = "%6.3f %6.3f %5d %5d %6.2f %6.2f %7.2f %7.2f %8.3f %8.1f %8.1f"
        y = yb - 3
        vars = (self.zb_BCG[k], self.zml_BCG[k], self.N1Mpc, self.N200,
                self.R200, self.r200 * 60, self.Mr_BCG[k], self.i_BCG[k],
                self.p_BCG[k], self.P_color, self.P_total)
        texto = format % vars
        ax2.text(x,
                 y,
                 texto,
                 family='monospace',
                 horizontalalignment='left',
                 fontsize=11)
        pylab.axis('off')
        pylab.ylim(ya, yb)

        pylab.savefig(
            os.path.join(self.outpath, "%s.png" % self.SCSname),
            dpi=300)
        pylab.close()
        # Clean up some files
        os.system("rm %s %s" % (self.jpeg_name, self.tiff_name))
        # Remove temporary fits files
        os.system("rm %s" % self.tmp_fits)
        return

    ##############################
    # Change the axes to arcmins
    ###############################
    def ax_to_arcmin(self, ds=1.0):  # ds in arcmin

        [xmin, xmax, ymin, ymax] = pylab.axis()
        dx = self.dx
        dy = self.dy
        scale = old_div(self.pixscale, 60.)  # in arcmin

        xo = old_div((xmin + xmax), 2.0)
        yo = old_div((ymin + ymax), 2.0)
        s1 = int((old_div(-dx, 2.0)) * scale)
        s2 = int((old_div(+dx, 2.0)) * scale)
        sx = numpy.arange(s1, s2 + 0.05, ds)

        xtext = []
        xtick = []
        for s in sx:
            x = xo + old_div(s, scale)  #+ ds/scale
            xtick.append(x)
            xtext.append("%.1f" % s)

        s1 = int((old_div(-dy, 2.0)) * scale)
        s2 = int((old_div(+dy, 2.0)) * scale)
        sy = numpy.arange(s1, s2 + 0.05, ds)

        ytext = []
        ytick = []
        for s in sy:
            y = yo + old_div(s, scale)  #+ ds/scale
            ytick.append(y)
            ytext.append("%.1f" % s)
        pylab.yticks(ytick, tuple(ytext))
        pylab.xticks(xtick, tuple(xtext))
        # Make sure we plot everithing
        pylab.xlim(xmin, xmax)
        pylab.ylim(ymin, ymax)
        return

    ########################################################
    # Modified/updated from find_clusters_ext_auto.py
    # Select galaxies around ID galaxy un redshift range
    ########################################################
    def select_members(self, i, dz=0.1, Mi_lim=-20.25):

        t0 = time.time()
        sout.write("# Selecting Cluster members... Ngal, N200, R200 ")
        # Get the relevant info for ith BCG
        zo = self.z_BCG[i]
        ra0 = self.ra_BCG[i]
        dec0 = self.dec_BCG[i]
        Mi_BCG = self.Mi_BCG[i]
        DM = self.DM_BCG[i]
        ID_BCG = self.id_BCG[i]

        # 1 - Select in position around ra0,dec0
        # Define 1h^-1 Mpc radius in degress @ zo
        R1Mpc = 1000 * 1.0 / self.h  # in kpc
        r1Mpc = old_div(
            astrometry.kpc2arc(zo, R1Mpc, self.cosmo), 3600.)  # in degrees.
        rcore = old_div(r1Mpc, 2.0)
        dist = astrometry.circle_distance(ra0,
                                          dec0,
                                          self.ra,
                                          self.dec,
                                          units='deg')
        mask_R1Mpc = numpy.where(dist <= r1Mpc, 1, 0)
        mask_rcore = numpy.where(dist <= rcore, 1, 0)
        arcmin2Mpc = old_div(
            astrometry.arc2kpc(zo, 60.0, self.cosmo),
            1000.0)  # scale between arcmin and Mpc

        # 2 - Select in redshift
        #dz = self.dz*(1 + zo)
        z1 = zo - dz
        z2 = zo + dz
        mask_z = numpy.where(land(self.z_ph >= z1, self.z_ph <= z2), 1, 0)

        # 3 - Select in brightness
        Mi_lim_zo = Mi_lim + self.evf['i'](zo) - self.evf['i'](0.1)
        mask_L1 = numpy.where(self.Mi <= Mi_lim_zo, 1, 0)  # Faint  cut > 0.4L*
        mask_L2 = numpy.where(self.Mi >= Mi_BCG, 1, 0)  # Bright cut < L_BCG

        # The final selection mask, position x redshift x Luminosity
        idx = numpy.where(mask_R1Mpc * mask_L1 * mask_L2 * mask_z == 1)
        idc = numpy.where(mask_rcore * mask_L1 * mask_L2 * mask_z == 1)

        # Shot versions handles
        gr = self.gr
        ri = self.ri
        iz = self.iz

        # Some simple 3-sigma clipping defined using r< rcore
        Nsigma = 3.0
        loop = 1
        converge = False
        while not converge:
            # The conditions to apply
            c1 = numpy.abs(gr[idc] - gr[idc].mean()) > Nsigma * numpy.std(
                gr[idc], ddof=1)
            c2 = numpy.abs(ri[idc] - ri[idc].mean()) > Nsigma * numpy.std(
                ri[idc], ddof=1)
            c3 = numpy.abs(iz[idc] - iz[idc].mean()) > Nsigma * numpy.std(
                iz[idc], ddof=1)
            iclip = numpy.where(lor(c1, c2,
                                    c3))  # where any of the conditions fails
            if len(iclip[0]) > 0:
                idc = numpy.delete(idc, iclip[0])  # Removed failed ones
                converge = False
            else:
                converge = True
            loop = loop + 1

        sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t0))

        # Or we can make a new mask where the condition's are true
        c1 = numpy.abs(self.gr - gr[idc].mean()) > Nsigma * numpy.std(gr[idc],
                                                                      ddof=1)
        c2 = numpy.abs(self.ri - ri[idc].mean()) > Nsigma * numpy.std(ri[idc],
                                                                      ddof=1)
        c3 = numpy.abs(self.iz - iz[idc].mean()) > Nsigma * numpy.std(iz[idc],
                                                                      ddof=1)
        mask_cm = numpy.where(lor(c1, c2, c3), 0, 1)  # where condition fails
        iR1Mpc = numpy.where(
            mask_R1Mpc * mask_L1 * mask_L2 * mask_z * mask_cm == 1)
        Ngal = len(iR1Mpc[0])
        sout.write("# Total: %s objects selected in 1h^-1Mpc around %s\n" %
                   (Ngal, self.id_BCG[i]))

        #############################################################################
        # We'll skip 200 measurement as they depend on the corrected values of Ngal
        # Now let's get R200 and N200
        R200 = 0.156 * (Ngal**0.6) / self.h  # In Mpc
        r200 = old_div(
            astrometry.kpc2arc(zo, R200 * 1000.0, self.cosmo),
            3600.)  # in degrees.
        mask_r200 = numpy.where(dist <= r200, 1, 0)
        i200 = numpy.where(
            mask_r200 * mask_L1 * mask_L2 * mask_z * mask_cm == 1)
        N200 = len(i200[0])
        self.i200 = i200
        self.N200 = N200
        self.R200 = R200
        self.r200 = r200
        self.L200 = self.Lr[i200].sum()
        ############################################################################

        # And the value for all galaxies up NxR1Mpc -- change if required.
        mask_R = numpy.where(dist <= 10 * r1Mpc, 1, 0)
        iR = numpy.where(mask_R * mask_L1 * mask_L2 * mask_z * mask_cm == 1)

        # Pass up
        self.iR = iR
        self.iR1Mpc = iR1Mpc
        self.N1Mpc = Ngal
        self.r1Mpc = r1Mpc  # in degress
        self.dist2BCG = dist
        self.arcmin2Mpc = arcmin2Mpc

        # Sort indices radially for galaxies < N*R1Mpc, will be used later
        i = numpy.argsort(self.dist2BCG[iR])
        self.ix_radial = iR[0][i]

        # We want to keep i200 and iR1Mpc to write out members.
        return Ngal, N200, R200  # iR1Mpc,i200

    ########################################################
    # Modified/updated from find_clusters_ext_auto.py
    # Select galaxies around ID galaxy un redshift range
    ########################################################
    def select_members_redshift(self,
                                i,
                                dz=0.05,
                                Mi_lim=-20.25,
                                zo=None,
                                radius=1000.0,
                                weight=None):

        t0 = time.time()
        sout.write("# Selecting Cluster members... Ngal, N200, R200 ")
        # Get the relevant info for ith BCG
        ra0 = self.ra_BCG[i]
        dec0 = self.dec_BCG[i]
        Mi_BCG = self.Mi_BCG[i]
        DM = self.DM_BCG[i]
        ID_BCG = self.id_BCG[i]
        if zo:
            print("Will use z:%.3f for cluster" % zo)
        else:
            zo = self.z_BCG[i]

        # 1 - Select in position around ra0,dec0
        # Define 1h^-1 Mpc radius in degress @ zo
        R1Mpc = radius * 1.0 / self.h  # in kpc
        r1Mpc = old_div(
            astrometry.kpc2arc(zo, R1Mpc, self.cosmo), 3600.)  # in degrees.
        rcore = old_div(r1Mpc, 2.0)
        dist = astrometry.circle_distance(ra0,
                                          dec0,
                                          self.ra,
                                          self.dec,
                                          units='deg')
        mask_R1Mpc = numpy.where(dist <= r1Mpc, 1, 0)
        mask_rcore = numpy.where(dist <= rcore, 1, 0)
        arcmin2Mpc = old_div(
            astrometry.arc2kpc(zo, 60.0, self.cosmo),
            1000.0)  # scale between arcmin and Mpc

        # 2 - Select in redshift
        #dz = self.dz*(1 + zo)
        z1 = zo - dz
        z2 = zo + dz
        mask_z = numpy.where(land(self.z_ph >= z1, self.z_ph <= z2), 1, 0)

        # 3 - Select in brightness
        Mi_lim_zo = Mi_lim + self.evf['i'](zo) - self.evf['i'](0.1)
        mask_L1 = numpy.where(self.Mi <= Mi_lim_zo, 1, 0)  # Faint  cut > 0.4L*
        mask_L2 = numpy.where(self.Mi >= Mi_BCG, 1, 0)  # Bright cut < L_BCG

        # The final selection mask, position x redshift x Luminosity
        idx = numpy.where(mask_R1Mpc * mask_L1 * mask_L2 * mask_z == 1)
        idc = numpy.where(mask_rcore * mask_L1 * mask_L2 * mask_z == 1)

        # Shot versions handles
        gr = self.gr
        ri = self.ri
        iz = self.iz

        # Some simple 3-sigma clipping defined using r< rcore
        Nsigma = 3.0
        loop = 1
        converge = False
        while not converge:
            # The conditions to apply
            c1 = numpy.abs(gr[idc] - gr[idc].mean()) > Nsigma * numpy.std(
                gr[idc], ddof=1)
            c2 = numpy.abs(ri[idc] - ri[idc].mean()) > Nsigma * numpy.std(
                ri[idc], ddof=1)
            c3 = numpy.abs(iz[idc] - iz[idc].mean()) > Nsigma * numpy.std(
                iz[idc], ddof=1)
            iclip = numpy.where(lor(c1, c2,
                                    c3))  # where any of the conditions fails
            if len(iclip[0]) > 0:
                idc = numpy.delete(idc, iclip[0])  # Removed failed ones
                converge = False
            else:
                converge = True
            loop = loop + 1

        # Put it back in case we missed the BCG
        #if self.idx_BCG[0][0] not in idc[0]:
        #    i0 = numpy.append(self.idx_BCG[0][0],idc[0])
        #    idc = (i0,)

        if weight:
            # Get the weighted-average redshift within the core:
            dz = 0.5 * numpy.abs(self.z2[idc] - self.z1[idc])
            z_cl, z_clrms = aux.statsw(self.z_ph[idc], weight=old_div(1.0, dz))
        else:
            # Get the average redshift within the core:
            z_cl = numpy.median(self.z_ph[idc])
            z_clrms = self.z_ph[idc].std()

        sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t0))

        # Or we can make a new mask where the condition's are true
        c1 = numpy.abs(self.gr - gr[idc].mean()) > Nsigma * numpy.std(gr[idc],
                                                                      ddof=1)
        c2 = numpy.abs(self.ri - ri[idc].mean()) > Nsigma * numpy.std(ri[idc],
                                                                      ddof=1)
        c3 = numpy.abs(self.iz - iz[idc].mean()) > Nsigma * numpy.std(iz[idc],
                                                                      ddof=1)
        mask_cm = numpy.where(lor(c1, c2, c3), 0, 1)  # where condition fails
        iR1Mpc = numpy.where(
            mask_R1Mpc * mask_L1 * mask_L2 * mask_z * mask_cm == 1)
        Ngal = len(iR1Mpc[0])
        sout.write("# Total: %s objects selected in %sh^-1Mpc around %s\n" %
                   (Ngal, old_div(radius, 1000.0), self.id_BCG[i]))

        #############################################################################
        # We'll skip 200 measurement as they depend on the corrected values of Ngal
        # Now let's get R200 and N200
        R200 = 0.156 * (Ngal**0.6) / self.h  # In Mpc
        r200 = old_div(
            astrometry.kpc2arc(zo, R200 * 1000.0, self.cosmo),
            3600.)  # in degrees.
        mask_r200 = numpy.where(dist <= r200, 1, 0)
        i200 = numpy.where(
            mask_r200 * mask_L1 * mask_L2 * mask_z * mask_cm == 1)
        N200 = len(i200[0])
        self.i200 = i200
        self.N200 = N200
        self.R200 = R200
        self.r200 = r200
        self.L200 = self.Lr[i200].sum()
        ############################################################################

        # And the value for all galaxies up NxR1Mpc -- change if required.
        mask_R = numpy.where(dist <= 10 * r1Mpc, 1, 0)
        iR = numpy.where(mask_R * mask_L1 * mask_L2 * mask_z * mask_cm == 1)

        # Pass up
        self.iR = iR
        self.iR1Mpc = iR1Mpc
        self.N1Mpc = Ngal
        self.r1Mpc = r1Mpc  # in degress
        self.dist2BCG = dist
        self.arcmin2Mpc = arcmin2Mpc

        # Sort indices radially for galaxies < N*R1Mpc, will be used later
        i = numpy.argsort(self.dist2BCG[iR])
        self.ix_radial = iR[0][i]

        # We want to keep i200 and iR1Mpc to write out members.
        return Ngal, N200, R200, z_cl, z_clrms

    ############################
    # Header for the info file
    ############################
    def write_head(self, k=0):
        header = ""
        header = header + "# catID  : %s \n" % self.id_BCG[k]
        header = header + "# RA,DEC : %s %s \n" % (self.RA, self.DEC)
        header = header + "# ra,dec : %s %s \n" % (self.ra_BCG[k],
                                                   self.dec_BCG[k])
        header = header + "# zo(BCG): %6.3f \n" % self.z_BCG[k]
        header = header + "# zcl    : %6.3f +/- %6.3f\n" % (self.zcl,
                                                            self.zcl_err)
        header = header + "# Ngal (1Mpc/h) : %3d   \n" % self.Ngal
        header = header + "# N200          : %3d   \n" % self.N200
        header = header + "# Cosmology     : %s, %s, %s \n" % self.cosmo
        header = header + "# R(1Mpc/h)     : %8.3f [arcmin] \n" % (self.r1Mpc *
                                                                   60.0)
        header = header + "# R200          : %8.3f [arcmin]\n" % (self.r200 *
                                                                  60.0)
        header = header + "# R200          : %8.3f [Mpc]\n" % self.R200
        header = header + "# 1 arcmin @ zo : %8.3f [Mpc] \n" % self.arcmin2Mpc
        header = header + "# \n"
        self.header = header
        return

    ########################################
    # Write the clusters members in a file
    ########################################
    def write_members(self, k=0):

        k = self.ix_radial  #= iR[0][i]

        # Do the header first
        self.write_head()

        # All E/S0 within N x R200
        iR = self.iR
        # The key describing the position
        self.keyR = self.dist2BCG * 0 - 1
        self.keyR[self.dist2BCG <= self.r200] = 0
        self.keyR[self.dist2BCG <= self.r1Mpc] = 1

        # Sort by distance
        filename = os.path.join(self.outpath, "%s.dat" % self.SCSname)
        header = self.header + "# %-22s %12s %12s %8s %8s %8s %8s %8s %8s %10s %10s %8s %5s" % (
            ' catID', 'RA', 'DEC', 'z_b', 'z_ml', 'm_r', 'm_i', 'Mr', 'Mi',
            'Lr[Mo]', 'Li[Mo]', 'd[arcmin]', 'R1Mpc')
        format = "%-24s %12.6f %12.6f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %10.3e %10.3e %8.3f %5d"
        vars = (self.id[k], self.ra[k], self.dec[k], self.z_b[k], self.z_ml[k],
                self.r[k], self.i[k], self.Mr[k], self.Mi[k], self.Lr[k],
                self.Li[k], self.dist2BCG[k] * 60, self.keyR[k])
        tableio.put_data(filename, vars, format=format, header=header)
        # The indexes redially sorted
        self.ix_radial = k
        return

        ##########################################
        # Compute the Background for the clusters
        ##########################################
    def background(self, k=0):

        ixr = self.ix_radial
        #zo  = self.z_BCG[k]
        zo = self.zcl

        # Store radially ordered
        r = self.dist2BCG[ixr] * 60.0  # in arcmin
        Lr = self.Lr[ixr]  # We do in the r-band as Reyes et al

        # Bin the Ngal/Lum data in log spacing
        n = 15
        rbin = mklogarray(0.0, r.max(), n)
        Nbin, rcenter = histo(r, rbin, center='yes')
        Lbin, rcenter = bin_data(r, Lr, rbin, center='yes')

        # Compute the area in each shell
        ir = numpy.indices(rbin.shape)[0]
        ir1 = ir[:-1]
        ir2 = ir[1:]
        r1 = rbin[ir1]
        r2 = rbin[ir2]
        abin = math.pi * (r2**2 - r1**2)
        PN = old_div(Nbin, abin)  # Number Surface density
        PL = old_div(Lbin, abin)  # Luminosity surface density

        # Compute the background median density both in Lum and Ngal
        # Between 4.0 - 9.0 r1Mpc
        R1 = 4.0 * self.r1Mpc * 60.0
        R2 = 9.0 * self.r1Mpc * 60.0

        if R2 >= r.max():
            R2 = r2.max()
            R1 = R2 - 2.0 * self.r1Mpc * 60.0

        PN_bgr = PN[land(rcenter > R1, rcenter < R2)]
        PL_bgr = PL[land(rcenter > R1, rcenter < R2)]
        r_bgr = rcenter[land(rcenter > R1, rcenter < R2)]

        # Get the mean values for the Ngal and Lr profiles, which will
        # be the correction per arcmin^2
        PN_mean = numpy.mean(PN_bgr)
        PL_mean = numpy.mean(PL_bgr)

        # Total number in area
        N_bgr = PN_bgr.sum()
        L_bgr = PL_bgr.sum()
        area_bgr = math.pi * (R2**2 - R1**2)

        # Get the correction for Number of galaxies and Luminosoty
        # For R200 we need to recompute R200 and N200 based on new
        # R200 value.
        area_r1Mpc = math.pi * (self.r1Mpc * 60.)**2  # in arcmin2
        self.Ngal_c = self.Ngal - PN_mean * area_r1Mpc
        if self.Ngal_c < 0:
            self.Ngal_c = 0.0
        self.R200_c = 0.156 * (self.Ngal_c**0.6) / self.h  # In Mpc
        self.r200_c = old_div((old_div(self.R200_c, self.arcmin2Mpc)), 60.0)
        area_r200_c = math.pi * (self.r200_c * 60.)**2  # in arcmin2
        area_r200 = math.pi * (self.r200 * 60.)**2  # in arcmin2
        self.i200_c = numpy.where(self.dist2BCG[ixr] <= self.r200_c)
        self.N200_c = len(self.i200_c[0]) - PN_mean * area_r200_c
        self.L200_c = Lr[self.i200_c].sum(
        ) - 0.3 * PL_mean * area_r200_c  # 0.3 factor from old code, why????

        #print self.Ngal
        #print PN
        #print r1
        #print rcenter
        #print R1,R2
        #print r.min(),r.max()
        #print "PN_mean",PN_mean
        #print PN_bgr
        #print area_r1Mpc
        #print self.Ngal_c
        #print self.r200_c
        #print self.R200_c

        # Errors for uncorrected valyes
        dL200 = self.Lr_err[self.i200].sum()
        self.d_Ngal = math.sqrt(self.Ngal)
        self.d_N200 = math.sqrt(self.N200)
        self.d_L200 = math.sqrt(dL200**2)

        # We estimate the errors
        dL200_c = self.Lr_err[self.i200_c].sum()

        self.d_Ngal_c2 = self.Ngal_c + (
            (old_div(area_r1Mpc, area_bgr))**2) * N_bgr
        self.d_N200_c2 = self.N200_c + (
            (old_div(area_r200_c, area_bgr))**2) * N_bgr
        self.d_L200_c2 = dL200_c**2 + (
            (old_div(area_r200_c, area_bgr))**2) * dL200_c**2

        # Avoid sqrt of negative number
        if self.d_Ngal_c2 < 0:
            self.d_Ngal_c = 0
        else:
            self.d_Ngal_c = math.sqrt(self.Ngal_c + ((old_div(
                area_r1Mpc, area_bgr))**2) * N_bgr)

        if self.d_N200_c2 < 0:
            self.d_N200_c = 0
        else:
            self.d_N200_c = math.sqrt(self.N200_c + ((old_div(
                area_r200_c, area_bgr))**2) * N_bgr)

        if self.d_L200_c2 < 0:
            self.d_L200_c = 0
        else:
            self.d_L200_c = math.sqrt(dL200_c**2 + ((old_div(
                area_r200_c, area_bgr))**2) * dL200_c**2)

            # Get the mass for corrected values
        (self.M_N200, self.M_L200) = Mass_calib(self.N200_c,
                                                self.L200_c,
                                                self.LBCG,
                                                zo,
                                                h=self.h)

        ####################################
        # Now plot the profiles + some info
        ####################################
        x_bg = [R1, R2]
        yN_bg = [PN_mean, PN_mean]
        yL_bg = [PL_mean, PL_mean]

        xx = [self.r1Mpc * 60., self.r1Mpc * 60.]
        rr = [self.r200_c * 60., self.r200_c * 60.]
        yy = [PN.min(), PN.max()]
        pylab.figure(1, figsize=(8, 8))
        pylab.subplot(2, 1, 1)
        pylab.plot(rcenter, PN, 'k-')
        pylab.plot(rcenter, PN, 'ko')
        pylab.plot(xx, yy, 'r-')
        pylab.plot(rr, yy, 'y-')
        pylab.plot(x_bg, yN_bg, 'g--')
        pylab.text(rcenter[0],
                   PN.min() * 1.2,
                   self.SCSname,
                   ha='left',
                   size=10)
        pylab.text(self.r1Mpc * 60.0 * 1.1,
                   PN.max() * 0.2,
                   "R1Mpc",
                   rotation='vertical',
                   ha='left',
                   size=10)
        pylab.text(self.r200_c * 60.0 * 1.1,
                   PN.max() * 0.2,
                   "R200",
                   rotation='vertical',
                   ha='left',
                   size=10)
        pylab.xlabel(r'$r {\rm (arcmin)}$', fontsize=14)
        pylab.ylabel(r'$N(r) {\rm arcmin}^{-2}$', fontsize=14)
        pylab.loglog()

        pylab.subplot(2, 1, 2)
        pylab.plot(rcenter, PL, 'k-')
        pylab.plot(rcenter, PL, 'ko')
        yy = [PL.min(), PL.max()]
        pylab.plot(xx, yy, 'r-')
        pylab.plot(rr, yy, 'y-')
        pylab.plot(x_bg, yL_bg, 'g--')
        pylab.text(rcenter[0],
                   PL.min() * 1.2,
                   self.SCSname,
                   ha='left',
                   size=10)
        pylab.text(self.r1Mpc * 60.0 * 1.1,
                   PL.max() * 0.2,
                   "R1Mpc",
                   rotation='vertical',
                   ha='left',
                   size=10)
        pylab.text(self.r200_c * 60.0 * 1.1,
                   PL.max() * 0.2,
                   "R200",
                   rotation='vertical',
                   ha='left',
                   size=10)
        pylab.xlabel(r'$r {\rm (arcmin)}$', fontsize=14)
        pylab.ylabel(r'$L(r) {\rm arcmin}^{-2}$', fontsize=14)
        outname = os.path.join(self.outpath, "%s_prof.png" % self.SCSname)
        pylab.loglog()
        try:
            pylab.savefig(outname)
            pylab.close()
        except:
            print("** ERROR: Could not write %s ***" % outname)
            pylab.close()
        return

    def write_info(self, k=0):

        header = ''
        header = header + "# %18s" % "name"
        header = header + "%14s %14s " % ('RA(deg)', 'DEC(deg)')
        #header = header + "%11s %11s " % ('RA','DEC')
        header = header + "%7s " * 2 % ('z_cl', 'err')
        header = header + "%7s " * 2 % ('zb', 'zml')
        header = header + "%7s " * 2 % ('R200_c', 'r200_c')
        header = header + "%7s " * 2 % ('Ngal_c', 'd_Ngal')
        header = header + "%7s " * 2 % ('N200_c', 'd_N200')
        header = header + "%9s " * 3 % ('L_BCG', 'L200_c', 'd_L200')
        header = header + "%9s " * 2 % ('M_N200', 'M_L200')
        header = header + "\n"

        # Corrected values
        cfile = os.path.join(self.outpath, "%s.cinfo" % self.SCSname)
        c = open(cfile, "w")
        c.write(header)
        c.write("%-20s" % self.SCSname)
        c.write("%14s %14s " % (self.ra_BCG[k], self.dec_BCG[k]))
        #c.write("%11s %11s " % (self.RA,self.DEC))
        c.write("%7.3f " * 2 % (self.zcl, self.zcl_err))
        c.write("%7.3f " * 2 % (self.zb_BCG[k], self.zml_BCG[k]))
        c.write("%7.2f " * 2 % (self.R200_c, self.r200_c * 60.0))
        c.write("%7.2f " * 2 % (self.Ngal_c, self.d_Ngal_c))
        c.write("%7.2f " * 2 % (self.N200_c, self.d_N200_c))
        c.write("%9.2e " * 3 % (self.LBCG, self.L200_c, self.d_L200_c))
        c.write("%9.2e " * 2 % (self.M_N200, self.M_L200))
        c.write("%4s " % self.zuse)
        c.write("\n")
        c.close()

        header = ''
        header = header + "# %18s" % " name"
        header = header + "%14s %14s " % ('RA(deg)', 'DEC(deg)')
        #header = header + "%11s %11s " % ('RA','DEC')
        header = header + "%7s " * 2 % ('z_cl', 'err')
        header = header + "%7s " * 2 % ('zb', 'zml')
        header = header + "%7s " * 2 % ('R200', 'r200')
        header = header + "%7s " * 2 % ('Ngal', 'd_Ngal')
        header = header + "%7s " * 2 % ('N200', 'd_N200')
        header = header + "%9s " * 3 % ('L_BCG', 'L200', 'd_L200')
        header = header + "%9s " * 2 % ('M_N200', 'M_L200')
        header = header + "\n"

        # UN-corrected values
        ofile = os.path.join(self.outpath, "%s.oinfo" % self.SCSname)
        o = open(ofile, "w")
        o.write(header)
        o.write("%-20s" % self.SCSname)
        o.write("%14s %14s " % (self.ra_BCG[k], self.dec_BCG[k]))
        #o.write("%11s %11s " % (self.RA,self.DEC))
        o.write("%7.3f " * 2 % (self.zcl, self.zcl_err))
        o.write("%7.3f " * 2 % (self.zb_BCG[k], self.zml_BCG[k]))
        o.write("%7.2f " * 2 % (self.R200, self.r200 * 60.0))
        o.write("%7.2f " * 2 % (self.Ngal, self.d_Ngal))
        o.write("%7.2f " * 2 % (self.N200, self.d_N200))
        o.write("%9.2e " * 3 % (self.LBCG, self.L200, self.d_L200))
        o.write("%9.2e " * 2 % Mass_calib(self.N200,
                                          self.L200,
                                          self.LBCG,
                                          self.z_BCG[0],
                                          h=self.h))
        o.write("%4s " % self.zuse)
        o.write("\n")
        o.close()
        return

        # #########################################################
        # # Compute Kcorr array to make linear interpolation later
        # #########################################################
        # def Kcorr_fit(sed='El_Benitez2003'):
        #     import scipy
        #     import scipy.interpolate
        #     k = {}
        #     zx = numpy.arange(0.0, 2.0, 0.005)
        #     kg = bpz_mix.Kcorr(zx,sed,'g_MOSAICII')
        #     kr = bpz_mix.Kcorr(zx,sed,'r_MOSAICII')
        #     ki = bpz_mix.Kcorr(zx,sed,'i_MOSAICII')
        #     kz = bpz_mix.Kcorr(zx,sed,'z_MOSAICII')
        #     k['g'] = scipy.interpolate.interp1d(zx,kg)
        #     k['r'] = scipy.interpolate.interp1d(zx,kr)
        #     k['i'] = scipy.interpolate.interp1d(zx,ki)
        #     k['z'] = scipy.interpolate.interp1d(zx,kz)
        #     return k


        ##################################################################
        # Read both kcorrection k(z) and evolution ev(z) from BC03 model
        ##################################################################
def KEfit(modelfile):

    import scipy
    import scipy.interpolate
    import tableio

    sout.write("# Getting K(z) and Ev(z) corrections from file:  %s\n" %
               modelfile)

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


############################################
# Read evolution ev(z) only from BC03 model
############################################
def evolfit(modelfile):

    import scipy
    import scipy.interpolate
    import tableio

    e = {}

    (z, e_g, e_r, e_i,
     e_z) = tableio.get_data(modelfile, cols=(0, 14, 15, 16, 17))

    e['g'] = scipy.interpolate.interp1d(z, e_g)
    e['r'] = scipy.interpolate.interp1d(z, e_r)
    e['i'] = scipy.interpolate.interp1d(z, e_i)
    e['z'] = scipy.interpolate.interp1d(z, e_z)

    return e


######################################
# BCG Probability function
# p = 0 for M dimmer than Mlimit
# p = 1 for M brighter than Mlimit
######################################
def p_BCG(M, Mlim, b=0.4, zp=0.5):
    x = M - Mlim
    return F_BCG(x, b, zp)


################################################################
# BCG priot aux function,
#################################################################
def F_BCG(x, b=0.4, zp=0.5):

    #print "will use zp:", zp
    #print "will use b:", b
    # Recenter at 50% (0.5) or at 68.2% (0.682)
    dx = old_div(math.log10(-math.log(zp)), b)
    u = x + dx
    phi = numpy.exp(-10**(b * u))
    return phi


#######################################################################
# Modified Schechter magnitude function from Postman et al (2001)
# uses alpha+2 rather than alpha+1 because of the extra 10^-0.4(m-m*)
# phi = (10^(-0.4(m-m*)))^(alpha+1) * exp[-10^(-0.4(m-m*))]
# PHI = phi*10^(-0.4(m-m*))
#######################################################################
def PHI(m, mstar, alpha):
    exp = numpy.exp
    a = 10**(-0.4 * (m - mstar))
    # Note (alpha+2) normally is just (alpha+1)
    phi = a**(alpha + 2) * exp(-a)
    return phi


###########################################
# Get m_star aparent mangnitude in i-band
###########################################
def mi_star(z, cosmo=(0.3, 0.7, 0.7)):

    # Set the cosmology
    c = cosmopy.set(cosmo)
    dlum = c.dlum(z)[0]

    Mb_star = -19.43 - 1.01 * z
    Mi_star = cosmology.reobs('El_Benitez2003',
                              m=Mb_star,
                              oldfilter="B_Johnson",
                              newfilter="i_MOSAICII")
    return Mi_star + 5.0 * math.log10(dlum) + 25


###################################################
# Read in the info file with candidates positions
###################################################
def read_candidates(file):

    ra = {}
    dec = {}
    sn = {}
    zo = {}
    zx = numpy.arange(0.1, 0.9, 0.1)

    for line in open(file).readlines():
        if line[0] == "#":
            continue

        vals = line.split()
        ID = vals[0]
        ra[ID] = float(vals[3])
        dec[ID] = float(vals[4])

        i = 0
        si_max = -99.
        for s in vals[5:]:
            try:
                si = float(s)
            except:
                si = 0.0

            if si > si_max:
                si_max = si
                zmax = zx[i]
            i = i + 1

        sn[ID] = si_max
        zo[ID] = zmax

    return ra, dec, zo, sn

####################################################
# Fake an ellipse using an N-sided polygon
#####################################################
import matplotlib.patches
import math
Polygon = matplotlib.patches.Polygon


def PEllipse(xxx_todo_changeme,
             xxx_todo_changeme1,
             resolution=100,
             angle=0.0,
             **kwargs):
    (xo, yo) = xxx_todo_changeme
    (A, B) = xxx_todo_changeme1
    pi = math.pi
    cos = math.cos
    sin = math.sin
    angle = -angle * pi / 180.  # hack to make it work, angle=-angle

    t = 2 * pi / resolution * numpy.arange(resolution)
    xtmp = A * numpy.cos(t)
    ytmp = B * numpy.sin(t)

    x = xtmp * cos(angle) - ytmp * sin(angle) + xo
    y = xtmp * sin(angle) + ytmp * cos(angle) + yo
    return Polygon(list(zip(x, y)), **kwargs)


##############################
# A circle as a polygon too
###############################
def PCircle(xxx_todo_changeme2, radius, resolution=100, **kwargs):
    (xo, yo) = xxx_todo_changeme2
    pi = math.pi
    cos = math.cos
    sin = math.sin
    t = 2 * pi / resolution * numpy.arange(resolution)
    xtmp = radius * numpy.cos(t)
    ytmp = radius * numpy.sin(t)
    x = xtmp + xo
    y = ytmp + yo
    return Polygon(list(zip(x, y)), **kwargs)


##########################################################
# Cuts fits file around (xo,yo)
# Uses getfits (faster) or imcopy(iraf) for bigger files
##########################################################
def cutfits(fitsin, xo, yo, dx, dy, fitsout):
    # Avoid crashing getfits as it fails for size ~ 1900 pix and above,
    # use iraf's imhead insted
    if dx > 1900 or dy > 1900:
        from .pyraf import iraf
        i1 = int(xo - old_div(dx, 2))
        i2 = int(xo + old_div(dx, 2))
        j1 = int(yo - old_div(dy, 2))
        j2 = int(yo + old_div(dy, 2))
        section = "[%s:%s,%s:%s]" % (i1, i2, j1, j2)
        iraf.imcopy("%s%s" % (fitsin, section), fitsout, verb=0)
    else:
        cmd = "getfits %s %s %s %s %s -o %s > /dev/null 2>&1 " % (
            fitsin, xo, yo, dx, dy, fitsout)
        os.system(cmd)
    return


#######################################
# make an array with power law growth
########################################
def mklogarray(x1, x2, n):

    if x1 > 0:
        i = numpy.indices((n + 1, ))[0] * 1.0
        x = x1 * (old_div(x2, x1))**(old_div(i, n))
        #dx = x1*( (x2/x1)**(i/n) - (x2/x1)**((i-1)/n))

    elif x1 == 0:
        i = numpy.indices((n, ))[0] * 1.0 + 1
        x = numpy.zeros((n + 1, ))
        #x[1:] = x2**(i/n)
        dx = (x2 + 1)**(old_div(i, n)) - (x2 + 1)**(old_div((i - 1), n))
        x[1:] = dx.cumsum()
    else:
        print("ERROR, x < 0")
        return

    return x


#################################################
# Make histogram using xbin, gives the same
# results as numpy.histogram
#################################################
def histo(x, xbin, center=None):

    n = len(xbin) - 1

    nbin = numpy.zeros(n).astype(int16)
    for i in range(n):
        if i == 0:
            nbin[i] = len(numpy.where(land(x >= xbin[i], x <= xbin[i + 1]))[0])
        else:
            nbin[i] = len(numpy.where(land(x > xbin[i], x <= xbin[i + 1]))[0])
    # Center and reduce to n-1
    if center:
        ix = numpy.indices(xbin.shape)[0]
        i1 = ix[:-1]
        i2 = ix[1:]
        dx = xbin[i2] - xbin[i1]
        xbin = xbin[:-1] + old_div(dx, 2.0)

    return nbin, xbin


################################################################
# Bin data in y(n) acoording to x(n) using bin spacing in xbin
###############################################################
def bin_data(x, y, xbin, center=None):

    n = len(xbin) - 1
    ybin = numpy.zeros(n).astype(float64)
    for i in range(n):
        if i == 0:
            idx = numpy.where(land(x >= xbin[i], x <= xbin[i + 1]))
        else:
            idx = numpy.where(land(x > xbin[i], x <= xbin[i + 1]))
        ybin[i] = y[idx].sum()
    # Center and reduce to n-1
    if center:
        ix = numpy.indices(xbin.shape)[0]
        i1 = ix[:-1]
        i2 = ix[1:]
        dx = xbin[i2] - xbin[i1]
        xbin = xbin[:-1] + old_div(dx, 2.0)

    return ybin, xbin


def get_R200(Ngal, h=0.7):
    return 0.156 * (Ngal**0.6) / h


def Mass123(N200, L200, LBCG, h=0.7):

    # Scale to convert Lum in units of 10^11*h^-2
    Lscale = old_div(1.e10, h**2)
    Mscale = old_div(1.e14, h)

    L200 = old_div(L200, Lscale)
    LBCG = old_div(LBCG, Lscale)

    M_N200 = Mscale * 1.43 * (old_div(N200, 20.))**1.20
    M_L200 = Mscale * 1.72 * (old_div(L200, 20.))**1.55
    M_LBCG = Mscale * 1.07 * (old_div(LBCG, 5.))**1.10

    return M_N200, M_L200, M_LBCG


def Mass_calib_old(N200, L200, LBCG, z, h=0.7):

    # The best fit parameters
    if z < 0.23:
        M0_N = 1.27
        M0_L = 1.81
        alphaN = 1.20
        alphaL = 1.27
        gammaN = 0.71
        gammaL = 0.40
        aN = 1.54
        aL = 7.77
        bN = 0.41
        bL = 0.67

    else:
        M0_N = 1.57
        M0_L = 1.76
        alphaN = 1.12
        alphaL = 1.30
        gammaN = 0.34
        gammaL = 0.26
        aN = 1.64
        aL = 7.92
        bN = 0.43
        bL = 0.66

    Lscale = old_div(1.e10, h**2)
    Mscale = old_div(1.e14, h)

    L200 = old_div(L200, Lscale)
    LBCG = old_div(LBCG, Lscale)

    LBCG_N = aN * N200**bN
    LBCG_L = aL * L200**bL

    M_N200 = Mscale * M0_N * (
        (old_div(N200, 20.0))**alphaN) * (old_div(LBCG, LBCG_N))**gammaN
    M_L200 = Mscale * M0_L * (
        (old_div(L200, 40.0))**alphaL) * (old_div(LBCG, LBCG_L))**gammaL

    return M_N200, M_L200


def Mass_calib(N200, L200, LBCG, z, h=0.7):

    # The best fit parameters
    if z < 0.23:
        M0_N = 1.27
        M0_L = 1.81
        alphaN = 1.20
        alphaL = 1.27
        gammaN = 0.71
        gammaL = 0.40
        aN = old_div(1.54, h**2)
        #aL     = 7.77/h**2 # bad value
        aL = old_div(0.61, h**2)
        bN = 0.41
        bL = 0.67

    else:
        M0_N = 1.57
        M0_L = 1.76
        alphaN = 1.12
        alphaL = 1.30
        gammaN = 0.34
        gammaL = 0.26
        aN = old_div(1.64, h**2)
        #aL     = 7.92/h**2 # bad value
        aL = old_div(0.58, h**2)
        bN = 0.43
        bL = 0.66

    L200 = L200 * (h**2) / 1.e10
    LBCG = LBCG * (h**2) / 1.e10

    #L200 = L200/1.e10
    #LBCG = LBCG/1.e10

    LBCG_N = aN * N200**bN
    LBCG_L = aL * L200**bL

    M_N200 = (old_div(1.e14, h)) * M0_N * (
        (old_div(N200, 20.0))**alphaN) * (old_div(LBCG, LBCG_N))**gammaN
    M_L200 = (old_div(1.e14, h)) * M0_L * (
        (old_div(L200, 40.0))**alphaL) * (old_div(LBCG, LBCG_L))**gammaL

    return M_N200, M_L200


def M_L200(L200, LBCG, z, h=0.7):

    # The best fit parameters
    if z < 0.23:
        M0_L = 1.81
        alphaL = 1.27
        gammaL = 0.40
        aL = old_div(7.77, h**2)
        bL = 0.67
    else:
        M0_L = 1.76
        alphaL = 1.30
        gammaL = 0.26
        aL = old_div(7.92, h**2)
        bL = 0.66

    L200 = L200 * (h**2) / 1.e10
    LBCG = LBCG * (h**2) / 1.e10
    M_L200 = (old_div(1.e14, h)) * M0_L * (
        (old_div(L200, 40.0))**alphaL) * (old_div(LBCG, LBCG_L))**gammaL

    return M_L200
