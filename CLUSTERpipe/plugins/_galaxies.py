import time
import sys
import os
import numpy
import math
from past.utils import old_div
from astropy.coordinates import SkyCoord

try:
    import extras
    import aux
except ImportError:
    sys.path.append(f'{os.environ["HOME"]}/Projects/'
                    'planckClusters/MOSAICpipe/pipe_utils')
    import extras
    import aux

# get the utils from the parent directory
try:
    from cluster_utils import (p_BCG, mklogarray, histo, bin_data, color)
except ImportError:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from cluster_utils import (p_BCG, mklogarray, histo, bin_data, color)

sout = sys.stderr
land = numpy.logical_and
lor = numpy.logical_or


##################################################################
# Define the sub-sample for BCGs candidates around a position
##################################################################
def get_BCG_candidates(self, Mr_limit=-22.71, p_lim=1e-4):

    t0 = time.time()
    sout.write("# Computing p_BCG probabilities...\n")

    # The Abs mag limit @ z=0.1 in the i-band
    _evo_model = self.evo_model.copy()
    _evo_model.set_normalization('r', 0.1, Mr_limit, vega=False)
    Mi_limit = _evo_model.get_absolute_mags(self.zf, filters='i', zs=0.1)

    # Evaluate the genertic mask for BCG only onece
    if not self.BCG_probs:
        # We get the limit at the z_ph of each candidate, corrected by z=0.1
        Mr_BCG_limit = (
            Mr_limit + self.ev_r -
            self.evo_model.get_ecorrects(self.zf, filters='r', zs=0.1)
        )  # + self.DM_factor
        Mi_BCG_limit = (
            Mi_limit + self.ev_i -
            self.evo_model.get_ecorrects(self.zf, filters='i', zs=0.1))

        # Evaluate the BCG Probability function, we
        # get the limit for each object
        self.p = p_BCG(self.Mr, Mr_BCG_limit)

        self.BCG_probs = True

        i_lim = 25.0
        star_lim = self.starlim
        p_lim = max(self.p) * 0.8
        sout.write("\tAvoiding BCG_prob < %.3f in BGCs\n" % p_lim)
        mask_p = numpy.where(self.p >= p_lim, 1, 0)
        mask_g = numpy.where(self.g < i_lim + 5, 1, 0)
        mask_r = numpy.where(self.r < i_lim + 2, 1, 0)
        mask_i = numpy.where(self.i < i_lim, 1, 0)
        mask_z = numpy.where(self.z < i_lim + 1, 1, 0)
        mask_t = numpy.where(self.type <= 2.0, 1, 0)

        # Avoid freakishly bright objects, 2.5 mags brighter than the
        # M_BCG_limit
        mask_br = numpy.where(self.Mr > Mr_BCG_limit - 2.5, 1, 0)
        mask_bi = numpy.where(self.Mi > Mi_BCG_limit - 2.5, 1, 0)

        # Put a more strict cut in class_star for bcg candidates
        sout.write("\tAvoiding CLASS_STAR > %s in BGCs\n" % star_lim)
        mask_star = numpy.where(self.class_star <= star_lim, 1, 0)

        # Construct the final mask now
        self.mask_BCG = (mask_t * mask_g * mask_r * mask_i * mask_z * mask_br *
                         mask_bi * mask_p * mask_star)

        self.BCG_masked = True

        # Model color only once
        #self.zx = numpy.arange(0.01, self.zlim, 0.01)
        # self.gr_model = cosmology.color_z(
        #     sed='El_Benitez2003',
        #     filter_new='g_MOSAICII',
        #     filter_old='r_MOSAICII',
        #     z=self.zx,
        #     calibration='AB')
        # self.ri_model = cosmology.color_z(
        #     sed='El_Benitez2003',
        #     filter_new='r_MOSAICII',
        #     filter_old='i_MOSAICII',
        #     z=self.zx,
        #     calibration='AB')
        # self.iz_model = cosmology.color_z(
        #     sed='El_Benitez2003',
        #     filter_new='i_MOSAICII',
        #     filter_old='z_MOSAICII',
        #     z=self.zx,
        #     calibration='AB')

        sout.write(" \tDone: %s\n" % extras.elapsed_time_str(t0))

    # Select the candidates now
    idx = numpy.where(self.mask_BCG == 1)

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

    # The r-band Luminosity of the BCGs
    self.LBCG = self.Lr[idx]

    # The distance to the candidate's position for each BCG, in arcmin
    if self.N_BCG:
        sout.write(color("\tFound %s BCG candidates\n" % self.N_BCG, 36, 1))
    else:
        sout.write(color("\tFound %s BCG candidates\n" % self.N_BCG, 31, 5))
    return


########################################################
# Modified/updated from find_clusters_ext_auto.py
# Select galaxies around ID galaxy un redshift range
########################################################
def select_members_radius(self, i, Mi_lim=-20.25, radius=500.0, zo=None):
    if zo:
        print("Will use z:%.3f for cluster" % zo)
    else:
        zo = self.z_ph[i]
        print("Will use z:%.3f for cluster" % zo)

    # Get the relevant info for ith BCG
    ra0 = self.ra[i]
    dec0 = self.dec[i]
    Mi_BCG = self.Mi[i]
    # DM = self.DM[i]
    ID_BCG = self.id[i]

    # Width of the redshift shell
    dz = self.dz

    t0 = time.time()
    sout.write("# Selecting Cluster members... Ngal, N200, R200 \n")

    # Calculate the M_star values
    Mstar = self.evo_model.get_absolute_mags(self.zf, filters='i', zs=zo)
    Mi_BCG = Mstar - 2.5 * numpy.log10(4.0)  # 4L* galaxy at z=zo

    ########
    ### 1 - Select in position around ra0,dec0

    # Define radius in degress @ zo
    R = radius  # in kpc
    r = self.cosmo.arcsec_per_kpc_proper(zo).value / 3600 * R  # in degrees

    rcore = r / 2.0

    pos0 = SkyCoord(ra0, dec0, unit='deg', frame='icrs')
    pos1 = SkyCoord(self.ra, self.dec, unit='deg', frame='icrs')
    dist = pos0.separation(pos1).value

    mask_R = numpy.where(dist <= r, 1, 0)
    mask_rcore = numpy.where(dist <= rcore, 1, 0)
    arcmin2Mpc = self.cosmo.kpc_proper_per_arcmin(
        zo).value / 1000  # scale between arcmin and Mpc

    # 2 - Select in redshift
    z1 = zo - dz
    z2 = zo + dz
    mask_z = numpy.where(land(self.z_ph >= z1, self.z_ph <= z2), 1, 0)

    # 3 - Select in brightness
    Mi_lim_zo = Mi_lim + self.evf['i'](zo) - self.evf['i'](0.1)
    mask_L1 = numpy.where(self.Mi <= Mi_lim_zo, 1, 0)  # Faint  cut > 0.4L*
    mask_L2 = numpy.where(self.Mi >= Mi_BCG, 1, 0)  # Bright cut < L_BCG

    # The final selection mask, position x redshift x Luminosity
    #idx = numpy.where(mask_R * mask_L1 * mask_L2 * mask_z == 1)[0]
    idc = numpy.where(mask_rcore * mask_L1 * mask_L2 * mask_z == 1)[0]

    # Shot versions handles
    gr = self.gr
    ri = self.ri

    # Some simple 3-sigma clipping defined using r< rcore
    Nsigma = 3.0
    loop = 1
    converge = False
    while not converge:
        # The conditions to apply
        c1 = numpy.abs(gr[idc] -
                       gr[idc].mean()) > Nsigma * numpy.std(gr[idc], ddof=1)
        c2 = numpy.abs(ri[idc] -
                       ri[idc].mean()) > Nsigma * numpy.std(ri[idc], ddof=1)
        iclip = numpy.where(lor(c1,
                                c2))[0]  # where any of the conditions fails
        if len(iclip) > 0:
            idc = numpy.delete(idc, iclip)  # Removed failed ones
            converge = False
        else:
            converge = True
        loop += 1

    #print(idc)
    #print(self.z_ph[idc])

    # Compute the weighted average and rms
    dz = 0.5 * numpy.abs(self.z2[idc] - self.z1[idc])
    # Fix zeros
    dz[dz == 0] = 1e-5
    z_cl, z_clrms = aux.statsw(self.z_ph[idc], weight=1.0 / dz)
    sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t0))

    # Or we can make a new mask where the condition's are true
    c1 = numpy.abs(self.gr -
                   gr[idc].mean()) > Nsigma * numpy.std(gr[idc], ddof=1)
    c2 = numpy.abs(self.ri -
                   ri[idc].mean()) > Nsigma * numpy.std(ri[idc], ddof=1)
    mask_cm = numpy.where(lor(c1, c2), 0, 1)  # where condition fails
    iRadius = numpy.where(mask_R * mask_L1 * mask_L2 * mask_z * mask_cm == 1)
    iRadius_all = numpy.where(mask_L1 * mask_L2 * mask_z * mask_cm == 1)
    Ngal = len(iRadius[0])
    try:
        sout.write(color(f"# Total: {Ngal} objects selected in "
                        f"{radius} [kpc] around {self.ID}\n", 36, 1))
    except AttributeError:
        sout.write(color(f"# Total: {Ngal} objects selected in "
                        f"{radius} [kpc] around {ra0:0.5f} {dec0:0.5f}\n", 36, 1))
    # Pass up
    self.iRadius = iRadius
    self.arcmin2Mpc = arcmin2Mpc
    self.dist2BCG = dist
    self.Lsum = self.Lr[iRadius].sum()
    self.Ngal = Ngal
    self.z_cl = z_cl
    self.z_clerr = z_clrms
    self.rdeg = r  # in degress
    self.r1Mpc = r  # naming fix for background estimates
    self.idc = idc  # galaxies used for mean redshift
    self.ID_BCG = ID_BCG

    # Sort indices radially for galaxies < N*R1Mpc, will be used later
    i = numpy.argsort(self.dist2BCG[iRadius_all])
    self.ix_radial = iRadius_all[0][i]

    return z_cl, z_clrms


##########################################
# Compute the Background for the clusters
##########################################
def background(self):
    ixr = self.ix_radial

    # No back substraction
    if self.Ngal <= 2:
        self.Ngal_c = self.Ngal
        print(
            color('Background -- Not enough galaxies found in cluster', 31, 5))
        return

    # Store radially ordered
    r = self.dist2BCG[ixr] * 60.0  # in arcmin
    Lr = self.Lr[ixr]  # We do in the r-band as Reyes et al

    # Bin the Ngal/Lum data in log spacing
    n = 10
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

    # Compute the background median density both in Lum and Ngal
    # Between 4.0 - 9.0 r1Mpc
    R1 = 4.0 * self.r1Mpc * 60.0
    R2 = 9.0 * self.r1Mpc * 60.0
    print("# Estimating Background between R1,R2 %.2f--%2.f[arcmin]" %
          (R1, R2))

    if R2 >= r.max():
        print(color('\tBackground R2 > image limits! -- recomputing', 31, 0))
        R2 = r2.max()
        R1 = R2 - 2.0 * self.r1Mpc * 60.0
    print("# Estimating Background between R1,R2 %.2f--%2.f[arcmin]" %
          (R1, R2))

    PN_bgr = PN[land(rcenter > R1, rcenter < R2)]

    # Get the mean values for the Ngal and Lr profiles, which will
    # be the correction per arcmin^2
    PN_mean = numpy.mean(PN_bgr)
    print('\tmean number of BG galaxies -- {}'.format(PN_mean))

    # Total number in area
    N_bgr = PN_bgr.sum()
    area_bgr = math.pi * (R2**2 - R1**2)

    # Get the correction for Number of galaxies and Luminosoty
    # For R200 we need to recompute R200 and N200 based on new
    # R200 value.
    area_r1Mpc = math.pi * (self.r1Mpc * 60.)**2  # in arcmin2
    self.Ngal_c = self.Ngal - PN_mean * area_r1Mpc
    if self.Ngal_c < 0:
        self.Ngal_c = 0.0

    print('---- test stuff -----')

    print(self.iclose)
    print(self.x_image[self.iclose], self.y_image[self.iclose])

    # print(self.Ngal)
    # print(PN)
    # print(r1)
    # print(rcenter)
    # print(R1,R2)
    # print(r.min(),r.max())
    # print("PN_mean",PN_mean)
    # print(PN_bgr)
    # print(area_r1Mpc)
    # print("Ngal ", self.Ngal)
    # print("Ngal_c", self.Ngal_c)
    #print("r200_c",self.r200_c)
    #print("R200_c",self.R200_c)

    self.d_Ngal_c2 = self.Ngal_c + ((old_div(area_r1Mpc, area_bgr))**2) * N_bgr

    # Avoid sqrt of negative number
    if self.d_Ngal_c2 < 0:
        self.d_Ngal_c = 0
    else:
        self.d_Ngal_c = math.sqrt(self.Ngal_c +
                                  ((old_div(area_r1Mpc, area_bgr))**2) * N_bgr)

    return


##########################################
# Compute the Background for the clusters
##########################################
def background_map(self):
    ixr = self.ix_radial

    # No back substraction
    if self.Ngal <= 2:
        self.Ngal_c = self.Ngal
        print(
            color('Background -- Not enough galaxies found in cluster', 31, 5))
        return

    # Store radially ordered
    r = self.dist2BCG[ixr] * 60.0  # in arcmin
    Lr = self.Lr[ixr]  # We do in the r-band as Reyes et al

    # Bin the Ngal/Lum data in log spacing
    n = 10
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

    # Compute the background median density both in Lum and Ngal
    # Here's is where we are going to make a couple of maps to compute the areas
    # for the background
    R1 = 3.0 * self.r1Mpc * 60.0
    R2 = r.max()  # go all the way out
    print("# Estimating Background @ r > 3mpc -- %.2f - %.2f [arcmin]" %
          (R1, R2))

    PN_bgr = PN[rcenter > R1]

    # Get the mean values for the Ngal and Lr profiles, which will
    # be the correction per arcmin^2
    PN_mean = numpy.mean(PN_bgr)
    print('\tmean number of BG galaxies -- {}'.format(PN_mean))

    # Total number in background area
    N_bgr = PN_bgr.sum()

    # get the area of the background. We'll make a 'blank' image the same size as
    # our input image and then sum over the pixels that are either in or out of
    # the cluster region.
    # cluster location
    a, b = round(self.x_image[self.iclose]), round(self.y_image[self.iclose])
    # size of the image
    n = self.jpg_array.shape[0]
    # cluster radius in arcseconds converted to pixels.
    r = R1 * 60 / self.pixscale

    # create pixel grid
    y, x = numpy.ogrid[-a:n - a, -b:n - b]
    # mask the cluster region
    mask = x * x + y * y <= r * r
    # create new 'bool' image
    img_array = numpy.ones((n, n), dtype='bool')
    # the cluster region becomes 'false' or zero
    img_array[mask] = False

    # sum the background region gives the number of pixels. Multiply by the pixel
    # scale to get the total area. Convert to arcminutes.
    area_bgr = img_array.sum() * self.pixscale / 60

    # Get the correction for Number of galaxies and Luminosoty
    # For R200 we need to recompute R200 and N200 based on new
    # R200 value.
    area_r1Mpc = math.pi * (self.r1Mpc * 60.)**2  # in arcmin2
    # use the inverse of the cluster mask to find the cluster area
    area_r1mpc = (n**2 - img_array.sum()) * self.pixscale / 60

    self.Ngal_c = self.Ngal - PN_mean * area_r1Mpc
    if self.Ngal_c < 0:
        self.Ngal_c = 0.0

    self.d_Ngal_c2 = self.Ngal_c + ((old_div(area_r1Mpc, area_bgr))**2) * N_bgr

    # Avoid sqrt of negative number
    if self.d_Ngal_c2 < 0:
        self.d_Ngal_c = 0
    else:
        self.d_Ngal_c = math.sqrt(self.Ngal_c +
                                  ((old_div(area_r1Mpc, area_bgr))**2) * N_bgr)

    return
