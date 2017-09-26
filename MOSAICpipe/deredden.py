#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from builtins import str
from past.utils import old_div
import os
import sys
from astropy.coordinates import SkyCoord
import numpy as np
from math import log10
import tableio


def compute_AEBV(filter='r_SDSS', sed='flat'):

    import bpz_tools
    """
    Return A(lambda)/E(B-V) for a given filter response and SED. Any
    SED can be used or a flat SED which is basically no SED
    convolved.
    Returns A(lambda)/A(V) and A(lambda)/E(B-v)
    """

    Vfilter = 'V_Bessell'
    Bfilter = 'B_Bessell'
    AV = 0.1

    # Get A(filter)
    if sed == 'flat':
        (xl, yf) = bpz_tools.get_filter(filter)
        ysed = xl * 0.0 + 1.0
    else:
        (xl, ysed, yf) = bpz_tools.get_sednfilter(sed, filter)
    xl = np.asarray(xl)
    yf = np.asarray(yf)
    ysed = np.array(ysed)

    Ax = AAV_ccm(xl)
    A_lambda = flux(xl, Ax, yf, units='f_l')
    mx_0 = -2.5 * log10(flux(xl, ysed, yf))  # Normal
    mx_1 = -2.5 * log10(flux(xl, ysed * 10**(-0.4 * Ax * AV), yf))  # Redder
    A_AV = old_div((mx_1 - mx_0), AV)

    # Get A(V)
    if sed == 'flat':
        (xl, yf) = bpz_tools.get_filter(Vfilter)
        ysed = xl * 0.0 + 1.0
    else:
        (xl, ysed, yf) = bpz_tools.get_sednfilter(sed, Vfilter)
    xl = np.asarray(xl)
    yf = np.asarray(yf)
    ysed = np.array(ysed) * 0.0 + 1.0
    Ax = AAV_ccm(xl)
    A_lambda = flux(xl, Ax, yf, units='f_l')
    mV_0 = -2.5 * log10(flux(xl, ysed, yf))  # Normal
    mV_1 = -2.5 * log10(flux(xl, ysed * 10**(-0.4 * Ax * AV), yf))  # Redder

    # Get A(B)
    if sed == 'flat':
        (xl, yf) = bpz_tools.get_filter(Bfilter)
        ysed = xl * 0.0 + 1.0
    else:
        (xl, ysed, yf) = bpz_tools.get_sednfilter(sed, Bfilter)
    xl = np.asarray(xl)
    yf = np.asarray(yf)
    ysed = np.array(ysed) * 0.0 + 1.0
    Ax = AAV_ccm(xl)
    A_lambda = flux(xl, Ax, yf, units='f_l')
    mB_0 = -2.5 * log10(flux(xl, ysed, yf))  # Normal
    mB_1 = -2.5 * log10(flux(xl, ysed * 10**(-0.4 * Ax * AV), yf))  # Redder

    # Compute A/E(B-V)
    A_EBV = old_div((mx_1 - mx_0), ((mB_1 - mV_1) - (mB_0 - mV_0)))

    return A_AV, A_EBV, old_div(A_EBV, A_AV)


def flux(xsr, ys, yr, ccd='yes', units='nu'):

    from scipy.integrate import trapz
    from math import sqrt
    """ Flux of spectrum ys observed through response yr,
        both defined on xsr
    Both f_nu and f_lambda have to be defined over lambda
    If units=nu, it gives f_nu as the output
    """
    clight_AHz = 2.99792458e18
    if ccd == 'yes': yr = yr * xsr
    norm = trapz(yr, xsr)
    f_l = old_div(trapz(ys * yr, xsr), norm)
    if units == 'nu':
        # Pivotal Wavelength
        lp = sqrt(old_div(norm, trapz(yr / xsr / xsr, xsr)))
        return f_l * lp**2 / clight_AHz
    else:
        return f_l


def getEBV_old(ra, dec):
    """ function recieves a coordinate pair, RA and DEC from the
        caller, converts to galactic coords and runs the dust_getval
        code, installed in the path. Returns an extinction correction
        in magnitudes and an error object (a list of strings) of
        possible reported problems with the region of the sky.
    """

    # convert ra and dec to l,b using astutil.
    # ra should be in hours (freaking iraf).
    # we emulate pipes to Stdin and Stdout so we don't write no stinking files
    #raanddec = [str(ra/15) +" "+str(dec)+" 2000"]
    raanddec = [str(ra) + " " + str(dec) + " 2000"]
    conversion = astutil.galactic(Stdin=raanddec, print_c="no", Stdout=1)
    # conversion is a list of strings, this has only one element:
    # eg. ['     227.5430   46.1912']
    # which is l and b
    gall = conversion[0].split()[0]
    galb = conversion[0].split()[1]

    # ok, we have l and b. now onto the extinction stuff.
    # build the dust_val command line
    cmd = "dust_getval " + gall + " " + galb + " interp=y verbose=n"
    output = _calc_ebv(cmd)
    # output is a list of strings, only one element in this case. looks like
    # [' 227.543  46.191      0.03452\n']
    # dust_getval returns the original coords and
    # the extinction correction in mags
    # which is the last piece of that string
    eBV = output[0].split()[2]
    # next run dust_getval with the mask option to look for anomolies in the maps
    cmd = "dust_getval " + gall + " " + galb + " map=mask verbose=n"
    mask = _calc_ebv(cmd)
    # quality is a string of data quality flags returned by dust_getval when
    # map=mask.
    # looks like
    # ' 227.543  46.191  3hcons OK      OK      OK      OK      OK      OK     \n'
    quality = mask[1]
    # return this with the extinction.
    return eBV, quality


def get_EBV(ra, dec):
    """ function recieves a coordinate pair, RA and DEC from the
        caller, converts to galactic coords and runs the dust_getval
        code, installed in the path. Returns an extinction correction
        in magnitudes and an error object (a list of strings) of
        possible reported problems with the region of the sky.
    """

    coords_eq = "/tmp/coords_radec.dat"
    coords_lb = "/tmp/coords_galac.dat"
    eBVdata = "/tmp/eBV.dat"

    # Check that they are arrays
    if isinstance(ra, float) or isinstance(dec, float):
        ra = np.asarray(ra)
        dec = np.asarray(dec)

    if len(ra) != len(dec):
        print("ERROR: RA,DEC must have same dimensions")
        return

    # convert ra and dec to l,b using astutil.  ra should be in hours
    # Write out a coords.dat file
    tableio.put_data(coords_eq, (ra, dec), format="%10.6f %10.6f 2000")
    # convert the ra/dec to l/b galactic coordinates
    c = SkyCoord(ra, dec, frame='icrs', unit='deg')
    tableio.put_data(coords_lb, (c.galactic.l.value, c.galactic.b.value),
                        format="%10.6f %10.6f 2000")

    #astutil.galactic(coords_eq, print_c="no", Stdout=coords_lb)

    # ok, we have l and b. now onto the extinction stuff. build the dust_val
    # command line
    cmd = "dust_getval infile=%s outfile=%s verbose=n interp=y" % (coords_lb,
                                                                   eBVdata)
    #_calc_ebv(cmd)
    os.system(cmd)
    # dust_getval returns the original coords and the extinction correction in
    # mags
    # [' 227.543  46.191      0.03452\n']
    # Get the array of eBV values for each coordinate position
    eBV = tableio.get_data(eBVdata, cols=(2, ))
    # Remove the temporary files
    os.system("rm %s %s %s" % (coords_eq, coords_lb, eBVdata))
    return eBV


def filterFactor(filter):
    """caller passes an ACS filter of the form "DET_FILTER" and function
    returns the extinction correction factor, a float, for that filter.  Now
    this function defines a dictionary of extinction correction factors
    directly, but this also exists as a file in $PIPELINE/maps.  It is
    anticipated that these values will not change often, if at all, hence, the
    dictionary is defined here rather than created on the fly from the file,
    but that could be changed if it is anticipated that these numbers might
    change a lot.

    eg.
    >>>filterFactor("HRC_F555W")
    '3.24695724147'
    """

    ffactors = {
        "g_MOSAICII": 3.88489537829,
        "r_MOSAICII": 2.78438802442,
        "i_MOSAICII": 2.06519949822,
        "z_MOSAICII": 1.39714057191,
        "g": 3.88489537829,
        "r": 2.78438802442,
        "i": 2.06519949822,
        "z": 1.39714057191,
        "HST_ACS_WFC_F435W": 4.11697363315,  # Values taken from APSIS
        "HST_ACS_WFC_F555W": 3.24230342654,
        "HST_ACS_WFC_F606W": 2.928500386,
        "HST_ACS_WFC_F814W": 1.84740717341,
        "HST_ACS_WFC_F475W": 3.74714182372,
        "HST_ACS_WFC_F625W": 2.67121669327,
        "HST_ACS_WFC_F775W": 2.01774028108,
        "HST_ACS_WFC_F850LP": 1.47335876958,
        "HST_ACS_WFC_F502N": 3.52366637215,
        "HST_ACS_WFC_F892N": 1.51713294198,
        "HST_ACS_WFC_F658N": 2.52193964747,
        "HST_ACS_WFC_F550M": 3.05672088958,
        "HST_ACS_HRC_F435W": 4.11227228078,
        "HST_ACS_HRC_F555W": 3.24695724147,
        "HST_ACS_HRC_F606W": 2.94741773243,
        "HST_ACS_HRC_F814W": 1.82249557542,
        "HST_ACS_HRC_F475W": 3.724544959,
        "HST_ACS_HRC_F625W": 2.67859346295,
        "HST_ACS_HRC_F775W": 2.02818977096,
        "HST_ACS_HRC_F850LP": 1.44407689634,
        "HST_ACS_HRC_F344N": 5.10086305785,
        "HST_ACS_HRC_F502N": 3.52410034736,
        "HST_ACS_HRC_F892N": 1.5170227107,
        "HST_ACS_HRC_F658N": 2.52245685895,
        "HST_ACS_HRC_F220W": 8.81083859281,
        "HST_ACS_HRC_F250W": 6.52297815722,
        "HST_ACS_HRC_F330W": 5.17376227866,
        "HST_ACS_HRC_F550M": 3.05823161415,
        # Geneated runnning:
        #deredden.compute_AEBV(filter='R_SPECIAL_FORS2',sed='flat')[1]
        #deredden.compute_AEBV(filter='I_BESS_FORS2',sed='flat')[1]
        #deredden.compute_AEBV(filter='z_GUNN_FORS2',sed='flat')[1]
        "R_SPECIAL": 2.6355966034891507,
        "I_BESS": 1.9785822369228967,
        "z_GUNN": 1.2740325163787274,
        # To handle the fact that there is no coverage for the IRAC CH1 and CH2
        # band
        "CH1": 0.0,
        "CH2": 0.0,
    }

    try:
        return ffactors[filter]
    except KeyError:
        return compute_AEBV(filter='K_KittPeak',sed='flat')[1]

#############################################################################################


def _calc_ebv(cmd):
    sproc = popen2.Popen3(cmd, 1)
    output = sproc.fromchild.readlines()
    errs = sproc.childerr.readlines()
    return output

# Compute A(lambda)/A(V) using the prescription from Cardelli, Clayton
# & Mathis (1998) uptaed in the optical-NIR using O'Donnell (1994).
def AAV_ccm(wavelength, rv=3.1):

    land = np.logical_and

    # Convert to inverse microns
    x = old_div(10000.0, wavelength)
    a = x * 0.0
    b = x * 0.0

    # Compute a(x) and b(x) for all cases
    # Case 1: x < 0.3
    ix = np.where(x < 0.3)
    Nsel = len(ix[0])
    if Nsel > 0:
        sys.exit("Wavelength out of range of extinction function")

    # Case 2: Infrared 0.3< x<1.1
    ix = np.where(land(x > 0.3, x <= 1.1))
    Nsel = len(ix[0])
    if Nsel > 0:
        y = x[ix]
        a[ix] = 0.574 * y**1.61
        b[ix] = -0.527 * y**1.61

        # Case 3; Optical/NIR  1.1 < x <= 3.3
    ix = np.where(land(x > 1.1, x <= 3.3))
    Nsel = len(ix[0])
    if Nsel > 0:
        y = x[ix] - 1.82
        # Carelli fit
        #a[ix] = 1.0 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 +
        #0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
        #b[ix] = 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 -
        #0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7
        # O'Donnell fit
        a[ix] = (1.0 + 0.104 * y - 0.609 * y**2 + 0.701 * y**3 + 1.137 *
                 y**4 - 1.718 * y**5 - 0.827 * y**6 + 1.647 * y**7 - 0.505 *
                 y**8)
        b[ix] = (1.952 * y + 2.908 * y**2 - 3.989 * y**3 - 7.985 * y**4 +
                 11.102 * y**5 + 5.491 * y**6 - 10.805 * y**7 + 3.347 * y**8)

        # Case 4: Mid-UV  3.3 < x <= 5.9
    ix = np.where(land(x > 3.3, x <= 5.9))
    Nsel = len(ix[0])
    if Nsel > 0:
        y = (x[ix] - 4.67)**2
        a[ix] = 1.752 - 0.316 * x[ix] - old_div(0.104, (y + 0.341))
        b[ix] = -3.090 + 1.825 * x[ix] + old_div(1.206, (y + 0.263))

        # Case 4: 5.9 < x < 8.0
    ix = np.where(land(x > 5.9, x <= 8.0))
    Nsel = len(ix[0])
    if Nsel > 0:
        y = (x[ix] - 4.67)**2
        a[ix] = 1.752 - 0.316 * x[ix] - old_div(0.104, (y + 0.341))
        b[ix] = -3.090 + 1.825 * x[ix] + old_div(1.206, (y + 0.263))
        y = x[ix] - 5.9
        a[ix] = a[ix] - 0.04473 * y**2 - 0.009779 * y**3
        b[ix] = b[ix] + 0.21300 * y**2 + 0.120700 * y**3

    # Case 5: 8 < x < 11 ; Far-UV
    ix = np.where(land(x > 8.0, x <= 11.0))
    Nsel = len(ix[0])
    if Nsel > 0:
        y = x[ix] - 8.0
        a[ix] = -1.072 - 0.628 * y + 0.137 * y**2 - 0.070 * y**3
        b[ix] = 13.670 + 4.257 * y - 0.420 * y**2 + 0.374 * y**3

        # Compute A(lambda)/A(V)
    AAV = a + old_div(b, rv)
    return AAV

# Compute values for the MOSAICII filters
def MOSAICII():
    filters = ('g_MOSAICII', 'r_MOSAICII', 'i_MOSAICII', 'z_MOSAICII')
    for filter in filters:
        print(filter, compute_AEBV(filter, 'flat')[1])
