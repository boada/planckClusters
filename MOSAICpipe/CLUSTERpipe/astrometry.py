from __future__ import print_function
from __future__ import division
# Based on Iraf's xy2rd task
# transforms (x,y) --> (RA,DEC)
from builtins import map
from builtins import range
from past.utils import old_div
import sys
import os

def xy2rd(x, y, fitsfile, units='hms', wrap=True):
    import numpy
    from pyfits import getheader
    from math import pi

    atan2 = numpy.arctan2
    sqrt = numpy.sqrt
    cos = numpy.cos
    sin = numpy.sin

    # Read the WCS information from the fitsfile
    header = getheader(fitsfile)

    CRPIX1 = header["CRPIX1"]
    CRPIX2 = header["CRPIX2"]
    CRVAL1 = header["CRVAL1"]
    CRVAL2 = header["CRVAL2"]
    CD1_1 = header["CD1_1"]
    CD1_2 = header["CD1_2"]
    CD2_1 = header["CD2_1"]
    CD2_2 = header["CD2_2"]

    # translate (x,y) to (ra, dec)
    xi = CD1_1 * (x - CRPIX1) + CD1_2 * (y - CRPIX2)
    eta = CD2_1 * (x - CRPIX1) + CD2_2 * (y - CRPIX2)

    # Transform into rads
    d2r = old_div(pi, 180.)
    xi = xi * d2r
    eta = eta * d2r
    ra0 = CRVAL1 * d2r
    dec0 = CRVAL2 * d2r

    ra = atan2(xi, cos(dec0) - eta * sin(dec0)) + ra0
    dec = atan2(eta * cos(dec0) + sin(dec0),
                sqrt((cos(dec0) - eta * sin(dec0))**2 + xi**2))

    ra = old_div(ra, d2r)
    dec = old_div(dec, d2r)
    if wrap:  # Wrap arund 360.0
        ra = numpy.where(ra < 0, ra + 360.0, ra)

    if units == 'hours':
        ra = dec2deg(old_div(ra, 15.))
        dec = dec2deg(dec)

    elif units == 'degrees':
        ra = dec2deg(ra)
        dec = dec2deg(dec)

    return ra, dec


def xy2rd_sdss(x, y, fitsfile, units='hms'):
    import numpy
    from pyfits import getheader
    from math import pi

    atan2 = numpy.arctan2
    sqrt = numpy.sqrt
    cos = numpy.cos
    sin = numpy.sin

    # Read the WCS information from the fitsfile
    header = getheader(fitsfile)

    CRPIX1 = header["CRPIX1"]
    CRPIX2 = header["CRPIX2"]
    CRVAL1 = header["CRVAL1"]
    CRVAL2 = header["CRVAL2"]
    CD1_1 = header["CD1_1"]
    CD1_2 = header["CD1_2"]
    CD2_1 = header["CD2_1"]
    CD2_2 = header["CD2_2"]

    # translate (x,y) to (ra, dec)
    xi = CD1_1 * (x - CRPIX1) + CD1_2 * (y - CRPIX2)
    eta = CD2_1 * (x - CRPIX1) + CD2_2 * (y - CRPIX2)

    # Transform into rads
    d2r = old_div(pi, 180.)
    xi = xi * d2r
    eta = eta * d2r
    ra0 = CRVAL1 * d2r
    dec0 = CRVAL2 * d2r

    ra = atan2(xi, cos(dec0) - eta * sin(dec0)) + ra0
    dec = atan2(eta * cos(dec0) + sin(dec0),
                sqrt((cos(dec0) - eta * sin(dec0))**2 + xi**2))

    ra = old_div(ra, d2r)
    dec = old_div(dec, d2r)
    #ra = numpy.where(ra<0,ra+360.0,ra)

    if units == 'hours':
        ra = dec2deg(old_div(ra, 15.))
        dec = dec2deg(dec)

    elif units == 'degrees':
        ra = dec2deg(ra)
        dec = dec2deg(dec)

    return dec, ra


# Based on Iraf's rd2xy task
# transforms (RA,DEC) --> (x,y)
def rd2xy(ra, dec, fitsfile, units='degrees'):
    import numpy
    from pyfits import getheader
    from math import pi

    atan2 = numpy.arctan2
    sqrt = numpy.sqrt
    cos = numpy.cos
    sin = numpy.sin

    # Read the WCS information from the fitsfile
    header = getheader(fitsfile)

    CRPIX1 = header["CRPIX1"]
    CRPIX2 = header["CRPIX2"]
    CRVAL1 = header["CRVAL1"]
    CRVAL2 = header["CRVAL2"]
    CD1_1 = header["CD1_1"]
    CD2_2 = header["CD2_2"]

    try:
        CD1_2 = header["CD1_2"]
    except:
        CD1_2 = 0.0

    try:
        CD2_1 = header["CD2_1"]
    except:
        CD2_1 = 0.0

    # invert the CD Matrix
    det = CD1_1 * CD2_2 - CD1_2 * CD2_1
    if det == 0.0:
        sys.exit("ERROR:singular CD matrix")

    ICD1_1 = old_div(CD2_2, det)
    ICD1_2 = old_div(-CD1_2, det)
    ICD2_1 = old_div(-CD2_1, det)
    ICD2_2 = old_div(CD1_1, det)

    # Transform into radians
    if units == 'hours':
        ra = ra * 15.0
    d2r = old_div(pi, 180.)
    ra0 = CRVAL1 * d2r
    dec0 = CRVAL2 * d2r
    ra = ra * d2r
    dec = dec * d2r

    bottom = sin(dec) * sin(dec0) + cos(dec) * cos(dec0) * cos(ra - ra0)
    #if bottom == 0.0:
    #    sys.exit("Unreasonable RA/Dec range")

    # translate (ra, dec) to (x, y)
    xi = cos(dec) * sin(ra - ra0) / bottom
    eta = old_div(
        (sin(dec) * cos(dec0) - cos(dec) * sin(dec0) * cos(ra - ra0)), bottom)
    xi = old_div(xi, d2r)
    eta = old_div(eta, d2r)

    x = ICD1_1 * xi + ICD1_2 * eta + CRPIX1
    y = ICD2_1 * xi + ICD2_2 * eta + CRPIX2

    return x, y


# Based on Iraf's rd2xy task
# transforms (RA,DEC) --> (x,y)
def rd2xy_sdss(dec, ra, fitsfile, units='degrees'):
    import numpy
    from pyfits import getheader
    from math import pi

    atan2 = numpy.arctan2
    sqrt = numpy.sqrt
    cos = numpy.cos
    sin = numpy.sin

    # Read the WCS information from the fitsfile
    header = getheader(fitsfile)

    CRPIX1 = header["CRPIX1"]
    CRPIX2 = header["CRPIX2"]
    CRVAL1 = header["CRVAL1"]
    CRVAL2 = header["CRVAL2"]
    CD1_1 = header["CD1_1"]
    CD2_2 = header["CD2_2"]

    try:
        CD1_2 = header["CD1_2"]
    except:
        CD1_2 = 0.0

    try:
        CD2_1 = header["CD2_1"]
    except:
        CD2_1 = 0.0

    # invert the CD Matrix
    det = CD1_1 * CD2_2 - CD1_2 * CD2_1
    if det == 0.0:
        sys.exit("ERROR:singular CD matrix")

    ICD1_1 = old_div(CD2_2, det)
    ICD1_2 = old_div(-CD1_2, det)
    ICD2_1 = old_div(-CD2_1, det)
    ICD2_2 = old_div(CD1_1, det)

    # Transform into radians
    if units == 'hours':
        ra = ra * 15.0
    d2r = old_div(pi, 180.)
    ra0 = CRVAL1 * d2r
    dec0 = CRVAL2 * d2r
    ra = ra * d2r
    dec = dec * d2r

    bottom = sin(dec) * sin(dec0) + cos(dec) * cos(dec0) * cos(ra - ra0)
    #if bottom == 0.0:
    #    sys.exit("Unreasonable RA/Dec range")

    # translate (ra, dec) to (x, y)
    xi = cos(dec) * sin(ra - ra0) / bottom
    eta = old_div(
        (sin(dec) * cos(dec0) - cos(dec) * sin(dec0) * cos(ra - ra0)), bottom)
    xi = old_div(xi, d2r)
    eta = old_div(eta, d2r)

    x = ICD1_1 * xi + ICD1_2 * eta + CRPIX1
    y = ICD2_1 * xi + ICD2_2 * eta + CRPIX2

    return x, y


# Calculates great-circle distances between the two points that is,
# the shortest distance over the earth's surface using the Haversine
# formula, see http://www.movable-type.co.uk/scripts/latlong.html
def circle_distance(ra1, dec1, ra2, dec2, units='deg'):

    import numpy
    from math import pi

    cos = numpy.cos
    sin = numpy.sin
    acos = numpy.arccos
    asin = numpy.arcsin

    if units == 'deg':
        ra1 = ra1 * pi / 180.
        ra2 = ra2 * pi / 180.
        dec1 = dec1 * pi / 180.
        dec2 = dec2 * pi / 180.

    x = sin(dec1) * sin(dec2) + cos(dec1) * cos(dec2) * cos(ra2 - ra1)
    x = numpy.where(x > 1.0, 1, x)  # Avoid x>1.0 values

    d = acos(x)
    #d = math.acos(math.sin(dec1)*math.sin(dec2) + math.cos(dec1)*math.cos(dec2) * math.cos(ra2-ra1))

    if units == 'deg':
        d = d * 180.0 / pi

    return d


def deg2dec(deg, sep=":"):
    vals = deg.split(sep)

    dd = float(vals[0])
    mm = old_div(float(vals[1]), 60.)
    ss = old_div(float(vals[2]), 3600.)

    if dd < 0:
        mm = -mm
        ss = -ss

    if vals[0] == '-00':
        mm = -mm
        ss = -ss

    return dd + mm + ss


# From decimal to degress, array or scalar
def dec2deg(dec, short=None, sep=":"):

    import numpy
    import sys

    dec = numpy.asarray(dec)
    # Keep the sign for later
    sig = numpy.where(dec < 0, -1, 1)

    dd = dec.astype("Int32")
    mm = (abs(dec - dd) * 60).astype("Int32")
    ss = (abs(dec - dd) * 60 - mm) * 60

    # If not a scalar
    if len(dec.shape) != 0:

        x = numpy.concatenate((sig, dd, mm, ss))
        ids = numpy.where(abs(ss - 60.) <= 1e-3)
        #print "IDS",ids
        ss[ids] = 0.0
        mm[ids] = mm[ids] + 1

        #n = len(dec)
        #x =  x.resize(3,len(dec)) # old numarray
        x = numpy.resize(x, (4, len(dec)))
        #x.swapaxes(0,1) # old numarray
        x = numpy.swapaxes(x, 0, 1)
        if short:
            return list(map(format_deg_short, x))
        else:
            return list(map(format_deg_long, x))

    else:
        if float(abs(ss - 60.)) < 1e-3:
            ss = 0.0
            mm = mm + 1
        return format_deg((sig, dd, mm, ss), short, sep)

    return


# From decimal to degress, array or scalar
def dec2deg_short(dec):

    import numpy

    dec = numpy.asarray(dec)

    dd = dec.astype("Int32")
    mm = (abs(dec - dd) * 60).astype("Int32")
    ss = (abs(dec - dd) * 60 - mm) * 60
    x = numpy.concatenate((dd, mm, ss))

    # If not a scalar
    if len(dec.shape) != 0:

        ids = numpy.where(abs(ss - 60.) <= 1e-3)
        print("IDS", ids)
        ss[ids] = 0.0
        mm[ids] = mm[ids] + 1

        #n = len(dec)
        x = x.resize(3, len(dec))
        x.swapaxes(0, 1)
        return list(map(format_deg_short, x))

    else:
        if float(abs(ss - 60.)) < 1e-3:
            ss = 0.0
            mm = mm + 1
        return format_deg_short((dd, mm, ss))

    return


# From decimal to degrees, only scalar
def dec2deg_simple(dec):

    dd = int(dec)
    mm = int(abs(dec - dd) * 60.)
    ss = (abs(dec - dd) * 60 - mm) * 60

    if abs(ss - 60) < 1e-5:
        ss = 0.0
        mm = mm + 1

    return dd, mm, ss


def format_deg(x, short, sep=":"):

    if x[0] < 0:
        sig = "-"
    else:
        sig = ""

    f1 = "%2d"

    if abs(float(x[1])) < 10:
        f1 = "0%1d"
    else:
        f1 = "%2d"

    if float(x[2]) < 10:
        f2 = "0%1d"
    else:
        f2 = "%2d"

    if float(x[3]) < 9.9:
        f3 = "0%.2f"
    else:
        f3 = "%.2f"

    if short == 'ra':
        format = sig + f1 + sep + f2 + ".%1d"
        return format % (abs(x[1]), x[2], int(old_div(x[3], 6)))

    if short:
        format = sig + f1 + sep + f2
        return format % (abs(x[1]), x[2])

    format = sig + f1 + sep + f2 + sep + f3
    return (format % (abs(x[1]), x[2], x[3]))[:-1]


# doesn't work for dec -00:01:23 ....
def format_deg_old(x, short, sep=":"):

    f1 = "%2d"

    if abs(float(x[0])) < 10:
        f1 = "0%1d"
    else:
        f1 = "%2d"

    if float(x[1]) < 10:
        f2 = "0%1d"
    else:
        f2 = "%2d"

    if float(x[2]) < 9.9:
        f3 = "0%.1f"
    else:
        f3 = "%.1f"

    if short:
        format = f1 + sep + f2
        return format % (x[0], x[1])

    format = f1 + sep + f2 + sep + f3
    return format % (x[0], x[1], x[2])


def format_deg_long(x):
    #return format_deg(x,short=None,sep=":")
    return format_deg(x, short=None, sep=":")


def format_deg_short(x):

    if x[0] < 0:
        sig = "-"
    else:
        sig = "+"

    if abs(float(x[1])) < 10:
        f1 = "0%1d"
    else:
        f1 = "%2d"

    if abs(float(x[2])) < 10:
        f2 = "0%1d"
    else:
        f2 = "%2d"

    if float(x[3]) < 10:
        f3 = "0%.1f"
    else:
        f3 = "%.1f"

    format = sig + f1 + ":" + f2
    return format % (abs(x[1]), x[2])


def arc2kpc(z, theta, cosmo):

    from math import pi

    # From sterad to arcsec to Mpc to kpc
    # We need to go from rad/Mpc
    # 1 arcsec is (pi/180) * 1/3600. rad x 1000 kpc
    # Da is in rad/Mpc
    scale = (old_div(pi, 180)) * 1 / 3600. * 1000
    Da = get_Da(z, cosmo)
    kpc = Da * theta * scale
    return kpc


def kpc2arc(z, kpc, cosmo):

    from math import pi
    # scale factor from Mpc -> kpc -> rad -> arcs
    scale = (old_div(180, pi)) * 3600. / 1000.
    Da = get_Da(z, cosmo)
    theta = kpc * scale / Da
    return theta


def get_Da(z, cosmo):

    import numpy

    try:
        from cosmopy import cosmopy
        #print "Will use cosmopy"
        c = cosmopy.set(cosmo)
        da = c.dang(z)

    except ImportError:
        import Cosmology
        import types
        #print "Will use Cosmology"
        c = Cosmology.Cosmology()
        c.Om, c.Ol, c.h = cosmo

        if type(z) is int or type(z) is float:
            da = (c.Da(z) / c.pc / 1e6)
            return da
        else:
            z = numpy.asarray(z)
            da = numpy.zeros(len(z)) * 0.0
            for i in range(len(z)):
                da[i] = (c.Da(z[i]) / c.pc / 1e6)

    if len(da) == 1:
        return da[0]

    return da


def get_Dl(z, cosmo):

    import numpy

    try:
        import cosmopy
        c = cosmopy.set(cosmo)
        dl = c.dlum(z)

    except:
        import Cosmology
        import types
        c = Cosmology.Cosmology()
        c.Om, c.Ol, c.h = cosmo

        if type(z) is int or type(z) is float:
            dl = (c.Dl(z) / c.pc / 1e6)
            return dl
        else:
            z = numpy.asarray(z)
            dl = numpy.zeros(len(z)) * 0.0
            for i in range(len(z)):
                dl[i] = (c.Dl(z[i]) / c.pc / 1e6)

    if len(dl) == 1:
        return dl[0]

    return dl


# Invert the dl(z) function to obtain z(dl), with dl in Mpc
def invert_dl(cosmo):
    import scipy
    import scipy.interpolate
    import numpy

    z = numpy.arange(0, 2, 0.001)
    dl = get_Dl(z, cosmo)
    return scipy.interpolate.interp1d(dl, z)


# wrapper for wcstool functions
def sky2xy(ra, dec, fitsfile):
    import os
    cmd = "sky2xy %s %s %s J2000" % (fitsfile, ra, dec)
    os.popen(cmd)
    sout = os.popen(cmd)
    vals = sout.readline().split()
    x = float(vals[4])
    y = float(vals[5])
    return x, y


def sky2xy_list(ra, dec, fitsfile):
    import os
    import tableio

    inlist = "/tmp/%s_list.sky2xy" % os.environ['USER']
    outlist = "/tmp/%s_list.sky2xy.out" % os.environ['USER']
    tableio.put_data(inlist, (ra, dec),
                     header='',
                     format="%s %s J2000",
                     append='no')
    cmd = "sky2xy %s @%s > %s" % (fitsfile, inlist, outlist)
    os.popen(cmd)
    (x, y) = tableio.get_data(outlist, cols=(4, 5))
    #print cmd
    os.system("rm %s" % inlist)
    os.system("rm %s" % outlist)
    return x, y


# wrapper for wcstool functions
def xy2sky(x, y, fitsfile, units='degrees'):
    import os
    import tableio

    if units == 'degrees':
        opts = ' -d '
    if units == 'hours':
        opts = ''

    cmd = "xy2sky %s %s %s %s" % (opts, fitsfile, x, y)

    sout = os.popen(cmd)
    vals = sout.readline().split()
    ra = vals[0]
    dec = vals[1]

    if units == 'degrees':
        ra = float(ra)
        dec = float(dec)

    return ra, dec


# Transform ix,iy grid using WCS info from images
# requires wcstools tasks xy2sky and sky2xy
def transform_grid(ACT_fits, optical_fits):

    import tableio
    from pyfits import getheader
    import numpy

    header = getheader(ACT_fits)
    nx = header['NAXIS1']
    ny = header['NAXIS2']
    (iy, ix) = numpy.indices((ny, nx))

    # flatten (x,y) indices arrays and add 1.0 and put in a tmp file
    x = ix.ravel() + 1
    y = iy.ravel() + 1
    tableio.put_data('/tmp/xy_file', (x, y), format="%7d %7d")

    # system call for xy2sky and sky2xy
    cmd1 = "xy2sky -d %s @/tmp/xy_file  > /tmp/radec_file" % ACT_fits
    cmd2 = "sky2xy  %s @/tmp/radec_file > /tmp/xynew_file" % optical_fits
    os.system(cmd1)
    os.system(cmd2)

    # Read in new grid and re-shape
    (ixnew, iynew) = tableio.get_data('/tmp/xynew_file', (4, 5))
    ix_new = ixnew.reshape(ix.shape) - 1.
    iy_new = iynew.reshape(iy.shape) - 1.

    # Clean up files
    os.remove('/tmp/xy_file')
    os.remove('/tmp/radec_file')
    os.remove('/tmp/xynew_file')

    return ix_new, iy_new
