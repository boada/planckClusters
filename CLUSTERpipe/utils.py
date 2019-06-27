import sys
import math
import numpy
import matplotlib.patches
from past.utils import old_div
import ezgal # BC03 model maker
import os

Polygon = matplotlib.patches.Polygon
sout = sys.stderr
int16 = numpy.int16
float64 = numpy.float64
land = numpy.logical_and

def KE(cosmo):

    # check to make sure we have defined the bpz filter path
    if not os.getenv('EZGAL_FILTERS'):
        os.environ['EZGAL_FILTERS'] = (f'{os.environ["HOME"]}/'
                                       'Projects/planckClusters/MOSAICpipe/'
                                       'bpz-1.99.3/FILTER/')

    model = ezgal.model('bc03_ssp_z_0.02_salp.model')
    model = model.make_exponential(1)
    model.set_cosmology(Om=cosmo.Om0, Ol=cosmo.Ode0, h=cosmo.h, w=cosmo.w(0))

    model.add_filter('g_MOSAICII.res', name='g')
    model.add_filter('r_MOSAICII.res', name='r')
    model.add_filter('i_MOSAICII.res', name='i')
    model.add_filter('z_MOSAICII.res', name='z')
    model.add_filter('K_KittPeak.res', name='K')

    # Blanton 2003 Normalization
    Mr_star = -20.44 + 5 * numpy.log10(cosmo.h) # abs mag.
    # set the normalization
    model.set_normalization('sloan_r', 0.1, Mr_star, vega=False)

    return model

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

    (z,
    k_g,
    k_r,
    k_i,
    k_z,
    e_g,
    e_r,
    e_i,
    e_z) = tableio.get_data(modelfile,
                            cols=(0, 12, 13, 14, 15, 17, 18, 19, 20))

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
    dx = math.log10(-math.log(zp)) / b
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
    Mi_star = cosmology.reobs(
        'El_Benitez2003',
        m=Mb_star,
        oldfilter="B_Johnson",
        newfilter="i_MOSAICII")
    return Mi_star + 5.0 * math.log10(dlum) + 25


####################################################
# Fake an ellipse using an N-sided polygon
#####################################################
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
    xtmp = (A + 5) * numpy.cos(t)
    ytmp = (B + 5) * numpy.sin(t)

    x = xtmp * cos(angle) - ytmp * sin(angle) + xo
    y = xtmp * sin(angle) + ytmp * cos(angle) + yo
    return Polygon(list(zip(x, y)), **kwargs)


##############################
# A circle as a polygon too
###############################
def PCircle(xxx_todo_changeme2, radius, resolution=100, **kwargs):
    (xo, yo) = xxx_todo_changeme2
    pi = math.pi
    t = 2 * pi / resolution * numpy.arange(resolution)
    xtmp = radius * numpy.cos(t)
    ytmp = radius * numpy.sin(t)
    x = xtmp + xo
    y = ytmp + yo
    return Polygon(list(zip(x, y)), **kwargs)


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


def Mass_calib(N200, L200, LBCG, z, h=1.0):

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

    LBCG_N = aN * N200**bN
    LBCG_L = aL * L200**bL

    M_N200 = (old_div(1.e14, h)) * M0_N * ((old_div(N200, 20.0))**alphaN) * (
        old_div(LBCG, LBCG_N))**gammaN
    M_L200 = (old_div(1.e14, h)) * M0_L * ((old_div(L200, 40.0))**alphaL) * (
        old_div(LBCG, LBCG_L))**gammaL
    M_LBCG = (old_div(1.e14, h)) * 1.07 * (old_div(LBCG, 5.0))**1.10

    return M_N200, M_L200, M_LBCG

def color(s, ncolor, nfont):

    colors = {'red': 31,
              'green': 32,
              'yellow': 33,
              'blue': 34,
              'purple': 35,
              'cyan': 36,
              'white': 37,
              'smoothgreen': '38;5;42',
              'magenta': '38;5;55',
              'turqoise': '38;5;50'}

    fonts = {'normal': 0,
             'bold': 1,
             'light': 2,
             'italic': 3,
             'underline': 4,
             'blink': 5}

    if isinstance(ncolor, str):
        try:
            ncolor = colors[ncolor]
        except KeyError:
            print('Color not understood')
            return s

    if isinstance(nfont, str):
        try:
            nfont = fonts[nfont]
        except KeyError:
            print('Font not understood')
            return s

    return "\033[{};{}m{}\033[0m".format(nfont, ncolor, s)
