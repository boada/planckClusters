from __future__ import print_function
from __future__ import division
from past.utils import old_div
import bpz_tools
from useful import match_resol
import numpy

# Hacked to use numpy and avoid import * commands
# FM
Float = numpy.float
less = numpy.less


def function(z, m, nt):
    """HDFN prior for the main six types of Benitez 2000
    Returns an numpy.array pi[z[:],:6]
    The input magnitude is F814W AB
    """
    mmax = 28.

    if nt != 6:
        print("Wrong number of template spectra!")
        sys.exit()

    global zt_at_a
    nz = len(z)
    momin_hdf = 20.
    if m > 32.: m = 32.
    if m < 20.: m = 20.
    a = numpy.array((2.465, 1.806, 1.806, 0.906, 0.906, 0.906))
    zo = numpy.array((0.431, 0.390, 0.390, 0.0626, 0.0626, 0.0626))
    km = numpy.array((0.0913, 0.0636, 0.0636, 0.123, 0.123, 0.123))
    fo_t = numpy.array((0.35, 0.25, 0.25))
    k_t = numpy.array((0.450, 0.147, 0.147))
    dm = m - momin_hdf
    zmt = numpy.clip(zo + km * dm, 0.01, 15.)
    zmt_at_a = zmt**(a)
    #We define z**a as global to keep it 
    #between function calls. That way it is 
    # estimated only once
    try:
        zt_at_a.shape
    except NameError:
        zt_at_a = numpy.power.outer(z, a)

        #Morphological fractions
    f_t = numpy.zeros((len(a), ), Float)
    f_t[:3] = fo_t * numpy.exp(-k_t * dm)
    f_t[3:] = old_div((1. - numpy.add.reduce(f_t[:3])), 3.)
    #Formula:
    #zm=zo+km*(m_m_min)
    #p(z|T,m)=(z**a)*numpy.exp(-(z/zm)**a)
    p_i = zt_at_a[:nz, :6] * numpy.exp(-numpy.clip(
        old_div(zt_at_a[:nz, :6], zmt_at_a[:6]), 0., 700.))
    #This eliminates the very low level tails of the priors
    norm = numpy.add.reduce(p_i[:nz, :6], 0)
    p_i[:nz, :6] = numpy.where(
        less(
            old_div(p_i[:nz, :6], norm[:6]), old_div(1e-2, float(nz))), 0.,
        old_div(p_i[:nz, :6], norm[:6]))
    norm = numpy.add.reduce(p_i[:nz, :6], 0)
    p_i[:nz, :6] = p_i[:nz, :6] / norm[:6] * f_t[:6]
    return p_i
