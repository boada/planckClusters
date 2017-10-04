from __future__ import print_function
from __future__ import division
from past.utils import old_div
from bpz_tools import *


def function(z, m, nt):
    """HDFN prior for the 11 types, including
    Manucci et al. El,S0,Sa,Sb,Sc and 6
    Kinney et al. starbursts
    Returns an array pi[z[:],:11]
    The input magnitude is F814W AB
    (Adapted from the hdfn prior of Benitez 2000)
    """
    mmax = 28.

    if nt != 11:
        print("Wrong number of template spectra!")
        sys.exit()

    global zt_at_a
    nz = len(z)
    momin_hdf = 20.
    if m > 32.: m = 32.
    if m < 20.: m = 20.
    a = array((2.465, 2.465, 1.806, 1.806, 1.806, 0.906, 0.906, 0.906, 0.906,
               0.906, 0.906))
    zo = array((0.431, 0.431, 0.390, 0.390, 0.390, 0.0626, 0.0626, 0.0626,
                0.0626, 0.0626, 0.0626))
    km = array((0.0913, 0.0913, 0.0636, 0.0636, 0.0636, 0.123, 0.123, 0.123,
                0.123, 0.123, 0.123))
    #fo_t=array((0.35,0.25,0.25))
    fo_t = array((0.175, 0.175, 0.166, 0.166, 0.166))
    k_t = array((0.450, 0.45, 0.147, 0.147, .147))
    dm = m - momin_hdf
    zmt = clip(zo + km * dm, 0.01, 15.)
    zmt_at_a = zmt**(a)
    #We define z**a as global to keep it 
    #between function calls. That way it is 
    # estimated only once
    try:
        zt_at_a.shape
    except NameError:
        zt_at_a = power.outer(z, a)

        #Morphological fractions
    f_t = zeros((len(a), ), Float)
    f_t[:5] = fo_t * exp(-k_t * dm)
    f_t[5:] = old_div((1. - add.reduce(f_t[:5])), 3.)
    #Formula:
    #zm=zo+km*(m_m_min)
    #p(z|T,m)=(z**a)*exp(-(z/zm)**a)
    p_i = zt_at_a[:nz, :11] * exp(-clip(old_div(zt_at_a[:nz, :11], zmt_at_a[:11]), 0.,
                                        700.))
    #This eliminates the very low level tails of the priors
    norm = add.reduce(p_i[:nz, :11], 0)
    p_i[:nz, :11] = where(
        less(old_div(p_i[:nz, :11], norm[:11]), old_div(1e-2, float(nz))), 0.,
        old_div(p_i[:nz, :11], norm[:11]))
    norm = add.reduce(p_i[:nz, :11], 0)
    p_i[:nz, :11] = p_i[:nz, :11] / norm[:11] * f_t[:11]
    return p_i
