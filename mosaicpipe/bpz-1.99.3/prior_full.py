import bpz_tools
from useful import match_resol
import numpy

# Hacked to use numpy and avoid import * commands
# FM
Float = numpy.float
less = numpy.less


def function(z, m, nt):
    """HDFN prior for the main six types of Benitez 2000
    Returns an array pi[z[:],:6]
    The input magnitude is F814W AB
    """
    mmax = 28.

    if nt != 6:
        print "Wrong number of template spectra!"
        sys.exit()

    global zt_at_a
    global zt_at_1p5
    global zt_at_2

    nz = len(z)
    momin_hdf = 20.

    if m <= 20.:
        xm = numpy.arange(12., 18.0)
        ft = numpy.array((0.55, 0.21, 0.21, .01, .01, .01))
        zm0 = numpy.array([0.021, 0.034, 0.056, 0.0845, 0.1155, 0.127]) * (2. /
                                                                           3.)

        if len(ft) != nt:
            print "Wrong number of templates!"
            sys.exit()

        nz = len(z)
        m = numpy.array([m])  #match_resol works with arrays
        m = numpy.clip(m, xm[0], xm[-1])
        zm = match_resol(xm, zm0, m)
        try:
            zt_2.shape
        except NameError:
            t2 = [2.] * nt
            zt_2 = numpy.power.outer(z, t2)
        try:
            zt_1p5.shape
        except NameError:
            t1p5 = [1.5] * nt
            zt_1p5 = numpy.power.outer(z, t1p5)

        zm_3 = numpy.power.outer(zm, 3)
        zm_1p5 = numpy.power.outer(zm, 1.5)
        p_i = 3. / 2. / zm_3 * zt_2[:, :] * numpy.exp(-numpy.clip(
            zt_1p5[:, :] / zm_1p5, 0., 700.))
        norm = numpy.add.reduce(p_i[:nz, :], 0)
        #Get rid of very low probability levels
        p_i[:nz, :] = numpy.where(
            numpy.less(p_i[:nz, :] / norm[:], 1e-5 / float(nz)), 0.,
            p_i[:nz, :] / norm[:])
        norm = numpy.add.reduce(p_i[:nz, :], 0)
        return p_i[:nz, :] / norm[:] * ft[:]

    else:

        m = numpy.minimum(numpy.maximum(20., m), 32)
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
        f_t[3:] = (1. - numpy.add.reduce(f_t[:3])) / 3.
        #Formula:
        #zm=zo+km*(m_m_min)
        #p(z|T,m)=(z**a)*numpy.exp(-(z/zm)**a)
        p_i = zt_at_a[:nz, :6] * numpy.exp(-numpy.clip(zt_at_a[:nz, :6] /
                                                       zmt_at_a[:6], 0., 700.))
        #This eliminates the very low level tails of the priors
        norm = numpy.add.reduce(p_i[:nz, :6], 0)
        p_i[:nz, :6] = numpy.where(
            less(p_i[:nz, :6] / norm[:6], 1e-2 / float(nz)), 0.,
            p_i[:nz, :6] / norm[:6])
        norm = numpy.add.reduce(p_i[:nz, :6], 0)
        p_i[:nz, :6] = p_i[:nz, :6] / norm[:6] * f_t[:6]
        return p_i
