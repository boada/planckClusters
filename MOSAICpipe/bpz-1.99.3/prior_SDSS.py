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
    """Prior based on the SDSS spectroscopic catalog
    This function defines a prior based only on the 
    magnitude of the objects. It assumes that the type 
    fraction does not depend on redshift
    It assumes a shape p(z)=z**2*exp(-(z/zm)**1.5)
    Returns an array pi[z[:],:6]
    The i-band should be  used as M_0
    """

    global zt_at_1p5
    global zt_at_2
    xm = numpy.arange(12., 18.0)
    ft = numpy.array((0.55, 0.21, 0.21, .01, .01, .01))
    zm0 = numpy.array([0.021, 0.034, 0.056, 0.0845, 0.1155, 0.127]) * (old_div(2., 3.))

    if len(ft) != nt:
        print("Wrong number of templates!")
        sys.exit()

    nz = len(z)
    m = numpy.array([m])  #match_resol works with arrays
    m = numpy.clip(m, xm[0], xm[-1])
    zm = match_resol(xm, zm0, m)
    #    zm=zm[0]
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
        old_div(zt_1p5[:, :], zm_1p5), 0., 700.))
    norm = numpy.add.reduce(p_i[:nz, :], 0)
    #Get rid of very low probability levels
    p_i[:nz, :] = numpy.where(
        less(old_div(p_i[:nz, :], norm[:]), old_div(1e-5, float(nz))), 0.,
        old_div(p_i[:nz, :], norm[:]))
    norm = numpy.add.reduce(p_i[:nz, :], 0)
    return p_i[:nz, :] / norm[:] * ft[:]
