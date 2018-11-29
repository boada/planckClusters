from __future__ import division
from __future__ import print_function
#This module contains several functions which calculate
#observational quantities affected by cosmology

from builtins import range
from builtins import object
from past.utils import old_div
import bpz_tools
import useful
import numpy
import glob
import os
import sys

f_z_sed = bpz_tools.f_z_sed
f_z_sed_AB = bpz_tools.f_z_sed_AB
equal = numpy.equal
log10 = numpy.log10
ABtoVega = bpz_tools.ABtoVega

cho = 2.99e3  # c/H_0 in Mpc
ht = 9.7776e9  # hubble time in h^-1 yr

#Get the ABflux files in stock
ab_db = []
ab_dir = bpz_tools.ab_dir
print("## AB_DIR: %s" % ab_dir, file=sys.stderr)
ab_db = glob.glob(ab_dir + '*.AB')

for i in range(len(ab_db)):
    ab_db[i] = os.path.basename(ab_db[i])
    ab_db[i] = ab_db[i][:-3]


#K-corrections and the like
def kcor(z, sed, filter):
    """K-correction in a giver filter for the spectrum SED at redshift z
    ( m=M+5log(D_L/10pc)+K(z)   )"""
    fo = f_z_sed(sed, filter)
    if type(z) == type(1.): z = numpy.array([z])
    k = 2.5 * log10((1. + z) * fo / f_z_sed(sed, filter, z))
    if len(k) == 1: return k[0]
    else: return k


def reobs(sed,
          m=0.,
          z_0=0.,
          oldfilter='I_LRIS',
          z_new=0.,
          newfilter='V_LRIS',
          cosmology=(0.3, 0.7, .7),
          madau='yes'):
    """Arguments: sed,m,z_0,oldfilter,z_new,newfilter,cosmology
    Takes a galaxy with m at redshift z_0 in oldfilter,
    SED=sed and produces its new magnitude in newfilter at z_new.
    Takes into account cosmological dimming and intergalactic Madau absorption
    The tuple cosmology=(omega,lambda,hubble_constant)
    """
    if sed[-4:] == '.sed': sed = sed[:-4]
    #single_z=type(z_new)==type(z_0)
    single_z = z_new.__class__.__name__[0:3] == z_0.__class__.__name__[0:3]

    if single_z:
        if z_0 == z_new and oldfilter == newfilter: return m
        z_new = numpy.array([z_new])

        #Calculate fnew
    model = '.'.join([sed, newfilter, 'AB'])
    model_path = os.path.join(ab_dir, model)

    #Check whether there are already AB files
    if madau == 'yes':
        if model[:-3] in ab_db:
            zo, f_mod_0 = useful.get_data(model_path, (0, 1))
            fnew = z_new * 0.
            for i in range(len(z_new)):
                fnew[i] = useful.match_resol(zo, f_mod_0, z_new[i])
        else:
            fnew = f_z_sed_AB(sed, newfilter, z_new, 'nu')
    else:
        fnew = f_z_sed(sed, newfilter, z_new, units='nu', madau=madau)

    fnew = numpy.where(
        equal(fnew, 0.), 99.,
        fnew)  # if the new flux is 0, returns 99. (code non-detection)

    #Calculate f_old
    model = '.'.join([sed, oldfilter, 'AB'])
    model_path = os.path.join(ab_dir, model)

    #Check whether there are already AB files
    if madau == 'yes':
        if model[:-3] in ab_db:
            zo, f_mod_0 = useful.get_data(model_path, (0, 1))
            f_old = useful.match_resol(zo, f_mod_0, z_0)
        else:
            f_old = f_z_sed_AB(sed, oldfilter, numpy.array([z_0]), units='nu')
    else:
        f_old = f_z_sed(sed, oldfilter, numpy.array([z_0]), units='nu')

    k = 2.5 * log10((old_div((1. + z_new), fnew)) * (old_div(f_old,
                                                             (1. + z_0))))

    if single_z and z_0 == z_new[0]:
        m_obs = m + k
        return m_obs[0]

        #Distance modulus
    dist = dist_mod(z_new, cosmology) - dist_mod(z_0, cosmology)
    m_obs = m + dist + k
    if single_z: return m_obs[0]
    else: return m_obs


def color_z(sed,
            filter_new,
            filter_old,
            z=numpy.arange(0., 1.5, 0.5),
            calibration='AB',
            file=None):
    """
    Calculates the color filter_new-filter_old at the redshift vector z.
    It can return the color in Vega or AB calibrations
    Usage:
    gr=color_z('El_cww','g_WFC','r_WFC',z=numpy.arange(0.,2.,0.1),'Vega')
    It also works with scalars, e.g.
    gr=color_z('El_cww','g_WFC','r_WFC',1.2)
    """
    try:
        n = len(z)
    except:
        z = numpy.array([z])
        n = 1
    color = z * 0.
    for i in range(len(z)):
        color[i] = reobs(sed, 0., z[i], filter_old, z[i], filter_new)

    if calibration == 'Vega':
        color += ABtoVega(0., filter_new) - ABtoVega(0., filter_old)

    if file == None:
        if n == 1: color = color[0]
        return color
    else:
        put_data(file, (z, color),
                 header='z     %s-%s(%s) ' %
                 (filter_new, filter_old, calibration))


def m_abs(m, z, sed, filter, cosmology=(0.3, .7, .7), filter2=None):
    """Arguments: m,z,sed,filter,cosmology,filter2
    If filter2 is used, returns the absolute magnitude
    in a different filter"""
    #print z,sed,filter
    mabs = m - dist_mod(z, cosmology) - kcor(z, sed, filter)

    if filter2 != None:
        #We add the color filter2-filter at z=0
        mabs = mabs + reobs(sed, 0., 0., filter, 0., filter2)
    return mabs

#Incluir algo que pase de magnitudes absolutas a luminosidades solares!!
#def luminosity(m,filter):
#    Basicamente hallar la magnitud absoluta del sol en el filtro que corresponda
#    m_sun y normalizar
#    return m-dist_mod(z,cosmo[0],cosmo[1],cosmo[2])

#COSMOLOGICAL DISTANCES


def dl_lambda(z, omega=.3, h=1.):
    """Aproximation for the luminosity distance
    for flat cosmologies with cosmological constant
    ApJSS, Ue Li Pen 120:4950, 1999"""

    if omega < 0.2:
        raise """omega less than 0.2: outside
    parameter range for the aproximation"""

    if h > 1. or h < .4:
        print("Wrong value for h", h)
        sys.exit()

    def eta(a, om):
        s = (old_div((1. - om), om))**(old_div(1., 3.))
        return 2.*numpy.sqrt(s**3+1.)*\
               (a**(-4)-0.1540*s/a**3+0.4304*s**2/a**2+
         0.19097*s**3/a+0.066941*s**4)**(old_div(-1.,8.))

    return 2.9979*1e5/(h*100.)*(1.+z)*\
        (eta(1.,omega)-eta(old_div(1.,(1.+z)),omega))


def dl_nolambda(z, omega=.3, h=.7):
    """Luminosity distance for a lambda=0
    universe"""
    cosa = sqrt(1. + omega * z)
    return 2.9979 * 1e5 / (h * 100.) * z * (1. + cosa + z) /\
            (1. + cosa + omega * z / 2.)


def dl(z, cosmology=(.3, .7, .7)):
    omega, l, h = cosmology
    if l > 0.:
        if l + omega != 1.: raise 'lambda>0 but no flat cosmology!'
        return dl_lambda(z, omega, h)
    if l == 0: return dl_nolambda(z, omega, h)
    if l < 0: raise 'lambda<0!!'


def da(z, cosmology=(.3, .7, .7)):
    return old_div(dl(z, cosmology), (1. + z)**2)

######New distance definitions. Hogg 1999, astro-ph/9905116############Not tested!


def dh(cosmology=(0.3, 0.7, 0.7)):
    return old_div(3000., cosmology[2])


def omega_k(cosmology=(0.3, 0.7, 0.7)):
    return 1. - cosmology[0] - cosmology[1]


def e(z, cosmology=(0.3, 0.7, 0.7)):
    o_k = omega_k(cosmology)
    o_m = cosmology[0]
    o_l = cosmology[1]
    return sqrt(o_m * (1. + z)**3 + o_k * (1. + z)**2 + o_l)


def lookback_time(z, cosmology):
    zp = numpy.arange(0., z, min([0.001, old_div(z, 100.)]))
    return 9.78e9 / cosmology[2] * trapz(1. / (1 + zp) / e(zp, cosmology), zp)

#def dc(z,cosmology=(0.3,0.7,0.7)):
#    dz=0.000001
#    xz=arange(0.,z+dz,dz)
#    ez=e(xz,cosmology)
#    return dh(cosmology)*trapz(ez,xz)

#def dm(z,cosmology=(0.3,0.7,0.7)):
#    o_k=omegak(cosmology)
#    if o_k>0:
#   return dh(cosmology)/sqrt(o_k)*sinh(sqrt(o_k)*dc(z,cosmology)/dh(cosmology))
#    elif o_k==0:
#   return dc(z,cosmology)
#    else:
#   o_k=abs(o_k)
#   return dh(cosmology)/sqrt(o_k)*sin(sqrt(o_k)*dc(z,cosmology)/dh(cosmology))

#def da_h(z,cosmology=(0.3,0.7,0.7)):
#    return dm(z,cosmology)/(1.+z)

#def dl_h(z,cosmology=(0.3,0.7,0.7)):
#    return dm(z,cosmology)*(1.+z)

######################## Hogg 1999 ########################################


class vc(object):
    def __init__(self,
                 z=0.57,
                 sed='Sbc_cww',
                 m=20.,
                 em=0.02,
                 filter='B_Johnson',
                 cosmo=(0.3, 0.7, 0.7),
                 vc_filter='B_Johnson'):
        """Generates velocity dispersion and error for a galaxy using TF or Faber--Jackson
    Inputs: redshift, magnitude, error, filter, spectral type, cosmo and filter for TF in the rest frame
    (FB Jackson always uses BJ as the rest frame filter)
    Usage:
    cosa=vc(0.55,'El_cww',20.,0.02,'I_Cousins',(0.3,0.7,0.7))
    cosa=vc(0.55,'Scd_cww',20.,0.02,'I_Cousins',(0.3,0.7,0.7),'H_Johnson')
    Uses the closest rest frame filter by default
    Assumes that the input magnitudes are AB
    """
        self.sed = sed

        #Everything has to be properly transformed from AB to Vega!!

        #Info about TF
        #Pierce and Tully 1992
        #vc=158.1*10.**(-(mabs+constant_TF)/slope_TF)
        #It actually only works for Sbc galaxies
        #For bluer stuff it is better to use the reddest filter, I_Cousins
        #H_Johnson gives weird results, it may be due to template problems

        self.filters_TF = ['B_Johnson', 'R_Cousins', 'I_Cousins', 'H_Johnson']
        self.centers_TF = [4477.8, 6648.33, 8086.4, 16509.64]
        self.slope_TF = [7.48, 8.23, 8.72, 9.50]
        self.constant_TF = [19.55, 20.46, 20.94, 21.67]
        self.error_TF = [0.14, 0.10, 0.10, 0.08]

        #Info about FB
        # Kochanek 1994, ApJ
        # sigma_*=(225+-22.5)*(L/L_*)**(.24+-0.03)
        # the error is approximate Kochanek 1996, magnitude is BJ
        # with M_B (BJ) = -19.9+5*log10(h)
        # Using L/L*=10.**[-0.4[M-M_*]]
        # sigma_*=225.*10.(-(mabs+19.9)/10.42)

        if sed == 'El_cww':
            self.m_abs = ABtoVega(m_abs(m, z, sed, filter, cosmo, 'BJ'), 'BJ')
            self.v_c = 225. * 10.**(.4 * (
                -self.m_abs - (19.9 - 5. * log10(cosmo[2]))) * 0.24)
            self.e_v_c = (old_div(25., 225.)) * self.v_c
            self.filter_v_c = 'BJ'
        else:
            if sed == 'Sbc_cww' or sed == 'Scd_cww':
                fc = filter_center(filter)
                #Look for the closest filter
                k = argmin(abs(numpy.array(self.centers_TF) - old_div(fc, (
                    1. + z))))
            elif sed == 'Im_cww' or sed == 'SB2_kin' or sed == 'SB3_kin':
                k = 2  #Use I_Cousins
            self.m_abs = ABtoVega(
                m_abs(m, z, sed, filter, cosmo, self.filters_TF[k]),
                self.filters_TF[k])
            self.v_c = 158.1 * 10.**(old_div(
                -(self.m_abs + self.constant_TF[k]), self.slope_TF[k]))
            self.e_v_c = self.v_c * 2.3 / self.slope_TF[k] * sqrt(
                em**2 + self.error_TF[k]**2)
            self.filter_v_c = self.filters_TF[k]


def dist_mod(z, cosmology=(0.3, .7, .7)):
    """Usage: dist_mod(z,cosmology)"""
    return 25. + 5. * log10(dl(z, cosmology))


def angular_size(length, z, cosmology=(0.3, 0.7, .7)):
    """Usage: angular_size(length,z,cosmology)"""
    return length / da(z, cosmology) / 3.141592654 * 180. * 3600.


def physical_size(angle, z, cosmology=(0.3, 0.7, .7)):
    """Usage: physical_size(angle,z,cosmology)
       Units: arcseconds, Mpc"""
    return angle / 360. / 60. / 60. * 2. * 3.141592654 * da(z, cosmology)


def lookback_time_open(z, omega=0.3, h=0.7):
    """Usage: lookback_time_open(z,omega,h)
       Units: Myr Approximation from Peacock"""
    omega_z = omega * (1. + z) / (1. + omega * z)
    h_z = h * (1. + z) * sqrt(1. + omega * z)
    t = h_z * (1. + old_div(omega_z**.6, 2.))
    return ht / t / 1e9


def test():
    test = 'reobs'
    Testing(test)
    z1, z2, f1, f2, t, c = 0.2, 0.8, 'V_LRIS', 'I_LRIS', 'El_cww', (0.3, 0.7,
                                                                    0.7)
    dr = reobs(t, 0., z1, f1, z2, f2, c)
    ds = (reobs(t, 0., 0., f1, 0., f2, c) + 5. * log10(old_div(
        dl(z2, c), dl(z1, c))) + (kcor(z2, t, f2) - kcor(z1, t, f1)))
    print('dr,ds')
    print(dr, ds)

    ask('More?')

    #The reobs function has been tested indirectly by using
    #bpz to estimate redshifts of objects whose colors were generated
    #by reobs, the agreement is perfect.
    #The rest of the cosmological functions can be tested
    #by comparing them with plots in the original references
    pass

    test = 'Distance modulus'
    Testing(test)
    print('Compare with Peebles Physical Cosmology, page 329')
    print('Values of Omega are 0.2,0.5,1.')
    z = numpy.arange(0.0001, 10., .01)
    omega = [0.2, 0.5, 1.]
    d = []
    dlambda = []
    p1 = FramedPlot()
    p1.title = 'Lambda = 0'
    p1.xrange = -0.2, 10.
    p1.yrange = 41., 52.
    p2 = FramedPlot()
    p2.title = 'Flat universes'
    p2.xrange = -0.2, 10.
    p2.yrange = 41., 52.
    for i in range(len(omega)):
        d.append(dist_mod(z, (omega[i], 0., 1.)))
        dlambda.append(dist_mod(z, (omega[i], 1. - omega[i], 1.)))
        p1.add(Curve(z, d[i]))
        p2.add(Curve(z, dlambda[i]))
    p1.show()
    p2.show()

    print()
    print()
    print()

    ask('More?')

    test = 'Cosmological distances'
    Testing(test)
    z = numpy.arange(0., 4., .01)
    da1 = old_div(da(z, (1., 0., 1.)), cho)
    da2 = old_div(da(z, (.3, 0., 1.)), cho)
    da3 = old_div(da(z, (.3, 0.7, 1.)), cho)
    p = FramedPlot()
    p.add(Curve(z, da1))
    p.add(Curve(z, da3))
    p.add(Curve(z, da2, style='dashed'))
    p.yrange = 0., 1.
    p.show()
    print("Compare with Cosmological Physics, page 93")

    print()
    print()
    print()

    ask('More?')

    test = 'K-corrections'
    Testing(test)
    print('Compare with Physical cosmology, page 331')

    z = numpy.arange(0., 1.5, .01)
    p = FramedPlot()
    p.xrange = 0., 1.5
    p.yrange = -1., 5.
    for tipo in ['El_cww', 'Sbc_cww', 'Scd_cww']:
        k = kcor(z, tipo, 'B_Johnson')
        p.add(Curve(z, k))
    p.show()

    print()
    print()
    print()


if __name__ == '__main__':
    test()
else:
    pass
