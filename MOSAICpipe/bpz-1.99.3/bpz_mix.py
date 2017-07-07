from bpz_tools import *
from cosmology import *
import sys
import numarray


# **********************************************
#        Auxiliar mixing functions
# ***********************************************

def reobs_mix(t,m=0.0,z_0=0.0,
              oldfilter='I_LRIS',
              z_new=0.0,
              newfilter='V_LRIS',
              cosmology=(0.3,0.7,.7),madau='yes',sed_lib='CWWSB_Benitez2003.list'):

    """ Returns the mixed reobs magnitude from a non-integer type """

    # Get the types and SED's names
    types = get_Types(t)
    seds  = get_str( sed_dir + sed_lib, 0)

    # non integer part of the T-type
    tfrac   = t % int(t)
    newmag = {}
    flux   = []

    for T in types:

        newmag[T] = reobs(seds[T],m=m,z_0=z_0,
                          oldfilter=oldfilter,
                          z_new=z_new,
                          newfilter=newfilter,
                          cosmology=cosmology,madau='yes')
        flux.append(10.**(-0.4*newmag[T]))

    # And compute the total flux using general shape
    if (tfrac == 0):
        return newmag[types[0]]
    else:
        flux_sum  = flux[0]*(1.0 - tfrac) + flux[1]*tfrac
        return -2.5*log10(flux_sum)



def M_abs_mix(m,z,t,filter,
              sed_lib='CWWSB_Benitez2003.list',
              cosmology=(0.3,.7,.7),
              filter2=None):

    """ Returns the absolute magnitude from a non-integer type """

    # Get the types and SED's names
    types = get_Types(t)
    seds  = get_str( sed_dir + sed_lib, 0)

    tfrac = t % int(t) # non integer part of the T-type
    abs  = {}
    Flux = []

    for T in types:
        abs[T]  = M_abs(m,z,
                        seds[T],
                        filter,
                        cosmology=cosmology,
                        filter2=filter2)
        Flux.append(10.**(-0.4*abs[T]))

    # And compute the total flux using general shape
    if (tfrac == 0):
        return abs[types[0]]
    else:
        Flux_sum  = Flux[0]*(1.0 - tfrac) + Flux[1]*tfrac
        return -2.5*log10(Flux_sum)


def Kcorr_mix(z,t,filter,sed_lib='CWWSB_Benitez2003.list'):

    """ Returns the Kcorrection from a non-integer type, using reobs """

    # Get the types and SED's names
    types = get_Types(t)
    seds  = get_str( sed_dir + sed_lib, 0)
    # non integer part of the T-type
    tfrac = t % int(t)
    #if tfrac == 0.0:
    #    return Kcorr(z,seds[types[0]],filter)

    if type(z)==type(1.):z=array([z])

    Fo = []
    Fz = []
    # For each type
    for T in types:

        sed = seds[T]
        if sed[-4:]=='.sed': sed=sed[:-4]
        # Get the models for each type
        model_name = join([sed,filter,'AB'],'.')
        model_path = os.path.join(ab_dir,model_name)

        # Search if the *.AB file exist
        if model_name[:-3] in ab_db:
            zo,f_mod_0 = get_data(model_path,(0,1))
            fo = f_mod_0[0]
            fz = z*0.
            for i in range(len(z)):
                fz[i] = match_resol(zo,f_mod_0,z[i])
        # if not, then compute
        else:

            fo = f_z_sed_AB(sed,filter,0.0,'nu')
            fz = f_z_sed_AB(sed,filter,z,'nu')

        Fo.append(fo)
        Fz.append(fz)

    if tfrac == 0.0:
        fo = Fo[0]
        fz = Fz[0]
    else:
        fo = Fo[0]*(1.0 - tfrac) + Fo[1]*tfrac
        fz = Fz[0]*(1.0 - tfrac) + Fz[1]*tfrac

    k = 2.5*log10((1.+z)*fo/fz)


    if len(k)==1:
        return k[0]
    else:
        return k

# *******************************************************
#    Extra functions to complement/modify bpz_tools
# *******************************************************
def Kcorr(z,sed,filter):

    """K-correction in a giver filter for the spectrum SED at redshift z
    ( m=M+5log(D_L/10pc)+K(z)   )

    This version tries to read first the values stored in the '*.AB'
    files in the BPZDIR/AB directory, in a similar way as the
    reobs routine from Txitxo. The results are very similar to the function
    kcor(z,sed,filter) in cosmology.py
    """

    if sed[-4:]=='.sed': sed=sed[:-4]
    if type(z)==type(1.):z=array([z])

    # Get the AB file name built
    model=join([sed,filter,'AB'],'.')
    model_path=os.path.join(ab_dir,model)

    # Search if the *.AB file exist
    if model[:-3] in ab_db:
        zo,f_mod_0 = get_data(model_path,(0,1))
        fo = f_mod_0[0]
        fz = z*0.
        for i in range(len(z)):
            fz[i] = match_resol(zo,f_mod_0,z[i])
    # if not, then compute
    else:
        fo = f_z_sed_AB(sed,filter,0.0,'nu')
        fz = f_z_sed_AB(sed,filter,z,'nu')

    k = 2.5*log10((1.+z)*fo/fz)

    #print (1.+z)*fo/fz

    if len(k)==1:
        return k[0]
    else:
        return k

def M_abs(m,z,sed,filter,cosmology=(0.3,.7,.7),filter2=None):

    """Arguments: m,z,sed,filter,cosmology,filter2
    If filter2 is used, returns the absolute magnitude
    in a different filter.
    Same as old m_abs, but now it uses Kcorr function instead
    of kcor, that reads AB files and is much faster
    """

    mabs = m - dist_mod(z,cosmology) - Kcorr(z,sed,filter)

    if filter2<>None:
	#We add the color filter2-filter at z=0
	mabs = mabs + reobs(sed,0.,0.,filter,0.,filter2)
    return mabs


def get_Types(type):

    """ Generic function to obtain the types'mix
    for a given non-integer galaxy type"""

    types = []
    if int(type) == type:
        types.append(int(type) - 1)
    else:
        types = [int(type) - 1, int(type)]

    return types



