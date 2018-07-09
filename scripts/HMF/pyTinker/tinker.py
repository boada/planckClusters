#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from builtins import range
from builtins import object
from past.utils import old_div
import numpy
import os
import sys
import tableio
from scipy import interpolate
#import scipy
import math
from cosmopy import cosmopy


class counts(object):
    """ A Class to computed dn/dM and dN/dz for the Universal Tinker mass function
    F.Menanteau, Rutgers, May 2012"""

    # Set up cosmology and location of files
    def __init__(
            self,
            Om=0.3,  # OMEGA_M
            OL=0.7,  # OMEGA_L
            H0=100.0,
            h0=None,
            Ob=0.04,  # OMEGA_B
            sigma8=0.9,
            spectral_index=0.961,
            TF_file=None,
            ITRANS=4):

        self.Ob = Ob
        self.Om = Om
        self.OL = OL
        self.h = old_div(H0, 100.)
        self.sigma8 = sigma8
        self.ITRANS = ITRANS
        self.spectral_index = spectral_index

        arch = check_arch()

        # PyTinkerPATH for location of aux files
        if arch['MACH'] == 'x86_64' and arch['OS'] == 'Linux':
            self.PyTikerPATH = os.path.join(os.environ['HOME'],
                            "Projects/planckClusters/scripts/HMF/pyTinker")
            self.tkbin = os.path.join(self.PyTikerPATH, "bin/massfunction")

        # Files to run MF_Code fortran code
        print("# Will use:%s " % (self.tkbin), file=sys.stderr)
        self.batFile = "/tmp/mf_PyTinker.bat"  # Bat file to run MF_Code

        # If ITRANS == 5 -- Eisenstein & Hu Transfer Function
        if h0 and self.ITRANS == 5:
            self.h0 = h0
            print(
                "# Will use ITRANS:%s and h0:%s\n" % (ITRANS, self.h0),
                file=sys.stderr)
        else:
            self.h0 = self.h

        # If ITRANS == 11 -- CAMB Transfer Function
        if TF_file and self.ITRANS == 11:
            self.TF_file = os.path.join(self.PyTikerPATH, "conf", TF_file)
            print(
                "# Will use ITRANS:%s with TF_FILE:%s\n" % (ITRANS,
                                                            self.TF_file),
                file=sys.stderr)
        else:
            self.TF_file = "dummy_camb.dat"

        return

    def get_dndz(self, Mlim):

        # 1 stdrad = (180/pi)**2 square-deg
        sterad2degsq = 1. / (180. / math.pi)**2

        dNdz = numpy.zeros_like(self.zx)
        Nsum = numpy.zeros_like(self.zx)

        # Set the distance object
        c = cosmopy.set((self.Om, self.OL, self.h))
        dV = c.dvol_comov(self.zx)  # Mpc^3/dz/strad

        if not isinstance(Mlim, (list, numpy.ndarray, tuple)):
            Mlim = numpy.ones_like(self.zx) * Mlim

        cumsum = 0.0
        print(Mlim)
        for i in range(self.zx.size):

            # Spline integration
            nc = spline_integral(self.M, self.dndM[i], x1=Mlim[i])
            if nc < 0:
                nc = 0
            dNdz[i] = dV[i] * nc * sterad2degsq
            cumsum += dNdz[i] * self.dz
            Nsum[i] = cumsum
            #Nsum[i] = dNdz[i] * self.dz

        return dNdz, Nsum

    # Get the dndM curves for all redshifts we want
    # We do this only once and
    def get_dndM(self, z1=0, z2=3.0, dz=0.01, delta=200):

        # Define the z-array and make it visible
        self.zx = numpy.arange(z1, z2 + dz, dz)
        self.dz = dz

        z = self.zx
        n = self.zx.size
        dndM = []
        for i in range(n):
            print("# Doing z:%.4f" % z[i], file=sys.stderr)
            x, y = self.run_dndM(z[i], delta=delta)
            dndM.append(y)

        # Make them visible
        self.M = numpy.array(x)
        self.dndM = numpy.array(dndM)
        return

    def run_dndM(self, z, delta=200):

        # Define overdensity first
        self.deltaHalo = delta

        # Write the input mf.bat file before running it.  We keep the
        # same name for the mf.bat and re-write it every time we call
        # the function.
        self.write_MF_input(z)

        # Make the call and read in the array
        cmd = "%s %s > /dev/null 2>&1" % (self.tkbin, self.batFile)
        #cmd = "%s %s" % (self.tkbin,self.batFile)
        os.system(cmd)

        # Read in the dndM array
        # M [Msun*h100], and dndM [h*100^4 counts/Mpc^3/Msun]
        (M, dndM) = self.read_dndM()
        return M, dndM

    # Read in the dndM array
    def read_dndM(self):
        (M, dndM) = tableio.get_data(self.dndMFile, cols=(0, 1))
        os.remove(self.dndMFile)
        return M, dndM

    def write_MF_input(self, z):

        self.dndMFile = "/tmp/tinker_paper_%.3f.dndM" % z
        self.rootFile = "/tmp/tinker_paper_%.3f" % z
        with open(self.batFile, 'w') as o:
            o.write(
                "%----------------------------------------------------------------------------\n"
            )
            o.write("%  Cosmological Parameters\n")
            o.write(
                "%----------------------------------------------------------------------------\n"
            )
            o.write("GAMMA           0.2\n")  # Only used for TF options 4
            o.write("OMEGA_M         %s\n" % self.Om)
            o.write("OMEGA_B         %s\n" % self.Ob)
            o.write("SIGMA_8         %s\n" % self.sigma8)
            o.write("RHO_CRIT        2.775E11		% h^2 M_sol/Mpc^3\n")
            o.write("SPECTRAL_INDX   %s\n" % self.spectral_index)
            o.write("HUBBLE          %s\n" % self.h0)  # Only used for TF option 5
            #o.write("HUBBLE          0.693\n")
            o.write("DELTA_CRIT      1.686\n")
            o.write("ITRANS          %s\n" %
                    self.ITRANS)  # We might want to change this one later
            o.write("TF_file         %s\n" % self.TF_file)
            o.write("DELTA_HALO      %s\n" % self.deltaHalo)
            o.write("REDSHIFT        %s\n" % z)
            o.write("root_filename   %s\n" % self.rootFile)
            o.write("%--------------------------\n")
            o.write("% TRANSFER FUNCTION OPTIONS\n")
            o.write("%--------------------------\n")
            o.write(
                "% ITRANS = 4 is for the Efstathiou, Bond, White transfer function. Uses GAMMA ONLY.\n"
            )
            o.write(
                "% ITRANS = 5 is Eisenstein & Hu, uses OMEGA_M, OMEGA_B, HUBBLE\n")
            o.write(
                "% ITRANS = 11 is an exterior file, must be in CMBFAST format, uses TF_file for filename\n"
            )
        return

def spline_integral(x, y, x1=None, x2=None):
    if not x1:  # Integration limits
        x1 = x.min()
    if not x2:
        x2 = x.max()
    tck = interpolate.splrep(x, y, s=0)  # Spline representation
    return interpolate.splint(x1, x2, tck)  # Integral

# Check the architecture of the system
def check_arch():

    name = {}

    cmd = 'uname -m'
    name['MACH'] = (os.popen(cmd).read()).strip()

    cmd = 'uname -s'
    name['OS'] = (os.popen(cmd).read()).strip()
    return name
