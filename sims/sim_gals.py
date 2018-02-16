#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from builtins import range
from builtins import object
import os
import time
import extras
import tableio
import astrometry
import cosmology
import math
import scipy
import scipy.interpolate
from astropy.io import fits as pyfits
#from pyfits import getheader
from pyraf import iraf
from iraf import artdata
import numpy

class simgal(object):
    def __init__(self,
                 fields,
                 Ngal=25,
                 filter='SLOAN-SDSS.r',
                 cosmo=(0.3, 0.7, 0.7),
                 pixscale=0.25):

        self.fields = fields
        self.filter = filter
        self.Ngal = Ngal
        self.pixscale = pixscale
        self.cosmo = cosmo

        self.pipeline = '/home/boada/Projects/planckClusters/MOSAICpipe'
        # The input datapath
        self.datapath = '/home/boada/Projects/planckClusters/data/proc2'

        # The output datapath
        self.outpath = '/home/boada/Projects/planckClusters/data/sims'

        # Check for environ vars
        # gotta set this for the SExtractor configuration files to work
        if not os.getenv('PIPE'):
            os.environ['PIPE'] = os.path.join(self.pipeline)

        # Get the z-magnitude relation -- only once
        self.get_zmag()

        # This are the values for SOAR
        self.NX = 1890
        self.NY = 1890

        return

    # Magnitude loop
    def mag_loop(self, m1, m2, dm=0.2, rh=3.0, Lstar=0.5, N=1):

        # The range in magnitudes to follow
        mx = numpy.arange(m1, m2, dm)
        for mag in mx:
            print("#\n# Starting loop for m=%.2f" % mag)
            self.iter_loop(mag, rh=rh, Lstar=Lstar, N=N)

        return

    # Iter loop for a give magnitude
    def iter_loop(self, mag, rh=3.0, Lstar=0.5, N=1):

        # Make the list of artificial galaxies, same for all fields
        self.rh = rh
        self.mag = mag

        # Loop repetirions
        for i in range(N):

            t1 = time.time()
            # pass up the iterarion number
            self.iter = i + 1
            self.make_gallist(mag, rh=rh, Lstar=Lstar, sseed=i + 1)

            # For every field do the magic
            for field in self.fields:
                # Create the fake image
                self.run_mkobject(field)
                # Run SExtractor on it
                self.SEx(field)
                # Match and store the catalog
                self.match_cats(field)

            # Merge them all in one file
            self.merge_matched()
            print("# %s " % extras.elapsed_time_str(t1))
            print("# ------\n")

        return

    # Merge the matched catalogs
    def merge_matched(self):

        # Set up names
        iter = format_iter(self.iter)
        MergedCat = os.path.join(self.outpath, "Catalogs",
                                 "iter_m%.2f_rh%.1f_%s.dat" %
                                 (self.mag, self.rh, iter))

        # Clean if exists
        if os.path.exists(MergedCat):
            os.remove(MergedCat)

        # Merge all field
        for field in self.fields:
            matchcat = os.path.join(self.outpath, "Catalogs", "%s.mch" % field)
            cmd = "cat %s >> %s" % (matchcat, MergedCat)
            os.system(cmd)

        print("# Merged results on %s:" % MergedCat)
        return

    # Match the sim and recovered catalogs
    def match_cats(self, field, dmax=4.0):

        matchcat = os.path.join(self.outpath, "Catalogs", "%s.mch" % field)
        mch = open(matchcat, "w")
        (x1, y1, m1) = tableio.get_data(self.GaList, cols=(0, 1, 2))
        (x2, y2, m2) = tableio.get_data(self.SExCat, cols=(1, 2, 3))
        idx = []

        ## Find matches for GaList on SEx one
        for i in range(len(x1)):
            d = numpy.sqrt((x1[i] - x2)**2 + (y1[i] - y2)**2)
            if d.min() < dmax:
                ix = numpy.argmin(d)
                idx.append(ix)
                mch.write("%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n" %
                          (x1[i], y1[i], m1[i], x2[ix], y2[ix], m2[ix], d[ix]))
        print("# Matched cat  on %s" % matchcat)
        idx = numpy.array(idx)
        print("# Matched %s/%s" % (len(idx), self.Ngal))
        return

    # Run SExtractor on the field
    def SEx(self, field):

        # The configuration file
        SExinpar = os.path.join(self.pipeline, 'confs', 'bcs_Catalog_SIM.inpar')

        #header = self.header
        filter = self.filter[-1]
        simu_fits = os.path.join(self.outpath, "Images", "%s%s_sim.fits" %
                                 (field, filter))
        out_cat = os.path.join(self.outpath, "Catalogs", "%s%s_sim.cat" %
                               (field, filter))
        zp_use = self.header['MAGZERO']

        opts = ''
        # Do the SEx
        print("# Will run SEx on %s" % field)
        cmd = "sex %s -CATALOG_NAME %s -MAG_ZEROPOINT %s -c %s %s > /dev/null" % (
            simu_fits, out_cat, zp_use, SExinpar, opts)
        print(cmd)
        os.system(cmd)

        # pass them up
        self.SExCat = out_cat
        return

    # Run artdata's mkobjects
    def run_mkobject(self, field):

        # Real and simulated data
        filter = self.filter[-1]
        real_fits = os.path.join(self.datapath, field, "%s%s.fits" %
                                 (field, filter))
        simu_fits = os.path.join(self.outpath, "Images", "%s%s_sim.fits" %
                                 (field, filter))

        # get header info
        self.header = pyfits.getheader(real_fits)

        # N exp and instrument values for SOAR/SOI
        Nexp = {}
        Nexp['g'] = 4
        Nexp['r'] = 3
        Nexp['i'] = 4
        seeing = 0.5 / self.pixscale
        #gain = 2.0  # e/ADU SOAR
        gain = self.header['GAIN']  # e/ADU
        #rdnoise = 4.4  # e SOAR
        rdnoise = 5.075  # e MOSAICII

        print("# Creating: %s" % simu_fits)
        print("# Based on: %s" % real_fits)

        # The header info
        #exptime = self.header['EXPTIME']
        exptime = 1
        zeropt = self.header['MAGZERO']
        # Default background artdata.mkobject.in ADU)
        artdata.mkobject.background = 0.0
        artdata.mkobject.title = "%s_%ssim" % (field, filter)  # Image title
        artdata.mkobject.objects = self.GaList
        artdata.mkobject.xoffset = 0.0
        artdata.mkobject.yoffset = 0.0
        artdata.mkobject.radius = seeing  # Seeing radius/scale artdata.mkobject)
        artdata.mkobject.exptime = exptime  # Exposure time
        artdata.mkobject.magzero = zeropt  # Magnitude zero point
        artdata.mkobject.gain = gain  # Gain artdata.mkobject.electrons/ADU
        artdata.mkobject.rdnoise = rdnoise  # Read noise artdata.mkobject.electrons
        artdata.mkobject.poisson = "no"  # Add Poisson noise?
        artdata.mkobject.seed = 1  # Random number seed
        artdata.mkobject.comments = 'yes'  # Add comments to image?

        # Run mkobjects
        print("# Creating list: %s" % self.GaList)
        artdata.mkobjects(input=real_fits,
                          output=simu_fits,
                          objects=self.GaList)
        return

    # Make the galaxy list for a given size, mag and Luminosity
    def make_gallist(self, m, rh=3.0, Lstar=0.5, sseed=1.0):

        self.GaList = '/tmp/gal_ell_m%.2f_rh%.1f.dist' % (m, rh)
        if os.path.exists(self.GaList):
            print("# Cleaning %s" % self.GaList)
            os.remove(self.GaList)

        # The size for the magnitude
        esize = self.size(m, rh=rh, Lstar=Lstar)

        print(esize)

        #esize=5
        # Transform into pixels
        esize = esize / self.pixscale

        #NX = self.NX - 40
        #NY = self.NY - 40
        NX = 512
        NY = 512

        # SPATIAL DISTRIBUTION
        artdata.gallist.interactive = "no"  # Interactive mode?
        # Spatial density function (uniform|hubble|file)
        artdata.gallist.spatial = "uniform"
        artdata.gallist.xmin = 40 # Minimum x coaordinate value
        artdata.gallist.xmax = NX  # Maximum x coordinate value
        artdata.gallist.ymin = 40  # Minimum y coordinate value
        artdata.gallist.ymax = NY  # Maximum y coordinate value
        # Seed for sampling the spatial probability function
        artdata.gallist.sseed = sseed

        # MAGNITUDE DISTRIBUTION
        # Luminosity function artdata.gallist.uniform|powlaw|schecter|file
        artdata.gallist.luminosity = "uniform"
        artdata.gallist.minmag = m  # Minimum magnitude
        artdata.gallist.maxmag = m  # Maximum magnitude

        # MORPHOLOGY DISTRIBUTION
        artdata.gallist.egalmix = 1.0  # Percentage of elliptical galaxies
        artdata.gallist.ar = 0.3  # Minimum elliptical galaxy axial ratio
        artdata.gallist.eradius = esize  # Maximum elliptical half flux radius
        artdata.gallist.sradius = 1.2  # Spiral/ellipical radius at same magnitude
        # Absorption in edge on spirals artdata.gallist.mag
        artdata.gallist.absorption = 1.2
        artdata.gallist.z = 0.1  # Minimum redshift

        artdata.gallist(self.GaList, self.Ngal)
        print("# Created gallist on: %s" % self.GaList)

        # Pass up
        self.rh = rh
        self.mag = m

        return

    # The redshift magnitude relation --  we do this only once
    def get_zmag(self):
        zx = numpy.arange(0.01, 10.0, 0.005)
        mstar = m_star(zx, filter_out=self.filter)
        # The z(mag) interpolated function for Mstar galaxy
        self.zm = scipy.interpolate.interp1d(mstar, zx)
        return

    # The redshift for a given magnitude, contained in size(m) function
    def zmag(self, m, Lstar=1.0):
        mx = m - 2.5 * math.log10(Lstar)
        return self.zm(mx)

    # The size[arc] magnitude relation for a give half-light radius rh in kpc
    def size(self, m, rh=3.0, Lstar=0.5):
        mx = m - 2.5 * math.log10(Lstar)
        print(mx)
        zx = self.zm(mx)
        print(zx)
        print(self.cosmo)
        print(rh)
        size = astrometry.kpc2arc(zx, rh, self.cosmo)
        print(size)
        return size


# observed mi_star as a function of redshift
def m_star(z, filter_out='r_SDSS', h=0.7, cosmo=(0.3, 0.7, 0.7)):
    dlum = cosmology.dl(z, cosmology=cosmo)
    # Blanton M*
    Mi_star = -21.22 - 5 * math.log10(h)  # + self.evf['i'](z)[0]
    M_star = cosmology.reobs('El_Benitez2003',
                             m=Mi_star,
                             oldfilter="SLOAN-SDSS.i",
                             newfilter=filter_out)
    return M_star + 5.0 * numpy.log10(dlum) + 25


# Simple match between catalogs
def simple_match(GaList, SExCat, dmax=5):

    matchcat = "/tmp/tmp.mch"
    mch = open(matchcat, "w")
    (x1, y1, m1) = tableio.get_data(GaList, cols=(0, 1, 2))
    (x2, y2, m2) = tableio.get_data(SExCat, cols=(1, 2, 3))
    idx = []
    ## Find matches for GaList on SEx one
    mch.write("# %6s %8s %8s %8s %8s %8s %8s\n" % (
        'x_sim', 'y_sim', 'mag_sim', 'x_SEx', 'y_SEx', 'MAG_AUTO', 'd[pix]'))
    for i in range(len(x1)):
        d = numpy.sqrt((x1[i] - x2)**2 + (y1[i] - y2)**2)
        if d.min() < dmax:
            ix = numpy.argmin(d)
            idx.append(ix)
            mch.write("%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n" %
                      (x1[i], y1[i], m1[i], x2[ix], y2[ix], m2[ix], d[ix]))
    print("Matched cat  on %s" % matchcat)
    idx = numpy.array(idx)
    print(len(idx))
    return


# Simple format
def format_iter(n):

    if n < 10:
        return "00%1d" % n
    elif n < 100:
        return "0%2d" % n
    else:
        return "%3d" % n


def main():

    t0 = time.time()

    # The SOAR fields to be used
    fields = ['PSZ1_G031.91+67.94', ]

    # Initialize the function
    m1 = 23.5
    m2 = 25.5
    dm = 0.25
    Lstar = 0.5
    rh = 0.5  # kpc

    # Initialize the class
    sim = simgal(fields)

    # Do the mag loop m1, m2
    sim.mag_loop(m1, m2, dm=dm, rh=rh, Lstar=Lstar, N=1)
    print(extras.elapsed_time_str(t0))
    return


if __name__ == "__main__":
    main()
