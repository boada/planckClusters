#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
from builtins import range
from builtins import object
import os
import time
import extras
import tableio
from astropy.io import fits as pyfits
from astropy.modeling.models import Sersic1D, Sersic2D
from astropy.modeling.models import Gaussian1D, Gaussian2D
import numpy
import multiprocessing
import traceback
import functools
from pyraf import iraf
from iraf import artdata
iraf.noao(_doprint=0)
iraf.artdata(_doprint=0)


def error(msg, *args):
    multiprocessing.log_to_stderr()
    return multiprocessing.get_logger().error(msg, *args)


def trace_unhandled_exceptions(func):
    @functools.wraps(func)
    def wrapped_func(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except Exception as e:
            error(traceback.format_exc())
            raise

    return wrapped_func


class AsyncFactory:
    def __init__(self, func, cb_func):
        self.func = func
        self.cb_func = cb_func
        self.pool = multiprocessing.Pool(maxtasksperchild=10,
                                         processes=multiprocessing.cpu_count())

    def call(self, *args, **kwargs):
        self.pool.apply_async(self.func, args, kwargs, self.cb_func)

    def wait(self):
        self.pool.close()
        self.pool.join()


class simgal(object):
    def __init__(self, filter, m1, m2, dm=0.2, Ngal=100, N=4, pixscale=0.25):

        self.filter = filter
        self.Ngal = Ngal
        self.pixscale = pixscale

        self.pipeline = '/home/boada/Projects/planckClusters/MOSAICpipe'
        # The input datapath
        self.datapath = '/home/boada/Projects/planckClusters/data/proc2'
        # The output datapath
        self.outpath = '/home/boada/Projects/planckClusters/data/sims'

        # Check for environ vars
        # gotta set this for the SExtractor configuration files to work
        if not os.getenv('PIPE'):
            os.environ['PIPE'] = os.path.join(self.pipeline)

        mx = numpy.arange(m1, m2, dm)
        for mag in mx:
            for i in range(N):
                self.iter = i + 1
                self.make_gallist(mag, sseed=i + 1)

        return

    # Magnitude loop
    def mag_loop(self, field, m1, m2, dm=0.2, N=1):

        # The range in magnitudes to follow
        mx = numpy.arange(m1, m2, dm)
        for mag in mx:
            print("#\n# Starting loop for m=%.2f" % mag)
            self.iter_loop(field, mag, N=N)
        return

    # Iter loop for a give magnitude
    @trace_unhandled_exceptions
    def iter_loop(self, field, mag, N=1):

        # Make the list of artificial galaxies, same for all fields
        self.mag = mag

        # Loop repetirions
        for i in range(N):
            t1 = time.time()
            # pass up the iterarion number
            self.iter = i + 1
            iter = format_iter(self.iter)
            self.GaList = '/tmp/gal_ell_m%.2f_rh%.1f_%s.dist' % (mag, 3.0,
                                                                 iter)
            #self.make_gallist(mag, sseed=i + 1)

            # Create the fake image
            self.mkobjects_homegrown(field)
            # Run SExtractor on it
            self.SEx(field)
            # Match and store the catalog
            self.match_cats(field)

            # Merge them all in one file
            self.merge_matched(field)
            print("# %s " % extras.elapsed_time_str(t1))
            print("# ------\n")
        return

    # Make the galaxy list for a given size, mag and Luminosity
    @trace_unhandled_exceptions
    def make_gallist(self, m, sseed=1.0):

        iter = format_iter(self.iter)
        GaList = '/tmp/gal_ell_m%.2f_rh%.1f_%s.dist' % (m, 3.0, iter)
        if os.path.exists(GaList):
            print("# Cleaning %s" % GaList)
            os.remove(GaList)

        # SPATIAL DISTRIBUTION
        artdata.gallist.interactive = "no"  # Interactive mode?
        # Spatial density function (uniform|hubble|file)
        artdata.gallist.spatial = "uniform"
        artdata.gallist.xmin = 1050  # Minimum x coaordinate value
        artdata.gallist.xmax = 1950  # Maximum x coordinate value
        artdata.gallist.ymin = 1050  # Minimum y coordinate value
        artdata.gallist.ymax = 1950  # Maximum y coordinate value
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
        #artdata.gallist.eradius = esize  # Maximum elliptical half flux radius
        artdata.gallist.sradius = 1.2  # Spiral/ellipical radius at same magnitude
        # Absorption in edge on spirals artdata.gallist.mag
        artdata.gallist.absorption = 1.2
        artdata.gallist.z = 0.1  # Minimum redshift

        artdata.gallist(GaList, self.Ngal)
        print("# Created gallist on: %s" % GaList)

        # Pass up
        self.mag = m

        return

    @trace_unhandled_exceptions
    def mkobjects_homegrown(self, field):

        filter = self.filter
        real_fits = os.path.join(self.datapath, field,
                                 "%s%s.fits" % (field, filter))
        simu_fits = os.path.join(self.outpath, "Images",
                                 "%s%s_sim.fits" % (field, filter))

        outfilename = simu_fits
        objectlistfilename = self.GaList
        seeing = 0.5 / self.pixscale

        radius = seeing

        with pyfits.open(real_fits) as infile:
            self.header = infile[0].header
            image = infile[0].data
            gain = infile[0].header['GAIN']
            magzero = infile[0].header['MAGZERO']

        # only get part of the total image
        # this will speed up the sextractor run
        limit = 1000 # make a 2000x2000 final image
        image = image[image.shape[1] // 2 - limit: image.shape[1] // 2 + limit,
                    image.shape[0] // 2 - limit: image.shape[0] // 2 + limit]

        xsize, ysize = image.shape

        # zero flux corresponding to magzero in ADU:
        zeroflux = 1.0 / gain

        # reads the IRAF GALLIST output.
        x, y, mag, size, ar, pa = numpy.loadtxt(objectlistfilename,
                                                usecols=(0, 1, 2, 4, 5, 6),
                                                unpack=True)

        flux = zeroflux * 10.0**(0.4 * (magzero - mag))
        galsize = numpy.sqrt(size**2 + radius**2)
        #galsize = numpy.ones_like(size) * radius
        # numpy.loadtxt() creates a view, but we need contiguous array
        x = numpy.copy(x)
        y = numpy.copy(y)

        image = self.make_galaxies_astropy(image, flux, galsize, x, y, ar, pa)
        #image = self.make_PS_astropy(image, flux, galsize, x, y, pa)

        # write out the image
        hdu = pyfits.PrimaryHDU(image)
        hdu.header['GAIN'] = (gain, 'the gain')
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto(outfilename, overwrite=True)

        return None

    def make_galaxies_astropy(self, image, flux, galsize, x, y, ar, pa, n=4):

        for f, s, xi, yi, pai, ari, in zip(flux, galsize, x, y, pa, ar):

            fluxlim = 0.0001 * f  # 0.1
            scale = 1  # arcsec/pixel
            r_e = s  # effective radius
            ellip = ari  # ellipticity
            theta = numpy.deg2rad(pai)  # position angle
            x_cent = 0  # x centroid
            y_cent = 0  # x centroid
            tot_flux = f  # total flux

            s1 = Sersic1D(amplitude=1, r_eff=r_e, n=n)
            r = numpy.arange(0, 1000, scale)
            s1_n = s1(r) / sum(s1(r))
            extent = numpy.where(s1_n * f > fluxlim)[0].max()

            if extent % 2 > 0:
                extent += 1

            ser_model = Sersic2D(r_eff=r_e,
                                 n=n,
                                 ellip=ellip,
                                 theta=theta,
                                 x_0=x_cent,
                                 y_0=y_cent)

            x = numpy.arange(-extent / 1., extent / 1., scale) + x_cent / scale
            y = numpy.arange(-extent / 1., extent / 1., scale) + y_cent / scale

            X, Y = numpy.meshgrid(x, y)

            img = ser_model(X, Y)
            img /= numpy.sum(img)
            img *= tot_flux

            xi, yi = int(xi), int(yi)
            # COLUMNS FIRST -- because FITS are silly
            image[yi - img.shape[1] // 2: yi + img.shape[1] // 2,
                  xi - img.shape[0] // 2: xi + img.shape[0] // 2] += img

        return image

    def make_PS_astropy(self, image, flux, galsize, x, y, pa):

        for f, s, xi, yi, pai, in zip(flux, galsize, x, y, pa):

            fluxlim = 0.0001 * f  # 0.1
            scale = 1  # arcsec/pixel
            theta = numpy.deg2rad(pai)  # position angle
            tot_flux = f  # total flux
            x_cent = 0
            y_cent = 0

            g1 = Gaussian1D(amplitude=1, stddev=s)
            r = numpy.arange(0, 100, scale)
            g1_n = g1(r) / sum(g1(r))
            extent = numpy.where(g1_n * f > fluxlim)[0].max()

            if extent % 2 > 0:
                extent += 1

            gaus_model = Gaussian2D(amplitude=1,
                                    x_stddev=s,
                                    y_stddev=s,
                                    theta=theta)

            x = numpy.arange(-extent / 1, extent / 1, scale) + x_cent / scale
            y = numpy.arange(-extent / 1, extent / 1, scale) + y_cent / scale

            X, Y = numpy.meshgrid(x, y)

            img = gaus_model(X, Y)
            img /= numpy.sum(img)
            img *= tot_flux

            xi, yi = int(xi), int(yi)
            # COLUMNS FIRST -- because FITS are silly
            image[yi - img.shape[1] // 2: yi + img.shape[1] // 2,
                  xi - img.shape[0] // 2: xi + img.shape[0] // 2] += img

        return image

    # Merge the matched catalogs
    @trace_unhandled_exceptions
    def merge_matched(self, field):

        # Set up names
        iter = format_iter(self.iter)
        MergedCat = os.path.join(self.outpath, "Catalogs", "%s_m%.2f_%s_%s.dat"
                                 % (field, self.mag, self.filter[-1], iter))

        # Clean if exists
        if os.path.exists(MergedCat):
            os.remove(MergedCat)

        # Merge all field
        matchcat = os.path.join(self.outpath, "Catalogs",
                                "{}{}.mch".format(field, self.filter[-1]))
        cmd = "cat %s >> %s" % (matchcat, MergedCat)
        os.system(cmd)

        print("# Merged results on %s:" % MergedCat)
        return

    # Match the sim and recovered catalogs
    @trace_unhandled_exceptions
    def match_cats(self, field, dmax=4.0):

        matchcat = os.path.join(self.outpath, "Catalogs",
                                "{}{}.mch".format(field, self.filter[-1]))
        with open(matchcat, 'w') as mch:
            (x1, y1, m1) = tableio.get_data(self.GaList, cols=(0, 1, 2))
            (x2, y2, m2) = tableio.get_data(self.SExCat, cols=(1, 2, 3))
            idx = []

            ## Find matches for GaList on SEx one
            for i in range(len(x1)):
                d = numpy.sqrt((x1[i] - x2)**2 + (y1[i] - y2)**2)
                if d.min() < dmax:
                    ix = numpy.argmin(d)
                    idx.append(ix)
                    mch.write("%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n" % (
                        x1[i], y1[i], m1[i], x2[ix], y2[ix], m2[ix], d[ix]))
        print("# Matched cat  on %s" % matchcat)
        idx = numpy.array(idx)
        print("# Matched %s/%s" % (len(idx), self.Ngal))
        return

    # Run SExtractor on the field
    def SEx(self, field):

        # The configuration file
        SExinpar = os.path.join(self.pipeline, 'confs',
                                'bcs_Catalog_SIM.inpar')

        #header = self.header
        filter = self.filter[-1]
        simu_fits = os.path.join(self.outpath, "Images",
                                 "%s%s_sim.fits" % (field, filter))
        #simu_fits = 'gals.fits'
        out_cat = os.path.join(self.outpath, "Catalogs",
                               "%s%s_sim.cat" % (field, filter))
        zp_use = self.header['MAGZERO']

        opts = ''
        # Do the SEx
        iter = format_iter(self.iter)
        print("# Will run SEx on %s -- iter %s" % (field, iter))
        cmd = "sex %s -CATALOG_NAME %s -MAG_ZEROPOINT %s -c %s %s > /dev/null" % (
            simu_fits, out_cat, zp_use, SExinpar, opts)
        print(cmd)
        os.system(cmd)

        # pass them up
        self.SExCat = out_cat

        # remove simulated fits
        #os.remove(simu_fits)

        return


# Simple format
def format_iter(n):

    if n < 10:
        return "00%1d" % n
    elif n < 100:
        return "0%2d" % n
    else:
        return "%3d" % n


def cb_func(self):
    print("PID: %d completed" % (os.getpid()))


def main():
    from glob import glob

    t0 = time.time()

    # Initialize the function
    m1 = 20
    m2 = 26
    dm = 0.2
    Ngal = 100
    N = 4
    filters = ['g', 'r', 'i', 'z', 'K']

    for filt in filters:

        # The fields to be used
        data_dir = '/home/boada/Projects/planckClusters/data/proc2/'
        files = glob('{}PSZ*/PSZ*{}.fits'.format(data_dir, filt), recursive=True)
        fields = [f.split('/')[-2] for f in files]

        #fields = ['PSZ1_G031.91+67.94', ]
        #fields = ['PSZ2_G125.55+32.72', ]
        #fields = ['PSZ2_G043.44-41.27', 'PSZ2_G029.66-47.63', ]
        #fields = ['PSZ2_G189.79-37.25']
        #fields = ['PSZ1_G183.26+12.25']

        # Initialize the class
        sim = simgal(filter=filt, m1=m1, m2=m2, dm=dm, Ngal=Ngal, N=N)

        # start the async factory
        async_worker = AsyncFactory(sim.mag_loop, cb_func)

        for field in fields:
            # Do the mag loop m1, m2
            async_worker.call(field, m1, m2, dm=dm, N=N)
            #sim.mag_loop(field, m1, m2, dm=dm, N=4)
        async_worker.wait()

    print(extras.elapsed_time_str(t0))
    return


if __name__ == "__main__":
    main()
