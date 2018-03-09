import math
import numpy as np
from astropy.io import fits
from pyraf import iraf
import os

# From Henry G. 10/2014

# mkobjects_iraf(): This is the original iraf task to plop down fake galaxies
# onto a real image.


def mkobjects_iraf(inimage, outimage, objfile, radius, magzero):

    from pyraf import iraf
    iraf.noao(_doprint=0)
    iraf.artdata(_doprint=0)

    if os.path.exists(
            outimage) == True:  # If the file exists already, IRAF doesn't tell you but just doesn't write the new file
        os.remove(outimage)
        print 'Deleted ' + outimage + ' first'

    iraf.mkobjects(input=inimage,
                   output=outimage,
                   objects=objfile,
                   radius=radius,
                   magzero=magzero)

    return None

# mkobjects_homegrown(): This is the main function. It is a replacement for the
# IRAF task "mkobjects()". It adds a bunch of fake galaxies to the input image.
#
#   - 'infilename' is the input image
#   - 'outfilename' is the output image
#   - 'objectlistfilename' is a file with a list of objects to create
#   - 'radius' is the minimum size of a fake galaxy image, given by the seeing
#   - 'magzero' is how deep the image is
#
# This function calls 'make_galaxies()', which can have many different
# implementations.


def mkobjects_homegrown(infilename, outfilename, objectlistfilename, radius,
                        magzero):
    infile = fits.open(infilename)
    image = infile[0].data
    gain = 1.0  #infile[0].header['GAIN']
    infile.close()
    xsize, ysize = image.shape

    # zero flux corresponding to magzero in ADU:
    zeroflux = 1.0 / gain

    x, y, mag, size = np.loadtxt(objectlistfilename,
                                 usecols=(0, 1, 2, 4),
                                 unpack=True)
    flux = zeroflux * 10.0**(0.4 * (magzero - mag))
    galsize = np.sqrt(size**2 + radius**2)
    x = np.copy(x)  # np.loadtxt() creates a view, but we need contiguous array
    y = np.copy(y)

    make_galaxies(image, flux, galsize, x, y)

    hdu = fits.PrimaryHDU(image)
    hdu.header['GAIN'] = (gain, 'the gain')
    hdulist = fits.HDUList([hdu])
    if os.path.exists(outfilename) == True:
        os.remove(outfilename)
        print 'Deleted ' + outfilename + ' first'  # This line is necessary else code throws error
    hdulist.writeto(outfilename)

    return None


# plain python implementation:
def make_galaxies_plain(image, flux, galsize, x, y):
    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            for f, s, xi, yi in zip(flux, galsize, x, y):
                R = math.sqrt((xi - i)**2 + (yi - j)**2)
                val = f * math.exp(-1.6783 * R / s)
                image[j, i] += val


# numpy implementation:
def make_galaxies_numpy(image, flux, galsize, x, y):
    """ Adds a fake object to 'image'. The object's flux is 'flux', and its
    size 'galsize' and positions 'x' and 'y' are given in pixels. """

    # in order to compare to the other implementations, run over the entire
    # image:
    ymax = image.shape[0]

    ii = np.arange(0, ymax, 1, dtype=int)
    for j in range(0, ymax):
        for f, s, xi, yi in zip([flux], [galsize], [x], [y]):
            R = np.sqrt((xi - ii)**2 + (yi - j)**2)
            val = f * np.exp(-1.6783 * R / s)
            image[j, :] += val

    return None


# numpy implementation:
def make_galaxies_numpy_fast(image, flux, galsize, x, y):
    """ Adds a fake object to 'image'. The object's flux is 'flux', and its
    size 'galsize' and positions 'x' and 'y' are given in pixels. """

    for f, s, xi, yi in zip([flux], [galsize], [x], [y]):
        # Select region to add object to.
        # Choosing the entire image is too inefficient, so choose the size that
        # is just large enough that the object's flux is 10% below the
        # zeropoint.
        fluxlim = 0.0001 * f  #0.1
        radius = np.log(fluxlim / f) * s / (-1.6783)

        xmin = int(np.floor(xi - radius))
        xmax = int(np.ceil(xi + radius))
        xmin = max(0, xmin)
        xmax = min(xmax, image.shape[1])

        ymin = int(np.floor(yi - radius))
        ymax = int(np.ceil(yi + radius))
        ymin = max(0, ymin)
        ymax = min(ymax, image.shape[0])

        ii = np.arange(xmin, xmax, 1, dtype=int)
        for j in range(ymin, ymax):
            R = np.sqrt((xi - ii)**2 + (yi - j)**2)
            val = f * np.exp(-1.6783 * R / s)
            image[j, xmin:xmax] += val

    return None

# configure which implementation to use:
#make_galaxies = make_galaxies_plain
make_galaxies = make_galaxies_numpy
#make_galaxies = make_galaxies_ufunc
#make_galaxies = make_galaxies_opencl
#make_galaxies = make_galaxies_numpy_fast

#mkobjects = mkobjects_homegrown
mkobjects = mkobjects_iraf

# some things need initialization:
if mkobjects == mkobjects_iraf:
    from pyraf import iraf
    iraf.noao(_doprint=0)
    iraf.artdata(_doprint=0)
