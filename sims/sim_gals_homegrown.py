import numpy as np
import os
from astropy.io import fits
from random import randint

def mkobjects_iraf(inimage, outimage, objfile, radius, magzero):

    from pyraf import iraf
    iraf.noao(_doprint=0)
    iraf.artdata(_doprint=0)

    if os.path.exists(outimage):
        os.remove(outimage)

    iraf.mkobjects(input=inimage,
                   output=outimage,
                   objects=objfile,
                   radius=radius,
                   magzero=magzero)

    return None

def mkobjects_homegrown(infilename, outfilename, objectlistfilename, radius,
                        magzero):
    infile = fits.open(infilename)
    image = infile[0].data
    gain = 1.0  # infile[0].header['GAIN']
    infile.close()
    xsize, ysize = image.shape

    # zero flux corresponding to magzero in ADU:
    zeroflux = 1.0 / gain

    # reads the IRAF GALLIST output.
    x, y, mag, size, ar, pa = np.loadtxt(objectlistfilename,
                                 usecols=(0, 1, 2, 4, 5, 6),
                                 unpack=True)
    flux = zeroflux * 10.0**(0.4 * (magzero - mag))
    galsize = np.sqrt(size**2 + radius**2)
    x = np.copy(x)  # np.loadtxt() creates a view, but we need contiguous array
    y = np.copy(y)

    make_galaxies(image, flux, galsize, x, y, ar, pa)

    hdu = fits.PrimaryHDU(image)
    hdu.header['GAIN'] = (gain, 'the gain')
    hdulist = fits.HDUList([hdu])
    if os.path.exists(outfilename):
        os.remove(outfilename)
    hdulist.writeto(outfilename)

    return None

# numpy implementation:
def make_galaxies(image, flux, galsize, x, y, pa, ar, profile='devauc'):
    """ Adds a fake object to 'image'. The object's flux is 'flux', and its
    size 'galsize' and positions 'x' and 'y' are given in pixels. """

    for f, s, xi, yi, pai, ari, in zip(flux, galsize, x, y, pa, ar):
        # Select region to add object to.
        # Choosing the entire image is too inefficient, so choose the size that
        # is just large enough that the object's flux is 10% below the
        # zeropoint.
        fluxlim = 0.0001 * f  # 0.1

        if profile == 'devauc':
            radius = s * np.log(fluxlim / f) ** 4 * 7.67**-4
        elif profile == 'expdisk':
            radius = s * np.log(fluxlim / f) * -1.6783**-1

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

            # from iraf.mkobjects manual
            dx = ii - xi
            dy = j - yi
            #dX = dx * np.cos(pai) + dy * np.sin(pai)
            #dY = (-dx * np.sin(pai) + dy * np.cos(pai)) / ari
            #R = np.sqrt(dX ** 2 + dY ** 2)
            R = np.sqrt(dx ** 2 + dy ** 2)

            print(R)
            if profile == 'devauc':
                val = f * np.exp(-7.67 * (R / s)**(1 / 4))
            elif profile == 'expdisk':
                val = f * np.exp(-1.6783 * R / s)
            image[j, xmin:xmax] += val

    return None

def make_galaxies_list(outfilename, numgals=25, minmag=20, maxmag=20, size=3,
                       xlim=(1, 512), ylim=(1, 512)):
    with open(outfilename, 'wt') as outfile:
        for i in range(numgals):
            x = randint(xlim[0], xlim[1])
            y = randint(ylim[0], ylim[1])
            outfile.write('{}\t{}\t{}\t{}\n'.format(x, y, minmag, size))

