#!/usr/bin/env python

import os
import sys
from cosmopy import cosmopy
import imp
import types

# float32 = numpy.float32
# float64 = numpy.float64
# int16 = numpy.int16
# nstr = numpy.char
# lge = numpy.greater_equal
# lle = numpy.less_equal

sout = sys.stderr

class PluginMeta(type):
    def __new__(cls, name, bases, dct):
        modules = [
            imp.load_source(filename, os.path.join(dct['plugindir'], filename))
            for filename in os.listdir(dct['plugindir'])
            if filename.endswith('.py')
        ]
        for module in modules:
            for name in dir(module):
                function = getattr(module, name)
                if isinstance(function, types.FunctionType):
                    dct[function.__name__] = function
        return type.__new__(cls, name, bases, dct)

class finder(metaclass=PluginMeta):

    plugindir = '/home/boada/Projects/planckClusters/CLUSTERpipe/plugins'

    def __init__(
            self,
            ctile,
            maglim=25.0,
            starlim=0.95,
            pixscale=0.25,
            zlim=1.8,
            zo=None,
            dz=0.05,
            radius=1000.0,  # Radius in kpc
            cosmo=(0.3, 0.7, 0.7),
            zuse="ZB",  # Use ZB (Bayesian) or ML (Max Like)
            outpath='plots',
            path='./',
            evolfile="0_1gyr_hr_m62_salp_Ks.color",
            p_lim=0.4,
            verb='yes'):

        # Check for environ vars
        self.home = os.environ['HOME']
        if not os.getenv('MOSAICpipe'):
            os.environ['MOSAICpipe'] = os.path.join(self.home, 'Projects',
                                                    'planckClusters',
                                                    'MOSAICpipe')
        self.MOSAICpipe = os.getenv('MOSAICpipe')

        self.zlim = zlim
        self.cosmo = cosmo
        self.evolfile = os.path.join(self.MOSAICpipe, "LIB/evol", evolfile)
        self.dz = dz
        self.zuse = zuse
        self.outpath = outpath
        self.path = path
        self.radius = radius
        self.zo = zo  # Input zo for the BCG

        # Set the cosmology now
        self.cset = cosmopy.set(self.cosmo)
        self.Om = cosmo[0]
        self.OL = cosmo[1]
        self.h = cosmo[2]
        self.Ho = self.h * 100.0

        self.ctile = ctile
        self.datapath = path
        self.catsfile = os.path.join(path, ctile, ctile + "_merged.cat")
        self.probsfile = os.path.join(path, ctile, ctile + "_probs.dat")

        # limits
        self.maglim = maglim
        self.starlim = starlim

        self.verb = verb
        self.ellipse = {}
        self.plot_around = None
        self.pixscale = pixscale

        self.read_cat()  # Read catalogs avoding, faint, high-z and 99 objects
        self.read_probs()  # Read probs function of objects from catalogs
        self.get_absmags()  # We compute Abs Mags for each object

        # Set the BCG masked to False, so we select BCGs on the firts run
        self.BCG_masked = False
        self.BCG_probs = False

        # Get the BCG candidates
        self.get_BCG_candidates(p_lim=p_lim)

        # Check if the destination folder exists
        if os.path.exists(self.outpath):
            sout.write("# Will put files to: %s\n" % self.outpath)
        else:
            sout.write("# Will create new folder: %s\n" % self.outpath)
            os.mkdir(self.outpath)

        self.rootname = os.path.join(self.datapath, self.ctile, self.ctile)

        return

def cmdline():

    from optparse import OptionParser

    # Read in the command line options
    USAGE = "usage:\t %prog <tilename> [options] \n"
    USAGE += "i.e.: %prog   --path ./PSZ2_G137.24+53.93/mosaic3/resampled/"
    parser = OptionParser(usage=USAGE)

    parser.add_option("--path", dest="path", default='./', help="Path to data")
    parser.add_option(
        "--radius", dest="radius", default=1000.0, help="Radius in kpc")
    parser.add_option("--zo", dest="zo", default=None, help="zo of cluster")
    parser.add_option("--dz", dest="dz", default=0.08, help="dz of shell")
    parser.add_option("--zuse", dest="zuse", default="ZB", help="use ZB or ML")
    parser.add_option("--dx", dest="dx", default=-1, help="cutout width")
    parser.add_option("--dy", dest="dy", default=-1, help="cutout heigth")
    parser.add_option(
        "--RA",
        dest="RA",
        default=None,
        help="Center Right Ascension - x pixel or SEX RA")
    parser.add_option(
        "--DEC",
        dest="DEC",
        default=None,
        help="Center Declination - y pixel or SEX DEC")

    # stellarity
    parser.add_option("--starlim", dest="starlim", default=0.95,
                      help="SExtractor star/gal limit: range 0:1")

    # add a bit to figure out Mosaic1/mosaic3
    parser.add_option(
        "--pixelscale",
        dest='pixelscale',
        default=-1,
        help='pixel scale of the instrument used MOS1:0.2666 or MOS3:0.25')

    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.error("Must supply at least one arguments required")

    return options, args

def main():

    opt, arg = cmdline()

    #print(opt)
    #print(arg)

    ctile = arg[0]
    radius = float(opt.radius)
    opt.dx = float(opt.dx)
    opt.dy = float(opt.dy)
    if opt.zo:
        zo = float(opt.zo)
    else:
        zo = None

    if float(opt.pixelscale) < 0:
        print('--pixelscale not defined.... assuming 0.25')
        pixelscale = 0.25
    else:
        pixelscale = float(opt.pixelscale)

    if opt.dx < 0 and opt.dy < 0:
        opt.dx = 5.1 * 60 / pixelscale
        opt.dy = 5.1 * 60 / pixelscale

    if not os.path.isfile('{}/{}irg.tiff'.format(opt.path, ctile)):
        print('{}/{}irg.tiff'.format(opt.path, ctile))
        print('RGB image does not exist -- probably only kband data')
        sys.exit()

    f = finder(
        ctile,
        maglim=26.0,
        starlim=float(opt.starlim),
        pixscale=pixelscale,
        zlim=1.8,
        zo=zo,
        dz=float(opt.dz),
        radius=radius,
        cosmo=(0.3, 0.7, 0.7),
        #zuse="ZB", # Use ZB (Bayesian) or ML (Max Like)
        zuse=opt.zuse,  # Use ZB (Bayesian) or ML (Max Like)
        outpath='plots',
        path=opt.path,
        evolfile="0_1gyr_hr_m62_salp_Ks.color",
        p_lim=0.4,
        verb='yes')

    f.jpg_read(dx=opt.dx, dy=opt.dy, RA=opt.RA, DEC=opt.DEC)
    f.jpg_display()

    return

if __name__ == '__main__':
    main()
