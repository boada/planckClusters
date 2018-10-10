#!/usr/bin/env python3

import os
import sys
import numpy as np
from math import log10
import time
import string


import types
import imp
# pipeline specific imports
from pipe_utils import tableio
from pipe_utils import deredden
from utils import elapsed_time

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


class combcat(metaclass=PluginMeta):
    ''' Combine, swarp and get catalogs '''

    plugindir = '/home/boada/Projects/planckClusters/MOSAICpipe/plugins'

    def __init__(self,
                 assocfile,
                 datapath='',
                 outpath='',
                 pixscale=0.2666,
                 dryrun=None,
                 noSWarp=False,
                 verb='yes'):

        self.assocfile = assocfile
        self.datapath = datapath
        self.tilename = os.path.basename(os.path.splitext(assocfile)[0])
        self.outpath = outpath
        self.dryrun = dryrun
        self.noSWarp = noSWarp
        self.verb = verb
        self.DetImage = None
        self.centered = None
        self.got_zeropt = False

        self.pipeline = '/home/boada/Projects/planckClusters/MOSAICpipe/'

        # Check for environ vars
        if not os.getenv('PIPE'):
            os.environ['PIPE'] = os.path.join(self.pipeline)

        if not os.getenv('BPZPATH'):
            os.environ['BPZPATH'] = os.path.join(self.pipeline + 'bpz-1.99.3')

        # Set the dust_directory
        os.environ['DUST_DIR'] = os.path.join(self.pipeline + 'LIB')

        # Set the pixel scale
        self.pixscale = pixscale

        # Read the association file
        self.read_assoc()

        # Initialize and set to zero the dust exctiction corrections
        self.XCorr = {}
        self.XCorrError = {}
        for filter in self.filters:
            self.XCorr[filter] = 0.0
            self.XCorrError[filter] = 0.0

        return


def cmdline():

    from optparse import OptionParser

    # Read in the command line options
    USAGE = '''usage:\t %prog <assoc> <in_datapath> <out_datapath> [options]
    i.e.: %prog ~/bcs-inventory/BCS0519-5448.assoc ~/PROC ~/BCS/PROC'''

    parser = OptionParser(usage=USAGE)

    parser.add_option(
        "--swarp",
        action="store_true",
        dest="SWarp",
        default=False,
        help="SWARP the mosaics")

    parser.add_option(
        "--swarpextras",
        action="store_true",
        dest="SWarpExtras",
        default=False,
        help="SWARP the mosaics to make multicolor images")

    parser.add_option(
        "--bpz",
        action="store_true",
        dest="BPZ",
        default=False,
        help="Calculate photometric redshifts")

    parser.add_option(
        "--sex", action="store_true", dest="SEx", default=False, help="No SEx")

    parser.add_option(
        "--dust",
        action="store_true",
        dest="Dust",
        default=False,
        help="Dust Extinction Correction")

    parser.add_option(
        "--WeightOnly",
        action="store_true",
        dest="WeightOnly",
        default=False,
        help="Create custom weights from images and exits")

    parser.add_option(
        "--Weight",
        action="store_true",
        dest="Weight",
        default=False,
        help="Create custom weights from images")

    parser.add_option(
        "--mask",
        action="store_true",
        dest="noMask",
        default=False,
        help="Do not generate the mask from weight")

    parser.add_option(
        "--useMask",
        action="store_true",
        dest="useMask",
        default=False,
        help="Use Trimming mask")

    parser.add_option(
        "--dryrun",
        action="store_true",
        dest="dryrun",
        default=False,
        help="Dry Run (only SExtractor remains)")

    parser.add_option(
        "--combtype",
        dest="combtype",
        default='MEDIAN',
        help="SWarp COMBINE TYPE")

    parser.add_option(
        "--noCleanUP",
        action='store_true',
        dest='noCleanUP',
        default=False,
        help='Whether or not to remove uncompressed files')

    parser.add_option(
        "--RGB",
        action='store_true',
        dest='RGB',
        default=False,
        help='Whether or not to create RGB images from mosaics')

    parser.add_option(
        "--astro",
        action='store_true',
        dest='Astro',
        default=False,
        help='Whether or not to astro calibrate the mosiacs.')

    parser.add_option(
        "--photo",
        action='store_true',
        dest='Photo',
        default=False,
        help='Whether or not to photo calibrate the mosaics.')

    parser.add_option(
        "--newfirm",
        action='store_true',
        dest='newfirm',
        default=False,
        help='Use the NEWFIRM imaging for BPZ and RGB.')

    parser.add_option(
        "--deblend",
        action='store_true',
        dest='deblend',
        default=False,
        help='Deblend mode in SExtractor.')

    (options, args) = parser.parse_args()

    if len(args) < 3:
        parser.error(USAGE + "\nMust supply at least one argument required")

    # Dry run turns off everything... almost
    if options.dryrun:
        options.SWarp = False
        options.Astro = False
        options.Photo = False
        options.SEx = False
        options.BPZ = False
        options.Mask = False
        options.Dust = True

    # Same for WeightOnly
    if options.WeightOnly:
        options.SWarp = False
        options.Astro = False
        options.Photo = False
        options.SEx = False
        options.BPZ = False
        options.Mask = False
        options.Dust = True

    return options, args


def main():

    # The start time
    tstart = time.time()

    # Get the command line options
    opt, arg = cmdline()

    assocfile = arg[0]
    inpath = arg[1]
    outpath = arg[2]

    # SWarp
    if not opt.SWarp:
        # Init the class
        c = combcat(
            assocfile,
            datapath=inpath,
            outpath=outpath,
            verb='yes',
            dryrun=opt.dryrun,
            noSWarp=not opt.SWarp)
        c.get_filenames()
    else:
        # Init the class
        c = combcat(
            assocfile,
            datapath=inpath,
            outpath=outpath,
            verb='yes',
            dryrun=opt.dryrun,
            noSWarp=not opt.SWarp)

        c.swarp_files(
            dryrun=not opt.SWarp,
            conf="SWarp-common.conf",
            combtype=opt.combtype)

    if opt.SWarpExtras:
        # Init the class
        c = combcat(
            assocfile,
            datapath=inpath,
            outpath=outpath,
            verb='yes',
            dryrun=opt.dryrun,
            noSWarp=not opt.SWarp)

        c.get_filenames()

        c.swarp_extras(
            dryrun=not opt.SWarpExtras,
            conf="SWarp-common.conf",
            combtype=opt.combtype,
            newfirm=opt.newfirm)

    if opt.Astro:
        c.get_astrometry(newfirm=opt.newfirm)

    if opt.Photo:
        c.get_zeropt(newfirm=opt.newfirm)

    # Make the detection image
    if opt.useMask:
        # Make the mask from the weights
        c.generate_masks(filters=('i', ), dryrun=opt.noMask)
        c.makeDetectionIma(filter='i')

    # SExtractor
    if opt.SEx:
        c.SEx(deblend=opt.deblend)

    # Dust Extinction Correction
    if opt.Dust:
        c.DustCorrection()

    # photometric redshifts
    if opt.BPZ:
        c.BuildColorCat(newfirm=opt.newfirm)
        c.runBPZ()

    # make RGB images (pngs)
    if opt.RGB:
        print('make rgb')
        c.make_RGB(newfirm=opt.newfirm)

    # cleanup
    if opt.noCleanUP or not opt.SWarp:
        if opt.noCleanUP or opt.SWarpExtras:
            pass
        else:
            print("CLEANUP!")
            c.cleanup_files()
    else:
        print("CLEANUP!")
        c.cleanup_files()

    elapsed_time(tstart, c.tilename)
    return


if __name__ == "__main__":
    main()
