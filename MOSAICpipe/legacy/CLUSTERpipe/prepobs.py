#!/usr/bin/env python

from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
import os, sys, time
import string, glob
import extras
import pyfits
from pyraf import iraf
from iraf import mscred
from iraf import images


class DataSet(object):
    ''' Define datasets, filters, bias, zero and object images'''

    #def __init__(self,obsnight,verb=None):
    def __init__(self,
                 rawnight,
                 procdir=None,
                 verb=None,
                 skyfilters=('i', 'z'),
                 Namps=4):

        self.modName = string.split(
            string.split(str(self))[0], '.')[0][1:]  # this is the module name
        self.verb = verb
        self.cwd = os.getcwd()
        self.tstart = time.time()
        self.Namps = Namps

        # Proc and Raw dirs
        self.rawdir = rawnight
        self.procdir = procdir

        # The filters for which we will generate Sflats
        self.skyfilters = skyfilters

        # Check for environ vars
        if not os.getenv('BCSPIPE'):
            os.environ['BCSPIPE'] = '/Users/felipe/BCSPIPE.coligue'
        self.BCSPIPE = os.getenv('BCSPIPE')
        print(self.BCSPIPE)

        # Add the the iraf env bcspipe
        iraf.set(bcspipe=self.BCSPIPE)

        return

    def copy_files(self):

        # Check if we need to copy files
        if self.procdir == None:
            print(" Will not copy files\n, Will analyze files on %s" %
                  self.rawdir,
                  file=sys.stderr)
            # The name of the reduce night
            self.obsnight = self.rawdir
            return

        # The name of the reduced night
        obsnight = os.path.basename(
            self.rawdir)  # the last part of /xxx/xxx/bcs051118/
        self.obsnight = os.path.join(self.procdir, obsnight)

        # Check if the destination folder exists
        if os.path.exists(self.procdir):
            print("Will copy files to: %s" % self.rawdir, file=sys.stderr)
        else:
            print("Will create new folder: %s" % self.rawdir, file=sys.stderr)
            os.mkdir(self.procdir)

        # And now we rsync the files
        print("Coping files from %s --> %s" % (self.rawdir, self.procdir),
              file=sys.stderr)
        cmd = "rsync -e ssh --delete --progress --exclude=\"*junk*\" -navL %s %s " % (
            self.rawdir, self.procdir)
        os.system(cmd)

        # Now for real
        cmd = "rsync -e ssh --delete --progress --exclude=\"*junk*\" -avL %s %s " % (
            self.rawdir, self.procdir)
        os.system(cmd)

        return

    def inventory(self, pattern="*.fits"):

        # Change to the nights directoty
        os.chdir(self.obsnight)
        print(" Changed from directory %s --> %s" % (self.cwd, self.obsnight),
              file=sys.stderr)
        print(" We are now in: %s " % os.getcwd(), file=sys.stderr)

        # Do the same in Iraf
        iraf.chdir(self.obsnight)
        print(" We are now in Iraf's: ", file=sys.stderr)
        iraf.pwd()

        # Get the all the fits files names and sort them out
        self.fitslist = glob.glob(
            pattern)  # Keep it relative, no absolute paths

        print(" Found %s fits files in %s " %
              (len(self.fitslist), self.obsnight),
              file=sys.stderr)
        print(" Will make data inventory now... this might take a while",
              file=sys.stderr)

        # Get the filternames
        self.getfilters()  # Returns self.filters

        # Make the images hash
        self.zeros = []
        self.dflats = {}
        self.object = {}
        self.forsky = {}  # BCS images only!
        for filter in self.filters:
            self.dflats[filter] = []
            self.object[filter] = []
            self.forsky[filter] = []

        # For through all files
        i = 0
        for file in self.fitslist:

            try:
                header = pyfits.getheader(file)
                filter1 = header['FILTER1']
                filter2 = header['FILTER2']
                obstype = header['OBSTYPE']
                exptime = header['EXPTIME']
                Namps = header['NAMPS']
            except:
                continue

            if filter1.split()[1] != 'Open':
                filter = filter1.split()[1]
            else:
                filter = filter2.split()[1]
            print("Found %s -- %s" % (filter, file))

            #try:
            #    filter = filter.split()[1]
            #except:
            #    filter = filter.split()[0]

            if i == 0:
                Namp0 = Namps

            if Namps != Namp0:
                print("ERROR: Namps mixed, %s - %s, image:%s" %
                      (Namo0, Namps, file))

            # Zero frames can have any filter they want
            if obstype == 'ZERO':
                self.zeros.append(file)
            # In case there a spurious filter
            if filter not in self.filters:
                continue
            elif exptime is None:
                print(" Skipping %s, no EXPTIME key" % file, file=sys.stderr)
                continue
            elif Namps is None:
                print(" Skipping %s, no NAMPS key" % file, file=sys.stderr)
                continue
            elif obstype == 'DFLAT' or obstype == 'FLAT':
                self.dflats[filter].append(file)
            elif obstype == 'OBJECT':
                self.object[filter].append(file)
                if float(exptime) >= 80:  # Select the program files only
                    self.forsky[filter].append(file)

                    # tweak the fiter name
            if filter in self.filters:
                tweak_filter(file)
            i = i + 1

        self.Namps = Namp0
        print(" Image inventory ready", file=sys.stderr)

    def setNames(self):

        # Make the image lists and files lists
        self.Zero = "Zero"
        self.Dflat = "Dflat"
        self.Sflat = "Sflat"
        self.Fringe = os.path.join(self.BCSPIPE, 'LIB/SOI/Fringe/Fringe')

        # Zero
        self.zerolist = extras.imfile(self.zeros)
        self.zerofiles = self.zeros

        # Dome Flats and science files
        self.flatfiles = []
        self.scifiles = []
        self.bckfiles = [
        ]  # Files used for background, with exptime > 100 sec!
        for filter in self.filters:
            self.flatfiles = self.flatfiles + self.dflats[filter]
            self.scifiles = self.scifiles + self.object[filter]
            self.bckfiles = self.bckfiles + self.forsky[filter]

        # Dome flats, science and bad-pixel lists, etc...
        #self.flatlist = extras.imlist(self.flatfiles)
        #self.scilist  = extras.imlist(self.scifiles)
        #self.msklist  = extras.imlist_append(self.scilist,"_msk")
        #self.skylist  = extras.imlist_append(self.scilist,"_sky")

        # Try with files instead of lists
        self.flatlist = extras.imfile(self.flatfiles)
        self.scilist = extras.imfile(self.scifiles)

        self.objlist = {}
        self.nfiles = {}
        self.obsfilters = []  # Filters that were actually observed
        for filter in self.filters:
            self.objlist[filter] = extras.imfile(self.object[filter])
            self.nfiles[filter] = len(self.object[filter])
            if self.nfiles[filter] > 0:
                self.obsfilters.append(filter)
            print(" %s files in filter %s" % (self.nfiles[filter], filter),
                  file=sys.stderr)

        return

    def getfilters(self, pattern="dflat*.fits"):

        # Get the list of dome flats to extract the filter names
        list = glob.glob(os.path.join(self.obsnight, pattern))
        list.sort()

        filters = []
        for file in list:

            try:
                header = pyfits.getheader(file)
                filter1 = header['FILTER1']
                filter2 = header['FILTER2']
                if filter1.split()[1] != 'Open':
                    filter = filter1.split()[1]
                else:
                    filter = filter2.split()[1]
                obstype = header['OBSTYPE']

                print("*** Found %s -- %s" % (filter, file))

            except:
                print("Skipping %s" % file)
                continue

            #try:
            #    filter = filter.split()[1]
            #except:
            #    filter = filter.split()[0]
            if obstype == 'DFLAT' and filter not in filters:

                filters.append(filter)

        filters.sort()
        print(" \tFound %s filters: %s (using dome flats)" %
              (len(filters), string.join(filters, ", ")),
              file=sys.stderr)
        self.filters = filters
        return

    # Change the header from TNX projection to TAN projection
    def fix_wcs_2012B(self):
        for file in self.scifiles:
            try:
                self.fix_wcs_file_2012B(file)
            except:
                print("Failed wcs fix for %s" % file)
        return

        # Do it....
    def fix_wcs_file_2012B(self, file):
        f = pyfits.open(file, "readonly")
        for i in range(4):
            n = i + 1
            header = f[n].header
            image = "%s[%s]" % (file, n)

            # Avoid if already corrcted
            try:
                WCSCOR = header['WCSCOR']
            except:
                WCSCOR = None
            # Skip
            if WCSCOR:
                print("Skipping WCS correction %s" % file)
                continue

            # Make it a tangential projection
            iraf.hedit(image,
                       'CTYPE1',
                       "RA---TAN",
                       add='no',
                       update='yes',
                       verify='no')
            iraf.hedit(image,
                       'CTYPE2',
                       "DEC--TAN",
                       add='no',
                       update='yes',
                       verify='no')
            iraf.hedit(image,
                       'EQUINOX',
                       2000.0,
                       delete='yes',
                       update='yes',
                       verify='no')
            iraf.hedit(image,
                       'EQUINOX',
                       2000.0,
                       add='yes',
                       update='yes',
                       verify='no')

        f.close()
        return

    # Fix the wcs header to value that work. Magically found in the
    # soar/MSU soi task soiwcs.cl. But only for science files to make it faster
    def fix_wcs(self):
        #for file in self.fitslist:
        for file in self.scifiles:
            try:
                self.fix_wcs_file(file)
            except:
                print("Failed wcs fix for %s" % file)
        return

    def fix_wcs_file(self, file):

        # Coefficients copied from soar.msu soiwcs.cl task
        crpix = [1049.003, 537.003, -24.997, -536.997]

        f = pyfits.open(file, "readonly")
        for i in range(4):
            n = i + 1
            header = f[n].header
            image = "%s[%s]" % (file, n)

            # Avoid if already corrcted
            try:
                WCSCOR = header['WCSCOR']
            except:
                WCSCOR = None
            # Skip
            if WCSCOR:
                print("Skipping WCS correction %s" % file)
                continue

            # Read in some useful info
            CRVAL1 = header['CRVAL1']
            CD2_2 = header['CD2_2']
            CD2_1 = header['CD2_1']
            CD1_2 = header['CD1_2']
            CDELT2 = header['CDELT2']
            DECPANGL = header['DECPANGL']
            header = None

            # Fix negative CD2_2 values
            if CD2_2 == -CDELT2:
                CD2_2 = CD2_2 * -1.0
                iraf.hedit(image,
                           'CD2_2',
                           CD2_2,
                           add='no',
                           update='yes',
                           verify='no')

            # Rotation angle
            CD2_2 = CD2_2 * math.cos(DECPANGL * math.pi / 180.0)
            iraf.hedit(image,
                       'CD2_2',
                       CD2_2,
                       add='no',
                       update='yes',
                       verify='no')

            # Make it a tangential projection
            iraf.hedit(image,
                       'CTYPE1',
                       "RA---TAN",
                       add='no',
                       update='yes',
                       verify='no')
            iraf.hedit(image,
                       'CTYPE2',
                       "DEC--TAN",
                       add='no',
                       update='yes',
                       verify='no')
            iraf.hedit(image,
                       'EQUINOX',
                       2000.0,
                       delete='yes',
                       update='yes',
                       verify='no')
            iraf.hedit(image,
                       'EQUINOX',
                       2000.0,
                       add='yes',
                       update='yes',
                       verify='no')
            # into degrees
            CRVAL1 = CRVAL1 * 15.0  # change it to degress
            iraf.hedit(image,
                       'CRVAL1',
                       CRVAL1,
                       add='no',
                       update='yes',
                       verify='no')
            iraf.hedit(image, 'CD1_2', 0, add='no', update='yes', verify='no')
            iraf.hedit(image, 'CD2_1', 0, add='no', update='yes', verify='no')
            # the magical values that work
            iraf.hedit(image,
                       "CRPIX1",
                       crpix[i],
                       verify='no',
                       add='no',
                       addonly='no',
                       show='yes',
                       delete='no',
                       update='yes')
            iraf.hedit(image,
                       "CRPIX2",
                       1024,
                       verify='no',
                       add='no',
                       addonly='no',
                       show='yes',
                       delete='no',
                       update='yes')
            # to mark that we did this
            iraf.hedit(image,
                       "WCSCOR",
                       'yes',
                       verify='no',
                       add='yes',
                       addonly='no',
                       show='no',
                       delete='no',
                       update='yes')

        f.close()
        return

    def ccdhedit(self):

        # DATASEC
        datasec = "[29:540,1:2048]"
        print("# updating the DATASEC %s " % datasec, file=sys.stderr)
        mscred.ccdhedit.extname = "*"
        mscred.ccdhedit.type = "string"
        mscred.ccdhedit("*.fits", "DATASEC", datasec)

        # TRIMSEC
        trimsec = "[29:540,1:2048]"
        print("# updating the TRIMSEC %s " % trimsec, file=sys.stderr)
        mscred.ccdhedit.extname = "*"
        mscred.ccdhedit.type = "string"
        mscred.ccdhedit("*.fits", "TRIMSEC", trimsec)

        # BIASSEC
        biassec = "[541:568,1:2048]"
        print("# updating the BIASSEC %s im1,im3" % biassec, file=sys.stderr)
        # amp#1
        mscred.ccdhedit.extname = "im1"
        mscred.ccdhedit.type = "string"
        mscred.ccdhedit("*.fits", "BIASSEC", biassec)
        # amp#3
        mscred.ccdhedit.extname = "im3"
        mscred.ccdhedit.type = "string"
        mscred.ccdhedit("*.fits", "BIASSEC", biassec)

        biassec = "[1:28,1:2048]"
        print("# updating the BIASSEC %s im2,im4" % biassec, file=sys.stderr)
        # amp#2
        mscred.ccdhedit.extname = "im2"
        mscred.ccdhedit.type = "string"
        mscred.ccdhedit("*.fits", "BIASSEC", biassec)
        # amp#4
        mscred.ccdhedit.extname = "im4"
        mscred.ccdhedit.type = "string"
        mscred.ccdhedit("*.fits", "BIASSEC", biassec)

        return

# ----------------------------------------------------
#         Functions to setIraf's parameters
# ----------------------------------------------------


def setIraf(obs, verb=None):
    # Set vars
    extras.cl_bye(verb)
    iraf.set(imclobber="yes")
    iraf.set(clobber="yes")
    iraf.set(Verbose='yes')

    # Initialize some tasks
    print(" setting up IRAF's initial params", file=sys.stderr)
    #set_instrument(verb)
    set_mscred(verb)
    set_ccdproc(obs, verb)
    return


def set_instrument(verb=None):

    print(" initializing mscred.setinstrument", file=sys.stderr)
    mscred.setinstrument.site = 'ctio'
    mscred.setinstrument.telescope = '4meter'
    mscred.setinstrument.instrument = 'Mosaic2'
    mscred.setinstrument.directory = 'mscdb$noao/'
    mscred.setinstrument.review = 'no'
    mscred.setinstrument.query_site = 'ctio'
    mscred.setinstrument.query_tel = '4meter'
    mscred.setinstrument.query_inst = 'Mosaic2'
    mscred.setinstrument.mode = 'al'
    return


def set_mscred(verb=None):

    # Check for environ vars
    if not os.getenv('BCSPIPE'):
        os.environ['BCSPIPE'] = '/Users/felipe/BCSPIPE.coligue'
    BCSPIPE = os.getenv('BCSPIPE')
    print(BCSPIPE)

    print(" initializing mscred", file=sys.stderr)
    mscred.pixeltype = "real real"
    mscred.verbose = 'yes'
    mscred.logfile = "logfile_%s" % time.strftime("%Y_%m_%d_%H:%M:%S_%Z")
    #mscred.plotfile   = ""
    #mscred.backup    = "once"
    mscred.backup = "none"
    mscred.bkuproot = "Raw/"
    #mscred.instrument = "mscdb$noao/ctio/4meter/Mosaic2.dat"
    #mscred.instrument = "/Users/felipe/iraf/tasks/soar/soi3.dat"
    mscred.instrument = os.path.join(BCSPIPE, "LIB/SOI/soi3.dat")
    #mscred.ampfile    = "amps"
    #mscred.ssfile     = "subsets"
    #mscred.im_bufsize = 4.0
    #mscred.graphics   = "stdgraph"
    #mscred.cursor     = ""
    #mscred.version    = "V4.8: May 11= 2004"
    #mscred.mode       = "ql"
    return


def set_ccdproc(obs, verb=None):
    '''Initialize ccdproc'''

    print(" initializing mscred.ccdproc", file=sys.stderr)

    mscred.ccdproc.images = ''
    mscred.ccdproc.output = ''
    mscred.ccdproc.bpmasks = ''
    mscred.ccdproc.ccdtype = ''
    mscred.ccdproc.noproc = 'no'

    #mscred.ccdproc.xtalkcor    = 'yes'
    #mscred.ccdproc.fixpix      = 'no'
    mscred.ccdproc.xtalkcor = 'no'
    mscred.ccdproc.fixpix = 'no'
    mscred.ccdproc.overscan = 'yes'
    mscred.ccdproc.trim = 'yes'
    mscred.ccdproc.zerocor = 'no'
    mscred.ccdproc.darkcor = 'no'
    mscred.ccdproc.flatcor = 'no'
    mscred.ccdproc.sflatcor = 'no'
    mscred.ccdproc.split = 'no'
    mscred.ccdproc.merge = 'no'
    #mscred.ccdproc.merge       = 'yes'

    #mscred.ccdproc.xtalkfile   = '!xtalkfil'
    #mscred.ccdproc.xtalkfile   = "bcspipe$LIB/MosaicII/xtalk/test2.dat"
    #mscred.ccdproc.xtalkfile   = os.path.join(obs.BCSPIPE,"LIB/MosaicII/xtalk/test2.dat")
    #mscred.ccdproc.fixfile     = 'BPM'
    mscred.ccdproc.saturation = 'INDEF'
    mscred.ccdproc.sgrow = 0
    mscred.ccdproc.bleed = 'INDEF'
    mscred.ccdproc.btrail = 20
    mscred.ccdproc.bgrow = 0
    mscred.ccdproc.biassec = '!biassec'
    mscred.ccdproc.trimsec = '!trimsec'
    mscred.ccdproc.zero = obs.Zero
    mscred.ccdproc.dark = 'Dark'
    mscred.ccdproc.flat = obs.Dflat + "*"
    mscred.ccdproc.sflat = obs.Sflat + "*"
    mscred.ccdproc.minreplace = 1.0

    mscred.ccdproc.interactive = 'no'
    mscred.ccdproc.function = 'legendre'
    mscred.ccdproc.order = 1
    mscred.ccdproc.sample = '*'
    mscred.ccdproc.naverage = 1
    mscred.ccdproc.niterate = 1
    mscred.ccdproc.low_reject = 3.0
    mscred.ccdproc.high_reject = 3.0
    mscred.ccdproc.grow = 0.0
    mscred.ccdproc.fd = ''
    mscred.ccdproc.fd2 = ''
    mscred.ccdproc.mode = 'h'
    extras.cl_bye(verb)

    return


def tweak_filter(file, Namps=4):
    import pyfits

    header = pyfits.getheader(file)
    filter1 = header['FILTER1']
    filter2 = header['FILTER2']

    if filter1.split()[1] != 'Open':
        filter = filter1.split()[1]
    else:
        filter = filter2.split()[1]

        #try:
        #    filter = filter1.split()[1]
        #except:
        #    filter = filter1.split()[0]

    f = pyfits.open(file, "readonly")
    for i in range(Namps):
        n = i + 1
        header = f[n].header
        image = "%s[%s]" % (file, n)
        filter1 = header['FILTER1']
        if filter != filter1:
            iraf.hedit(image,
                       'FILTER1',
                       filter,
                       add='no',
                       update='yes',
                       verify='no')
            iraf.hedit(image,
                       'FILTER',
                       filter,
                       add='yes',
                       update='yes',
                       verify='no')
        else:
            print("Skipping filter %s %s" % (filter, image))
    f.close()
    return
