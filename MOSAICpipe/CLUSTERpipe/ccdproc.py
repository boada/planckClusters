from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
import os, sys
import string, glob
import extras
import pyfits
import numpy
from pyraf import iraf
from iraf import mscred
from iraf import proto


class ccdproc(object):
    ''' Passes of ccdproc '''

    def __init__(self, obs):

        self.verb = obs.verb
        self.scilist = obs.scilist
        #self.bplist   = obs.bplist

        self.scifiles = obs.scifiles
        self.filters = obs.filters
        self.Namps = obs.Namps

        self.Zero = obs.Zero
        self.Dflat = obs.Dflat
        self.Sflat = obs.Sflat
        self.object = obs.object

        self.forsky = obs.forsky

        self.obsfilters = obs.obsfilters  # the filters actually observed
        self.nfiles = obs.nfiles
        self.BCSPIPE = obs.BCSPIPE
        self.frgfilters = ['i', 'z']
        self.skyfilters = obs.skyfilters  # filters for which we generate Sflats
        self.Fringe = obs.Fringe

        return

    def pass_one(self):
        ''' First pass of ccdproc '''
        extras.cl_bye(self.verb)

        print(" Making 1st pass on images: %s " % len(self.scifiles),
              file=sys.stderr)
        print(" \t\t\toverscan, zero, trim", file=sys.stderr)
        print(" Reading file from %s" % self.scilist, file=sys.stderr)

        # Set more params
        mscred.ccdproc.ccdtype = 'object'
        mscred.ccdproc.noproc = 'no'

        mscred.ccdproc.xtalkcor = 'no'
        mscred.ccdproc.fixpix = 'no'
        mscred.ccdproc.overscan = 'yes'
        mscred.ccdproc.trim = 'yes'
        mscred.ccdproc.zerocor = 'yes'
        mscred.ccdproc.darkcor = 'no'
        mscred.ccdproc.flatcor = 'no'
        mscred.ccdproc.sflatcor = 'no'
        mscred.ccdproc.split = 'no'
        mscred.ccdproc.merge = 'no'
        #mscred.ccdproc.merge    = 'yes'

        #mscred.ccdproc.xtalkfile   = '!xtalkfil'
        #mscred.ccdproc.fixfile     = 'BPM'
        mscred.ccdproc.saturation = '!saturate'
        mscred.ccdproc.sgrow = 1
        mscred.ccdproc.bleed = 20000
        mscred.ccdproc.btrail = 20
        mscred.ccdproc.bgrow = 0
        mscred.ccdproc.biassec = '!biassec'
        mscred.ccdproc.trimsec = '!trimsec'
        mscred.ccdproc.zero = self.Zero
        mscred.ccdproc.dark = 'Dark'
        mscred.ccdproc.flat = self.Dflat + "*"
        mscred.ccdproc.sflat = self.Sflat + "*"
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

        #mscred.ccdproc.images      = self.scilist
        mscred.ccdproc.output = ''
        #mscred.ccdproc.bpmasks     = self.bplist

        print("# Will use scilist %s" % self.scilist)

        mscred.ccdproc(self.scilist)
        extras.cl_bye(self.verb)

        return

    def pass_two(self):
        ''' Second pass of ccdproc '''
        extras.cl_bye(self.verb)

        print(" Making 2nd pass on images: %s " % len(self.scifiles),
              file=sys.stderr)
        print(" \t\t\tdome flats", file=sys.stderr)

        # Set more params
        mscred.ccdproc.ccdtype = 'object'
        mscred.ccdproc.noproc = 'no'

        mscred.ccdproc.xtalkcor = 'no'
        mscred.ccdproc.fixpix = 'no'
        mscred.ccdproc.overscan = 'yes'
        mscred.ccdproc.trim = 'yes'
        mscred.ccdproc.zerocor = 'yes'
        mscred.ccdproc.darkcor = 'no'
        mscred.ccdproc.flatcor = 'yes'
        mscred.ccdproc.sflatcor = 'no'
        mscred.ccdproc.split = 'no'
        mscred.ccdproc.merge = 'no'

        #mscred.ccdproc.xtalkfile   = '!xtalkfil'
        #mscred.ccdproc.fixfile     = 'BPM'
        mscred.ccdproc.saturation = '!saturate'
        mscred.ccdproc.sgrow = 1
        mscred.ccdproc.bleed = 20000
        mscred.ccdproc.btrail = 20
        mscred.ccdproc.bgrow = 0
        mscred.ccdproc.biassec = '!biassec'
        mscred.ccdproc.trimsec = '!trimsec'
        mscred.ccdproc.zero = self.Zero
        mscred.ccdproc.dark = 'Dark'
        mscred.ccdproc.flat = self.Dflat + "*"
        mscred.ccdproc.sflat = self.Sflat + "*"
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

        mscred.ccdproc.images = self.scilist
        mscred.ccdproc.output = ''
        #mscred.ccdproc.bpmasks     = self.bplist

        # Run it
        mscred.ccdproc()

        extras.cl_bye(self.verb)
        return

    def pass_three(self, skyfilters=('i', 'z')):

        self.skyfilters = skyfilters
        ''' Thrid pass of ccdproc '''
        extras.cl_bye(self.verb)

        print(" Making 3rd pass on images: %s " % len(self.scifiles),
              file=sys.stderr)
        print(" \t\t\tsky flats", file=sys.stderr)

        # Set more params
        mscred.ccdproc.ccdtype = 'object'
        mscred.ccdproc.noproc = 'no'

        mscred.ccdproc.xtalkcor = 'no'
        mscred.ccdproc.fixpix = 'no'
        mscred.ccdproc.overscan = 'no'
        mscred.ccdproc.trim = 'no'
        mscred.ccdproc.zerocor = 'no'
        mscred.ccdproc.darkcor = 'no'
        mscred.ccdproc.flatcor = 'no'
        mscred.ccdproc.sflatcor = 'yes'
        mscred.ccdproc.split = 'no'
        mscred.ccdproc.merge = 'no'

        #mscred.ccdproc.xtalkfile   = '!xtalkfil'
        #mscred.ccdproc.fixfile     = 'BPM'
        mscred.ccdproc.saturation = '!saturate'
        mscred.ccdproc.sgrow = 1
        mscred.ccdproc.bleed = 20000
        mscred.ccdproc.btrail = 20
        mscred.ccdproc.bgrow = 0
        mscred.ccdproc.biassec = '!biassec'
        mscred.ccdproc.trimsec = '!trimsec'
        mscred.ccdproc.zero = self.Zero
        mscred.ccdproc.dark = 'Dark'
        mscred.ccdproc.flat = self.Dflat + "*"
        mscred.ccdproc.sflat = self.Sflat + "*"
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

        #mscred.ccdproc.images      = self.scilist
        #mscred.ccdproc.output      = ''
        #mscred.ccdproc.bpmasks     = self.bplist

        # Do it filter by filter to avoid non Sflat images
        for filter in self.obsfilters:

            if len(self.forsky[filter]) < 1:
                print("# Skipping %s filter, no sky images available" % filter,
                      file=sys.stderr)
                continue

            if filter not in self.skyfilters:
                print("# Skipping %s filter, no need sky image for " % filter,
                      file=sys.stderr)
                continue

            mscred.ccdproc.images = extras.imfile(self.object[filter])
            mscred.ccdproc.output = ''
            mscred.ccdproc.bpmasks = extras.imfile_append(self.object[filter],
                                                          "_BPM")

            # Run it
            mscred.ccdproc()
            extras.cl_bye(self.verb)

        return

    # interpolate flats with bpm mask
    def fixpix(self):

        N = self.Namps
        for file in self.scifiles:

            print("Fixpix %s" % file)

            for i in range(N):
                image = "%s[%s]" % (file, i + 1)
                mask = os.path.join(self.BCSPIPE,
                                    "LIB/MosaicII/BPM/bpmSum%s_0511.pl" %
                                    str(i + 1))
                proto.fixpix(image,
                             mask,
                             linterp="INDEF",
                             cinterp="INDEF",
                             verbose="yes",
                             pixels="no")
                extras.cl_bye()

        return

    # Correct science images for sky bpm
    def fixpix_sky(self):

        N = self.Namps
        for filter in self.obsfilters:

            flat = self.Sflat + filter
            bpmdir = "%s_BPM" % flat

            print("Correcting for bpm for %s" % flat)

            for file in self.forsky[filter]:

                for i in range(N):
                    mask = "%s/mask_%s.pl" % (bpmdir, i + 1)
                    scifile = "%s[%s]" % (file, i + 1)
                    proto.fixpix(scifile,
                                 mask,
                                 linterp="INDEF",
                                 cinterp="INDEF",
                                 verbose="yes",
                                 pixels="no")
                    extras.cl_bye()

        return

    # Apply fringe correction
    def correct_fringe(self):

        print(" Correcting images with fringe map", file=sys.stderr)

        # Run it for each filter with fringe
        for filter in self.frgfilters:

            print(" Correcting Fringe for filter %s " % filter,
                  file=sys.stderr)
            # Set up names
            fringe = self.Fringe + filter

            # Avoid if not observed
            if filter not in self.obsfilters:
                print(" No %s-band observations... skipping" % filter,
                      file=sys.stderr)
                continue

            for ima in self.forsky[filter]:
                print("\t%s" % ima, file=sys.stderr)
            print("\t\t\t\t with fringe: %s" % fringe, file=sys.stderr)

            #print images, fringe, extras.imlist_append(images,"_msk"), extras.imlist_append(images,"_sky")

            extras.cl_bye(self.verb)
            mscred.rmfringe.input = extras.imfile(self.forsky[filter])
            mscred.rmfringe.output = ''
            mscred.rmfringe.fringe = fringe
            #mscred.rmfringe.masks       = extras.imfile_append(self.forsky[filter],"_msk")
            #mscred.rmfringe.fringemasks = ""
            #mscred.rmfringe.background  = extras.imfile_append(self.forsky[filter],"_sky")
            mscred.rmfringe.ncblk = "5"
            mscred.rmfringe.nlblk = "5"
            mscred.rmfringe.extfit = ""
            mscred.rmfringe.logfile = ""
            mscred.rmfringe.verbose = "yes"
            mscred.rmfringe.mode = "h"
            mscred.rmfringe()
            extras.cl_bye(self.verb)

        extras.cl_bye(self.verb)
        return

    def fringe_hack(self):

        # Finge hack for MEF file

        # For each filter
        for filter in self.frgfilters:
            fringe_pattern = self.Fringe + "%s.fits" % filter

            for file in self.forsky[filter]:
                # Read in the fringe pattern
                frg = pyfits.open(fringe_pattern, "readonly")
                ima = pyfits.open(file, "update")
                for i in range(4):
                    n = i + 1
                    # Get the scaling
                    b = ima[n].data
                    scale = numpy.median(b) * 0.003
                    new = ima[n].data - scale * frg[n].data
                    ima[n].data = new

                ima.verify('fix')
                ima.flush()
                ima.close()
                frg.close()
                print("Fringe Corrected %s file" % file)
        return

    # MOSAIC the science frames together
    def mos_frames(self):

        iraf.soar(_doprint=0)
        iraf.soi(_doprint=0)

        for file in self.scifiles:
            input = file
            output = "MOS_" + file
            print("MOSAICING %s -->  %s " % (input, output), file=sys.stderr)
            iraf.soimosaic(input, output, verbose=None)

        return
