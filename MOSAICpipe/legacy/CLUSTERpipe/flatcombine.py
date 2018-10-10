from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
import os, sys
import string, glob
import extras
from pyraf import iraf
from iraf import mscred
from iraf import images
from iraf import proto


class flatcombine(object):
    ''' Creates the dome flats image'''

    def __init__(self, obs, process='yes'):

        self.process = process
        self.flatfiles = obs.flatfiles
        self.flatlist = obs.flatlist
        self.verb = obs.verb
        self.Dflat = obs.Dflat
        self.Zero = obs.Zero
        self.filters = obs.filters
        self.Namps = obs.Namps

        self.obsfilters = obs.obsfilters  # the filters actually observed
        self.nfiles = obs.nfiles
        self.BCSPIPE = obs.BCSPIPE

        return

    def combine(self):

        extras.cl_bye(self.verb)

        print(" Making superbias image with %s files:" % len(self.flatfiles),
              file=sys.stderr)
        for i in self.flatfiles:
            print("\t%s" % i, file=sys.stderr)
        for filter in self.filters:
            print("\t\t\t--> %s_%s.fits" % (self.Dflat, filter),
                  file=sys.stderr)

        # Modify the params in ccdproc that changed
        mscred.ccdproc.ccdtype = 'flat'
        mscred.ccdproc.noproc = 'no'
        mscred.ccdproc.zerocor = 'yes'

        # And the params for flatcombine
        mscred.flatcombine.combine = "average"
        mscred.flatcombine.reject = "avsigclip"
        mscred.flatcombine.ccdtype = "DFLAT"
        mscred.flatcombine.process = self.process
        mscred.flatcombine.subsets = 'yes'
        mscred.flatcombine.delete = 'no'
        mscred.flatcombine.scale = "mode"
        mscred.flatcombine.statsec = ""
        mscred.flatcombine.nlow = 1
        mscred.flatcombine.nhigh = 1
        mscred.flatcombine.nkeep = 1
        mscred.flatcombine.mclip = 'yes'
        mscred.flatcombine.lsigma = 3.0
        mscred.flatcombine.hsigma = 3.0
        mscred.flatcombine.rdnoise = "0."
        mscred.flatcombine.gain = "1."
        mscred.flatcombine.snoise = "0."
        mscred.flatcombine.pclip = -0.5
        mscred.flatcombine.blank = 1.
        mscred.flatcombine.mode = "hl"
        # Make the call

        print("# Will use flatlist %s" % self.flatlist)

        mscred.flatcombine(self.flatlist, output=self.Dflat)

        print("# flatcombine is ready")

        extras.cl_bye(self.verb)
        return

    # interpolate flats with bpm mask
    def fixpix(self):

        N = self.Namps
        for filter in self.obsfilters:

            flat = self.Dflat + filter

            print("Fixpix %s" % flat)

            for i in range(N):
                image = "%s[%s]" % (flat, i + 1)
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
