from __future__ import print_function
from builtins import object
import os, sys
import shutil
import string, glob
import extras
from pyraf import iraf
from iraf import mscred


class zerocombine(object):
    ''' Creates the super bias image'''

    def __init__(self, obs, process='no'):

        self.zerofiles = obs.zerofiles
        self.zerolist = obs.zerolist
        self.verb = obs.verb
        self.Zero = obs.Zero
        self.process = process

        self.obsfilters = obs.obsfilters  # the filters actually observed
        self.nfiles = obs.nfiles

        return

    def combine(self):

        print("\n Making superbias image with %s files:" % len(self.zerofiles),
              file=sys.stderr)
        if self.verb:
            for i in self.zerofiles:
                print("\t%s" % i, file=sys.stderr)
            print("\t\t\t--> %s" % self.Zero + ".fits", file=sys.stderr)

        # Avoid clobbering -- way around imclobber = yes not working
        if os.path.isfile(self.Zero + ".fits"):
            print(" Moving of %s.fits to %s.fits.old" % (self.Zero, self.Zero),
                  file=sys.stderr)
            shutil.move(self.Zero + ".fits", self.Zero + ".fits.old")

        # And the params for zerocombine
        iraf.mscred.zerocombine.combine = 'average'
        iraf.mscred.zerocombine.reject = 'minmax'
        iraf.mscred.zerocombine.ccdtype = 'zero'
        iraf.mscred.zerocombine.process = self.process
        iraf.mscred.zerocombine.delete = 'no'
        iraf.mscred.zerocombine.scale = 'none'
        iraf.mscred.zerocombine.statsec = ''
        iraf.mscred.zerocombine.nlow = 0
        iraf.mscred.zerocombine.nhigh = 1
        iraf.mscred.zerocombine.nkeep = 1
        iraf.mscred.zerocombine.mclip = 'yes'
        iraf.mscred.zerocombine.lsigma = 3.0
        iraf.mscred.zerocombine.hsigma = 3.0
        iraf.mscred.zerocombine.rdnoise = '0.'
        iraf.mscred.zerocombine.gain = '1.'
        iraf.mscred.zerocombine.snoise = '0.'
        iraf.mscred.zerocombine.pclip = -0.5
        iraf.mscred.zerocombine.blank = 0.0
        iraf.mscred.zerocombine.mode = 'h'
        # Make the call
        iraf.mscred.zerocombine(self.zerolist, output=self.Zero)
        extras.cl_bye(self.verb)
        return
