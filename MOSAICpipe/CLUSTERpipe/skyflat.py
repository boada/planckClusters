from __future__ import print_function
from builtins import str
from builtins import range
from builtins import object
import os, sys
import string, glob
import extras
from pyraf import iraf
from iraf import mscred
from iraf import nproto
from iraf import proto
from iraf import images


class skyflat(object):
    ''' Perform the several steps required to make the super sky-flats'''

    def __init__(self, obs):

        self.verb = obs.verb
        self.scilist = obs.scilist
        self.bplist = obs.bplist

        self.scifiles = obs.scifiles
        self.filters = obs.filters
        self.Namps = obs.Namps

        self.Zero = obs.Zero
        self.Dflat = obs.Dflat
        self.Sflat = obs.Sflat
        self.Fringe = obs.Fringe
        self.object = obs.object
        self.objlist = obs.objlist

        self.forsky = obs.forsky

        self.obsfilters = obs.obsfilters  # the filters actually observed
        self.nfiles = obs.nfiles

        self.msklist = obs.msklist
        self.skylist = obs.skylist
        self.frgfilters = ['i', 'z']
        self.BCSPIPE = obs.BCSPIPE

        self.skyfilters = obs.skyfilters  # The filters for which we generate sky flats and fringe too

        return

    def objmasks(self):

        extras.cl_bye(self.verb)
        # Fix objmasks1
        nproto.objmasks1.exps = ""
        nproto.objmasks1.gains = ""
        nproto.objmasks1.catalogs = ""
        nproto.objmasks1.catdefs = ""
        nproto.objmasks1.dodetect = "yes"
        nproto.objmasks1.dosplit = "no"
        nproto.objmasks1.dogrow = "yes"
        nproto.objmasks1.doevaluate = "no"
        nproto.objmasks1.skytype = "block"
        nproto.objmasks1.fitstep = "10"
        nproto.objmasks1.fitblk1d = "10"
        nproto.objmasks1.fithclip = "2.0"
        nproto.objmasks1.fitlclip = "3.0"
        nproto.objmasks1.fitxorder = "1"
        nproto.objmasks1.fityorder = "1"
        nproto.objmasks1.fitxterms = "half"
        nproto.objmasks1.blknsubblks = "2"
        nproto.objmasks1.updatesky = "yes"
        nproto.objmasks1.sigavg = "4.0"
        nproto.objmasks1.sigmax = "4.0"
        nproto.objmasks1.bpval = "INDEF"
        nproto.objmasks1.splitmax = "INDEF"
        nproto.objmasks1.splitstep = "0.4"
        nproto.objmasks1.splitthresh = "5.0"
        nproto.objmasks1.sminpix = "8"
        nproto.objmasks1.ssigavg = "10.0"
        nproto.objmasks1.ssigmax = "5.0"
        nproto.objmasks1.magzero = "INDEF"
        nproto.objmasks1.mode = "ql"

        # One filter at a time
        for filter in self.filters:

            print(" running object masks for filter %s on images" % filter,
                  file=sys.stderr)
            for i in self.forsky[filter]:
                print("\t%s" % i, file=sys.stderr)

            # The old way --  all at one
            #nproto.objmasks.images    = self.scilist
            #nproto.objmasks.objmasks  = self.msklist
            #nproto.objmasks.skys      = self.skylist
            nproto.objmasks.images = extras.imfile(self.forsky[filter])
            nproto.objmasks.objmasks = extras.imfile_append(
                self.forsky[filter], "_msk")
            nproto.objmasks.omtype = "numbers"
            nproto.objmasks.skys = extras.imfile_append(self.forsky[filter],
                                                        "_sky")
            nproto.objmasks.sigmas = ""
            nproto.objmasks.masks = "!BPM"
            nproto.objmasks.extnames = ""
            nproto.objmasks.logfiles = "STDOUT"
            nproto.objmasks.blkstep = "1"
            nproto.objmasks.blksize = "-10"
            nproto.objmasks.convolve = "block 3 3"
            nproto.objmasks.hsigma = "3.0"
            nproto.objmasks.lsigma = "10.0"
            nproto.objmasks.hdetect = "yes"
            nproto.objmasks.ldetect = "no"
            nproto.objmasks.neighbors = "8"
            nproto.objmasks.minpix = "6"
            nproto.objmasks.ngrow = "2"
            nproto.objmasks.agrow = "2.0"
            nproto.objmasks.mode = "h"  # to avoid prompt
            # Run
            nproto.objmasks()
            extras.cl_bye(self.verb)

        return

    def super_skyflat(self, ver=""):

        self.output = self.Sflat + ver
        print(" Creating super sky-flat with %s images" % len(self.scifiles),
              file=sys.stderr)

        extras.cl_bye()
        #mscred.sflatcombine.input       = self.scilist
        #mscred.sflatcombine.output      = self.output
        mscred.sflatcombine.combine = "average"
        #mscred.sflatcombine.combine     = "median"
        mscred.sflatcombine.reject = "ccdclip"
        mscred.sflatcombine.ccdtype = "object"
        mscred.sflatcombine.subsets = "yes"
        mscred.sflatcombine.masktype = "!objmask"
        mscred.sflatcombine.maskvalue = "0.0"
        mscred.sflatcombine.scale = "mode"
        mscred.sflatcombine.statsec = ""
        mscred.sflatcombine.nkeep = "1"
        mscred.sflatcombine.nlow = "1"
        mscred.sflatcombine.nhigh = "1"
        mscred.sflatcombine.mclip = "yes"
        mscred.sflatcombine.lsigma = "6.0"
        mscred.sflatcombine.hsigma = "3.0"
        mscred.sflatcombine.rdnoise = "rdnoise"
        mscred.sflatcombine.gain = "gain"
        mscred.sflatcombine.snoise = "0."
        mscred.sflatcombine.pclip = "-0.5"
        mscred.sflatcombine.blank = "1.0"
        mscred.sflatcombine.grow = "3.0"
        mscred.sflatcombine.fd = ""
        mscred.sflatcombine.mode = "h"
        #mscred.sflatcombine()

        # Go filter by filter to avoid iraf's crash
        for filter in self.obsfilters:

            if filter not in self.skyfilters:
                print("# Skipping %s filter, no need for Sflat" % filter,
                      file=sys.stderr)
                continue

            if len(self.forsky[filter]) < 4:
                print("# Skipping %s filter, no sky images available" % filter,
                      file=sys.stderr)
                continue

            for l in self.forsky[filter]:
                print("\t%s" % l, file=sys.stderr)
            print("\t\t\t\t --> %s %s" % (self.output, filter),
                  file=sys.stderr)

            #mscred.sflatcombine.input       = self.objlist[filter]
            mscred.sflatcombine.input = extras.imfile(self.forsky[filter])
            mscred.sflatcombine.output = self.output
            mscred.sflatcombine()
            extras.cl_bye(self.verb)

        extras.cl_bye(self.verb)
        return

    def make_fringe(self, ver=""):

        extras.cl_bye()
        #mscred.mscmedian.input        = skyflat
        #mscred.mscmedian.output       = tmpmedian
        mscred.mscmedian.xwindow = 129
        mscred.mscmedian.ywindow = 129
        mscred.mscmedian.outtype = "median"
        mscred.mscmedian.zloreject = -20000.0
        mscred.mscmedian.zhireject = 30000.0
        mscred.mscmedian.verbose = "yes"
        mscred.mscmedian.fmedian = "yes"
        mscred.mscmedian.hmin = -20000.0
        mscred.mscmedian.hmax = 30000.0
        mscred.mscmedian.zmin = -20000.0
        mscred.mscmedian.zmax = 30000.0
        mscred.mscmedian.fd = ""
        mscred.mscmedian.mode = "h"

        # Run it for each filter with fringe
        for filter in self.frgfilters:

            # Avoid if not observed
            if filter not in self.obsfilters:
                print(" No %s-band observations... skipping" % filter,
                      file=sys.stderr)
                continue

            # Avoid if Sflat was not created
            if len(self.forsky[filter]) < 4:
                print("# Skipping %s filter, no sky images available" % filter,
                      file=sys.stderr)
                continue

            print(" Making Fringe for filter %s " % filter, file=sys.stderr)

            skyflat = self.Sflat + ver + filter
            tmpmedian = "tmpMedian_" + filter
            print("\trunning mscmedian...\t %s --> %s " % (skyflat, tmpmedian),
                  file=sys.stderr)
            if os.path.isfile(tmpmedian + ".fits"):
                os.remove(tmpmedian + ".fits")
                print(" removed file: %s" % (tmpmedian + ".fits"),
                      file=sys.stderr)

                # 1 - Make the median
            mscred.mscmedian(skyflat, tmpmedian, 129, 129)
            extras.cl_bye()

            # 2 - run mscarith on the sky image
            Fringe = self.Fringe + filter
            print("\trunning mscarith...\t %s - %s = %s " %
                  (skyflat, tmpmedian, Fringe),
                  file=sys.stderr)

            mscred.mscarith(skyflat, "-", tmpmedian, Fringe)
            extras.cl_bye()

        print(" Fringe images ready", file=sys.stderr)
        return

    def correct_fringe(self):

        print(" Correcting images with fringe map", file=sys.stderr)

        # Run it for each filter with fringe
        for filter in self.frgfilters:

            # Avoid if not observed
            if filter not in self.obsfilters:
                print(" No %s-band observations... skipping" % filter,
                      file=sys.stderr)
                continue

            if filter not in self.obsfilters:
                print(" No %s-band observations... skipping" % filter,
                      file=sys.stderr)
                continue

            print(" Correcting Fringe for filter %s " % filter,
                  file=sys.stderr)

            # Set up names
            fringe = self.Fringe + filter

            # Avoid if not observed
            if filter not in self.obsfilters:
                print(" No %s-band observations... skipping" % filter,
                      file=sys.stderr)
                continue

            # Avoid if Sflat was not created
            if len(self.forsky[filter]) < 4:
                print("# Skipping %s filter, no Fringe images available" %
                      filter,
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
            mscred.rmfringe.masks = extras.imfile_append(self.forsky[filter],
                                                         "_msk")
            mscred.rmfringe.fringemasks = ""
            mscred.rmfringe.background = extras.imfile_append(
                self.forsky[filter], "_sky")
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

    # interpolate flats with bpm mask
    def fixpix(self):

        N = self.Namps
        for filter in self.obsfilters:

            # Avoid if Sflat was not created
            if len(self.forsky[filter]) < 4:
                print("# Skipping %s filter, no sky images available" % filter,
                      file=sys.stderr)
                continue

            flat = self.Sflat + filter
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

        return

    # Fix badpixel in the Sflat images, that were transported from
    # science images. These are usually present at the bottom edge of
    # some i,z band images, which subsequentaly appear in the Fringe
    # images and really mess up the combined image. The solution is to
    # interpolate all the pixels with values of <=1 with fixpix and
    # imexpr.
    def badpix(self, fix=None):

        N = self.Namps

        express = "a > 1 ? 0 : 1"

        for filter in self.obsfilters:

            # Avoid if Sflat was not created
            if len(self.forsky[filter]) < 4:
                print("# Skipping %s filter, no sky images available" % filter,
                      file=sys.stderr)
                continue

            flat = self.Sflat + filter
            bpmdir = "%s_BPM" % flat

            print("Making bpm for %s" % flat)

            # Check if the destination folder exists
            if os.path.exists(bpmdir):
                print("Will put mask files to: %s" % bpmdir, file=sys.stderr)
            else:
                print("Will create new folder: %s" % bpmdir, file=sys.stderr)
                os.mkdir(bpmdir)

            for i in range(N):

                image = "%s[%s]" % (flat, i + 1)
                mask = "%s/mask_%s.pl" % (bpmdir, i + 1)

                # make bpm
                images.imutil.imexpr(express, mask, image)
                extras.cl_bye()

                # Correct the Sky flat
                if fix:
                    # Interpolate bad pixels
                    proto.fixpix(image,
                                 mask,
                                 linterp="INDEF",
                                 cinterp="INDEF",
                                 verbose="yes",
                                 pixels="no")
                    extras.cl_bye()

        return
