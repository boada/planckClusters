#!/usr/bin/env python

from __future__ import print_function
from builtins import object
from pyraf import iraf
from iraf import mscred
import sys
import os
import extras


class astrometry(object):
    ''' Fix the WCS system with new plate solution '''

    def __init__(self, obs):

        self.scifiles = obs.scifiles
        self.scilist = obs.scilist
        self.verb = obs.verb
        self.object = obs.object
        self.filters = obs.filters
        self.BCSPIPE = obs.BCSPIPE

        #self.update_wcs()
        #self.msccmatch()

        return

    def update_wcs(self, wcs="MosaicII.db"):
        ''' Update - fix the wcs solution '''

        wcs = os.path.join(self.BCSPIPE, "LIB/MosaicII/wcs", wcs)
        print("\t Updating wcs solution for", file=sys.stderr)
        for file in self.scifiles:
            print(" \t\t%s" % file, file=sys.stderr)
        print("\t\t\twith %s" % wcs, file=sys.stderr)

        mscred.mscsetwcs(self.scilist, wcs)
        extras.cl_bye(self.verb)

        return

    def msccmatch(self):

        #mscred.msccmatch.coords         = "!mscgetcat $I $C magmin=14.0 magmax=19.0"
        mscred.msccmatch.coords = "!mscgetcat $I $C magmin=12.0 magmax=21.0 catalog = NOAO:USNO-A2"
        mscred.msccmatch.coords = "!mscgetcat $I $C magmin=12.0 magmax=21.0 catalog = CADC:USNO-A2"
        mscred.msccmatch.usebpm = "yes"
        mscred.msccmatch.verbose = "yes"
        mscred.msccmatch.nsearch = "60"
        mscred.msccmatch.search = "60.0"
        mscred.msccmatch.rsearch = "2.0"

        mscred.msccmatch.cbox = "11"
        mscred.msccmatch.maxshift = "3.0"
        mscred.msccmatch.csig = "0.1"
        mscred.msccmatch.cfrac = "0.5"
        mscred.msccmatch.listcoords = "no"

        mscred.msccmatch.nfit = "60"
        mscred.msccmatch.rms = "1.0"
        mscred.msccmatch.fitgeometry = "general"
        mscred.msccmatch.reject = "2.5"
        mscred.msccmatch.update = "yes"
        mscred.msccmatch.interactive = "no"
        mscred.msccmatch.fit = "yes"
        mscred.msccmatch.graphics = "stdgraph"
        mscred.msccmatch.cursor = ""
        mscred.msccmatch.accept = "yes"
        mscred.msccmatch.mode = "h"

        errorfile = "error.log"

        filenames = []
        for line in open(self.scilist[1:]).readlines(
        ):  # Trick to avoid the '@' in filenames

            file = line.split()[0]
            filenames.append(file)

            print("Doing file: %s" % file, file=sys.stderr)

            error = None
            #mscred.msccmatch.coords = "!mscgetcat $I $C magmin=9.0 magmax=17.0"
            #mscred.msccmatch.nfit   = "60"

            print("\t Running msccmatch on: %s" % file, file=sys.stderr)
            mscred.msccmatch(file, Stderr=errorfile, Stdout="STDOUT")
            extras.cl_bye(self.verb)

            for line in open(errorfile).readlines():
                #print "line:",line
                vals = line.split()

                if "ERROR:" in vals:
                    print("msccmatch problem, retrying... ", file=sys.stderr)
                    error = "yes"
                    continue

            if error:
                # Use brigher lower limit
                print("\t Re-Running msccmatch on: %s" % file, file=sys.stderr)
                mscred.msccmatch.coords = "!mscgetcat $I $C magmin=8.0 magmax=16.0"
                mscred.msccmatch.nfit = "20"
                mscred.msccmatch(file, nfit=20)
                extras.cl_bye(self.verb)

                # Again
                mscred.msccmatch(file, nfit=20)
                extras.cl_bye(self.verb)

        return
