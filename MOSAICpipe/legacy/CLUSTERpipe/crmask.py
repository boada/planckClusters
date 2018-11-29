from __future__ import print_function
from builtins import range
from builtins import object
import glob, sys, os
import extras
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import crutil
from iraf import craverage
from iraf import images
from iraf import proto


class crmask(object):
    ''' CR rejection and combined mask creation '''

    def __init__(self, obs):

        self.scifiles = obs.scifiles
        self.scilist = obs.scilist
        self.verb = obs.verb
        self.object = obs.object
        self.filters = obs.filters
        self.Namps = obs.Namps

        return

    def craverage(self):
        ''' Run IRAF craverage '''

        N = self.Namps

        craverage.input = "List"
        craverage.output = ""
        craverage.crmask = ""
        craverage.average = ""
        craverage.sigma = ""

        craverage.navg = "7"
        craverage.nrej = "3"
        craverage.nbkg = "5"
        craverage.nsig = "50"
        craverage.var0 = "0.0"
        craverage.var1 = "0.0"
        craverage.var2 = "0.0"

        craverage.crval = "1"
        craverage.lcrsig = "10.0"
        craverage.hcrsig = "3.75"
        craverage.crgrow = "0.0"

        craverage.objval = "0"
        craverage.lobjsig = "10.0"
        craverage.hobjsig = "5.5"
        craverage.objgrow = "0.0"
        craverage.mode = "al"

        filenames = []
        for line in open(self.scilist[1:]).readlines(
        ):  # Trick to avoid the '@' in filenames

            file = line.split()[0]
            filenames.append(file)

            print("Doing file: %s" % file, file=sys.stderr)

            for i in range(N):
                print("craverage %s[%s]" % (file, i + 1), file=sys.stderr)
                input = "%s[%s]" % (file, i + 1)
                crmask = "%s_crmask_%s" % (file, i + 1)
                craverage(input, output="", crmask=crmask)
                extras.cl_bye()

        # And now we do the mask combination and fixpix
        express = "max(a,b)"
        for ima in filenames:

            # Merging masks
            print(" Merging masks for %s" % ima, file=sys.stderr)
            for i in range(N):

                n = i + 1

                hfile = "%s[%s]" % (ima, n)
                mask1 = "%s_crmask_%s.fits[1]" % (ima, n)
                mask2 = "%s_BPM/bpm_im%s.pl" % (ima, n)
                maskout = "%s_BPM/bpmcr_im%s.fits" % (ima, n)

                images.imutil.imexpr(express, maskout, mask1, mask2)

                print("\t\tUpdating header information for %s" % hfile,
                      file=sys.stderr)

                images.imutil.hedit(hfile,
                                    "BPM",
                                    maskout,
                                    add='no',
                                    delete='no',
                                    ver='no',
                                    show='yes',
                                    update='yes')
                images.imutil.hedit(maskout,
                                    "BPM",
                                    maskout,
                                    add='no',
                                    delete='no',
                                    ver='no',
                                    show='yes',
                                    update='yes')
                extras.cl_bye()

            # Fixing bad pixels
            print(" Fixpix Interpolating bad pixels %s" % ima, file=sys.stderr)
            for i in range(N):

                n = i + 1
                image = "%s[%s]" % (ima, n)
                mask = "%s_BPM/bpmcr_im%s.fits" % (ima, n)

                proto.fixpix(image,
                             mask,
                             linterp="INDEF",
                             cinterp="INDEF",
                             verbose="yes",
                             pixels="no")
                extras.cl_bye()
