from __future__ import print_function


# ----------------------------------------------------
#         Functions to setIraf's parameters
# ----------------------------------------------------
def setIraf(verb=None):

    extras.cl_bye(verb)
    # Set vars
    iraf.set(imclobber="yes")
    iraf.set(clobber="yes")
    iraf.set(Verbose='yes')
    # Initialize some tasks
    set_instrument()
    set_mscred()
    set_ccdproc()
    return


def set_instrument():

    print(" setting up instrument params", file=sys.stderr)
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


def set_mscred():

    print("initializing mscred params", file=sys.stderr)
    mscred.pixeltype = "real real"
    mscred.verbose = 'yes'
    mscred.logfile = "logfile"
    mscred.plotfile = ""
    #mscred.backup    = "once"
    mscred.backup = "none"
    mscred.bkuproot = "Raw/"
    mscred.instrument = "mscdb$noao/ctio/4meter/Mosaic2.dat"
    mscred.ampfile = "amps"
    mscred.ssfile = "subsets"
    mscred.im_bufsize = 4.0
    mscred.graphics = "stdgraph"
    mscred.cursor = ""
    mscred.version = "V4.8: May 11= 2004"
    mscred.mode = "ql"
    return


def set_ccdproc():
    '''Initialize ccdproc'''

    print(" initializing mscred.ccdproc params", file=sys.stderr)

    mscred.ccdproc.images = ''
    mscred.ccdproc.output = ''
    mscred.ccdproc.bpmasks = ''
    mscred.ccdproc.ccdtype = ''
    mscred.ccdproc.noproc = 'no'

    mscred.ccdproc.xtalkcor = 'yes'
    mscred.ccdproc.fixpix = 'no'
    mscred.ccdproc.overscan = 'yes'
    mscred.ccdproc.trim = 'yes'
    mscred.ccdproc.zerocor = 'no'
    mscred.ccdproc.darkcor = 'no'
    mscred.ccdproc.flatcor = 'no'
    mscred.ccdproc.sflatcor = 'no'
    mscred.ccdproc.split = 'no'
    mscred.ccdproc.merge = 'no'

    mscred.ccdproc.xtalkfile = '!xtalkfil'
    mscred.ccdproc.fixfile = 'BPM'
    mscred.ccdproc.saturation = 'INDEF'
    mscred.ccdproc.sgrow = 0
    mscred.ccdproc.bleed = 'INDEF'
    mscred.ccdproc.btrail = 20
    mscred.ccdproc.bgrow = 0
    mscred.ccdproc.biassec = '!biassec'
    mscred.ccdproc.trimsec = '!trimsec'
    mscred.ccdproc.zero = 'Zero'
    mscred.ccdproc.dark = 'Dark'
    mscred.ccdproc.flat = 'Dflat*'
    mscred.ccdproc.sflat = 'Sflat*'
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
