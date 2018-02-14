from __future__ import division
from __future__ import print_function
# Some extra functions for the BCS pipeline
from builtins import object
from past.utils import old_div
import string
from astropy.io import fits as pyfits
import os
import sys
import random
import numpy
from numpy.ma import logical_and


##################################
# String and list manipulation
##################################
def list_prepend(list, string):
    new = []
    for l in list:
        new.append(string + l)
    return new

def list_append(list, string, sep='.'):
    new = []
    for l in list:

        try:
            (l1, l2) = l.split(sep)
            #new.append(l1+string+sep+l2)
            new.append(l1 + string)
        except:
            new.append(l + string)
    return new


def list_trick(array, trick="-"):
    n = len(array)
    return list(trick * n).join(",")

def imlist(list, sep=".fits"):

    new = []
    list.sort()
    # Remove the .fits extension
    for l in list:
        try:
            (l1, l2) = l.split(sep)
            new.append(l1)
        except:
            new.append(l)

    return new.join(", ")

def imlist_append(list, s=""):

    new = []
    for l in list.split(","):
        new.append(l + s)

    return new.join(",")

def imlist_prepend(list, s=""):

    new = []
    for l in list.split(","):
        new.append(s + l)

    return new.join(",")

def imfile(list, sep=".fits"):

    tmp = "tmp%s.list" % int(random.uniform(1000, 9999))
    new = open(tmp, 'w')

    list.sort()
    # Remove the .fits extension
    for l in list:
        try:
            (l1, l2) = l.split(sep)
            new.write("%s\n" % l1)
        except:
            new.write("%s\n" % l)
    new.close()

    return "@" + tmp


def imfile_append(list, s="", sep=".fits"):

    # Make the tmpfile name
    tmp = "tmp%s.list" % int(random.uniform(1000, 9999))
    new = open(tmp, 'w')

    # Remove the .fits extension
    list.sort()
    for l in list:
        try:
            (l1, l2) = l.split(sep)
            new.write("%s%s\n" % (l1, s))
        except:
            new.write("%s%s\n" % (l, s))
    new.close()

    return "@" + tmp

########################################
# Read and write fits simple functions
#######################################


# To read in simple fits file
def rfits(filename):
    import pyfits

    ff = pyfits.open(filename, "readonly")
    data = ff[0].data
    hdr = ff[0].header
    ff.close()
    return data, hdr


# To write in a simple file
def writefits(array, filename):
    import pyfits
    import os
    # Writes fits files of a given array
    newfits = pyfits.HDUList()
    hdu = pyfits.PrimaryHDU()
    hdu.data = array
    newfits.append(hdu)

    # Remove old version of file before
    if os.path.isfile(filename):
        os.remove(filename)

    newfits.writeto(filename)
    newfits.close
    return


wfits = writefits


def fitsInfo_pyfits(fitsfile, ext=0, verb=None):
    """returns a the filter name and object type of the BCS fits file"""

    if verb:
        print("Getting header info for: %s" % fitsfile, file=sys.stderr)

    f = pyfits.open(fitsfile)
    filter = f[ext].header.get('FILTER')
    obstype = f[ext].header.get('OBSTYPE')
    f.close()
    del f

    if not filter:
        return None, obstype

    if len(filter.split()) > 1:
        filter = filter.split()[0]

    return filter, obstype


def fitsInfo(fitsfile, verb=None):
    """returns a the filter name and object type of the BCS fits file"""

    if verb:
        print("Getting header info for: %s" % fitsfile, file=sys.stderr)

    header = fitsheader(fitsfile)

    try:
        filter = header['FILTER']
    except:
        filter = None

    try:
        obstype = header['OBSTYPE']
    except:
        obstype = None

    try:
        exptime = header['EXPTIME']
    except:
        exptime = None

    try:
        Namps = header['NAMPS']
    except:
        Namps = None

    return filter, obstype, exptime, Namps


# To avoid problems with Iraf
def cl_bye(verb=None):
    from pyraf import iraf
    if verb:
        print(" Clearing iraf --> iraf.flpr()", file=sys.stderr)
    iraf.flpr()
    iraf.flpr()
    return


# Fits header information, faster and dirtier than pyfits
def fitsheader(fitsfile, comments=None):
    ''' Return fits header, uses fitshead'''
    # Check if we have fitshead

    check_exe('fitshead', verb=None)

    cmd = 'fitshead %s' % fitsfile
    (stdout) = os.popen(cmd)

    header = {}
    comment = {}

    for line in stdout.readlines():

        # Get the key value
        key = string.strip(line[0:8])

        # Ignore comments and blank spaces
        if key == "" or key == "COMMENT" or key == "END":
            continue

        # Split and strip value from comment
        vals = line[10:].split(' / ')
        strip_val = string.strip(vals[0].replace("'", ""))

        # Check if comments exist for the keyword
        if len(vals) > 1:
            comment[key] = string.strip(vals[1])
        else:
            comment[key] = ''

        # Case 1, if string starts with "'"
        if vals[0][0] == "'" or strip_val[0].isalpha():
            value = strip_val

        # Case 2, float or integer
        else:
            try:
                value = int(vals[0])
            except:
                try:
                    value = float(vals[0])
                except:
                    value = vals[0]  # just a string

        header[key] = value

    # If we require comments, give them back
    if comments:
        return header, comment
    else:
        return header


# Checks if program is in PATH
def check_exe(exe, verb="yes"):

    path = os.environ['PATH'].split(':')
    for p in path:

        f = os.path.join(p, exe)
        if os.path.isfile(f):
            if verb:
                print("# Found %s in %s" % (exe, f), file=sys.stderr)
            return f
    print("# ERROR: Couldn't find %s" % exe, file=sys.stderr)
    return


def elapsed_time(t1, text=''):
    import time
    t2 = time.time()
    hh = int(old_div((t2 - t1), 3600.))
    mm = int(old_div((t2 - t1), 60) - hh * 60)
    ss = (t2 - t1) - 60 * mm - 3600 * hh
    print("Elapsed time: %dh %dm %2.2fs %s" % (hh, mm, ss, text))
    print("Elapsed time: %dh %dm %2.2fs %s" % (hh, mm, ss, text),
          file=sys.stderr)
    #print >>sys.stderr,"Elapsed time: %dm %2.2fs" % ( int( (t2-t1)/60.), (t2-t1) - 60*int((t2-t1)/60.))
    return


def elapsed_time_str(t1):
    import time
    t2 = time.time()
    hh = int(old_div((t2 - t1), 3600.))
    mm = int(old_div((t2 - t1), 60) - hh * 60)
    ss = (t2 - t1) - 60 * mm - 3600 * hh
    return "Elapsed time: %dh %dm %2.4fs" % (hh, mm, ss)


# Propagate some fits headers keyword into children files
def copy_headkeys(parent, child, keys, after='DATE', verb=None, n=0):

    import time
    ''' Propagate list of fits keywords using pyfits and fitsheaderlibrary '''

    import pyfits

    # Read in the parent information
    p = pyfits.open(parent)
    cards = p[n].header.ascardlist()

    # Open thr child info, update mode
    f = pyfits.open(child, mode="update")  # open a FITS file
    hdr = f[n].header  # the primary HDU header

    for key in keys:
        cards[key].verify('fix')
        hdr.update(key, cards[key].value, cards[key].comment, after=after)

    hdr.add_comment('-------------------------------------', after=after)
    hdr.add_comment('Date: %s' % time.ctime(), after=after)
    hdr.add_comment('Keywords propagated from %s' % parent, after=after)
    hdr.add_comment('BCS-Rutgers pipeline by F. Menanteau', after=after)
    hdr.add_comment('-------------------------------------', after=after)
    #f.verify('fix')

    # Close the file
    f.close()
    return


# Transform form Dec to Sexagesimal hh:mm:ss format
def dec2sex(hms):

    hh = int(hms)
    mm = int((hms - hh) * 60)
    ss = ((hms - hh) * 60 - mm) * 60

    if abs(mm) < 10:
        mm = "0%d" % abs(mm)
    else:
        mm = "%2d" % abs(mm)

    if abs(ss) < 10:
        ss = "0%.2f" % abs(ss)
    else:
        ss = "%5.2f" % abs(ss)

    #return "%2d:%d:%.2f" % (hh,abs(mm),abs(ss))
    return "%s:%s:%s" % (hh, mm, ss)


# Parses the header of a SExtractor catalog and stores the column
# number (in python' order starting from zero) for every keyword columns
# It returs a dictionary with the column number for evry keyword.
def SEx_head(catalog, verb='yes'):

    if verb:
        print("\r Parsing SEx head for:", catalog, file=sys.stderr)

    # Dictionary with column numbers
    SExcols = {}

    # Read the SExtractor catalog
    for line in open(catalog).readlines():

        if line[0] != '#':
            break

        if line[:2] == "##":
            continue

        try:
            line = line.strip()
            vals = line.split()
            col = vals[1]
            key = vals[2]
            SExcols[key] = int(col) - 1
            if verb:
                print("# %-20s %s" % (key, SExcols[key] + 1), file=sys.stderr)
        except:
            continue

    return SExcols

def deNAN(a, value=0.0):
    nans = numpy.logical_not(numpy.less(a, 0.0) + numpy.greater_equal(
        a, 0.0))
    return numpy.where(nans, value, a)

# Sigma clipping
def sigclip(x, Nsig=2.0, eps=1e-6, ids=None):

    xo = x.mean()
    xlo = x.mean() - Nsig * x.stddev()
    xhi = x.mean() + Nsig * x.stddev()
    #x = numpy.clip(x,xlo,xhi)

    idx = numpy.where(logical_and(x > xlo, x < xhi))
    xn = x[idx]

    if x.stddev() == 0.0:
        if ids:
            idx = numpy.indices(x.shape)[0]
            print("idx", idx)
            return idx
        else:
            return x

    i = 0
    while abs(1 - old_div(xo, xn.mean())) > eps:

        xo = xn.mean()
        xlo = xo - Nsig * xn.stddev()
        xhi = xo + Nsig * xn.stddev()

        idx = numpy.where(logical_and(x > xlo, x < xhi))
        xn = x[idx]

        i = i + 1

    if ids:
        return idx
    else:
        #return x
        return xn

# Sigma clipping


def sigclip_w(x, wt=None, Nsig=2.0, eps=1e-6, ids=None):

    i = 1
    s = stats(x, wt)
    xo = s.mean
    std = s.std
    xlo = s.mean - Nsig * std
    xhi = s.mean + Nsig * std

    idx = numpy.where(logical_and(x > xlo, x < xhi))
    xn = x[idx]
    wn = wt[idx]

    i = 0
    #print i,xo,stats(xn,wn).mean,stats(xn,wn).std, xn.mean(), xn.stddev(), len(xn),xlo,xhi
    while abs(1 - old_div(xo, stats(xn, wn).mean)) > eps:

        #print i,xo,stats(xn,wn).mean,stats(xn,wn).std, xn.mean(), xn.stddev(), len(xn),xlo,xhi

        s = stats(xn, wn)
        xo = s.mean
        xlo = xo - Nsig * s.std
        xhi = xo + Nsig * s.std

        # Clip it
        idx = numpy.where(logical_and(x > xlo, x < xhi))
        xn = x[idx]
        wn = wt[idx]
        i = i + 1

    if ids:
        return idx
    else:
        return xn, wn


def log2(x):
    from math import log
    return old_div(log(x), log(2))


def isflag(x, flag):
    exp = []
    while x > 0:
        p = int(log2(x))
        x = x - 2**p
        exp.append(p)

    if flag in exp:
        return True
    else:
        return False

    return

###########################
#  Weighted statitsics
###########################

class stats(object):
    ''' a class to get the mean and std deviation fron a
    numpy array for a given (optional) weight'''

    def __init__(self, x, weight=None):

        self.x = numpy.asarray(x)
        self.w = weight
        if self.w is None:
            self.mean = x.mean()
            self.std = x.stddev()
            self.var = self.std**2
        else:
            self.w = numpy.asarray(self.w)
            self._stats_w()
        return

    def _stats_w(self):
        self.mean = old_div((self.w * self.x).sum(), self.w.sum())
        self.var = self.w * (self.x - self.mean)**2
        self.var = old_div(self.var.sum(), self.w.sum())
        #self.var  = 1.0 / self.w.sum()
        #self.var  = 1.0 / (self.w**2).sum()
        self.std = numpy.sqrt(self.var)
        return


def asnumpy(x):
    ''' Convert into numpy most of types'''

    if type(x) is type(numpy.array([1])):
        return x
    # Check if list
    if type(x) is list:
        return numpy.array(x)
    if type(x) is int or type(x) is float:
        return numpy.array([x])

    try:
        import numpy
        if type(x) is type(numpy.array([1])):
            return numpy.asarray(x)
    except:
        print("# Numpy not found", file=sys.stderr)

    print("# Unknown type cannot convert to numpy", file=sys.stderr)
    return x

