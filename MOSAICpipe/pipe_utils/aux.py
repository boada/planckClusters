from __future__ import print_function
from __future__ import division
from builtins import map
from builtins import range
from builtins import object
from past.utils import old_div
import sys
import numpy
from numpy.ma import logical_and

def n2N(m, type='f'):
    x = numpy.resize(numpy.ravel(m), m.shape)
    return x.astype(type)


# Converts numpy.arrays to Numeric arrays
def N2n(m, type='Float'):
    x = numpy.resize(numpy.ravel(m), m.shape)
    return x.astype(type)


# Converts scalar value to Numeric arrays of 1-element
def s2N(m):
    return numpy.array([m])

def map2numpy(m):
    return numpy.array(list(map(float, list(m))))


def asnumpy(x):
    ''' Convert into numpy.most of types'''

    if type(x) is type(numpy.array([1])):
        return x
    # Check if list
    if type(x) is list:
        return numpy.array(x)
    if type(x) is int or type(x) is float:
        return numpy.array([x])

    try:
        if type(x) is type(numpy.array([1])):
            return numpy.asarray(x)
    except:
        print("# Numpy not found", file=sys.stderr)

    print("# Unknown type cannot convert to numpy.", file=sys.stderr)
    return x

###########################
#  Weighted statitsics
###########################

class stats(object):
    ''' a class to get the mean and std deviation fron a
    numpy.array for a given (optional) weight'''

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
        #self.var  = self.w*(self.x - self.mean)**2
        #self.var  = self.var.sum() / self.w.sum()
        self.var = old_div(1.0, self.w.sum())
        self.std = numpy.sqrt(self.var)
        return


def statsw(x, weight=None):
    ''' Same as the above class, but a function rather than class,
    it returns the mean and std, for an (optional) weight'''

    x = numpy.asarray(x)
    if weight is None:
        return x.mean(), x.stddev()

    elif (weight.sum() == 0):
        print("# warning: statsw weight sum is zero", file=sys.stderr)
        return x.mean(), x.stddev()
    else:
        weight = numpy.asarray(weight)
        mean = old_div((weight * x).sum(), weight.sum())
        #var  = ( weight*(x - mean)**2 ).sum() / weight.sum()
        var = old_div(1.0, weight.sum())
        std = numpy.sqrt(var)
        return mean, std


# Get the rms and variance of an numpy.array
def rms(x):
    ''' Get the rms and variance of an numpy.array'''
    var = (old_div((x - x.mean())**2, len(x))).sum()
    return numpy.sqrt(var)

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

# Sorting dictionary by values, Daniel Schult, 2004/01/23 Here's some
# code to sort the keys of a dictionary by their associated values.
# http://aspn.activestate.com/ASPN/Python/Cookbook/Recipe/52306


def sort_by_value(d):
    """ Returns the keys of dictionary d sorted by their values """
    items = list(d.items())
    backitems = [[v[1], v[0]] for v in items]
    backitems.sort()
    return [backitems[i][1] for i in range(0, len(backitems))]


# To sort string list numerically, rather than alphabetycally
# x = ['1','1000','3','15','120']
# x.sort(numeric_string)
# x = ['1', '3', '15', '120', '1000']
def numeric_string(x, y):
    return int(x) - int(y)

###################################################
# Routines to bin data and make histogram of data
###################################################


def bin_data(x, y, x1, x2, dx, center=None):
    ''' to bin data in y acording to x, bewteen x1 and x2
    at a step dx'''

    # Check if we want to use the center of the bin
    if center:
        dc = old_div(dx, 2.0)
    else:
        dc = 0
    ibin = (old_div((x - x1), dx)).astype(numpy.int)
    xbin = numpy.arange(x1, x2, dx) + dc
    ybin = numpy.zeros(len(xbin)).astype(numpy.float)

    # Fill in the bins according to ibin one by one
    # slow, but not if nbin is small
    for i in range(len(xbin)):
        ybin[i] = y[ibin == i].sum()
    return xbin, ybin


def bin_data_ave(x, y, x1, x2, dx, center=None):
    ''' to bin data in y acording to x, bewteen x1 and x2
    at a step dx'''

    # Check if we want to use the center of the bin
    if center:
        dc = old_div(dx, 2.0)
    else:
        dc = 0
    ibin = (old_div((x - x1), dx)).astype(numpy.int)
    xbin = numpy.arange(x1, x2, dx) + dc
    ybin = numpy.zeros(len(xbin)).astype(numpy.float)

    # Fill in the bins according to ibin one by one
    # slow, but not if nbin is small
    for i in range(len(xbin)):
        ybin[i] = old_div(y[ibin == i].sum(), len(y[ibin == i]))
    return xbin, ybin


def hist(x, x1, x2, dx, center=None):
    ''' Makes a histogram of the array x(n) between x1, x2 and at a
    step of with dx'''

    # Check if we want to use the center of the bin
    if center == 'yes':
        dc = old_div(dx, 2.0)
    elif center == 'left':
        dc = 0
    elif center == 'right':
        dc = dx
    else:
        dc = 0

    # trick to get the indx for the right
    ibin = (old_div((x - x1), dx)).astype(numpy.int)
    xbin = numpy.arange(x1, x2, dx) + dc
    nbin = numpy.zeros(len(xbin)).astype(numpy.float)

    # Fill in the bins according to ibin one by one
    # slow, but not if nbin is small
    for i in range(len(xbin)):
        nbin[i] = len(x[ibin == i])
    return xbin, nbin


def bin_by_data(x, y, x1, x2, dx, center=None):
    ''' to bin data in y acording to x, bewteen x1 and x2
    at a step dx --- NO SUM IS DONE!!!'''

    # Check if we want to use the center of the bin
    if center:
        dc = old_div(dx, 2.0)
    else:
        dc = 0

    #ibin = ( (x - x1)/dx).astype('Int')
    ibin = (old_div((x - x1), dx)).astype(numpy.int)
    xbin = numpy.arange(x1, x2, dx) + dc
    #ybin = numpy.zeros(len(xbin)).astype('Float')
    ybin = []

    # Fill in the bins according to ibin one by one
    # slow, but not if nbin is small
    for i in range(len(xbin)):
        in_bin = y[ibin == i]
        ybin.append(in_bin)
        #ybin[i] = y[ibin == i]
    return xbin, ybin


def FilterName(fitsfile):
    import os
    ''' Return filters name for ACS images, uses fitshead'''

    cmd1 = 'fitshead %s | grep FILTER1' % fitsfile
    cmd2 = 'fitshead %s | grep FILTER2' % fitsfile

    (stdout) = os.popen(cmd1)
    vals = stdout.readline().split("'")
    f1 = vals[1].strip()

    (stdout) = os.popen(cmd2)
    vals = stdout.readline().split("'")
    f2 = vals[1].strip()

    # Select the filter name
    if (f1 == "CLEAR1L"):
        FILTER = f2
    else:
        FILTER = f1
    return FILTER


def dec2sex(dec):

    dd = int(dec)
    mm = int(abs(dec - dd) * 60.)
    ss = (abs(dec - dd) * 60 - mm) * 60
    return (dd, mm, ss)


def log2(x):
    from numpy import log as ln
    return old_div(ln(x), ln(2))


def inpath(program, verb=None):
    """ Checks if program is in the user's path """
    import os.path
    for path in os.environ['PATH'].split(':'):
        if os.path.exists(os.path.join(path, program)):
            if verb:
                print("# program: %s found in: %s" % (program, os.path.join(
                    path, program)))
            return 1
    if verb: print("# program: %s NOT found in user's path " % program)
    return 0

# Simple Sigma clipping

def sigclip(x, Nsig=3.0, eps=1e-6, ids=None):

    xo = x.mean()
    xlo = x.mean() - Nsig * x.std(ddof=1)
    xhi = x.mean() + Nsig * x.std(ddof=1)
    #x = numpy.clip(x,xlo,xhi)

    idx = numpy.where(logical_and(x > xlo, x < xhi))
    xn = x[idx]

    if x.std(ddof=1) == 0.0:
        if ids:
            idx = numpy.indices(x.shape)[0]
            print("idx", idx)
            return idx
        else:
            return x

    i = 0
    while abs(1 - old_div(xo, xn.mean())) > eps:

        xo = xn.mean()
        xlo = xo - Nsig * xn.std(ddof=1)
        xhi = xo + Nsig * xn.std(ddof=1)

        idx = numpy.where(logical_and(x > xlo, x < xhi))
        xn = x[idx]

        i = i + 1

    if ids:
        return idx
    else:
        #return x
        return xn
