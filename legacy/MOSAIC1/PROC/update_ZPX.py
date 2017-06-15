#!/usr/bin/env python

import pyfits
import os,sys
import glob


newproj ="TPV"

fnames = sys.argv[1:]

for fname in fnames:

    hdulist = pyfits.open(fname, mode='update')
    Nima = len(hdulist)-1
    for n in range(Nima):
        k = n+1
        hdulist[k].header['CTYPE1'] = 'RA---%s' % newproj 
        hdulist[k].header['CTYPE2'] = 'DEC--%s' % newproj
        
    hdulist.flush()  # changes are written
    hdulist.close()
    print "# Updated %s to %s" % (fname, newproj)
