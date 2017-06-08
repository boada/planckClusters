#!/usr/bin/env python

import os,sys
import numpy
import tableio
import pylab


seds = sys.argv[1:]



pylab.figure(1)
fmax = 0.0
for sed in seds:

    (w,f) = tableio.get_data(sed,cols=(0,1))
    #fmax = max(f.max(),fmax)
    f = f/f.max()
    #print fmax
    pylab.semilogx(w,f)


pylab.xlim(1000,5e4)
pylab.show()
pylab.close()
    
