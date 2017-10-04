#!/usr/bin/env python

from __future__ import division
from past.utils import old_div
import os, sys
import tableio
import pylab

seds = sys.argv[1:]

pylab.figure(1)
fmax = 0.0
for sed in seds:

    (w, f) = tableio.get_data(sed, cols=(0, 1))
    #fmax = max(f.max(),fmax)
    f = old_div(f, f.max())
    #print fmax
    pylab.semilogx(w, f, label=sed)

pylab.xlim(1000, 5e4)
pylab.legend()
pylab.show()
pylab.close()
