#!/usr/bin/env python

import os, sys
import time

# Count time
tstart = time.time()

if len(sys.argv) < 2:
    print("usage: %s <tilename> <path>" % sys.argv[0])
    sys.exit(2)

tilename = sys.argv[1]
try:
    path = sys.argv[2]
except:
    path = os.path.join(os.environ['HOME'], "KittPeak-data/MOSAIC1/COMB")

# Order is important
filters = ('Red', 'Green', 'Blue')

cmd = 'stiff '
for filter in filters:
    file = "%s%s.fits " % (tilename, filter)
    fitsname = os.path.join(path, tilename, file)
    cmd = cmd + fitsname

confdir = os.path.join(path, tilename, 'stiff_MOS1.conf')
if os.path.exists(confdir):
    conf = confdir
else:
    MOSAIC1pipe = os.environ['MOSAIC1pipe']
    conf = os.path.join(MOSAIC1pipe, 'LIB/pars/stiff_MOS1.conf')

print("Will use conf: %s" % conf)
outtiff = "%s_RGB.tiff " % os.path.join(path, tilename, tilename)
cmd = cmd + "-c %s -OUTFILE_NAME %s" % (conf, outtiff)
print(cmd)
os.system(cmd)
print("# Done %s" % tilename)

#outjpeg = "%s.jpg "  % os.path.join(path,tilename,tilename)
#cmd = "convert -quality 100 %s %s" % (outtiff,outjpeg)
#print cmd
#os.system(cmd)
