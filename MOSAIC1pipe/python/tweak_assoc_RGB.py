#!/usr/bin/env python

import os,sys

filter   = {}
filters  = []
filelist = []
infiles  = []
exptime  = {}
airmass  = {}

# Make image list per filter
files        = {}
exptimes = {}

try:
    assocfile    = sys.argv[1]
    newassocfile = sys.argv[2]
except:
    print "%s <assocfile> <newassocfile>" % sys.argv[0]
    sys.exit()
    

print >>sys.stderr,"# Will read %s" % assocfile


# Read in the assoc file
for line in open(assocfile).readlines():
    
    if line[0] == "#":
        continue
    
    vals = line.split()
    fname           = os.path.basename(vals[0])
    infile          = vals[0]
    filter[fname]   = vals[1]

    exptime[fname]  = float(vals[2])
    airmass[fname]  = float(vals[3])
    
    #print fname, filter[fname], airmass[fname]
    
    # A list of the files
    filelist.append(fname)
    infiles.append(infile)

    if vals[1] not in filters:
        filters.append(vals[1])


# Re-pack
newfilters  = ['Blue','Green','Red']
newfilelist = {}

for f in newfilters:
    newfilelist[f] = []

for fname in filelist:

    f = filter[fname]

    if (f == 'g' or f == 'r') and fname not in newfilelist['Blue']:
        newfilelist['Blue'].append(fname)

    if (f == 'r' or f == 'i') and fname not in newfilelist['Green']:
        newfilelist['Green'].append(fname)

    if (f == 'i' or f == 'z') and fname not in newfilelist['Red']:
        newfilelist['Red'].append(fname)
        

# Now print all of the new information

print "# Will write new assocfile to %s" % (newassocfile)
onew = open(newassocfile,'w')


for filtername in newfilters:
    for fname in newfilelist[filtername]:
        onew.write("%s %-6s %5.1f %.3f\n" % (fname, filtername, exptime[fname], airmass[fname]))

onew.close()
