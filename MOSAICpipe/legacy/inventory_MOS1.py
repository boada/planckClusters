#!/usr/bin/env python3

import glob
import os
import sys
from astropy.io.fits import getheader

USAGE = "usage:  " + os.path.split(sys.argv[0])[1] + " <rawdir> <outdir>"
try:
    rawdir = sys.argv[1]
    outdir = sys.argv[2]
except:
    sys.exit(USAGE)

pattern = "*.fits"

# get the cwd
cwd = os.getcwd()

# Go to the directory
os.chdir(rawdir)

# Make the fits files data inventory
full_list = glob.glob(pattern)
if not len(full_list) > 0:
    pattern = '*.fz'
    full_list = glob.glob(pattern)

filters = ('z', 'i', 'r', 'g')
objects = []
imalist = {}
FILTER = {}
EXPTIME = {}
AIRMASS = {}
tiles = []

print("Found %s files, please wait... this will take a while" % len(full_list),
      file=sys.stderr)

counter = 1
for file in full_list:

    print("Reading %-45s ... (%4s/%4s)" % (file, counter, len(full_list)),
          file=sys.stderr)

    header = getheader(file)
    counter = counter + 1
    # Weed out the calibration files we don't care about
    try:
        if header['OBSTYPE'] != 'object':
            continue
    except KeyError:
        print("# OBSTYPE not present for %s" % file)
        print("# Trying next extention")
        try:
            header = getheader(file, ext=1)
            if header['OBSTYPE'] != 'object':
                continue
        except KeyError:
            print("# OBSTYPE not present for %s" % file)
        except IndexError:
            print("# %s only has a single extention" % file)

    # Keep the values
    try:
        OBJECT = header['OBJECT']  # [:-1]
    except:
        print("OBJECT KEY NOT FOUND FOR:", file)
        OBJECT = header['FILENAME'][0:11]

    TILE = OBJECT[0:-1]

    try:
        FNAME = header['FILTER'][0]
        FILTER[file] = header['FILTER'][0]
        AIRMASS[file] = header['AIRMASS']
    except:
        continue

    try:
        EXPTIME[file] = header['EXPTIME']
    except:
        EXPTIME[file] = "undef"

    # OBJECT list per filter
    #if OBJECT not in objects:
    #    objects.append(OBJECT)
    #    imalist[OBJECT] = {}
    #    for filter in filters:
    #        imalist[OBJECT][filter] = []

    #imalist[OBJECT][FNAME].append(file)

    # TILE list per filter
    if TILE not in tiles:
        tiles.append(TILE)
        imalist[TILE] = {}
        for filter in filters:
            imalist[TILE][filter] = []

    imalist[TILE][FNAME].append(file)
    header = None

# Write them out
print(" Will write results to: %s" % outdir, file=sys.stderr)
os.chdir(outdir)
objects.sort()
for TILE in tiles:

    o = open(TILE + ".assoc", "w")
    print("# %s" % TILE)
    filters = list(imalist[TILE].keys())
    filters.sort()
    for filter in filters:
        for file in imalist[TILE][filter]:
            print(file, filter, EXPTIME[file], AIRMASS[file])
            o.write("%s %s %s %s\n" %
                    (file, filter, EXPTIME[file], AIRMASS[file]))
    o.close()

    #for OBJECT,filenames in imalist.items():
    #    print OBJECT
    #
    #    for filter FILTER
    #    for file in filenames:
    #        print "\t %20s %s %s" % (file,FILTER[file],EXPTIME[file])
