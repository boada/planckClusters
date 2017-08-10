#!/usr/bin/env python3

import glob
import os
import sys
from astropy.io.fits import getheader

# this is a specific file that is meant to be run at the top level of the proc
# directory to make the associations for all of the objects in the lower
# directories. I decided to do it this way because there are a lot of objects
# and I couldn't be bothered trying to run it for each one. I also wanted to
# make sure that it only worked on the mosaic data.

# made a couple small changes that might let it work with the `find` command
# from the command line.

USAGE = "usage:  " + os.path.split(sys.argv[0])[1] + " <rawdir> <outdir>"

dirs = [dirs for _, dirs, _ in os.walk('./')][0] # only want top level
# get the cwd
cwd = os.getcwd()
for d in dirs:
    os.chdir(cwd)
    rawdir = './{}'.format(d)
    outdir = '.'

    pattern = "*.fits"

# Go to the directory
    try:
        os.chdir(rawdir)
    except FileNotFoundError:
        continue

# Make the fits files data inventory
    full_list = glob.glob(pattern)
    if not len(full_list) > 0:
        pattern = '*.fz'
        full_list = glob.glob(pattern)

    filters = ('g', 'r', 'i', 'z', 'I', 'K')
    objects = []
    imalist = {}
    FILTER = {}
    EXPTIME = {}
    AIRMASS = {}
    tiles = []

    print("Found %s files, please wait... this will take a while" %
          len(full_list), file=sys.stderr)

    counter = 1
    for file in full_list:

        # make sure we are only working on mosaic data
        if not 'mosaic' and 'resampled' in file:
            continue

        print("Reading %-45s ... (%4s/%4s)" % (file, counter, len(full_list)),
            file=sys.stderr)
        try:
            header = getheader(file)
        except FileNotFoundError:
            counter += 1
            continue

        counter += 1
        # Weed out the calibration files we don't care about
        try:
            if header['prodtype'] != 'image':
                continue
        except KeyError:
            print("# OBSTYPE not present for %s" % file)
            print("# Trying next extention")
            try:
                header = getheader(file, ext=1)
                if header['OBSTYPE'] != 'object':
                    continue
                if header['PRODTYPE'] != 'image':
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

        TILE = OBJECT[0:]

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

