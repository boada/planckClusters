#!/usr/bin/env python3

from __future__ import print_function
from builtins import str
from builtins import map
from builtins import object
import os
import sys
import glob
from SEx_reader import SEx_reader, SEx_reader_multi
import tableio
import re
import astrometry
import numpy as np


class bcs_catalogs(object):
    """Reads in the BCS catalogs and stores the information"""

    def __init__(self,
                 tile,
                 path,
                 det_filter='i',
                 filters=('g', 'r', 'i', 'z'),
                 verb=None,
                 probs='yes',
                 bpz=True):

        if verb:
            print("\r Initializing class ", file=sys.stderr)

        # Set up the standard apsis output file names
        self.path = path
        self.tile = tile
        self.verb = verb
        self.det_filter = det_filter
        self.filters = filters

        self._get_names()
        self._read_SExcats()

        if os.path.exists(self.bpzcat):
            self._read_bpz()
        else:
            print("No BPZ file to read")

        if os.path.exists(self.colorcat):
            self._read_multicolor()
        else:
            print("No ColorCat file to read")

        if probs:
            try:
                self.read_probs()
                print("\tReading ", self.bpzprobs, file=sys.stderr)
            except FileNotFoundError:
                print("No BPZ probs to read")

            try:
                self.read_probs_flat()
                print("\tReading ", self.bpzprobs_flat, file=sys.stderr)
            except FileNotFoundError:
                print("No BPZ flat probs to read")

        if os.path.exists(self.dustcat):
            self.read_dust()
        else:
            print("No DUST file to read")

        self.xyID = None
        self.rdID = None
        return

    def _get_names(self):

        self.catname = {}
        self.imaname = {}

#        self.bpzcat = os.path.join(self.path, self.tile, self.tile + ".bpz")
#        self.bpzprobs = os.path.join(self.path, self.tile,
#                                     self.tile + ".probs")
#        self.colorcat = os.path.join(self.path, self.tile,
#                                     self.tile + ".color")
#        self.dustcat = os.path.join(self.path, self.tile, self.tile + ".dust")
#        self.detcat = os.path.join(self.path, self.tile,
#                                   self.tile + self.det_filter + ".cat")
#        self.detIma = os.path.join(self.path, self.tile,
#                                   self.tile + self.det_filter + ".fits")

#        self.sci_images = glob.glob(os.path.join(self.path, self.tile,
#                                                 self.tile + "[g,r,i,z].fits"))
#        self.wgt_images = glob.glob(os.path.join(
#            self.path, self.tile, self.tile + "[g,r,i,z]_weight.fits"))

        self.bpzcat = os.path.join(self.path, self.tile + ".bpz")
        self.bpzprobs = os.path.join(self.path, self.tile + ".probs")
        self.colorcat = os.path.join(self.path, self.tile + ".color")
        self.dustcat = os.path.join(self.path, self.tile + ".dust")
        self.detcat = os.path.join(self.path, self.tile + self.det_filter + ".cat")
        self.detIma = os.path.join(self.path, self.tile + self.det_filter + ".fits")

        self.sci_images = glob.glob(os.path.join(self.path, self.tile +
                                                 "[g,r,i,z].fits"))
        self.wgt_images = glob.glob(os.path.join(self.path, self.tile +
                                                 "[g,r,i,z]_weight.fits"))
 
        for filter in self.filters:
            self.catname[filter] = os.path.join(self.path, self.tile + filter +
                                                "_cal.cat")
            self.imaname[filter] = os.path.join(self.path, self.tile + filter +
                                                ".fits")

        # And for alt probs flat prior == ML
        self.bpzprobs_flat = os.path.join(self.path, self.tile + ".probs_flat")
#            self.catname[filter] = os.path.join(self.path, self.tile,
#                                                self.tile + filter + ".cat")
#            self.imaname[filter] = os.path.join(self.path, self.tile,
#                                                self.tile + filter + ".fits")

        # And for alt probs flat prior == ML
#        self.bpzprobs_flat = os.path.join(self.path, self.tile,
#                                          self.tile + ".probs_flat")
        return

    def _read_SExcats(self):

        self.cat = {}
        self.SExcols = {}

        for filter in self.filters:

            if self.verb:
                print("\r Reading:", self.catname[filter], file=sys.stderr)
            c = None
            try:
                c = SEx_reader(self.catname[filter], verb=self.verb)
                self.cat[filter] = c.cat
                self.SExcols[filter] = c.SExcols
            except FileNotFoundError:
                continue

        return

    def _read_bpz(self):

        if self.verb:
            print("\r Reading:", self.bpzcat, file=sys.stderr)

        self.Z_B = {}
        self.Z_B_MIN = {}
        self.Z_B_MAX = {}
        self.T_B = {}
        self.ODDS = {}
        self.Z_ML = {}
        self.T_ML = {}
        self.CHI2 = {}

        for line in open(self.bpzcat).readlines():

            vals = line.split()

            if vals[0][0] == "#":
                continue

            ID = str(vals[0])
            self.Z_B[ID] = float(vals[1])
            self.Z_B_MIN[ID] = float(vals[2])
            self.Z_B_MAX[ID] = float(vals[3])
            self.T_B[ID] = float(vals[4])
            self.ODDS[ID] = float(vals[5])
            self.Z_ML[ID] = float(vals[6])
            self.T_ML[ID] = float(vals[7])
            self.CHI2[ID] = float(vals[8])

        return

    def read_dust(self, filters=('g', 'r', 'i', 'z')):

        if self.verb:
            print("\r Reading:", self.dustcat, file=sys.stderr)

        self.Xcorr = {}
        self.XcorrErr = {}
        for filter in filters:
            self.Xcorr[filter] = {}
            self.XcorrErr[filter] = {}

        for line in open(self.dustcat).readlines():

            vals = line.split()

            if vals[0][0] == "#":
                continue

            ID = str(vals[0])
            self.Xcorr['g'][ID] = float(vals[1])
            self.XcorrErr['g'][ID] = float(vals[2])
            self.Xcorr['i'][ID] = float(vals[3])
            self.XcorrErr['i'][ID] = float(vals[4])
            self.Xcorr['r'][ID] = float(vals[5])
            self.XcorrErr['r'][ID] = float(vals[6])
            try:
                self.Xcorr['z'][ID] = float(vals[7])
                self.XcorrErr['z'][ID] = float(vals[8])
            except KeyError:
                pass

        return

    def read_probs(self):

        # The reg expresion to compile
        regexp_point = re.compile(r"arange\("
                                  r"(?P<z1>[0-9]+.[0-9]+),"
                                  r"(?P<z2>[0-9]+.[0-9]+),"
                                  r"(?P<dz>[0-9]+.[0-9]+)\)")

        # Dictionary of probability arrays
        self.PROBS = {}
        for line in open(self.bpzprobs).readlines():

            point = regexp_point.search(line)

            # Extract the information if a point was selected
            if point:
                z1 = float(point.group('z1'))
                z2 = float(point.group('z2'))
                dz = float(point.group('dz'))
                self.zprobs = np.arange(z1, z2, dz)
                self.probs_header = '# ID  p_bayes(z)  where z=arange(%.4f,%.4f,%.4f)' % (
                    z1, z2, dz)

            fields = line.split()
            if fields[0][0] == "#":
                continue
            ID = fields[0]
            self.PROBS[ID] = np.asarray(list(map(float, fields[1:])))

        return

    def read_probs_flat(self):

        # The reg expresion to compile
        regexp_point = re.compile(r"arange\("
                                  r"(?P<z1>[0-9]+.[0-9]+),"
                                  r"(?P<z2>[0-9]+.[0-9]+),"
                                  r"(?P<dz>[0-9]+.[0-9]+)\)")

        # Dictionary of probability arrays
        self.PROBS_FLAT = {}
        for line in open(self.bpzprobs_flat).readlines():

            point = regexp_point.search(line)

            # Extract the information if a point was selected
            if point:
                z1 = float(point.group('z1'))
                z2 = float(point.group('z2'))
                dz = float(point.group('dz'))
                self.zprobs_flat = np.arange(z1, z2, dz)
                self.probs_header = '# ID  p_bayes(z)  where z=arange(%.4f,%.4f,%.4f)' % (
                    z1, z2, dz)

            fields = line.split()
            if fields[0][0] == "#":
                continue
            ID = fields[0]
            self.PROBS_FLAT[ID] = np.asarray(list(map(float, fields[1:])))

        return

    def _read_multicolor(self):

        if self.verb:
            print("\r Reading:", self.colorcat, file=sys.stderr)

        # Parse the header
        self._Color_head()

        self.MAG_BPZ = {}
        self.MAG_BPZERR = {}

        for filter in self.filters:
            self.MAG_BPZ[filter] = {}
            self.MAG_BPZERR[filter] = {}

        for line in self.colorstr:

            vals = line.split()
            if vals[0][0] == "#":
                continue

            ID = str(vals[0])

            for filter in self.filters:

                try:
                    key1 = filter + "_MOSAICII_MAG_ISO"
                    key2 = filter + "_MOSAICII_MAGERR_ISO"
                    self.MAG_BPZ[filter][ID] = float(vals[self.Ccols[key1]])
                    self.MAG_BPZERR[filter][ID] = float(vals[self.Ccols[key2]])
                except KeyError:
                    key1 = filter + "_KittPeak_MAG_ISO"
                    key2 = filter + "_KittPeak_MAGERR_ISO"
                    self.MAG_BPZ[filter][ID] = float(vals[self.Ccols[key1]])
                    self.MAG_BPZERR[filter][ID] = float(vals[self.Ccols[key2]])

#                try:
#                    key1 = filter + "_SDSS_MAG_ISO"
#                    key2 = filter + "_SDSS_MAGERR_ISO"
#                    self.MAG_BPZ[filter][ID] = float(vals[self.Ccols[key1]])
#                    self.MAG_BPZERR[filter][ID] = float(vals[self.Ccols[key2]])
#                except KeyError:
#                    key1 = filter + "_MOSAICII_MAG_ISO"
#                    key2 = filter + "_MOSAICII_MAGERR_ISO"
#                    self.MAG_BPZ[filter][ID] = float(vals[self.Ccols[key1]])
#                    self.MAG_BPZERR[filter][ID] = float(vals[self.Ccols[key2]])

                #print filter,float(vals[self.Ccols[key1]]),float(vals[self.Ccols[key2]])

        return

    def _Color_head(self):

        if self.verb:
            print("\r Parsing Header for:", self.colorcat, file=sys.stderr)

        self.Ccols = {}

        # Read the detection catalog
        self.colorstr = open(self.colorcat).readlines()

        for line in self.colorstr:

            if line[0] != '#':
                break

            if line[:2] == "##":
                continue

            try:
                line = line.strip()
                vals = line.split()
                col = vals[1]
                key = vals[2]
                self.Ccols[key] = int(col) - 1
                #if self.verb:
                #    print >>sys.stderr, "# %-20s %s" % (key,self.Ccols[key]+1)
            except:
                continue
        if self.verb:
            print("# Read %s columns" % len(list(self.Ccols.values())),
                  file=sys.stderr)
        return

    def get_xyID(self, filter='i'):

        # Columns to read in as array
        col_keys = ['NUMBER', 'X_IMAGE', 'Y_IMAGE']
        col_list = []

        # Tuple with required columns
        for key in col_keys:
            col_list.append(self.SExcols[filter][key])

        col_tuple = tuple(col_list)

        (self.NUMBER, self.X_IMAGE,
         self.Y_IMAGE) = tableio.get_data(self.detcat, col_tuple)

        self.xyID = 1
        return

    def get_rdID(self, filter='i'):

        # Columns to read in as array
        col_keys = ['NUMBER', 'X_WORLD', 'Y_WORLD']
        col_list = []

        # Tuple with required columns
        for key in col_keys:
            col_list.append(self.SExcols[filter][key])

        col_tuple = tuple(col_list)

        (self.NUMBER, self.X_WORLD,
         self.Y_WORLD) = tableio.get_data(self.detcat, col_tuple)

        self.rdID = 1
        return

    def nearest(self, x, y):

        if not self.xyID:
            self.get_xyID()

        # Make the distance array
        D = np.sqrt((x - self.X_IMAGE)**2 + (y - self.Y_IMAGE)**2)
        idx = D.argmin()
        ID = int(self.NUMBER[idx])
        return str(ID)

    def nearest_rd(self, x, y):

        if not self.rdID:
            self.get_rdID()

        # Make the distance array
        dmin = 1e20
        X = self.X_WORLD
        Y = self.Y_WORLD
        N = self.NUMBER
        D = np.sqrt((x - X)**2 + (y - Y)**2)
        idx = D.argmin()
        if D[idx] < dmin:
            ID = str(int(N[idx]))
            dmin = min(dmin, D[idx])

        return str(ID)

    def near_rd(self, x, y, dmin, out_d=None):

        # dmin in arcmin
        if not self.rdID:
            self.get_rdID()

        IDsel = None
        dsel = None

        IDsel = []
        dsel = {}  # distance to seld
        # Make the distance array

        X = self.X_WORLD
        Y = self.Y_WORLD
        N = self.NUMBER
        #D = 60.0*np.sqrt( (x-X)**2 + (y-Y)**2)
        D = 60.0 * astrometry.circle_distance(x, y, X, Y,
                                              units='deg')  # in arcmins
        idx, = np.where(D < dmin)

        for i in idx:
            ID = str(int(N[i]))
            IDsel.append(ID)
            dsel[ID] = D[i]

        if out_d:
            return IDsel, dsel

        return IDsel


class bcs_catalogs_tiles(object):
    """Reads in the BCS catalogs and stores the information"""

    def __init__(self,
                 tiles,
                 path,
                 det_filter='i',
                 filters=('g', 'r', 'i', 'z'),
                 verb=None,
                 probs='yes'):

        if verb:
            print("\r Initializing class ", file=sys.stderr)

        # Set up the standard apsis output file names
        self.path = path
        self.tiles = tiles
        self.verb = verb
        self.det_filter = det_filter
        self.filters = filters

        self.seenBPZ = None
        self.seenColor = None
        self.seenProbs = None
        self.seenDust = None

        # Read in all the catalogs
        self.read_SExcats()

        for self.tile in self.tiles:

            self.get_names()
            self.read_bpz()

            if probs:
                try:
                    print("Reading ", self.bpzprobs)
                    self.read_probs()
                except FileNotFoundError:
                    print("No BPZ probs to read")

            self.read_multicolor()

            if os.path.exists(self.dustcat):
                self.read_dust()
            else:
                print("No DUST file to read")

        self.xyID = None
        self.rdID = None
        return

    def get_names(self):

        self.catname = {}
        self.imaname = {}

        self.bpzcat = os.path.join(self.path, self.tile, self.tile + ".bpz")
        self.bpzprobs = os.path.join(self.path, self.tile,
                                     self.tile + ".probs")
        self.colorcat = os.path.join(self.path, self.tile,
                                     self.tile + ".color")
        self.dustcat = os.path.join(self.path, self.tile, self.tile + ".dust")
        self.detcat = os.path.join(self.path, self.tile,
                                   self.tile + self.det_filter + ".cat")
        self.detIma = os.path.join(self.path, self.tile,
                                   self.tile + self.det_filter + ".fits")

        self.sci_images = glob.glob(os.path.join(self.path, self.tile,
                                                 self.tile + "[g,r,i,z].fits"))
        self.wgt_images = glob.glob(os.path.join(
            self.path, self.tile, self.tile + "[g,r,i,z]_weight.fits"))

        return

    def read_SExcats(self):

        self.cat = {}
        self.SExcols = {}
        for filter in self.filters:

            # Initialize one per filter
            c = None
            c = SEx_reader_multi(verb=self.verb)
            for tile in self.tiles:

                catalog = os.path.join(self.path, tile, tile + filter + ".cat")
                if self.verb:
                    print("\r Reading:", catalog, file=sys.stderr)

                c.read_catalog(catalog, preID=tile)

            self.cat[filter] = c.cat
            self.SExcols[filter] = c.SExcols

        return

    def read_bpz(self):

        if self.verb:
            print("\r Reading:", self.bpzcat, file=sys.stderr)

        if not self.seenBPZ:
            self.Z_B = {}
            self.Z_B_MIN = {}
            self.Z_B_MAX = {}
            self.T_B = {}
            self.ODDS = {}
            self.Z_ML = {}
            self.T_ML = {}
            self.CHI2 = {}
            self.seenBPZ = 1

        for line in open(self.bpzcat).readlines():

            vals = line.split()

            if vals[0][0] == "#":
                continue

            ID = self.tile + "_" + str(vals[0])
            self.Z_B[ID] = float(vals[1])
            self.Z_B_MIN[ID] = float(vals[2])
            self.Z_B_MAX[ID] = float(vals[3])
            self.T_B[ID] = float(vals[4])
            self.ODDS[ID] = float(vals[5])
            self.Z_ML[ID] = float(vals[6])
            self.T_ML[ID] = float(vals[7])
            self.CHI2[ID] = float(vals[8])

        return

    def read_dust(self, filters=('g', 'r', 'i', 'z')):

        if self.verb:
            print("\r Reading:", self.dustcat, file=sys.stderr)

        # Initialize only once
        if not self.seenDust:
            self.Xcorr = {}
            self.XcorrErr = {}
            for filter in filters:
                self.Xcorr[filter] = {}
                self.XcorrErr[filter] = {}
            self.seenDust = 1

        for line in open(self.dustcat).readlines():

            vals = line.split()

            if vals[0][0] == "#":
                continue

            ID = self.tile + "_" + str(vals[0])
            self.Xcorr['g'][ID] = float(vals[1])
            self.XcorrErr['g'][ID] = float(vals[2])
            self.Xcorr['i'][ID] = float(vals[3])
            self.XcorrErr['i'][ID] = float(vals[4])
            self.Xcorr['r'][ID] = float(vals[5])
            self.XcorrErr['r'][ID] = float(vals[6])
            self.Xcorr['z'][ID] = float(vals[7])
            self.XcorrErr['z'][ID] = float(vals[8])

        return

    def read_multicolor(self):

        if self.verb:
            print("\r Reading:", self.colorcat, file=sys.stderr)

        if not self.seenColor:

            self.MAG_BPZ = {}
            self.MAG_BPZERR = {}

            for filter in self.filters:
                self.MAG_BPZ[filter] = {}
                self.MAG_BPZERR[filter] = {}
            self.seenColor = 1

        # Parse the header
        self.Color_head()

        for line in self.colorstr:

            vals = line.split()
            if vals[0][0] == "#":
                continue

            ID = self.tile + "_" + str(vals[0])

            for filter in self.filters:
                try:
                    key1 = filter + "_MOSAICII_MAG_ISO"
                    key2 = filter + "_MOSAICII_MAGERR_ISO"
                    self.MAG_BPZ[filter][ID] = float(vals[self.Ccols[key1]])
                    self.MAG_BPZERR[filter][ID] = float(vals[self.Ccols[key2]])
                except KeyError:
                    key1 = filter + "_KittPeak_MAG_ISO"
                    key2 = filter + "_KittPeak_MAGERR_ISO"
                    self.MAG_BPZ[filter][ID] = float(vals[self.Ccols[key1]])
                    self.MAG_BPZERR[filter][ID] = float(vals[self.Ccols[key2]])

#                try:
#                    key1 = filter + "_SDSS_MAG_ISO"
#                    key2 = filter + "_SDSS_MAGERR_ISO"
#                    self.MAG_BPZ[filter][ID] = float(vals[self.Ccols[key1]])
#                    self.MAG_BPZERR[filter][ID] = float(vals[self.Ccols[key2]])
#                except KeyError:
#                    key1 = filter + "_MOSAICII_MAG_ISO"
#                    key2 = filter + "_MOSAICII_MAGERR_ISO"
#                    self.MAG_BPZ[filter][ID] = float(vals[self.Ccols[key1]])
#                    self.MAG_BPZERR[filter][ID] = float(vals[self.Ccols[key2]])

                #print filter,float(vals[self.Ccols[key1]]),float(vals[self.Ccols[key2]])

        return

    def Color_head(self):

        if self.verb:
            print("\r Parsing Header for:", self.colorcat, file=sys.stderr)

        self.Ccols = {}

        # Read the detection catalog
        self.colorstr = open(self.colorcat).readlines()

        for line in self.colorstr:

            if line[0] != '#':
                break

            if line[:2] == "##":
                continue

            try:
                line = line.strip()
                vals = line.split()
                col = vals[1]
                key = vals[2]
                self.Ccols[key] = int(col) - 1
                #if self.verb:
                #    print >>sys.stderr, "# %-20s %s" % (key,self.Ccols[key]+1)
            except:
                continue
        if self.verb:
            print("# Read %s columns" % len(list(self.Ccols.values())),
                  file=sys.stderr)
        return

    def get_rdID(self, filter='i'):

        self.NUMBER = {}
        self.X_WORLD = {}
        self.Y_WORLD = {}

        # Columns to read in as array
        col_keys = ['NUMBER', 'X_WORLD', 'Y_WORLD']
        col_list = []

        # Tuple with required columns
        for key in col_keys:
            col_list.append(self.SExcols[filter][key])

        col_tuple = tuple(col_list)

        # For each tile read as arrays
        for tile in self.tiles:

            detcat = os.path.join(self.path, tile,
                                  tile + self.det_filter + ".cat")
            (self.NUMBER[tile], self.X_WORLD[tile],
             self.Y_WORLD[tile]) = tableio.get_data(detcat, col_tuple)

        self.rdID = 1
        return

    def nearest_rd(self, x, y):

        if not self.rdID:
            self.get_rdID()

        # Make the distance array
        dmin = 1e20
        for tile in self.tiles:

            X = self.X_WORLD[tile]
            Y = self.Y_WORLD[tile]
            N = self.NUMBER[tile]
            D = np.sqrt((x - X)**2 + (y - Y)**2)
            idx = D.argmin()

            if D[idx] < dmin:
                ID = tile + "_" + str(int(N[idx]))

            dmin = min(dmin, D[idx])

        return str(ID)

    def near_rd(self, x, y, dmin, out_d=None):

        # dmin in arcmin
        if not self.rdID:
            self.get_rdID()

        IDsel = None
        dsel = None

        IDsel = []
        dsel = {}  # distance to seld
        # Make the distance array
        for tile in self.tiles:

            X = self.X_WORLD[tile]
            Y = self.Y_WORLD[tile]
            N = self.NUMBER[tile]
            #D  = np.sqrt( (x-X)**2 + (y-Y)**2)
            D = 60.0 * astrometry.circle_distance(x, y, X,
                                                  Y, units='deg')  # in arcmins
            idx, = np.where(D < dmin)

            for i in idx:
                ID = tile + "_" + str(int(N[i]))
                IDsel.append(ID)
                dsel[ID] = D[i]

        if out_d:
            return IDsel, dsel

        return IDsel

    def get_xyID(self, filter='i'):

        # Columns to read in as array
        col_keys = ['NUMBER', 'X_IMAGE', 'Y_IMAGE']
        col_list = []

        # Tuple with required columns
        for key in col_keys:
            col_list.append(self.SExcols[filter][key])

        col_tuple = tuple(col_list)

        (self.NUMBER, self.X_IMAGE,
         self.Y_IMAGE) = tableio.get_data(self.detcat, col_tuple)

        self.xyID = 1
        return

    def nearest_xy(self, x, y):

        if not self.xyID:
            self.get_xyID()

        # Make the distance array
        D = np.sqrt((x - self.X_IMAGE)**2 + (y - self.Y_IMAGE)**2)
        idx = D.argmin()
        ID = self.tile + "_" + str(int(self.NUMBER[idx]))
        return str(ID)

    def read_probs(self):

        # The reg expresion to compile
        regexp_point = re.compile(r"arange\("
                                  r"(?P<z1>[0-9]+.[0-9]+),"
                                  r"(?P<z2>[0-9]+.[0-9]+),"
                                  r"(?P<dz>[0-9]+.[0-9]+)\)")

        # Dictionary of probability arrays
        self.PROBS = {}
        for line in open(self.bpzprobs).readlines():

            point = regexp_point.search(line)

            # Extract the information if a point was selected
            if point:
                z1 = float(point.group('z1'))
                z2 = float(point.group('z2'))
                dz = float(point.group('dz'))
                self.zprobs = np.arange(z1, z2, dz)
                self.probs_header = '# ID  p_bayes(z)  where z=arange(%.4f,%.4f,%.4f)' % (
                    z1, z2, dz)

            fields = line.split()
            if fields[0][0] == "#":
                continue
            ID = self.tile + "_" + fields[0]
            self.PROBS[ID] = np.asarray(list(map(float, fields[1:])))

        return
