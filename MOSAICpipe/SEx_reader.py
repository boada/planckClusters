#!/usr/bin/env python

from __future__ import print_function
from builtins import str
from builtins import object
import sys


class SEx_reader(object):
    """Reads in the SExtractor catalogs and stores the information in a dictionary"""

    def __init__(self, catalog, preID=None, verb=None):

        self.verb = verb
        self.catalog = catalog

        if preID:
            self.preID = preID + '_'
        else:
            self.preID = ''

        self._read_catalog()

        return

    def _SEx_head(self):

        if self.verb:
            print("\r Parsing SEx head for:", self.catalog, file=sys.stderr)

        self.SExcols = {}

        # Read the detection catalog
        self.catstr = open(self.catalog).readlines()

        for line in self.catstr:

            if line[0] != '#':
                break

            if line[:2] == "##":
                continue

            try:
                line = line.strip()
                vals = line.split()
                col = vals[1]
                key = vals[2]
                self.SExcols[key] = int(col) - 1
                if self.verb:
                    print("# %-20s %s" % (key, self.SExcols[key] + 1),
                          file=sys.stderr)
            except:
                continue

        self.ncols = len(list(self.SExcols.values()))
        self.colnames = list(self.SExcols.keys())

        #print >>sys.stderr,"# Read %s columns" % self.ncols

        return

    def _read_catalog(self):

        if self.verb:
            print("\r Reading:", self.catalog, file=sys.stderr)

        self._SEx_head()

        self.cat = {}
        self.IDs = []

        # Make the dictionary
        self.cat['SHAPE'] = {}
        for key in list(self.SExcols.keys()):
            self.cat[key] = {}

        for line in self.catstr:

            vals = line.split()

            if vals[0][0] == "#":
                continue

            ID = self.preID + str(vals[0])
            self.IDs.append(ID)
            for key in list(self.SExcols.keys()):

                if key == 'NUMBER' or key == 'FLAGS':
                    self.cat[key][ID] = int(vals[self.SExcols[key]])
                else:
                    self.cat[key][ID] = float(vals[self.SExcols[key]])

                # Replace MAGERR == 0 to by something more sensible
                if key[0:6] == 'MAGERR' and self.cat[key][ID] == 0:
                    self.cat[key][ID] = 1.0e-5

            # Add the SHAPE param
            self.cat['SHAPE'][ID] = (
                self.cat['X_IMAGE'][ID], self.cat['Y_IMAGE'][ID],
                self.cat['A_IMAGE'][ID] * self.cat['KRON_RADIUS'][ID],
                self.cat['B_IMAGE'][ID] * self.cat['KRON_RADIUS'][ID],
                self.cat['THETA_IMAGE'][ID])

        self.nrows = len(list(self.cat['NUMBER'].values()))

        return


class SEx_reader_multi(object):
    """Reads in the SExtractor catalogs and stores the information in a dictionary"""

    def __init__(self, verb=None):

        self.verb = verb
        self.cat = {}
        self.IDs = []
        self.readcat = None
        return

    # Function to read the catalog
    def read_catalog(self, catalog, preID=None):

        self.catalog = catalog
        self.preID = preID

        if self.verb:
            print("\r Reading:", self.catalog, file=sys.stderr)

        self._SEx_head()

        # Make the dictionaries only once
        if not self.readcat:
            self.cat['SHAPE'] = {}
            for key in list(self.SExcols.keys()):
                self.cat[key] = {}
            self.readcat = 'yes'

        # Go thru the whole as string
        for line in self.catstr:

            vals = line.split()

            if vals[0][0] == "#":
                continue

            ID = self.preID + "_" + str(vals[0])
            self.IDs.append(ID)
            for key in list(self.SExcols.keys()):

                if key == 'NUMBER' or key == 'FLAGS':
                    self.cat[key][ID] = int(vals[self.SExcols[key]])
                else:
                    self.cat[key][ID] = float(vals[self.SExcols[key]])

                # Replace MAGERR == 0 to by something more sensible
                if key[0:6] == 'MAGERR' and self.cat[key][ID] == 0:
                    self.cat[key][ID] = 1.0e-5

            # Add the SHAPE param
            self.cat['SHAPE'][ID] = (
                self.cat['X_IMAGE'][ID], self.cat['Y_IMAGE'][ID],
                self.cat['A_IMAGE'][ID] * self.cat['KRON_RADIUS'][ID],
                self.cat['B_IMAGE'][ID] * self.cat['KRON_RADIUS'][ID],
                self.cat['THETA_IMAGE'][ID])

        #self.nrows = len( self.cat['NUMBER'].values())
        return

    def _SEx_head(self):

        if self.verb:
            print("\r Parsing SEx head for:", self.catalog, file=sys.stderr)

        self.SExcols = {}

        # Read the detection catalog
        self.catstr = open(self.catalog).readlines()

        for line in self.catstr:

            if line[0] != '#':
                break

            if line[:2] == "##":
                continue

            try:
                line = line.strip()
                vals = line.split()
                col = vals[1]
                key = vals[2]
                self.SExcols[key] = int(col) - 1
                #if self.verb:
                #    print >>sys.stderr, "# %-20s %s" % (key,self.SExcols[key]+1)
            except:
                continue

        self.ncols = len(list(self.SExcols.values()))
        self.colnames = list(self.SExcols.keys())

        #print >>sys.stderr,"# Read %s columns" % self.ncols

        return

#catalog = "/home/felipe/BCS/PROC/BCS2322-5417/BCS2322-5417i.cat"
#c = SEx_reader(catalog,verb=1)

#print c.IDs
#print c.colnames
#print c.nrows
#print c.ncols
