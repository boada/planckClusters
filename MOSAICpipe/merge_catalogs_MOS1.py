#!/usr/bin/env python

from bcs_catalogs import *
import os
import sys


def main(tile, path, filters):
    # check that the data for the filters is really there
    filters_real = []
    for f in filters:
        if os.path.isfile('{}{}_cal.cat'.format(tile, f)):
            filters_real.append(f)

    filters = filters_real

    # Make the header
    make_header(tile, path)

    # Go for each tile, read and print
    header = 1
    print("Doing tile: %s" % tile, file=sys.stderr)
    read_cats(tile, path, filters_real, det_filter='i', header=header)
    return

def read_cats(tilename, path, filters, det_filter='i', header=None):

    mcatfile = os.path.join(path, tilename, "%s_merged.cat" % tilename)
    probfile = os.path.join(path, tilename, "%s_probs.dat" % tilename)
    mcat = open(mcatfile, "a")  # append mode
    prob = open(probfile, "a")  # append mode

    # Read in the catalogs for that tile
    c = bcs_catalogs(tilename, path, verb=1, filters=filters)

    # The header for the prob file
    if header:
        prob.write("%s  \n" % c.probs_header)

    IDsel = []
    for ID in list(c.cat['i']['SHAPE'].keys()):

        IDname = "%s_%s" % (tilename, ID)

        goodID = 1
        cat = c.cat[det_filter]

        # mask out the stars
        if cat['CLASS_STAR'][ID] > 1:
            goodID = 0

            # If OK add them to the list
        if goodID:
            IDsel.append(ID)

            # Now write the whole thing tile by tile
            mcat.write("%-20s " % IDname)
            mcat.write("%15.10s " % cat['X_WORLD'][ID])
            mcat.write("%15.10s " % cat['Y_WORLD'][ID])

            for filter in filters:

                # If dust, use corrected values
                # There isn't actually an Xcorr -- so this is going to fail
                try:
                    mcat.write("%8.3f " % (
                        c.cat[filter]['MAG_AUTO'][ID] - c.Xcorr[filter][ID]))
                    mcat.write("%8.3f " % (c.cat[filter]['MAGERR_AUTO'][ID] +
                                           c.XcorrErr[filter][ID]))
                except AttributeError:
                    mcat.write("%8.3f " % c.cat[filter]['MAG_AUTO'][ID])
                    mcat.write("%8.3f " % c.cat[filter]['MAGERR_AUTO'][ID])
                #mcat.write("%8.3f " % c.cat[filter]['MAG_AUTO'][ID])
                #mcat.write("%8.3f " % c.cat[filter]['MAGERR_AUTO'][ID])

                # Write the SN, but avoid overflow
                if (c.cat[filter]['MAGERR_AUTO'][ID] > 100):
                    SN = 0.0
                else:
                    SN = 1. / (10 **
                               (0.4 * c.cat[filter]['MAGERR_AUTO'][ID]) - 1)
                mcat.write("%8.3f " % SN)

            for filter in filters:
                mcat.write("%8.3f " % c.MAG_BPZ[filter][ID])
                mcat.write("%8.3f " % c.MAG_BPZERR[filter][ID])

            mcat.write("%5.2f " % c.Z_B[ID])
            mcat.write("%5.2f " % c.Z_B_MIN[ID])
            mcat.write("%5.2f " % c.Z_B_MAX[ID])
            mcat.write("%5.2f " % c.T_B[ID])
            mcat.write("%6.4f " % c.ODDS[ID])
            mcat.write("%5.2f " % c.Z_ML[ID])
            mcat.write("%5.2f " % c.T_ML[ID])
            mcat.write("%8.2f " % c.CHI2[ID])
            mcat.write("%5.2f " % cat['CLASS_STAR'][ID])
            mcat.write("%8.2f " %
                       (cat['A_IMAGE'][ID] * cat['KRON_RADIUS'][ID]))
            mcat.write("%8.2f " %
                       (cat['B_IMAGE'][ID] * cat['KRON_RADIUS'][ID]))
            mcat.write("%8.2f " % cat['THETA_IMAGE'][ID])
            mcat.write("%8.2f " % cat['X_IMAGE'][ID])
            mcat.write("%8.2f " % cat['Y_IMAGE'][ID])
            mcat.write("\n")

            # and the probs file
            format = len(c.PROBS[ID]) * '%.3e '
            prob.write("%-20s " % IDname)
            prob.write(format % tuple(c.PROBS[ID]))
            prob.write("\n")

    return

def make_header(tilename, path):

    s = sys.stdout
    mcatfile = os.path.join(path, tilename, "%s_merged.cat" % tilename)
    probfile = os.path.join(path, tilename, "%s_probs.dat" % tilename)
    mcat = open(mcatfile, "w")  # create a new one
    prob = open(probfile, "w")  # create a new one

    mcat.write("#   0  IDNAME           \n")
    mcat.write("#   1  RA (degrees)     \n")
    mcat.write("#   2  DEC(degress)     \n")

    i = 3
    for filter in filters:
        mcat.write("# %3d  %s MAG_AUTO     \n" % (i, filter))
        mcat.write("# %3d  %s MAGERR_AUTO  \n" % (i + 1, filter))
        mcat.write("# %3d  %s S/N  \n" % (i + 2, filter))
        i += 3

    for filter in filters:
        mcat.write("# %3d  %s MAG_BPZ       \n" % (i, filter))
        mcat.write("# %3d  %s MAG_BPZERR    \n" % (i + 1, filter))
        i += 2

    extra_keys = ['Z_B      Bayesian photo-z', 'Z_B_MIN  Min photo-z',
                  'Z_B_MAX  Max photo-z', 'T_B      Bayesian Spectral Type',
                  'ODDS     Odds', 'Z_ML     Likelihood photo-z',
                  'T_ML     Likelihood Spectral Type', 'CHI-SQUARED',
                  'CLASS_STAR (1=star, 0=galaxy)',
                  'A_IMAGE x KRON_RADIUS (pixels)',
                  'B_IMAGE x KRON_RADIUS (pixels)', 'THETA_IMAGE (degrees)',
                  'X_IMAGE (pix)', 'Y_IMAGE (pix)']

    for key in extra_keys:
        mcat.write("# %3d  %s   \n" % (i, key))
        i += 1

    return


if __name__ == "__main__":
    tile = sys.argv[1]
    path = sys.argv[2]

    filters = ('g', 'r', 'i', 'z')
    main(tile, path, filters)
