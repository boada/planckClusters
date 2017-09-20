import tableio
from glob import glob
import os
import sys
import numpy as np

def SEx_head(catalog, verb='yes'):

    if verb:
        print("\r Parsing SEx head for:", catalog, file=sys.stderr)

    # Dictionary with column numbers
    SExcols = {}

    # Read the SExtractor catalog
    for line in open(catalog).readlines():

        if line[0] != '#':
            break

        if line[:2] == "##":
            continue

        try:
            line = line.strip()
            vals = line.split()
            col = vals[1]
            key = vals[2]
            SExcols[key] = int(col) - 1
            if verb:
                print("# %-20s %s" % (key, SExcols[key] + 1), file=sys.stderr)
        except:
            continue

    return SExcols

# Taken from APSIS fUtil
def deNAN(a, value=0.0):
    nans = np.logical_not(np.less(a, 0.0) + np.greater_equal(
        a, 0.0))
    return np.where(nans, value, a)

# Build th color catalog to use when computing the photo-z
# Adapted from JHU APSIS pipeline
def BuildColorCat(tilename, combcat, filters=['g', 'r', 'i', 'z', 'K'],
                  newfirm=True):

    # The default output names
    colorCat = tilename + "_complete.catalog"

    print('Processing catalogs... for: ', tilename, file=sys.stderr)

    flux = {}
    fluxerr = {}

    m = {}
    em = {}

    # Get the detection catalog required columns
    outColumns = ['NUMBER', 'X_IMAGE', 'Y_IMAGE']
    detCatalog = combcat['i']
    detcols = SEx_head(detCatalog, verb=None)
    detectionList = []
    for key in outColumns:
        detectionList.append(detcols[key])
    # the get_data function requires a tuple
    detectionColumns = tuple(detectionList)
    detection_variables = tableio.get_data(detCatalog,
                                           detectionColumns)

    # Read in the MAG_ISO and MAG_ISOERR from each catalog
    for filter in filters:
        if not newfirm and filter == 'K':
            continue
        # get the zeropoint Info
        tmp = np.genfromtxt('photometry_control_star_{}.dat'.format(
                            filter), names=True, dtype=None)
        zpoint = tmp['ZP']

        # Get the columns
        sexcols = SEx_head(combcat[filter], verb=None)

        ## Info for flux columns
        fluxList = []
        fluxList.append(sexcols['FLUX_ISO'])
        fluxList.append(sexcols['FLUXERR_ISO'])
        fluxColumns = tuple(
            fluxList)  # the get_data function interface requires a tuple

        # Get the array using tableio
        flux[filter], fluxerr[filter] = tableio.get_data(combcat[filter],
                                                                fluxColumns)
        m[filter] = flux[filter] * 0.0
        em[filter] = flux[filter] * 0.0

        # Fix the NAN values
        flux[filter] = deNAN(flux[filter])

        # Those objects with flux equal or less than 0 are assigned a
        # magnitude of 99 and a limiting magnitude equal to their
        # SExtractor photometric error. This is interpreted by BPZ as a
        # nondetection with zero flux and 1-sigma error equal to the
        # limiting magnitude

        #nondetected = np.less_equal(flux[filter], 0.0) * \
        #              np.greater(fluxerr[filter], 0.0)

        # update: There are a lot of really small positive values. I am
        # going to modify this to look for values really close to zero.

        nondetected = (flux[filter] < 1E-3) & (fluxerr[filter] > 0.0)

        # Those objects with error flux and flux equal to 0 are assigned a
        # magnitude of -99
        # and a flux of 0, which is interpreted by SExtractor as a
        # non-observed object

        nonobserved = np.less_equal(fluxerr[filter], 0.0)

        # When flux error > 100*(flux), mark as nonobserved (Benitez,
        # 24-Oct-03).

        nonobserved = np.where(fluxerr[filter] > 100 *
                                    (abs(flux[filter])), True,
                               nonobserved)

        detected = np.logical_not(nonobserved + nondetected)

        print(filter, zpoint)

        flux[filter] = np.clip(flux[filter], 1e-100, 1e100)
        m[filter] = np.where(detected,
                                  -2.5 * np.log10(abs(flux[filter])) +
                                  zpoint, m[filter])
        m[filter] = np.where(nondetected, 99.0, m[filter])
        m[filter] = np.where(nonobserved, -99.0, m[filter])

        # clip values from being too small or large, i.e. 0 or inf.
        fluxerr[filter] = np.clip(fluxerr[filter], 1e-100, 1e100)
        em[filter] = np.where(
            detected,
            2.5 * np.log10(1.0 + abs(fluxerr[filter] / flux[filter])), em[filter])
        em[filter] = np.where(
            nondetected,
            2.5 * np.log10(abs(fluxerr[filter])) - zpoint, em[filter])
        em[filter] = np.where(nonobserved, 0.0, em[filter])

        if filter == 'K':
            outColumns.append(filter + '_KittPeak_MAG_ISO')
            outColumns.append(filter + '_KittPeak_MAGERR_ISO')
        else:
            outColumns.append(filter + '_MOSAICII_MAG_ISO')
            outColumns.append(filter + '_MOSAICII_MAGERR_ISO')

    # Prepare the header
    header = \
           '## ' + '\n' + \
           '## BPZ Catalog file for Observation: ' + tilename + \
            '\n' + \
            '## (This file was generated automatically by' + \
            'the BCS Rutgers pipeline)\n##\n'
    for i in range(len(outColumns)):
        header = header + '# ' + str(i + 1) + '\t' + outColumns[i] + '\n'

        # Prepare the data
    vars = list(detection_variables)
    for filter in filters:
        if not newfirm and filter == 'K':
            continue
        vars.append(m[filter])
        vars.append(em[filter])

    variables = tuple(vars)
    format = '%i\t %10.2f %10.2f' + '%10.4f  ' * (len(variables) - 3)
    print('Writing data to multicolor catalog...', file=sys.stderr)
    tableio.put_data(colorCat,
                     variables,
                     header=header,
                     format=format,
                     append='no')
    print('Multicolor catalog complete.', file=sys.stderr)

    return


if __name__ == '__main__':
    dirs = [dirs for _, dirs, _ in os.walk('./')][0] # only want top level
    cwd = os.getcwd()
    for d in dirs:
        if 'PSZ' not in d:
            continue
        os.chdir(d)
        catalogs = glob('*_cal.cat')
        if len(catalogs) == 0:
            os.chdir(cwd)
            continue
        else:
            print(d, catalogs)

        # make the combined catalog Dictionary and filter information
        combcat = {}
        filters = []
        for c in catalogs:
            filter = c.split('_cal.cat')[0][-1]
            combcat[filter] = c
            filters.append(filter)
        tilename = c.split('_cal.cat')[0][:-1]
        try:
            BuildColorCat(tilename, combcat, filters)
        except KeyError:
            pass
        os.chdir(cwd)
