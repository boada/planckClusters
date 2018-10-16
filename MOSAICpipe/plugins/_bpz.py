import os
import sys
import time
import numpy as np
from astropy.io.fits import getheader
import subprocess
import shlex

# get the utils from the parent directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import SEx_head, deNAN
from pipe_utils import tableio
from add_catalogs import match_SEx, add_Speczs

# Build th color catalog to use when computing the photo-z
# Adapted from JHU APSIS pipeline
def BuildColorCat(self, newfirm=False):
    ''' Builds a photometry catalog that is going to be given to BPZ for the
    photometric redshift calculation. This function also cleans the original
    sextractor catalogs of the objects that have really low detections or
    are otherwise questionable detections. I don't think the cleaning is
    perfect, but is a lot better than it was before.

    '''

    print()
    # The default output names
    self.colorCat = self.tilename + ".color"
    self.columnsFile = self.tilename + ".columns"

    print('# Processing catalogs... for: ', self.tilename, file=sys.stderr)

    flux = {}
    fluxerr = {}
    mag = {}
    magerr = {}

    m = {}
    em = {}

    # Get the detection catalog required columns
    outColumns = ['NUMBER', 'X_IMAGE', 'Y_IMAGE']
    detCatalog = self.combcat['i']
    detcols = SEx_head(detCatalog, verb=None)
    detectionList = []
    for key in outColumns:
        detectionList.append(detcols[key])
    # the get_data function requires a tuple
    detectionColumns = tuple(detectionList)
    detection_variables = tableio.get_data(detCatalog,
                                           detectionColumns)

    # Read in the MAG_ISO and MAG_ISOERR from each catalog
    for filter in self.filters:
        if not newfirm and filter == 'K':
            continue
        # get the zeropoint Info

        input = self.combima[filter] + ".fits"
        hdr = getheader(input)
        zpoint = hdr['MAGZERO']
        zp_error = hdr['MAGZSIG']

        # Get the columns
        sexcols = SEx_head(self.combcat[filter], verb=None)

        ## Info for flux columns
        fluxList = []
        fluxList.append(sexcols['FLUX_ISO'])
        fluxList.append(sexcols['FLUXERR_ISO'])
        fluxList.append(sexcols['MAG_AUTO'])
        fluxList.append(sexcols['MAGERR_AUTO'])
        fluxColumns = tuple(
            fluxList)  # the get_data function interface requires a tuple

        # Get the array using tableio
        (flux[filter],
         fluxerr[filter],
         mag[filter],
         magerr[filter]) = tableio.get_data(self.combcat[filter],
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

        #nonobserved = np.where(fluxerr[filter] > 100 *
        big_errors = np.where(fluxerr[filter] > 100 *
                                    (abs(flux[filter])), True,
                               nonobserved)

        detected = np.logical_not(nonobserved + nondetected)

        flux[filter] = np.clip(flux[filter], 1e-100, 1e100)
        m[filter] = np.where(detected,
                                  -2.5 * np.log10(abs(flux[filter])) +
                                  zpoint - self.XCorr[filter],
                             m[filter])
        m[filter] = np.where(nondetected, 99.0, m[filter])
        m[filter] = np.where(nonobserved, -99.0, m[filter])

        # use the mag_auto where we have big errors
        m[filter] = np.where(big_errors, mag[filter] - self.XCorr[filter],
                             m[filter])

        # clip values from being too small or large, i.e. 0 or inf.
        fluxerr[filter] = np.clip(fluxerr[filter], 1e-100, 1e100)

        # detected
        em[filter] = np.where(
            detected,
            2.5 * np.log10(1.0 + abs(fluxerr[filter] / flux[filter])) +
            self.XCorrError[filter], em[filter])

        # non-detected
        em[filter] = np.where(
            nondetected,
            2.5 * np.log10(abs(fluxerr[filter])) - zpoint, em[filter])
        em[filter] = np.where(nonobserved, 0.0, em[filter])

        # use the magerr_auto where we have big ERRORS
        em[filter] = np.where(big_errors, magerr[filter], em[filter])

        #outColumns.append(filter +'_SDSS_MAG_ISO')
        #outColumns.append(filter +'_SDSS_MAGERR_ISO')
        if filter == 'K':
            outColumns.append(filter + '_KittPeak_MAG_ISO')
            outColumns.append(filter + '_KittPeak_MAGERR_ISO')
        else:
            outColumns.append(filter + '_MOSAICII_MAG_ISO')
            outColumns.append(filter + '_MOSAICII_MAGERR_ISO')

    # Prepare the header
    header = \
           '## ' + time.ctime() + '\n' + \
           '## BPZ Catalog file for Observation: ' + self.tilename + \
            '\n' + \
            '## This file was generated automatically by' + \
            'the BCS Rutgers pipeline.\n'

    try:
        if not self.XCorr[filter].all() == 0:
            header += '## This file HAS been dust corrected.\n##\n'
        else:
            header += '## This file HAS NOT been dust corrected.\n##\n'
    except AttributeError:
        header += '## This file HAS NOT been dust corrected.\n##\n'

    for i in range(len(outColumns)):
        header += '# ' + str(i + 1) + '\t' + outColumns[i] + '\n'

        # Prepare the data
    vars = list(detection_variables)
    for filter in self.filters:
        if not newfirm and filter == 'K':
            continue
        vars.append(m[filter])
        vars.append(em[filter])

    variables = tuple(vars)
    format = '%i\t %10.2f %10.2f' + '%10.4f  ' * (len(variables) - 3)
    print('# Writing data to multicolor catalog...', file=sys.stderr)
    tableio.put_data(self.colorCat,
                     variables,
                     header=header,
                     format=format,
                     append='no')
    print('# Multicolor catalog complete.', file=sys.stderr)

    # And now write .columns file
    with open(self.columnsFile, 'w') as cfile:
        cfile.write('## ' + time.ctime() + '\n')
        cfile.write('## BPZ .columns file for Observation: {}\n'.format(
            self.tilename))
        cfile.write('## (This file was generated by the '
                'BCS Rutgers pipeline)\n##\n')

        i = len(detection_variables)
        for filter in self.filters:
            if not newfirm and filter == 'K':
                continue
            # Get the zeropoint information
            input = self.combima[filter] + ".fits"
            hdr = getheader(input)
            zpoint = hdr['MAGZERO']
            zp_error = hdr['MAGZSIG']
            cal_catalog = hdr['MAGZCAT']

            # this says whether or not we should use the photometric catalog
            # calibrating filter or whether we should use the filters from
            # kpno to calculate the redshifts.
            # currently we are using the photometric catalog filter.
            kpno = False
            photocat = True

            if kpno:
                if filter == 'i':
                    n_mo = str(i + 1)

                if filter == 'K':
                    cfile.write('%s_KittPeak\t %s,%s\t AB\t %.2f\t 0.0\n' %
                            (filter, i + 1, i + 2, zp_error))
                else:
                    cfile.write('%s_MOSAICII\t %s,%s\t AB\t %.2f\t 0.0\n' %
                            (filter, i + 1, i + 2, zp_error))

            if photocat:
                if 'sdss' in cal_catalog.lower():
                    cal_system = 'SLOAN-SDSS'
                elif 'panstarrs' in cal_catalog.lower():
                    cal_system = 'PAN-STARRS-PS1'
                elif '2mass' in cal_catalog.lower():
                    cal_system = '2MASS-2MASS'
                else:
                    print('Catalog system not understood!')

                if filter == 'i':
                    n_mo = str(i + 1)

                if filter == 'K':
                    cfile.write('%s.%ss\t %s,%s\t AB\t %.2f\t 0.0\n' %
                        (cal_system, filter, i + 1, i + 2, zp_error))
                else:
                    cfile.write('%s.%s\t %s,%s\t AB\t %.2f\t 0.0\n' %
                        (cal_system, filter, i + 1, i + 2, zp_error))

            i += 2

        cfile.write('M_0\t%s\n' % n_mo)
    return

def BuildColorCat_auto(self, newfirm=False):
    ''' Builds a photometry catalog that is going to be given to BPZ for the
    photometric redshift calculation. This function also cleans the original
    sextractor catalogs of the objects that have really low detections or
    are otherwise questionable detections. I don't think the cleaning is
    perfect, but is a lot better than it was before.

    '''

    print()
    # The default output names
    self.colorCat = self.tilename + ".color"
    self.columnsFile = self.tilename + ".columns"

    print('# Processing catalogs... for: ', self.tilename, file=sys.stderr)

    mag = {}
    magerr = {}

    m = {}
    em = {}

    # Get the detection catalog required columns
    outColumns = ['NUMBER', 'X_IMAGE', 'Y_IMAGE']
    detCatalog = self.combcat['i']
    detcols = SEx_head(detCatalog, verb=None)
    detectionList = []
    for key in outColumns:
        detectionList.append(detcols[key])
    # the get_data function requires a tuple
    detectionColumns = tuple(detectionList)
    detection_variables = tableio.get_data(detCatalog,
                                           detectionColumns)

    # Read in the MAG_ISO and MAG_ISOERR from each catalog
    for filter in self.filters:
        if not newfirm and filter == 'K':
            continue
        # get the zeropoint Info
        tmp = np.genfromtxt('photometry_control_star_{}.dat'.format(
                            filter), names=True, dtype=None)
        zpoint = tmp['ZP']
        zp_error = tmp['ZP_sig']

        # Get the columns
        sexcols = SEx_head(self.combcat[filter], verb=None)

        ## Info for flux columns
        magList = []
        magList.append(sexcols['MAG_AUTO'])
        magList.append(sexcols['MAGERR_AUTO'])
        magColumns = tuple(
            magList)  # the get_data function interface requires a tuple

        # Get the array using tableio
        mag[filter], magerr[filter] = tableio.get_data(
            self.combcat[filter], magColumns)
        m[filter] = mag[filter] - self.XCorr[filter]
        em[filter] = magerr[filter] + self.XCorrError[filter]

        print(filter, zpoint)

        #outColumns.append(filter +'_SDSS_MAG_ISO')
        #outColumns.append(filter +'_SDSS_MAGERR_ISO')
        if filter == 'K':
            outColumns.append(filter + '_KittPeak_MAG_ISO')
            outColumns.append(filter + '_KittPeak_MAGERR_ISO')
        else:
            outColumns.append(filter + '_MOSAICII_MAG_ISO')
            outColumns.append(filter + '_MOSAICII_MAGERR_ISO')

    # Prepare the header
    header = \
           '## ' + time.ctime() + '\n' + \
           '## BPZ Catalog file for Observation: ' + self.tilename + \
            '\n' + \
            '## This file was generated automatically by' + \
            'the BCS Rutgers pipeline.\n'

    try:
        if not self.XCorr[filter].all() == 0:
            header += '## This file HAS been dust corrected.\n##\n'
        else:
            header += '## This file HAS NOT been dust corrected.\n##\n'
    except AttributeError:
        header += '## This file HAS NOT been dust corrected.\n##\n'

    for i in range(len(outColumns)):
        header += '# ' + str(i + 1) + '\t' + outColumns[i] + '\n'

        # Prepare the data
    vars = list(detection_variables)
    for filter in self.filters:
        if not newfirm and filter == 'K':
            continue
        vars.append(m[filter])
        vars.append(em[filter])

    variables = tuple(vars)
    format = '%i\t %10.2f %10.2f' + '%10.4f  ' * (len(variables) - 3)
    print('# Writing data to multicolor catalog...', file=sys.stderr)
    tableio.put_data(self.colorCat,
                     variables,
                     header=header,
                     format=format,
                     append='no')
    print('# Multicolor catalog complete.', file=sys.stderr)

    # And now write .columns file
    with open(self.columnsFile, 'w') as cfile:
        cfile.write('## ' + time.ctime() + '\n')
        cfile.write('## BPZ .columns file for Observation: {}\n'.format(
            self.tilename))
        cfile.write('## (This file was generated by the '
                'BCS Rutgers pipeline)\n##\n')

        i = len(detection_variables)
        for filter in self.filters:
            if not newfirm and filter == 'K':
                continue
            # Get the zeropoint information
            tmp = ascii.read('photometry_control_star_{}.dat'.format(
                                filter))
            zpoint = tmp['ZP'].data[0]
            zp_error = tmp['ZP_sig'].data[0]
            cal_catalog = tmp['[6]'].data[0]

            # this says whether or not we should use the photometric catalog
            # calibrating filter or whether we should use the filters from
            # kpno to calculate the redshifts.
            # currently we are using the photometric catalog filter.
            kpno = False
            photocat = True

            if kpno:
                if filter == 'i':
                    n_mo = str(i + 1)

                if filter == 'K':
                    cfile.write('%s_KittPeak\t %s,%s\t AB\t %.2f\t 0.0\n' %
                            (filter, i + 1, i + 2, zp_error))
                else:
                    cfile.write('%s_MOSAICII\t %s,%s\t AB\t %.2f\t 0.0\n' %
                            (filter, i + 1, i + 2, zp_error))

            if photocat:
                if 'sdss' in cal_catalog.lower():
                    cal_system = 'SLOAN-SDSS'
                elif 'panstarrs' in cal_catalog.lower():
                    cal_system = 'PAN-STARRS-PS1'
                elif '2mass' in cal_catalog.lower():
                    cal_system = '2MASS-2MASS'
                else:
                    print('Catalog system not understood!')

                if filter == 'i':
                    n_mo = str(i + 1)

                if filter == 'K':
                    cfile.write('%s.%ss\t %s,%s\t AB\t %.2f\t 0.0\n' %
                        (cal_system, filter, i + 1, i + 2, zp_error))
                else:
                    cfile.write('%s.%s\t %s,%s\t AB\t %.2f\t 0.0\n' %
                        (cal_system, filter, i + 1, i + 2, zp_error))

            i += 2

        cfile.write('M_0\t%s\n' % n_mo)
    return


# Run Benitez BPZ
def runBPZ(self, Specz=True, newfirm=False):
    """Runs BPZ on the multicolor catalog file using the .columns """

    print()
    # first we update with Specz's if we want to
    print('# Match Catalogs -- Add spec-zs')
    match_SEx(self.tilename, self.filters, newfirm=False)
    try:
        if not self.XCorr['i'].all() == 0:
            add_Speczs(self.tilename, dust=True, newfirm=False)
        else:
            add_Speczs(self.tilename, dust=False, newfirm=False)
    except AttributeError:
        add_Speczs(self.tilename, dust=False, newfirm=False)
    #add_extra_photometry(self.tilename)

    print('# Starting photometric redshift determination...',
          file=sys.stderr)
    bpz = os.path.join(os.environ['BPZPATH'], 'bpz.py ')
    #bpzcat = self.tilename + ".bpz"
    bpzprobs = self.tilename + ".probs"

    cmd = 'python ' + bpz + self.colorCat + \
        ' -ZMAX 1.8 -VERBOSE no -INTERP 2 -DZ 0.01 -SPECTRA ' + \
        'CWWSB_Benitez2003.list -PRIOR full -PROBS_LITE ' + bpzprobs

    if not self.dryrun:
        print(cmd)
        print("Running full prior", file=sys.stderr)
        p = subprocess.Popen(shlex.split(cmd))
        #p.wait(timeout=600) # this prevents really long running. for testing
        p.wait()
        print("Photo-z ready", file=sys.stderr)
    else:
        print(cmd)
    return
