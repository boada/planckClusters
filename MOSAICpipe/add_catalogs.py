def match_SEx(tilename, filters, newfirm=False):
    from astropy.io import ascii
    from astropy import wcs
    from astropy.table import Column
    from astropy.coordinates import SkyCoord
    from astropy import units as u

    # get the wcs info from the i-band because it is the dectection image
    if os.path.isfile('{}i.fits'.format(tilename)):
        w = wcs.WCS('{}i.fits'.format(tilename))
    else:
        print('# No i-band detection image aborting...')
        return

    # read in the sdss catalog. We are doing this first because we only need to
    # do it once
    try:
        sdss_cat = ascii.read('/home/boada/Projects/'
                         'planckClusters/data/extern/SDSS/{}/'
                         '{}_SDSS_catalog.csv'.format(tilename, tilename))
        if len(sdss_cat) < 2:
            print('# SDSS TOO SHORT!')
            sdss = False
        else:
            sdss = True
    except FileNotFoundError:
        print('# SDSS CATALOG NOT FOUND!')
        sdss = False

    try:
        ps1_cat = ascii.read('/home/boada/Projects/'
                         'planckClusters/data/extern/PS1/{}/'
                         '{}_PS1_catalog.csv'.format(tilename, tilename))
        if len(ps1_cat) < 2:
            print('# PS1 TOO SHORT!')
        else:
            ps1 = True
    except FileNotFoundError:
        print('# PS1 CATALOG NOT FOUND!')
        ps1 = False

    try:
        twoMASS_cat = ascii.read('/home/boada/Projects/'
                         'planckClusters/data/extern/2MASS/{}/'
                         '{}_2MASS_catalog.csv'.format(tilename, tilename))
        if len(twoMASS_cat) < 2:
            print('# 2MASS TOO SHORT!')
            twoMASS = False
        else:
            twoMASS = True
    except FileNotFoundError:
        print('# 2MASS CATALOG NOT FOUND!')
        twoMASS = False

    # need these coordinates for the matching
    if sdss:
        s_coord = SkyCoord(ra=sdss_cat['ra'] * u.degree, dec=sdss_cat['dec'] *
                       u.degree)
    if ps1:
        p_coord = SkyCoord(ra=ps1_cat['ramean'] * u.degree,
                        dec=ps1_cat['decmean'] * u.degree)
    if twoMASS:
        tm_coord = SkyCoord(ra=twoMASS_cat['ra'] * u.degree,
                    dec=twoMASS_cat['dec'] * u.degree)

    for filter in filters:
        if not newfirm and filter == 'K':
            continue
        # only do the Kband when we get there
        if filter == 'K':
            if os.path.isfile('{}{}_cal.cat'.format(tilename, filter)):
                kband = True
            elif os.path.isfile('{}{}.cat'.format(tilename, filter)):
                kband = True
            else:
                kband = False
                continue
        else:
            kband = False

        try:
            cat = ascii.read('{}{}_cal.cat'.format(tilename, filter))
            cal = True
        except FileNotFoundError:
            cat = ascii.read('{}{}.cat'.format(tilename, filter))
            cal = False

        # calulate the RA/DEC of all of the detections
        ra, dec = w.all_pix2world(cat['X_IMAGE'], cat['Y_IMAGE'], 0)
        ra = Column(ra, name='RA')
        dec = Column(dec, name='DEC')
        try:
            cat.add_column(ra)
        except ValueError:
            pass
        try:
            cat.add_column(dec)
        except ValueError:
            pass

        # need these coordinates for the matching
        c_coord = SkyCoord(ra=cat['RA'] * u.degree, dec=cat['DEC'] * u.degree)

        if sdss and not kband:
            # match the two catalogs -- idxc for cat, idxs for sdss
            idxc, idxs, d2d, d3d = s_coord.search_around_sky(c_coord, 1 * u.arcsec)

            # make some data to catch
            d = np.ones(len(cat)) * 99.0 # 99 is the non-detection value in SEx...

            # add some extra info from SDSS -- Specz, Photoz
            cols = []
            cols.append(Column(d, name='sdss_{}'.format(filter)))
            cols.append(Column(d, name='sdss_{}_err'.format(filter)))
            cols.append(Column(d, name='sdss_petro_{}'.format(filter)))
            cols.append(Column(d, name='sdss_petro_{}_err'.format(filter)))
            cols.append(Column(d, name='sdss_objid', dtype=np.long))
            cols.append(Column(d, name='sdss_specz'))
            cols.append(Column(d, name='sdss_specz_err'))
            cols.append(Column(d, name='sdss_photoz'))
            cols.append(Column(d, name='sdss_photoz_err'))
            cols.append(Column(d, name='sdss_type'))
            for col in cols:
                try:
                    cat.add_column(col)
                except ValueError:
                    pass

            # merge the matches
            cat['sdss_objid'][idxc] = sdss_cat['objid'][idxs]
            cat['sdss_{}'.format(filter)][idxc] = sdss_cat[
                                            'fiberMag_{}'.format(filter)][idxs]
            cat['sdss_{}_err'.format(filter)][idxc] = sdss_cat[
                                            'fiberMagErr_{}'.format(filter)][idxs]
            cat['sdss_petro_{}'.format(filter)][idxc] = \
                                    sdss_cat['petroRad_{}'.format(filter)][idxs]
            cat['sdss_petro_{}_err'.format(filter)][idxc] = \
                                    sdss_cat['petroRadErr_{}'.format(filter)][idxs]
            cat['sdss_specz'][idxc] = sdss_cat['specz'][idxs]
            cat['sdss_specz_err'][idxc] = sdss_cat['specz_err'][idxs]
            cat['sdss_photoz'][idxc] = sdss_cat['photoz'][idxs]
            cat['sdss_photoz_err'][idxc] = sdss_cat['photoz_err'][idxs]
            cat['sdss_type'][idxc] = sdss_cat['type'][idxs]

        #### NOW THE PANSTARRS DATA ####
        if ps1 and not kband:
            idxc, idxp, d2d, d3d = p_coord.search_around_sky(c_coord, 1 *
                                                             u.arcsec)
            # make some data to catch
            d = np.ones(len(cat)) * 99.0 # 99 is the non-detection value in SEx...
            col = Column(d, name='ps1_{}'.format(filter))
            col_err = Column(d, name='ps1_{}_err'.format(filter))
            try:
                cat.add_column(col)
                cat.add_column(col_err)
            except ValueError:
                pass
            # merge the matches
            cat['ps1_{}'.format(filter)][idxc] = ps1_cat[
                                            '{}meanpsfmag'.format(filter)][idxp]
            cat['ps1_{}_err'.format(filter)][idxc] = \
                                ps1_cat['{}meanpsfmagerr'.format(filter)][idxp]

        #### NOW THE 2MASS DATA -- IF KBAND ####
        if kband and twoMASS:
            idxc, idxp, d2d, d3d = tm_coord.search_around_sky(c_coord,
                                                              1 * u.arcsec)
            # make some data to catch
            # 99 is the non-detection value in Sextractor
            d = np.ones(len(cat)) * 99.0
            col = Column(d, name='2MASS_{}'.format(filter))
            col_err = Column(d, name='2MASS_{}_err'.format(filter))
            try:
                cat.add_column(col)
                cat.add_column(col_err)
            except ValueError:
                pass
            # merge the matches
            cat['2MASS_{}'.format(filter)][idxc] = twoMASS_cat[
                                            '{}mag'.format(filter)][idxp]
            cat['2MASS_{}_err'.format(filter)][idxc] = twoMASS_cat[
                                            'e_{}mag'.format(filter)][idxp]

        # write out the results
        cat.write('tmp.color', format='ascii.commented_header', overwrite=True)

        # mv the old color catalog
        if cal:
            os.rename('{}{}_cal.cat'.format(tilename, filter),
                      '{}{}_cal.cat.orig'.format(tilename, filter))
        else:
            os.rename('{}{}.cat'.format(tilename, filter),
                      '{}{}.cat.orig'.format(tilename, filter))

        with open('tmp.color') as f:
            line = f.readline()
            line = line.split()
            if cal:
                file = '{}{}_cal.cat'.format(tilename, filter)
            else:
                file = '{}{}.cat'.format(tilename, filter)
            with open(file, 'wt') as f2:
                f2.write('## ' + time.ctime() + '\n')
                f2.write('## Matched Catalog file for Observation: '
                            '{}\n'.format(tilename))
                f2.write('## (This file was generated by the '
                        'BCS Rutgers pipeline)\n##\n')
                for i, l in enumerate(line[1:]):
                    f2.write('{} {} {}\n'.format(line[0], i + 1, l))
                line = f.readline()
                while line:
                    f2.write(line)
                    line = f.readline()

        os.remove('tmp.color')

def add_Speczs(tilename, dust=False, newfirm=False):
    ''' This function adds the spec-z's to the color catalog. The previous
    function adds all of the extra information to the sextractor catalogs. So
    you can run either one, or both. Because they don't do the same thing even
    though it looks like it.

    '''

    from astropy.io import ascii
    from astropy import wcs
    from astropy.table import Column
    from astropy.coordinates import SkyCoord
    from astropy import units as u

    w = wcs.WCS('{}i.fits'.format(tilename))
    cat = ascii.read('{}.color'.format(tilename))
    ra, dec = w.all_pix2world(cat['X_IMAGE'], cat['Y_IMAGE'], 0)
    ra = Column(ra, name='RA')
    dec = Column(dec, name='DEC')
    cat.add_column(ra)
    cat.add_column(dec)

    # match the catalogs
    try:
        sdss_cat = ascii.read('/home/boada/Projects/'
                         'planckClusters/data/extern/SDSS/{}/'
                         '{}_SDSS_catalog.csv'.format(tilename, tilename))
        sdss = True
        if len(sdss_cat) < 2:
            sdss = False
    except FileNotFoundError:
        print('SDSS CATALOG NOT FOUND!')
        sdss = False
    if sdss:
        c_coord = SkyCoord(ra=cat['RA'] * u.degree, dec=cat['DEC'] * u.degree)
        s_coord = SkyCoord(ra=sdss_cat['ra'] * u.degree, dec=sdss_cat['dec'] *
                       u.degree)
        idxc, idxs, d2d, d3d = s_coord.search_around_sky(c_coord, 1 * u.arcsec)

        # set the fill value
        sdss_cat['specz'].fill_value = -99.0

    # make a new column to catch the results
    zspec = np.ones(len(cat)) * -99.0
    if sdss:
        zspec[idxc] = sdss_cat['specz'].filled()[idxs]
    zspec = Column(zspec, name='Z_S')
    cat.add_column(zspec)
    cat.write('tmp.color', format='ascii.commented_header', overwrite=True)

    # mv the old color catalog
    os.rename('{}.color'.format(tilename), '{}.color.orig'.format(tilename))

    with open('tmp.color') as f:
        line = f.readline()
        line = line.split()
        with open('{}.color'.format(tilename), 'wt') as f2:
            f2.write('## ' + time.ctime() + '\n')
            f2.write('## BPZ Catalog file for Observation: '
                        '{}\n'.format(tilename))
            f2.write('## This file was generated by the '
                      'BCS Rutgers pipeline.\n')

            if dust:
                f2.write('## This file HAS been dust corrected.\n##\n')
            else:
                f2.write('## This file HAS NOT been dust corrected.\n##\n')

            for i, l in enumerate(line[1:]):
                f2.write('{} {} {}\n'.format(line[0], i + 1, l))
            line = f.readline()
            while line:
                f2.write(line)
                line = f.readline()

    # fix the columns file
    with open('{}.columns'.format(tilename), 'at') as f:
        f.write('Z_S {}'.format(i + 1))

    os.remove('tmp.color')

    return

def add_extra_photometry(tilename, filters=['u']):
    ''' This function adds the spec-z's to the color catalog. The previous
    function adds all of the extra information to the sextractor catalogs. So
    you can run either one, or both. Because they don't do the same thing even
    though it looks like it.

    '''

    from astropy.io import ascii
    from astropy.table import Column
    from astropy.coordinates import SkyCoord
    from astropy import units as u

    # read in the sdss catalog. We are doing this first because we only need to
    # do it once
    try:
        cat = ascii.read('/home/boada/Projects/'
                         'planckClusters/data/extern/SDSS/{}/'
                         '{}_SDSS_catalog.csv'.format(tilename, tilename))
        sdss = True
        ps1 = False
        if len(cat) < 2:
            print('# SDSS CATALOG TOO SHORT! -- TRYING PS1')
            sdss = False
    except FileNotFoundError:
        sdss = False
        print('# SDSS CATALOG NOT FOUND! -- TRYING PS1')
    if not sdss:
        try:
            cat = ascii.read('/home/boada/Projects/'
                             'planckClusters/data/extern/PS1/{}/'
                             '{}_PS1_catalog.csv'.format(tilename, tilename))
            ps1 = True
            if len(cat) < 2:
                print('# PS1 CATALOG TOO SHORT!')
                return
        except FileNotFoundError:
            print('# PS1 CATALOG NOT FOUND!')
            return
    # need these coordinates for the matching
    if sdss:
        s_coord = SkyCoord(ra=cat['ra'] * u.degree, dec=cat['dec'] * u.degree)
    elif ps1:
        s_coord = SkyCoord(ra=cat['ramean'] * u.degree, dec=cat['decmean'] *
                           u.degree)
    else:
        return

    color = ascii.read('{}.color'.format(tilename))
    c_coord = SkyCoord(ra=color['RA'] * u.degree, dec=color['DEC'] * u.degree)

    # match the catalogs
    idxc, idxs, d2d, d3d = s_coord.search_around_sky(c_coord, 1 * u.arcsec)

    for filter in filters:
        if sdss:
            # set the fill value
            cat[filter].fill_value = -99.0
            cat['{}_err'.format(filter)].fill_value = 0.

            # make a new column to catch the results
            newPhot = np.ones(len(color)) * -99.0
            newPhot_err = np.zeros(len(color))
            newPhot[idxc] = cat[filter].filled()[idxs]
            newPhot_err[idxc] = cat['{}_err'.format(filter)].filled()[idxs]

            newPhot = Column(newPhot, name='sdss_{}'.format(filter))
            newPhot_err = Column(newPhot_err, name='sdss_{}_err'.format(filter))
            color.add_column(newPhot)
            color.add_column(newPhot_err)
        elif ps1:
            cat[filter].fill_value = -99.0
            cat['{}_err'.format(filter)].fill_value = 0.

            # make a new column to catch the results
            newPhot = np.ones(len(color)) * -99.0
            newPhot_err = np.zeros(len(color))
            newPhot[idxc] = cat['{}meanpsfmag'.format(filter)].filled()[idxs]
            newPhot_err[idxc] = cat['{}meanpsfmagerr'.format(
                                                    filter)].filled()[idxs]

            newPhot = Column(newPhot, name='ps1_{}'.format(filter))
            newPhot_err = Column(newPhot, name='ps1_{}_err'.format(filter))
            color.add_column(newPhot)
            color.add_column(newPhot_err)
        else:
            return

    color.write('tmp.color', format='ascii.commented_header', overwrite=True)

    with open('tmp.color') as f:
        line = f.readline()
        line = line.split()
        with open('{}.color'.format(tilename), 'wt') as f2:
            f2.write('## ' + time.ctime() + '\n')
            f2.write('## BPZ .columns file for Observation: '
                        '{}\n'.format(tilename))
            f2.write('## (This file was generated by the '
                      'BCS Rutgers pipeline)\n##\n')
            for i, l in enumerate(line[1:]):
                f2.write('{} {} {}\n'.format(line[0], i + 1, l))
            line = f.readline()
            while line:
                f2.write(line)
                line = f.readline()

    # fix the columns file
    with open('{}.columns'.format(tilename)) as f:
        line = f.readline()
        with open('tmp.columns', 'w') as f2:
            while line:
                if 'M_0' in line:
                    f2.write('SLOAN-SDSS.{}\t {},{}'.format(filter, i, i + 1))
                    f2.write('\t AB 0.05\t 0.0\n')
                f2.write(line)
                line = f.readline()

    # mv the old color catalog
    os.rename('tmp.columns'.format(tilename), '{}.columns'.format(tilename))
    os.remove('tmp.color')

    return
