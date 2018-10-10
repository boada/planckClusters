# Find the eBV dust correction for each source in the catalogs
def DustCorrection(self):
    ''' This figures out the dust extinction and corrects the sextractor
    photometry that has been cleaned by the BuildColorCat function. It also
    puts the dust corrections into a series of dictions that are used by
    BuildColorCat. So if we don't run this function it doesn't include the
    dust correction. This is even true after it writes a dust file. I think
    the dust file is really just there for us to inspect for funny stuff.

    '''

    print()
    self.DustCat = self.tilename + ".dust"

    # Get RA,DEC from the detection catalog
    detCatalog = self.combcat['i']
    detcols = SEx_head(detCatalog, verb=None)
    cols = (detcols['NUMBER'], detcols['X_WORLD'], detcols['Y_WORLD'])
    (id, ra, dec) = tableio.get_data(detCatalog, cols)
    outColumns = ['ID', ]

    # Get e(B-V) for every source in the detection catalog
    print("# Computing e(B-V) for all %s ra,dec" % len(ra), file=sys.stderr)
    self.eBV = deredden.get_EBV(ra, dec)
    print("# Done...", file=sys.stderr)

    # Prepare the header for the output file
    header = '## {}\n'.format(time.ctime()) + \
             '## Dust correction extinction ' +\
                'for each object/filter in: {}\n'.format(self.tilename) +\
             '## This file was generated automatically by the BCS ' +\
                'Rutgers pipeline\n' +\
             '## These must be subtracted from the SExtractor ' +\
                'magnitudes \n' +\
             '## Dust Correction e(B-V), mean, min, max: ' +\
             '{0:.4f}, {0:.4f}, {0:.4f}\n'.format(self.eBV.mean(),
                                            self.eBV.min(), self.eBV.max())
    VarsOut = [id]

    # Get the dust extinction correction for each filter
    for filter in self.filters:
        self.XCorr[filter] = deredden.filterFactor(filter) * self.eBV
        self.XCorrError[filter] = self.XCorr[filter] * 0.16
        # Some more work on the header
        header += "## Dust Correction %s, mean, min, max:  %.4f %.4f, %.4f mags\n" % (
            filter, self.XCorr[filter].mean(), self.XCorr[filter].min(),
            self.XCorr[filter].max())
        outColumns.append(filter + '_MOSAICII Dust Correction')
        outColumns.append(filter + '_MOSAICII Dust Correction Error')
        VarsOut.append(self.XCorr[filter])
        VarsOut.append(self.XCorrError[filter])

    #print outColumns
    i = 0
    header += '# ' + str(i + 1) + '\t' + outColumns[i] + '\n'
    for filter in self.filters:
        header += '# {}\t{}\n'.format(str(i + 2), outColumns[i + 1])
        header += '# {}\t{}\n'.format(str(i + 3), outColumns[i + 2])
        i += 2

    vars = tuple(VarsOut)
    format = '%8i' + '%10.5f  ' * (len(vars) - 1)
    print('# Writing Dust Extinction Catalog...', file=sys.stderr)
    tableio.put_data(self.DustCat,
                     vars,
                     header=header,
                     format=format,
                     append='no')
    print('# Dust file complete.', file=sys.stderr)

    return
