def replace_vals_image(infits, outfits, repval=1):
    # Replace all values in MOSAIC multi-extension fits file by
    # repvalue

    hdulist = fits.open(infits)
    Nima = len(hdulist)
    print("Replacing with n=%s on %s images, %s --> %s" % (
        repval, Nima - 1, infits, outfits), file=sys.stderr)
    for i in range(Nima)[1:]:
        hdulist[i].data = (hdulist[i].data * 0 + 1).astype('UInt8')
    hdulist.verify('fix')
    hdulist.writeto(outfits)
    hdulist.close()
    return

def put_headerKeywords(origFile, newFile, keywords, ra, dec):
    ''' This takes all of the keywords from the original fits image and makes
    sure they are present in the new fits image. This is important for the
    mosaic3 data. For that data it doesn't seem to propigate the header
    keywords correctly. This function makes sure those keywords are present.

    '''

    print(origFile, newFile)
    with fits.open(origFile) as oF:
        with fits.open(newFile, mode='update') as nF:
            oHDR = oF[0].header # original header
            nHDR = nF[0].header # new header

            print('# Updating {} with new keywords'.format(newFile))
            for kw in keywords:
                try:
                    nHDR[kw]
                except KeyError: # only update if it doesn't exist
                    try:
                        if kw == 'RA':
                            nHDR[kw] = ra
                        elif kw == 'DEC':
                            nHDR[kw] = dec
                        else:
                            nHDR[kw] = oHDR[kw]
                    except KeyError: # sometimes the original keys aren't there
                        print('#{} Not Found'.format(kw))
    return

def put_airmass(filename, airmass):

    # Open the file, update mode
    with fits.open(filename, mode="update") as f: # open a FITS file
        hdr = f[0].header  # the primary HDU header

        print("# Updating %s with AIRMASS=%s after SWarp" % (
            filename, airmass), file=sys.stderr)
        hdr['AIRMASS'] = (airmass, 'After SWarp equivalent airmass')

        # check the file
        f.verify('fix')

    return

def put_exptime(filename, exptime):

    # Open the file, update mode
    with fits.open(filename, mode="update") as f: # open a FITS file
        hdr = f[0].header  # the primary HDU header

        print("# Updating %s with EXPTIME=%s after SWarp" % (
            filename, exptime), file=sys.stderr)
        hdr['EXPTIME'] = (exptime, 'After SWarp equivalent exptime (sec)')

        # check the file
        f.verify('fix')
    return

def chtype_replace_nonzero(infits, repval=1, out=None, verb=None):

    # Replace all values in !=0 in fits file to 1
    # certain value file by repvalue and changes the type to UInt8

    if verb:
        print("\tChange type UInt8 on %s -- Replacing pixels not eq 0 --> %s"
              % (infits, repval), file=sys.stderr)

    hdulist = fits.open(infits, mode='update')
    Nima = len(hdulist)
    if Nima == 1:
        ima = hdulist[0].data
        ima = np.where(ima != 0, repval, ima).astype('UInt8')
        hdulist[0].data = ima
    else:
        for i in range(Nima)[1:]:
            ima = hdulist[i].data
            ima = np.where(ima != 0, repval, ima).astype('UInt8')
            hdulist[i].data = ima
    hdulist.verify('silentfix')
    hdulist.flush()
    hdulist.close()
    return


def chtype_fits(infits, type='UInt8', out=None, verb=None):
    ''' Don't really know what this is supposed to do. '''

    if verb:
        print("\tChange type on %s to ---> %s" % (infits, type),
              file=sys.stderr)

    hdulist = fits.open(infits, mode='update')
    Nima = len(hdulist)
    if Nima == 1:
        ima = (hdulist[0].data).astype(type)
        hdulist[0].data = ima
    else:
        for i in range(Nima)[1:]:
            ima = (hdulist[i].data).astype(type)
            hdulist[i].data = ima
    hdulist.verify('silentfix')
    hdulist.flush()
    hdulist.close()
    return

def check_exe(exe, verb="yes"):
    ''' Checks to make sure we have the appropriate system command available.
    If we don't it raises an exception.

    '''

    path = os.environ['PATH'].split(':')
    for p in path:
        f = os.path.join(p, exe)
        if os.path.isfile(f):
            if verb:
                print("# Found %s in %s" % (exe, f), file=sys.stderr)
            return True
    # it wasn't found
    print("# ERROR: Couldn't find %s" % exe, file=sys.stderr)
    raise FileNotFoundError(exe)
    return False

def elapsed_time(t1, text=''):
    import time
    t2 = time.time()
    hh = int((t2 - t1) / 3600.)
    mm = int((t2 - t1) / 60 - hh * 60)
    ss = (t2 - t1) - 60 * mm - 3600 * hh
    print("Elapsed time: %dh %dm %2.2fs %s" % (hh, mm, ss, text))
    print("Elapsed time: %dh %dm %2.2fs %s" % (hh, mm, ss, text),
          file=sys.stderr)
    #print >>sys.stderr,"Elapsed time: %dm %2.2fs" % ( int( (t2-t1)/60.),
    #(t2-t1) - 60*int((t2-t1)/60.))
    return


def elapsed_time_str(t1):
    import time
    t2 = time.time()
    hh = int((t2 - t1) / 3600.)
    mm = int((t2 - t1) / 60 - hh * 60)
    ss = (t2 - t1) - 60 * mm - 3600 * hh
    return "Elapsed time: %dh %dm %2.4fs" % (hh, mm, ss)


# Parses the header of a SExtractor catalog and stores the column
# number (in python' order starting from zero) for every keyword columns
# It returs a dictionary with the column number for evry keyword.
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
        except IndexError:
            continue

    return SExcols

def write_ds9(d):
    with open('sdss_specz2.reg', 'wt') as f:
        f.write('# Region file formart: DS9 version 4.1\n')
        f.write('# Filename: sdss_specz.reg\n')
        mask = np.where(d['Z_S'] > -1)[0]
        for ra, dec, z in zip(d['RA'][mask], d['DEC'][mask], d['Z_S'][mask]):
            f.write('fk5;circle({},{},2")'.format(ra, dec))
            f.write('# width=2 text={' + str(z) + '} \n')

########################################
# Read and write fits simple functions
#######################################


# To read in simple fits file
def rfits(filename, hdu=0):

    ff = fits.open(filename, "readonly")
    data = ff[hdu].data
    hdr = ff[hdu].header
    ff.close()
    return data, hdr


# To write in a simple file
def writefits(array, filename):
    import os
    # Writes fits files of a given array
    newfits = fits.HDUList()
    hdu = fits.PrimaryHDU()
    hdu.data = array
    newfits.append(hdu)

    # Remove old version of file before
    if os.path.isfile(filename):
        os.remove(filename)

    newfits.writeto(filename)
    newfits.close
    return


wfits = writefits

# Transform form Dec to Sexagesimal hh:mm:ss format
def dec2sex(hms):

    hh = int(hms)
    mm = int((hms - hh) * 60)
    ss = ((hms - hh) * 60 - mm) * 60

    if abs(mm) < 10:
        mm = "0%d" % abs(mm)
    else:
        mm = "%2d" % abs(mm)

    if abs(ss) < 10:
        ss = "0%.2f" % abs(ss)
    else:
        ss = "%5.2f" % abs(ss)

        #return "%2d:%d:%.2f" % (hh,abs(mm),abs(ss))
    return "%s:%s:%s" % (hh, mm, ss)

# Taken from APSIS fUtil
def deNAN(a, value=0.0):
    nans = np.logical_not(np.less(a, 0.0) + np.greater_equal(
        a, 0.0))
    return np.where(nans, value, a)

# Create a mask file from the weights
def mask_from_weight(infile, outfile, value=0):

    (data, hdr) = rfits(infile)

    newdata = np.where(data > 0, 1, 0)

    newfits = fits.HDUList()
    hdu = fits.PrimaryHDU()
    #hdu.data   = newdata.astype("Int16") # Just as integer
    hdu.data = newdata.astype("UInt8")  # Just as ushort integer
    hdu.header = hdr

    # Remove old version of file before
    if os.path.isfile(outfile):
        os.remove(outfile)

    newfits.append(hdu)
    newfits.writeto(outfile)
    newfits.close
    return

# Create a mask file from the weights
def weight_from_dqfile(infile, outfile, clobber=False):
    # Remove old version of file before
    if os.path.isfile(outfile):
        if clobber:
            os.remove(outfile)
        else:
            print("# Skipping creation, %s image exists" % outfile)
            return

    with fits.open(infile, mode='readonly') as inHDU:
        if len(inHDU) > 2:
            newhdu = fits.HDUList()
            newhdu.append(inHDU[0])
            for hdu in inHDU[1:]:
                data = hdu.data
                newdata = np.where(data > 0, 0, 1)
                hdu.data = newdata.astype("UInt8")  # Just as ushort integer
                newhdu.append(hdu)
        else:
            data = inHDU[-1].data
            newdata = np.where(data > 0, 0, 1)
            inHDU[-1].data = newdata.astype("UInt8")  # Just as ushort integer
            newhdu = fits.HDUList()
            newhdu.append(inHDU[-1])

        newhdu.writeto(outfile, clobber=clobber)
        newhdu.close()
    return
