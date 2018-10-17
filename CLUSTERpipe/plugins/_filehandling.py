import time
import sys
import numpy
import re
import os
import scipy.misc as sci_misc
try:
    import tableio
    import extras
except ImportError:
    sys.path.append('/home/boada/Projects/planckClusters/MOSAICpipe/pipe_utils')
    import tableio
    import extras

sout = sys.stderr
lor = numpy.logical_or
nstr = numpy.char

#########################################
# Read in the big catalog of photometry
#########################################
def read_cat(self):
    t1 = time.time()

    cols = (1, 2, 23, 27, 26, 28, 29, 30, 3, 4, 6, 7, 9, 10, 12, 13, 15,
            16, 17, 18, 19, 20, 21, 22, 31, 32, 33, 34, 35, 36)

    sout.write("# Reading cols:%s\n# Reading cats from: %s... \n" %
               (cols, self.catsfile))
    (ra,
    dec,
    z_b,
    odds,
    t_b,
    z_ml,
    t_ml,
    chi,
    g,
    g_err,
    r,
    r_err,
    i,
    i_err,
    z,
    z_err,
    g_bpz,
    g_berr,
    r_bpz,
    r_berr,
    i_bpz,
    i_berr,
    z_bpz,
    z_berr,
    class_star,
    a_image,
    b_image,
    theta,
    x_image,
    y_image) = tableio.get_data(self.catsfile, cols=cols)

    (id) = tableio.get_str(self.catsfile, cols=(0, ))

    ############################################
    # Choose the photo-z to use, ml or bayesian
    ############################################
    sout.write("# Will use %s redshifts\n" % self.zuse)
    if self.zuse == "ML":
        z_ph = z_ml
        # t = t_ml
    elif self.zuse == "ZB":
        z_ph = z_b
        # t = t_b

    i_lim = self.maglim
    odds_lim = 0.80  # not currently used
    star_lim = self.starlim

    # Clean up according to BPZ
    sout.write("# Avoiding magnitudes -99 and 99 in BPZ \n")
    g_mask = numpy.where(lor(g_bpz == 99, g_bpz == -99), 0, 1)
    r_mask = numpy.where(lor(r_bpz == 99, r_bpz == -99), 0, 1)
    i_mask = numpy.where(lor(i_bpz == 99, i_bpz == -99), 0, 1)
    z_mask = numpy.where(lor(z_bpz == 99, z_bpz == -99), 0, 1)
    bpz_mask = g_mask * r_mask * i_mask * z_mask

    # Clean up to avoid 99 values and very faint i_mag values
    sout.write("# Avoiding magnitudes 99 in MAG_AUTO \n")
    #g_mask = numpy.where( g >= 99,    0 , 1)
    #r_mask = numpy.where( r >= 99,    0 , 1)
    #i_mask = numpy.where( i >= i_lim, 0 , 1)
    #z_mask = numpy.where( z >= 99,    0 , 1)
    sout.write("# Avoiding magnitudes i > %s in MAG_AUTO \n" % i_lim)

    # Clean by class_star
    sout.write("# Avoiding CLASS_STAR > %s \n" % star_lim)
    mask_star = numpy.where(class_star > star_lim, 0, 1)

    # Clean up by odds
    #sout.write( "# Avoiding ODDS < %s in BPZ \n" % odds_lim)
    odds_mask = numpy.where(odds > odds_lim, 1, 0)
    odds_mask = 1

    # Avoid z> zlim objects too.
    #sout.write( "# Avoiding objects with z > %s " % self.zlim)
    zp_mask = numpy.where(z_ph > self.zlim, 0, 1)
    zp_mask = 1

    # Clean up by BPZ type
    # sout.write('# Avoiding objects with type > %s' % t)
    tp_mask = 1

    # The final 'good' mask
    mask_good = bpz_mask * odds_mask * mask_star * zp_mask * tp_mask
    idx = numpy.where(mask_good == 1)

    # Make ids a Char String in numarray
    self.id = nstr.array(id)[idx]

    # Only keep the 'good' one, avoid -99 and 99 values in BPZ mags
    self.ra = ra[idx]
    self.dec = dec[idx]
    self.z_b = z_b[idx]
    self.odds = odds[idx]

    self.z_ml = z_ml[idx]
    self.t_ml = t_ml[idx]
    self.t_b = t_b[idx]
    self.t_ml = t_ml[idx]
    self.chi = chi[idx]

    ############################################
    # Choose the photo-z to use, ml or bayesian
    ############################################
    if self.zuse == "ML":
        self.z_ph = self.z_ml
        self.type = self.t_ml
    elif self.zuse == "ZB":
        self.z_ph = self.z_b
        self.type = self.t_b

    self.g = g[idx]
    self.r = r[idx]
    self.i = i[idx]
    self.z = z[idx]
    self.g_err = g_err[idx]
    self.r_err = r_err[idx]
    self.i_err = i_err[idx]
    self.z_err = z_err[idx]

    self.g_bpz = g_bpz[idx]
    self.r_bpz = r_bpz[idx]
    self.i_bpz = i_bpz[idx]
    self.z_bpz = z_bpz[idx]
    self.g_berr = g_berr[idx]
    self.r_berr = r_berr[idx]
    self.i_berr = i_berr[idx]
    self.z_berr = z_berr[idx]

    self.class_star = class_star[idx]
    self.a_image = a_image[idx]
    self.b_image = b_image[idx]
    self.theta = theta[idx]
    self.x_image = x_image[idx]
    self.y_image = y_image[idx]

    # Color of selected galaxies
    self.gr = self.g_bpz - self.r_bpz
    self.ri = self.r_bpz - self.i_bpz
    self.iz = self.i_bpz - self.z_bpz

    # Min and and max values in RA/DEC
    self.ramin = self.ra.min()
    self.ramax = self.ra.max()
    self.decmin = self.dec.min()
    self.decmax = self.dec.max()

    self.idx_cat = idx

    sout.write(" \tDone: %s\n" % extras.elapsed_time_str(t1))
    return

####################################
# Read in the the probabilty file
####################################
def read_probs(self):

    # The reg expresion to compile
    regexp_point = re.compile(r"arange\("
                              r"(?P<z1>[0-9]+.[0-9]+),"
                              r"(?P<z2>[0-9]+.[0-9]+),"
                              r"(?P<dz>[0-9]+.[0-9]+)\)")
    t0 = time.time()
    sout.write("# Reading probs from :%s... \n" % self.probsfile)

    # probability arrays
    probs = []
    for line in open(self.probsfile).readlines():

        fields = line.split()
        if fields[0][0] == "#":
            point = regexp_point.search(line)
            # Extract the information if a point was selected
            if point:
                z1 = float(point.group('z1'))
                z2 = float(point.group('z2'))
                dz = float(point.group('dz'))
                zx = numpy.arange(z1, z2, dz)
            continue
        probs.append(numpy.asarray(list(map(float, fields[1:]))))

    # Transform the list into an N array
    p_z = numpy.asarray(probs)

    # select same galaxies as in catalogs we just read
    self.p_z = p_z[self.idx_cat][:]
    self.zx = zx
    sout.write(" \tDone: %s\n" % extras.elapsed_time_str(t0))

    t1 = time.time()
    # Get the 1-sigma z1, z2 limits for each galaxy
    # Cumulatibe P(<z) function for each selected galaxy
    self.Psum = numpy.cumsum(self.p_z, axis=1)
    sout.write("# Getting +/- 1sigma (z1,z2) limits for each galaxy \n")
    self.z1 = self.ra * 0.0
    self.z2 = self.ra * 0.0

    # One by one in the list
    for i in range(len(self.ra)):
        i1 = numpy.where(self.Psum[i, :] >= 0.159)[0][0]
        i2 = numpy.where(self.Psum[i, :] > 0.842)[0][0]
        self.z1[i] = self.zx[i1]
        self.z2[i] = self.zx[i2]

    sout.write(" \tDone: %s\n" % extras.elapsed_time_str(t1))
    return

##################################################
# Read in the jpg file and corresponding fitsfile
##################################################
def jpg_read(self, dx=1200, dy=1200, RA=None, DEC=None):

    # The fitsfile with the wcs information
    self.fitsfile = os.path.join(self.datapath, self.ctile + 'i.fits')
    if os.path.isfile(self.datapath + self.ctile + 'irg.tiff'):
        self.jpgfile = os.path.join(self.datapath, self.ctile + 'irg.tiff')
    else:
        self.jpgfile = os.path.join(self.datapath, self.ctile + '.tiff')
    t0 = time.time()
    print("# Reading %s" % self.jpgfile, file=sys.stderr)
    self.jpg_array = sci_misc.imread(self.jpgfile)
    print("\tDone in %.3f sec." % (time.time() - t0), file=sys.stderr)

    # Get the shape of the array
    print('# Orig. Image Size: %s, %s, %s' % self.jpg_array.shape)
    (self.ny, self.nx, self.nz) = self.jpg_array.shape

    # pass up to class
    self.RA = RA
    self.DEC = DEC

    if float(dx) < 0 or float(dy) < 0:
        self.dx = self.nx / 2.0
        self.dy = self.ny / 2.0
    elif dx == 0 or dy == 0:
        self.dx = float(self.nx)
        self.dy = float(self.ny)
    else:
        self.dx = float(dx)
        self.dy = float(dy)

    if isinstance(RA, str) and isinstance(DEC, str):
        if ':' in RA and ':' in DEC:
            RA = astrometry.hms2dec(RA)
            DEC = astrometry.deg2dec(DEC)
            self.xo, self.yo = astrometry.rd2xy(RA, DEC, self.fitsfile)
            print(self.xo, self.yo)
            yo_tmp = self.yo
        else:
            RA = float(RA)
            DEC = float(DEC)
    elif isinstance(RA, float) and isinstance(DEC, float):
        self.xo = RA
        self.yo = DEC
        yo_tmp = self.yo
    elif RA is None and DEC is None:
        self.xo = self.nx / 2.0
        self.yo = self.ny / 2.0
        yo_tmp = self.yo
    else:
        print('Center not understood')
        sys.exit()

    #print(self.xo, self.yo)

    # a little fix when not centered -- it has something to do with the way
    # the image is being displayed. Without this fix, the the catalog is in
    # the correct non-centered position, but the image was off. This fixes
    # that bug, but I didn't think hard enough to understand why.
    self.yo = self.ny - yo_tmp

    # Select the limits of the image to display
    x1 = int(self.xo - self.dx)
    x2 = int(self.xo + self.dx)
    y1 = int(self.yo - self.dy)
    y2 = int(self.yo + self.dy)

    # this is part of the little fix described above.
    self.yo = yo_tmp

    # Get the region to use for plotting
    self.jpg_region = self.jpg_array[y1:y2, x1:x2, :]
    (self.ny, self.nx, self.nz) = self.jpg_region.shape

    # print the cropped region's size
    print('# New Image Size: %s, %s, %s' % self.jpg_region.shape)
    print('\txo : %s' % self.xo)
    print('\tyo : %s' % self.yo)
    print('\tdx : %s' % self.dx)
    print('\tdy : %s' % self.dy)

    return

def write_info(self, blank=False):
    ''' This writes out the info for the cluster we found. If we didn't
    find a cluster it still writes out a file with nothing but zeros. That
    will make sure we always have a file even for the fields where we
    didn't find anything.

    '''

    if blank:
        filename = os.path.join(self.datapath, self.ctile,
                                self.ctile + ".info")

        print("# Will write info to %s" % filename)
        s = open(filename, "w")
        head = ("# %-18s %12s %12s %7s %7s %7s %5s %5s %10s %10s %8s %8s "
                "%8s %8s %8s %8s %8s %11s\n" %
                ('ID_BCG', 'RA', 'DEC', 'zBCG', 'z_cl', 'z_clerr', 'Ngal',
                 'Ngal_c', 'L_i', 'L_iBCG', 'Mr', 'Mi', 'r', 'i', 'p_BCG',
                 'R[kpc]', 'area[%]', 'Confidence'))
        format = ("%20s %12s %12s %7.3f %7.3f %7.3f %5d %5d %10.3e %10.3e"
                  "%8.2f %8.2f %8.2f %8.2f %8.3f %8.1f %8.2f %2d\n")

        vars = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        s.write(head)
        s.write(format % vars)
        s.close()
        return

    filename = os.path.join(self.datapath, self.ctile,
                            self.ctile + ".info")
    print(" Will write info to %s" % filename)

    i = self.iBCG

    RA = astrometry.dec2deg(self.ra[i] / 15.)
    DEC = astrometry.dec2deg(self.dec[i])

    s = open(filename, "w")
    head = ("# %-18s %12s %12s %7s %7s %7s %5s %5s %10s %10s %8s %8s "
            "%8s %8s %8s %8s %8s %11s\n" %
            ('ID_BCG', 'RA', 'DEC', 'zBCG', 'z_cl', 'z_clerr', 'Ngal',
             'Ngal_c', 'L_i', 'L_iBCG', 'Mr', 'Mi', 'r', 'i', 'p_BCG',
             'R[kpc]', 'area[%]', 'Confidence'))
    format = ("%20s %12s %12s %7.3f %7.3f %7.3f %5d %5d %10.3e %10.3e"
              "%8.2f %8.2f %8.2f %8.2f %8.3f %8.1f %8.2f %2d\n")

    vars = (self.ID_BCG, RA, DEC, self.z_ph[i], self.z_cl, self.z_clerr,
            self.Ngal, self.Ngal_c, self.Lsum, self.Lr[i], self.Mr[i],
            self.Mi[i], self.r[i], self.i[i], self.p[i], self.radius,
            self.area_fraction, self.confidence)
    s.write(head)
    s.write(format % vars)
    s.close()
    return

# Write down all relevant redshift information
def write_redshift(self):

    filename = os.path.join(self.datapath, self.ctile,
                            self.ctile + ".redshifts")
    print("Will write redshits to %s" % filename)

    # The BCG redshift interval
    iBCG = self.iBCG
    zx = self.zx
    dz = zx[1] - zx[0]
    zo = self.z_ph[iBCG]
    tck = scipy.interpolate.splrep(zx, self.p_z[iBCG], s=0)
    ds = 0.001
    x = numpy.arange(zx[2], zx[-3], ds)
    y = scipy.interpolate.splev(x, tck, der=0)
    y = numpy.where(y < 0, 0, y)

    # Get the 1 sigma error bars (68.2%)
    Psum = numpy.cumsum(y * ds / dz)
    i1 = numpy.where(Psum >= 0.159)[0][0]
    i2 = numpy.where(Psum > 0.841)[0][0]
    z1_68 = x[i1]
    z2_68 = x[i2]

    # And the 2-sigma (95.4%)
    i1 = numpy.where(Psum >= 0.023)[0][0]
    z1_95 = x[i1]
    try:
        i2 = numpy.where(Psum > 0.977)[0][0]
    except IndexError:
        z2_95 = -1
    z2_95 = x[i2]

    s = open(filename, "w")
    head = ("# %-18s " + "%8s " * 7 + "\n") % ('ID_BCG', 'z_cl', 'err',
                                               'zBCG', 'z1', 'z2', 'z1',
                                               'z2')
    format = ("%20s " + "%8.3f " * 7 + "\n")
    vars = (self.ID_BCG, self.z_cl, self.z_clerr, zo, z1_68, z2_68, z1_95,
            z2_95)
    s.write(head)
    s.write(format % vars)
    s.close()
    return

# Write the members information out
def write_members(self):

    filename = os.path.join(self.datapath, self.ctile,
                            self.ctile + ".members")
    m = open(filename, "w")
    print("Will write members to %s" % filename)

    head = ("# %-23s %15s %15s  %6s %6s %6s %6s %6s  %6s  %6s  "
            "%6s  %6s  %6s \n" %
            ("ID", "RA", "DEC", "ZB", "TB", "ZML", "TML", "g_mag", "g_err",
             "r_mag", "r_err", "i_mag", "i_err"))
    m.write(head)
    for i in self.iRadius[0]:
        format = "%-25s %15f %15f  %6.3f %6.3f %6.3f %6.2f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n"
        vars = (self.id[i], self.ra[i], self.dec[i], self.z_b[i],
                self.t_b[i], self.z_ml[i], self.t_ml[i], self.g[i],
                self.g_err[i], self.r[i], self.r_err[i], self.i[i],
                self.i_err[i])

        m.write(format % vars)
    m.close()
    return
