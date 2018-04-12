#!/usr/bin/env python

import os
import sys
import numpy
import math
import re
import time
import scipy
import scipy.misc as sci_misc
import pylab
from cosmopy import cosmopy
import matplotlib.patches
from past.utils import old_div

try:
    import cosmology
except ImportError:
    sys.path.append('/home/boada/Projects/MOSAICpipe')
    import cosmology
import aux
import tableio
import extras
import astrometry
# fix large image error
import PIL
PIL.Image.MAX_IMAGE_PIXELS = None

Polygon = matplotlib.patches.Polygon

float32 = numpy.float32
float64 = numpy.float64
int16 = numpy.int16
nstr = numpy.char
land = numpy.logical_and
lge = numpy.greater_equal
lle = numpy.less_equal
lor = numpy.logical_or

sout = sys.stderr

class finder:

    def __init__(self, ctile, maglim=25.0,
                 pixscale=0.25,
                 zlim=1.8,
                 zo=None,
                 dz=0.05,
                 radius=1000.0, # Radius in kpc
                 cosmo=(0.3, 0.7, 0.7),
                 zuse="ZB", # Use ZB (Bayesian) or ML (Max Like)
                 outpath='plots',
                 path='./',
                 evolfile="0_1gyr_hr_m62_salp_Ks.color",
                 p_lim=0.4,
                 verb='yes'):

        # Check for environ vars
        self.home = os.environ['HOME']
        if not os.getenv('MOSAICpipe'):
            os.environ['MOSAICpipe'] = os.path.join(self.home, 'Projects',
                                                    'planckClusters',
                                                    'MOSAICpipe')
        self.MOSAICpipe = os.getenv('MOSAICpipe')

        self.zlim = zlim
        self.cosmo = cosmo
        self.evolfile = os.path.join(self.MOSAICpipe, "LIB/evol", evolfile)
        self.dz = dz
        self.zuse = zuse
        self.outpath = outpath
        self.path = path
        self.radius = radius
        self.zo = zo  # Input zo for the BCG

        # Set the cosmology now
        self.cset = cosmopy.set(self.cosmo)
        self.Om = cosmo[0]
        self.OL = cosmo[1]
        self.h = cosmo[2]
        self.Ho = self.h * 100.0

        self.ctile = ctile
        self.datapath = path
        self.catsfile = os.path.join(path, ctile, ctile + "_merged.cat")
        self.probsfile = os.path.join(path, ctile, ctile + "_probs.dat")
        self.maglim = maglim

        self.verb = verb
        self.ellipse = {}
        self.plot_around = None
        self.pixscale = pixscale

        self.read_cat()  # Read catalogs avoding, faint, high-z and 99 objects
        self.read_probs()  # Read probs function of objects from catalogs
        self.get_absmags()  # We compute Abs Mags for each object

        # Set the BCG masked to False, so we select BCGs on the firts run
        self.BCG_masked = False
        self.BCG_probs = False

        # Get the BCG candidates
        self.get_BCG_candidates(p_lim=p_lim)

        # Check if the destination folder exists
        if os.path.exists(self.outpath):
            sout.write("# Will put files to: %s\n" % self.outpath)
        else:
            sout.write("# Will create new folder: %s\n" % self.outpath)
            os.mkdir(self.outpath)

        self.rootname = os.path.join(self.datapath, self.ctile, self.ctile)

        return

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
        odds_lim = 0.80 # not currently used
        star_lim = 0.80

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

        sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t1))
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
        sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t0))

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

        sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t1))
        return

    ################################################
    # Get the absolute magnitudes for each object
    ################################################
    def get_absmags(self):

        # Distance modulus, dlum and dangular
        self.dlum = self.cset.dlum(self.z_ph)
        self.dang = self.cset.dang(self.z_ph)
        self.DM = 25.0 + 5.0 * numpy.log10(self.dlum)

        t0 = time.time()
        # Get the absolute magnitudes, *not including evolution*, only Kcorr
        # We use a BPZ E's template for Kcorr
        #sout.write("# Computing absolute magnitudes interpolating Kcorr ")
        #k = Kcorr_fit(sed='El_Benitez2003')

        # Alternatibely we can get both the kcorr and the evol from
        # the *.color file from BC03 *.ised file
        sout.write("# Computing absolute magnitudes interpolating "
                   "konly from BC03 model \n")
        k, ev = KEfit(self.evolfile)

        sout.write("# Computing evolution ev(z) for each galaxy \n")
        self.ev_g = ev['g'](self.z_ph)
        self.ev_r = ev['r'](self.z_ph)
        self.ev_i = ev['i'](self.z_ph)
        self.ev_z = ev['z'](self.z_ph)

        # Also get the luminosities in Msun
        # taken from http://www.ucolick.org/~cnaw/sun.html
        self.Msun = {}
        self.Msun['g'] = 5.11
        self.Msun['r'] = 4.65
        self.Msun['i'] = 4.54
        self.Msun['z'] = 4.52
        self.Msun['Ks'] = 5.14
        # ^^ http://www.astronomy.ohio-state.edu/~martini/usefuldata.html

        # Mags k-corrected
        self.Mg = self.g - self.DM - k['g'](self.z_ph)
        self.Mr = self.r - self.DM - k['r'](self.z_ph)
        self.Mi = self.i - self.DM - k['i'](self.z_ph)
        self.Mz = self.z - self.DM - k['z'](self.z_ph)

        self.Lg = 10.0**(-0.4 * (self.Mg - self.Msun['g']))
        self.Lr = 10.0**(-0.4 * (self.Mr - self.Msun['r']))
        self.Li = 10.0**(-0.4 * (self.Mi - self.Msun['i']))
        self.Lz = 10.0**(-0.4 * (self.Mz - self.Msun['z']))
        self.Lg_err = self.Lg * self.g_err / 1.0857
        self.Lr_err = self.Lr * self.r_err / 1.0857
        self.Li_err = self.Li * self.i_err / 1.0857
        self.Lz_err = self.Lz * self.z_err / 1.0857

        # Pass it up to the class
        self.kcorr = k
        self.evf = ev
        sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t0))
        return

    ##################################################################
    # Define the sub-sample for BCGs candidates around a position
    ##################################################################
    def get_BCG_candidates(self, Mr_limit=-22.71, p_lim=1e-4):

        t0 = time.time()
        sout.write("# Computing p_BCG probabilities... ")

        # The Abs mag limit @ z=0.1 in the i-band
        Mi_limit = cosmology.reobs('El_Benitez2003',
                                   m=Mr_limit,
                                   oldfilter="r_MOSAICII",
                                   newfilter="i_MOSAICII")

        # Evaluate the genertic mask for BCG only onece
        if not self.BCG_probs:

            # We get the limit at the z_ph of each candidate, corrected by z=0.1
            Mr_BCG_limit = Mr_limit + self.ev_r - self.evf['r'](
                0.1)  # + self.DM_factor
            Mi_BCG_limit = Mi_limit + self.ev_i - self.evf['i'](
                0.1)  # + self.DM_factor

            # Evaluate the BCG Probability function, we
            # get the limit for each object
            self.p = p_BCG(self.Mr, Mr_BCG_limit)

            self.BCG_probs = True

            i_lim = 25.0
            star_lim = 0.8
            p_lim = max(self.p) * 0.8
            sout.write("# Avoiding BCG_prob < %.3f in BGCs\n" % p_lim)
            mask_p = numpy.where(self.p >= p_lim, 1, 0)
            mask_g = numpy.where(self.g < i_lim + 5, 1, 0)
            mask_r = numpy.where(self.r < i_lim + 2, 1, 0)
            mask_i = numpy.where(self.i < i_lim, 1, 0)
            mask_z = numpy.where(self.z < i_lim + 1, 1, 0)
            mask_t = numpy.where(self.type <= 2.0, 1, 0)

            # Avoid freakishly bright objects, 2.5 mags brighter than the
            # M_BCG_limit
            mask_br = numpy.where(self.Mr > Mr_BCG_limit - 2.5, 1, 0)
            mask_bi = numpy.where(self.Mi > Mi_BCG_limit - 2.5, 1, 0)

            # Put a more strict cut in class_star for bcg candidates
            sout.write("# Avoiding CLASS_STAR > %s in BGCs\n" % star_lim)
            mask_star = numpy.where(self.class_star <= star_lim, 1, 0)

            # Construct the final mask now
            self.mask_BCG = (mask_t * mask_g * mask_r * mask_i * mask_z *
                                mask_br * mask_bi * mask_p * mask_star)

            self.BCG_masked = True

            # Model color only once
            #self.zx = numpy.arange(0.01, self.zlim, 0.01)
            self.gr_model = cosmology.color_z(sed='El_Benitez2003',
                                              filter_new='g_MOSAICII',
                                              filter_old='r_MOSAICII',
                                              z=self.zx,
                                              calibration='AB')
            self.ri_model = cosmology.color_z(sed='El_Benitez2003',
                                              filter_new='r_MOSAICII',
                                              filter_old='i_MOSAICII',
                                              z=self.zx,
                                              calibration='AB')
            self.iz_model = cosmology.color_z(sed='El_Benitez2003',
                                              filter_new='i_MOSAICII',
                                              filter_old='z_MOSAICII',
                                              z=self.zx,
                                              calibration='AB')

            sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t0))

        # Select the candidates now
        idx = numpy.where(self.mask_BCG == 1)

        # And pass up to to class
        self.idx_BCG = idx
        self.id_BCG = self.id[idx]
        self.ra_BCG = self.ra[idx]
        self.dec_BCG = self.dec[idx]
        self.p_BCG = self.p[idx]
        self.z_BCG = self.z_ph[idx]
        self.t_BCG = self.type[idx]
        self.N_BCG = len(idx[0])
        self.Mi_BCG = self.Mi[idx]
        self.Mr_BCG = self.Mr[idx]
        self.DM_BCG = self.DM[idx]  # distance modulus
        self.dang_BCG = self.dang[idx]  # distance modulus

        self.zml_BCG = self.z_ml[idx]
        self.tml_BCG = self.t_ml[idx]
        self.zb_BCG = self.z_b[idx]
        self.tb_BCG = self.t_b[idx]
        self.class_BCG = self.class_star[idx]
        self.a_BCG = self.a_image[idx]
        self.b_BCG = self.b_image[idx]
        self.theta_BCG = self.theta[idx]

        # r,i-band stuff
        self.r_BCG = self.r[idx]
        self.i_BCG = self.i[idx]

        # Get the 1-sigma intervals
        self.z1_BCG = self.z1[idx]
        self.z2_BCG = self.z2[idx]

        # The r-band Luminosity of the BCGs
        self.LBCG = self.Lr[idx]

        # The distance to the candidate's position for each BCG, in arcmin
        sout.write("# Found %s BCG candidates\n" % self.N_BCG)

        return

    ##################################################
    # Read in the jpg file and corresponding fitsfile
    ##################################################
    def jpg_read(self, dx=1200, dy=1200, RA=None, DEC=None):

        # The fitsfile with the wcs information
        self.fitsfile = os.path.join(self.datapath,
                                     self.ctile + 'Detec.fits')
        self.jpgfile = os.path.join(self.datapath, self.ctile + 'irg.tiff')
        t0 = time.time()
        print("Reading %s" % self.jpgfile, file=sys.stderr)
        self.jpg_array = sci_misc.imread(self.jpgfile)
        print("Done in %.3f sec." % (time.time() - t0), file=sys.stderr)

        # Get the shape of the array
        print('Orig. Image Size: %s, %s, %s' % self.jpg_array.shape)
        (self.ny, self.nx, self.nz) = self.jpg_array.shape

        # pass up to class
        self.RA = RA
        self.DEC = DEC

        if float(dx) < 0 or float(dy) < 0:
            self.dx = self.nx / 2.0
            self.dy = self.ny / 2.0

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

        print(self.xo, self.yo)

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
        print('New Image Size: %s, %s, %s' % self.jpg_region.shape)
        print('xo : %s' % self.xo)
        print('yo : %s' % self.yo)
        print('dx : %s' % self.dx)
        print('dy : %s' % self.dy)

        return

    ##############################
    # Change the axes to arcmins
    ###############################
    def ax_to_arcmin(self, ds=1.0):  # ds in arcmin

        [xmin, xmax, ymin, ymax] = pylab.axis()
        dx = 2 * self.dx
        dy = 2 * self.dy
        scale = self.pixscale / 60.  # in arcmin

        xo = (xmin + xmax) / 2.0
        yo = (ymin + ymax) / 2.0
        s1 = int((-dx / 2.0) * scale)
        s2 = int((+dx / 2.0) * scale)

        sx = numpy.arange(s1, s2 + 0.05, ds)

        xtext = []
        xtick = []
        for s in sx:
            x = xo + s / scale  # + ds/scale
            xtick.append(x)
            xtext.append("%.1f" % s)

        s1 = int((-dy / 2.0) * scale)
        s2 = int((+dy / 2.0) * scale)
        sy = numpy.arange(s1, s2 + 0.05, ds)

        ytext = []
        ytick = []
        for s in sy:
            y = yo + s / scale  # + ds/scale
            ytick.append(y)
            ytext.append("%.1f" % s)
        pylab.yticks(ytick, tuple(ytext))
        pylab.xticks(xtick, tuple(xtext))
        # Make sure we plot everithing
        pylab.xlim(xmin, xmax)
        pylab.ylim(ymin, ymax)
        return

    ####################################
    # The loop click-n-display routine
    ####################################
    def jpg_display(self):

        # Measure time
        pylab.close('all')
        t0 = time.time()
        print("Displaying... be patient", file=sys.stderr)
        # Display
        self.ax = pylab.figure(1, figsize=(10, 10))
        pylab.imshow(self.jpg_region, origin='upper', interpolation='bicubic')
        # Change ax to arcmin
        self.ax_to_arcmin(ds=2.0)
        pylab.xlabel("x[arcmin]")
        pylab.ylabel("y[arcmin]")
        pylab.title(self.ctile)
        nameA = self.ctile + "_figA.png"
        nameB = self.ctile + "_figB.png"
        self.figAname = os.path.join(self.datapath, self.ctile, nameA)
        self.figBname = os.path.join(self.datapath, self.ctile, nameB)
        pylab.savefig(self.figAname,
                      transparent=True,
                      dpi=100,
                      bbox_inches='tight')
        print("Done in %s sec." % (time.time() - t0), file=sys.stderr)
        print("Wrote Fig A, %s " % (self.figAname), file=sys.stderr)
        # Draw the search zone
        #self.draw_zone(n=8)
        self.draw_zone2()
        # register this function with the event handler
        pylab.connect('button_press_event', self.get_object)
        pylab.show()
        return

    #####################################
    # Draw the zone where to select
    ######################################

    def draw_zone2(self):
        ''' This draws a 5' and 2' circle where we think the clusters will be. This
        should be the region we select clusters from.

        '''

        (nx, ny, nz) = self.jpg_region.shape
        center = nx / 2, ny / 2
        r1_pixels = 5 * 60.0 / self.pixscale
        r2_pixels = 2 * 60.0 / self.pixscale
        ax = pylab.gca()
        Cc1 = PCircle(center,
                    r1_pixels,
                    resolution=80,
                    fill=0,
                    edgecolor="white",
                    linestyle='dashed',
                    linewidth=0.5)
        Cc2 = PCircle(center,
                    r2_pixels,
                    resolution=80,
                    fill=0,
                    edgecolor="white",
                    linestyle='dashed',
                    linewidth=0.5)
        ax.add_patch(Cc1)
        ax.add_patch(Cc2)
        return

    def draw_zone(self, n=8):
        # Draw the region to select from
        n = float(n)
        (nx, ny, nz) = self.jpg_region.shape
        dx = nx / n
        dy = ny / n
        xx = [dx, (n - 1) * dx, (n - 1) * dx, dx, dx]
        yy = [dy, dy, (n - 1) * dy, (n - 1) * dy, dy]
        pylab.plot(xx, yy, 'w--', linewidth=0.5)
        return

    #######################################
    # This is the loop routine interactive
    #######################################
    def get_object(self, event):

        # Print the commands:
        print('Hold one of the following keys and click the image to issue '
              'associated command')
        print('q:\t quit\n'
              'b:\t identify BCG candidates\n'
              'c:\t plot potential cluster members\n'
              'j:\t plot the probabilities\n'
              'v:\t write info onto the figure\n'
              '1-3:\t write confidence info onto the figure. 1: High 3: Low\n'
              'w:\t write out the result\n'
              'h:\t recenter the figure onto the clicked location\n'
              'i:\t toggle between optical and Optical + IR image')
        print('You used:\t %s' % event.key)

        if event.key == 'q' or event.key == 'Q':
            sys.exit()
            pylab.close('all')
            return

        # Remap to right positions
        #ximage = event.xdata
        #yimage = self.ny - event.ydata
        # Remap to original coordinate system
        ximage = event.xdata + (self.xo - self.dx)
        yimage = (self.ny - event.ydata) + (self.yo - self.dy)

        #print('clicked location: %s, %s' % (event.xdata, self.ny - event.ydata))
        #print('interp. location: %s, %s' % (ximage, yimage))

        # Fins the closest one
        self.get_nearest(ximage, yimage)
        iclose = self.iclose

        # Plot BCG candidates
        if event.key == 'b':
            self.ellipse_BCGs()
            return

        # Plot around ID's redshift
        if event.key == 'c':
            if self.zo:
                self.z_ph[iclose] = self.zo
            # Select members using distance and 3sigma clipping
            z_cl, z_err = self.select_members_radius(self.iclose,
                                                     radius=self.radius,
                                                     zo=self.zo)
            z_cl, z_err = self.select_members_radius(self.iclose,
                                                     radius=self.radius,
                                                     zo=z_cl)
            self.iBCG = self.iclose
            self.ellipse_members()
            self.background()
            print("\t Ngal: %d" % self.Ngal)
            print("\t Ngal_c: %d" % self.Ngal_c)
            print("\t z_cl: %.3f +/- %.3f" % (self.z_cl, self.z_clerr))
            print("\t L   : %.3e [Lsun]" % self.Lsum)
            print("\t Mi  : %6.2f " % self.Mi[iclose])
            print("\t Mr  : %6.2f " % self.Mr[iclose])

            return

        # Plot probs
        if event.key == 'j':
            self.plot_probs()
            return
        if event.key == '1' or event.key == '2' or event.key == '3':
            try:
                self.conf_back.remove()
            except AttributeError:
                pass
            try:
                self.conf_front.remove()
            except AttributeError:
                pass

            if event.key == '1':
                conf = 'High'
            elif event.key == '2':
                conf = 'Medium'
            elif event.key == '3':
                conf = 'Low'
            else:
                # this should never happen
                return
            text = '{} Confidence'.format(conf)
            self.confidence = int(event.key)

            xo = 2 * self.dx - 80
            yo = 80
            self.conf_back = pylab.text(xo + 2,
                       self.ny - yo + 2,
                       text,
                       color='black',
                       fontsize=18, ha='right')
            self.conf_front = pylab.text(xo, self.ny - yo, text, color='white',
                                    fontsize=18, ha='right')
            pylab.draw()
            #pylab.show()
            return

        if event.key == 'v':
            try:
                self.txt_back.remove()
            except AttributeError:
                pass
            try:
                self.txt_front.remove()
            except AttributeError:
                pass

            iclose = self.iclose
            text = "z$_{cl}$ = %.3f +/- %.3f\n" \
                    "z$_{BCG}$ = %.3f\n" \
                    "N$_{galc}$ = %d (%d)\n" \
                    "R = %d[kpc]" % (self.z_cl, self.z_clerr, self.z_ph[iclose],
                                                           self.Ngal_c,
                                                           self.Ngal,
                                                           self.radius)
            xo = 80
            yo = 80
            self.txt_back = pylab.text(xo + 2,
                       self.ny - yo + 2,
                       text,
                       color='black',
                       fontsize=18)
            self.txt_front = pylab.text(xo, self.ny - yo, text, color='white',
                                    fontsize=18)
            pylab.draw()
            #pylab.show()
            return

        if event.key == 'w':
            self.write_info()
            self.write_redshift()
            self.write_members()
            pylab.savefig(self.figBname,
                          transparent=False,
                          dpi=100,
                          bbox_inches='tight')
            pylab.close('all')
            sys.exit()
            return

        if event.key == 'h':
            xloc, yloc = self.click(event)
            self.jpg_read(dx=self.dx, dy=self.dy, RA=xloc, DEC=yloc)
            self.jpg_display()

        if event.key == 'i':
            try:
                if self.toggled:
                    print('Reading Op/IR Image...')
                    self.toggled = False
            except AttributeError:
                print('Reading Optical Image...')
                self.toggled = True
            print(self.RA, self.DEC)
            self.jpg_read(dx=self.dx, dy=self.dy, RA=self.RA, DEC=self.DEC,
                          toggle=self.toggled)
            self.jpg_display()

        # Print info
        self.click(event)
        self.print_info()
        self.handle_ellipses(event)

        pylab.draw()
        pylab.show()

        return

    ########################################################
    # Modified/updated from find_clusters_ext_auto.py
    # Select galaxies around ID galaxy un redshift range
    ########################################################
    def select_members_radius(self, i, Mi_lim=-20.25, radius=500.0, zo=None):

        # Width of the redshift shell
        dz = self.dz

        t0 = time.time()
        sout.write("# Selecting Cluster members... Ngal, N200, R200 \n")
        # Get the relevant info for ith BCG
        ra0 = self.ra[i]
        dec0 = self.dec[i]
        Mi_BCG = self.Mi[i]
        #DM = self.DM[i]
        ID_BCG = self.id[i]
        if zo:
            print("Will use z:%.3f for cluster" % zo)
        else:
            zo = self.z_ph[i]
        # 1 - Select in position around ra0,dec0
        # Define radius in degress @ zo
        R = radius  # in kpc
        r = astrometry.kpc2arc(zo, R, self.cosmo) / 3600.  # in degrees.
        rcore = r / 2.0
        dist = astrometry.circle_distance(ra0,
                                          dec0,
                                          self.ra,
                                          self.dec,
                                          units='deg')
        mask_R = numpy.where(dist <= r, 1, 0)
        mask_rcore = numpy.where(dist <= rcore, 1, 0)
        arcmin2Mpc = astrometry.arc2kpc(
            zo, 60.0, self.cosmo) / 1000.0  # scale between arcmin and Mpc

        # 2 - Select in redshift
        z1 = zo - dz
        z2 = zo + dz
        mask_z = numpy.where(land(self.z_ph >= z1, self.z_ph <= z2), 1, 0)

        # 3 - Select in brightness
        Mi_lim_zo = Mi_lim + self.evf['i'](zo) - self.evf['i'](0.1)
        mask_L1 = numpy.where(self.Mi <= Mi_lim_zo, 1, 0)  # Faint  cut > 0.4L*
        mask_L2 = numpy.where(self.Mi >= Mi_BCG, 1, 0)  # Bright cut < L_BCG

        # The final selection mask, position x redshift x Luminosity
        #idx = numpy.where(mask_R * mask_L1 * mask_L2 * mask_z == 1)[0]
        idc = numpy.where(mask_rcore * mask_L1 * mask_L2 * mask_z == 1)[0]

        # Shot versions handles
        gr = self.gr
        ri = self.ri

        # Some simple 3-sigma clipping defined using r< rcore
        Nsigma = 3.0
        loop = 1
        converge = False
        while not converge:
            # The conditions to apply
            c1 = numpy.abs(gr[idc] - gr[idc].mean()) > Nsigma * numpy.std(
                gr[idc], ddof=1)
            c2 = numpy.abs(ri[idc] - ri[idc].mean()) > Nsigma * numpy.std(
                ri[idc], ddof=1)
            iclip = numpy.where(lor(c1, c2))[
                0]  # where any of the conditions fails
            if len(iclip) > 0:
                idc = numpy.delete(idc, iclip)  # Removed failed ones
                converge = False
            else:
                converge = True
            loop += 1

            # Get the average redshift within the core:
            #z_cl    = self.z_ph[idc].mean()
            #z_clrms = self.z_ph[idc].std()

        #print(idc)
        #print(self.z_ph[idc])

        # Compute the weighted average and rms
        dz = 0.5 * numpy.abs(self.z2[idc] - self.z1[idc])
        # Fix zeros
        dz[dz == 0] = 1e-5
        z_cl, z_clrms = aux.statsw(self.z_ph[idc], weight=1.0 / dz)
        sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t0))

        # Or we can make a new mask where the condition's are true
        c1 = numpy.abs(self.gr - gr[idc].mean()) > Nsigma * numpy.std(gr[idc],
                                                                      ddof=1)
        c2 = numpy.abs(self.ri - ri[idc].mean()) > Nsigma * numpy.std(ri[idc],
                                                                      ddof=1)
        mask_cm = numpy.where(lor(c1, c2), 0, 1)  # where condition fails
        iRadius = numpy.where(
                    mask_R * mask_L1 * mask_L2 * mask_z * mask_cm == 1)
        iRadius_all = numpy.where(
                    mask_L1 * mask_L2 * mask_z * mask_cm == 1)
        Ngal = len(iRadius[0])
        sout.write("# Total: %s objects selected in %s [kpc] around %s\n" %
                   (Ngal, radius, self.ID))

        # Pass up
        self.iRadius = iRadius
        self.arcmin2Mpc = arcmin2Mpc
        self.dist2BCG = dist
        self.Lsum = self.Lr[iRadius].sum()
        self.Ngal = Ngal
        self.z_cl = z_cl
        self.z_clerr = z_clrms
        self.rdeg = r  # in degress
        self.r1Mpc = r # naming fix for background estimates
        self.idc = idc  # galaxies used for mean redshift
        self.ID_BCG = ID_BCG

        # Sort indices radially for galaxies < N*R1Mpc, will be used later
        i = numpy.argsort(self.dist2BCG[iRadius_all])
        self.ix_radial = iRadius_all[0][i]

        return z_cl, z_clrms

    ##########################################
    # Compute the Background for the clusters
    ##########################################
    def background(self, k=0):
        ixr = self.ix_radial

        # No back substraction
        if self.Ngal <= 2:
            self.Ngal_c = self.Ngal
            print('Background -- Not enough galaxies found in cluster')
            return

        # Store radially ordered
        r = self.dist2BCG[ixr] * 60.0  # in arcmin
        Lr = self.Lr[ixr]  # We do in the r-band as Reyes et al

        # Bin the Ngal/Lum data in log spacing
        n = 10
        rbin = mklogarray(0.0, r.max(), n)
        Nbin, rcenter = histo(r, rbin, center='yes')
        Lbin, rcenter = bin_data(r, Lr, rbin, center='yes')

        # Compute the area in each shell
        ir = numpy.indices(rbin.shape)[0]
        ir1 = ir[:-1]
        ir2 = ir[1:]
        r1 = rbin[ir1]
        r2 = rbin[ir2]
        abin = math.pi * (r2**2 - r1**2)
        PN = old_div(Nbin, abin)  # Number Surface density

        # Compute the background median density both in Lum and Ngal
        # Between 4.0 - 9.0 r1Mpc
        R1 = 4.0 * self.r1Mpc * 60.0
        R2 = 9.0 * self.r1Mpc * 60.0
        print("# Estimating Background between R1,R2 %.2f--%2.f[arcmin]" %
              (R1, R2))

        if R2 >= r.max():
            R2 = r2.max()
            R1 = R2 - 2.0 * self.r1Mpc * 60.0
        print("# Estimating Background between R1,R2 %.2f--%2.f[arcmin]" %
              (R1, R2))

        PN_bgr = PN[land(rcenter > R1, rcenter < R2)]

        # Get the mean values for the Ngal and Lr profiles, which will
        # be the correction per arcmin^2
        PN_mean = numpy.mean(PN_bgr)
        print('mean number of BG galaxies -- {}'.format(PN_mean))

        # Total number in area
        N_bgr = PN_bgr.sum()
        area_bgr = math.pi * (R2**2 - R1**2)

        # Get the correction for Number of galaxies and Luminosoty
        # For R200 we need to recompute R200 and N200 based on new
        # R200 value.
        area_r1Mpc = math.pi * (self.r1Mpc * 60.)**2  # in arcmin2
        self.Ngal_c = self.Ngal - PN_mean * area_r1Mpc
        if self.Ngal_c < 0:
            self.Ngal_c = 0.0

        # print self.Ngal
        # print PN
        # print r1
        # print rcenter
        # print R1,R2
        # print r.min(),r.max()
        # print "PN_mean",PN_mean
        # print PN_bgr
        # print area_r1Mpc
        print("Ngal ", self.Ngal)
        print("Ngal_c", self.Ngal_c)
        # print "r200_c",self.r200_c
        # print "R200_c",self.R200_c

        self.d_Ngal_c2 = self.Ngal_c + (
            (old_div(area_r1Mpc, area_bgr))**2) * N_bgr

        # Avoid sqrt of negative number
        if self.d_Ngal_c2 < 0:
            self.d_Ngal_c = 0
        else:
            self.d_Ngal_c = math.sqrt(self.Ngal_c + ((old_div(
                area_r1Mpc, area_bgr))**2) * N_bgr)

        return

    #####################################
    # Draw the BGCs and cluster members
    #####################################
    def ellipse_members(self, k=0):

        iclose = self.iclose

        ax = pylab.gca()
        # Delete all patches, reset ellipses before redraw
        del ax.patches[2:]
        self.ellipses = {}
        #zo = self.z_ph[iclose]
        pylab.title("%s" % (self.ctile))

        # construct the ellipses for each members
        for i in self.iRadius[0]:
            #ra = self.ra[i]
            #dec = self.dec[i]
            a = self.a_image[i]
            b = self.b_image[i]
            theta = self.theta[i]  # *math.pi/180.0

            # move to cropped reference frame
            xgal = self.x_image[i] - (self.xo - self.dx)
            ygal = self.y_image[i] - (self.yo - self.dy)
            # Change the referece pixel to reflect jpg standards where the
            # origin is at (0,ny), is the upper left corner
            ygal = self.ny - ygal

            if i == self.iclose:
                ec = 'yellow'
            else:
                ec = 'red'
            E = PEllipse((xgal, ygal), (a, b),
                         resolution=80,
                         angle=theta,
                         fill=0,
                         edgecolor=ec,
                         linewidth=1.0)
            self.ellipse[i] = E
            ax.add_patch(E)

        ## And a circle of [kpc] in radius
        Xo = self.x_image[iclose] - (self.xo - self.dx)
        Yo = self.y_image[iclose] - (self.yo - self.dy)
        # Change the referece pixel to reflect jpg standards where the
        # origin is at (0,ny), is the upper left corner
        Yo = self.ny - Yo

        r_pixels = self.rdeg * 3600.0 / self.pixscale
        C = PCircle((Xo, Yo),
                    r_pixels,
                    resolution=80,
                    fill=0,
                    edgecolor="white",
                    linestyle='solid',
                    linewidth=0.5)
        ax.add_patch(C)

        # Get the area coverage
        self.area_in_circle(Xo, Yo, r_pixels)

        pylab.draw()
        #pylab.show()
        return

    # Get the area inside circle
    def area_in_circle(self, xo, yo, r_pixels):
        (ix, iy) = numpy.indices((self.nx, self.ny))
        d = numpy.sqrt((xo - ix)**2 + (yo - iy)**2)
        pix_in = numpy.where(d <= r_pixels)
        area_in = float(len(pix_in[0]))
        area_tot = math.pi * r_pixels**2
        self.area_fraction = area_in / area_tot
        return

        #####################################
        # Draw the BGCs candidates
        #####################################

    def ellipse_BCGs(self):

        ax = pylab.gca()
        # Delete all patches, reset ellipses before redraw
        del ax.patches[2:]
        self.ellipses = {}
        # construct the ellipses for each members
        for i in self.idx_BCG[0]:
            #ra = self.ra[i]
            #dec = self.dec[i]
            a = self.a_image[i]
            b = self.b_image[i]
            theta = self.theta[i]  # *math.pi/180.0

            # move to cropped reference frame
            xgal = self.x_image[i] - (self.xo - self.dx)
            ygal = self.y_image[i] - (self.yo - self.dy)
            # Change the referece pixel to reflect jpg standards where the
            # origin is at (0,ny), is the upper left corner
            ygal = self.ny - ygal

            ec = 'cyan'
            E = PEllipse((xgal, ygal), (a, b),
                         resolution=80,
                         angle=theta,
                         fill=0,
                         edgecolor=ec,
                         linewidth=1)
            self.ellipse[i] = E
            ax.add_patch(E)

        pylab.draw()
        pylab.show()
        return

    # Get the nearest object in catalog
    def get_nearest(self, x, y):
        distance = numpy.sqrt((self.x_image - x)**2 + (self.y_image - y)**2)
        self.iclose = numpy.argmin(distance)
        self.ID = self.id[self.iclose]
        return

    # Plot the redshift distribution of the members
    def redshift_members(self):
        pylab.figure(3)
        pylab.hist(self.z_ph[self.iRadius])
        pylab.draw()
        pylab.show()
        return

    # Print basic info to the screen
    def print_info(self):

        i = self.iclose
        gr = self.g_bpz[i] - self.r_bpz[i]
        ri = self.r_bpz[i] - self.i_bpz[i]

        ra = astrometry.dec2deg(self.ra[i] / 15.)
        dec = astrometry.dec2deg(self.dec[i])
        print("-------------------------------")
        print(" Object %s" % self.ID)
        print(" RA,DEC: %s %s" % (ra, dec))
        print(" X,Y:\t%s %s" % (self.x_image[i], self.y_image[i]))
        print(" T_B:\t%6.3f " % self.t_b[i])
        print(" Z_B:\t%6.3f (%.3f)" % (self.z_b[i], self.odds[i]))
        print(" Z_ML:\t%6.3f (%.3f)" % (self.z_ml[i], self.chi[i]))
        print(" Mi:\t%6.2f " % self.Mi[i])
        print(" Mr:\t%6.2f " % self.Mr[i])
        print(" p_BCG:\t%6.3f " % self.p[i])
        print(" g :\t%6.2f (%.3f) " % (self.g[i], self.g_err[i]))
        print(" r :\t%6.2f (%.3f) " % (self.r[i], self.r_err[i]))
        print(" i :\t%6.2f (%.3f) " % (self.i[i], self.i_err[i]))
        print(" z :\t%6.2f (%.3f) " % (self.z[i], self.z_err[i]))
        print(" g-r:\t%6.2f " % gr)
        print(" r-i:\t%6.2f " % ri)
        print(" Stellarity:\t%.2f " % self.class_star[i])
        print("--------------------------------")
        return

    def write_info(self):

        filename = os.path.join(self.datapath, self.ctile,
                                self.ctile + ".info")
        print("Will write info to %s" % filename)

        i = self.iBCG

        RA = astrometry.dec2deg(self.ra[i] / 15.)
        DEC = astrometry.dec2deg(self.dec[i])

        s = open(filename, "w")
        head = ("# %-18s %12s %12s %7s %7s %5s %10s %10s %8s %8s "
                "%8s %8s %8s %8s %8s %11s\n" % ('ID_BCG', 'RA', 'DEC', 'zBCG',
                                           'z_cl', 'Ngal', 'L_i', 'L_iBCG',
                                           'Mr', 'Mi', 'r', 'i', 'p_BCG',
                                           'R[kpc]', 'area[%]', 'Confidence'))
        format = ("%20s %12s %12s %7.3f %7.3f %5d %10.3e %10.3e %8.2f %8.2f "
                  "%8.2f %8.2f %8.3f %8.1f %8.2f %2d\n")
        vars = (self.ID_BCG, RA, DEC, self.z_ph[i], self.z_cl, self.Ngal,
                self.Lsum, self.Lr[i], self.Mr[i], self.Mi[i], self.r[i],
                self.i[i], self.p[i], self.radius, self.area_fraction,
                self.confidence)
        s.write(head)
        s.write(format % vars)
        s.close()
        return

    # Click for testing x,y recovery
    def click(self, event):

        xo = self.xo
        yo = self.yo

        ximage = event.xdata + (self.xo - self.dx)
        yimage = (self.ny - event.ydata) + (self.yo - self.dy)

        print('you clicked', event.xdata, self.ny - event.ydata)
        print('xo,yo      ', xo, yo)
        print('ximage,yimage', ximage, yimage)

        ra, dec = astrometry.xy2rd(ximage, yimage, self.fitsfile)
        RA = astrometry.dec2deg(ra / 15)
        DEC = astrometry.dec2deg(dec)
        print("ra,dec,filename", RA, DEC, self.fitsfile)
        return ximage, yimage

    def handle_ellipses(self, event, figure=1):

        i = self.iclose
        # Grab the plot
        pylab.figure(figure)
        ax = pylab.gca()

        # Delete all patches, before redraw
        del ax.patches[2:]

        # construct the ellipse for the current display
        a = self.a_image[i]
        b = self.b_image[i]
        theta = self.theta[i]  # *math.pi/180.0
        # move to cropped reference frame
        xgal = self.x_image[i] - (self.xo - self.dx)
        ygal = self.y_image[i] - (self.yo - self.dy)
        # Change the referece pixel to reflect jpg standards where the
        # origin is at (0,ny), is the upper left corner
        ygal = self.ny - ygal

        ec = 'white'
        E = PEllipse((xgal, ygal), (a, b),
                     resolution=80,
                     angle=theta,
                     fill=0,
                     edgecolor=ec,
                     linewidth=1.0)
        ax.add_patch(E)
        return

        # Erase current ID from list and don't draw anything else
        if event.key == 'd' or event.key == 'D':
            try:
                del self.ellipse[i]
            except FileNotFoundError:
                print('PROBLEM! Line 1173')
                #print("Ellipse for ID:%s is not defined" % ID)
            current_ell = None
        else:
            self.ellipse[i] = PEllipse((xgal, ygal), (a, b),
                                       resolution=100,
                                       angle=theta,
                                       fill=0,
                                       edgecolor="yellow")
            current_ell = PEllipse((xgal, ygal), (a, b),
                                   resolution=100,
                                   angle=theta,
                                   fill=0,
                                   edgecolor="white")

        # Draw all ellipses in list
        for id in list(self.ellipse.keys()):
            if id == i:
                continue
            ax.add_patch(self.ellipse[id])

        # And add current selection if defined
        if current_ell:
            ax.add_patch(current_ell)

        pylab.draw()
        # pylab.show()
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
        #dz1 = zo - x[i1]
        #dz2 = x[i2] - zo

        # And the 2-sigma (95.4%)
        i1 = numpy.where(Psum >= 0.023)[0][0]
        i2 = numpy.where(Psum > 0.977)[0][0]
        z1_95 = x[i1]
        z2_95 = x[i2]
        #dz1 = zo - x[i1]
        #dz2 = x[i2] - zo

        s = open(filename, "w")
        head = ("# %-18s " + "%8s " * 7 + "\n") % (
            'ID_BCG', 'z_cl', 'err', 'zBCG', 'z1', 'z2', 'z1', 'z2')
        format = ("%20s " + "%8.3f " * 7 + "\n")
        vars = (self.ID_BCG, self.z_cl, self.z_clerr, zo, z1_68, z2_68, z1_95,
                z2_95)
        s.write(head)
        s.write(format % vars)
        s.close()
        return

    # Handle the probability functions for the BCG and the core members
    def plot_probs(self):

        # The BCG
        iBCG = self.iBCG
        zBCG = self.z_ph[iBCG]
        dz1 = zBCG - self.z1[iBCG]
        dz2 = self.z2[iBCG] - zBCG

        # The probs
        zx = self.zx
        p_zBCG = self.p_z[iBCG]

        pylab.close(2)
        pylab.figure(2)
        pylab.plot(zx, p_zBCG, 'k-', alpha=0.8)
        for k in self.idc:
            if k == iBCG:
                pass
            else:
                pylab.plot(zx, self.p_z[k], 'r-', lw=0.5, alpha=0.8)
        pylab.xlabel("Redshift")
        pylab.ylabel("p(z)")

        n = 15.0
        pos = 0.25
        size = 12
        (x1, x2, y1, y2) = pylab.axis()
        dx = abs(x2 - x1)
        dy = abs(y2 - y1)
        xtxt = dx - dx * pos / 3.0
        ytxt = dy - 1.5 * dy / float(n)
        text = "zBCG = %.3f +%.3f/-%.3f)\n" % (zBCG, dz2, dz1)
        text = text + "z_cl = %.3f +/- %.3f" % (self.z_cl, self.z_clerr)
        pylab.text(xtxt,
                   ytxt,
                   text,
                   color='k',
                   size=size,
                   horizontalalignment='right')

        pylab.tight_layout()
        pylab.show()

        return

    # Write the members information out
    def write_members(self):

        filename = os.path.join(self.datapath, self.ctile,
                                self.ctile + ".members")
        m = open(filename, "w")
        print("Will write members to %s" % filename)

        head = ("# %-23s %15s %15s  %6s %6s %6s %6s %6s  %6s  %6s  "
               "%6s  %6s  %6s \n" % ("ID", "RA", "DEC", "ZB", "TB", "ZML",
                                     "TML", "g_mag", "g_err", "r_mag", "r_err",
                                     "i_mag", "i_err"))
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


##################################################################
# Read both kcorrection k(z) and evolution ev(z) from BC03 model
##################################################################
def KEfit(modelfile):

    import scipy
    import scipy.interpolate
    import tableio

    sout.write("# Getting K(z) and Ev(z) corrections from file:  %s\n" %
               modelfile)

    e = {}
    k = {}

    (z,
    k_g,
    k_r,
    k_i,
    k_z,
    e_g,
    e_r,
    e_i,
    e_z) = tableio.get_data(modelfile,
                            cols=(0, 12, 13, 14, 15, 17, 18, 19, 20))

    # K-only correction at each age SED,
    k['g'] = scipy.interpolate.interp1d(z, k_g)
    k['r'] = scipy.interpolate.interp1d(z, k_r)
    k['i'] = scipy.interpolate.interp1d(z, k_i)
    k['z'] = scipy.interpolate.interp1d(z, k_z)

    # Evolution term alone
    e['g'] = scipy.interpolate.interp1d(z, e_g)
    e['r'] = scipy.interpolate.interp1d(z, e_r)
    e['i'] = scipy.interpolate.interp1d(z, e_i)
    e['z'] = scipy.interpolate.interp1d(z, e_z)

    return k, e


############################################
# Read evolution ev(z) only from BC03 model
############################################
def evolfit(modelfile):

    import scipy
    import scipy.interpolate
    import tableio

    e = {}

    (z, e_g, e_r, e_i,
     e_z) = tableio.get_data(modelfile, cols=(0, 14, 15, 16, 17))

    e['g'] = scipy.interpolate.interp1d(z, e_g)
    e['r'] = scipy.interpolate.interp1d(z, e_r)
    e['i'] = scipy.interpolate.interp1d(z, e_i)
    e['z'] = scipy.interpolate.interp1d(z, e_z)

    return e


######################################
# BCG Probability function
# p = 0 for M dimmer than Mlimit
# p = 1 for M brighter than Mlimit
######################################
def p_BCG(M, Mlim, b=0.4, zp=0.5):
    x = M - Mlim
    return F_BCG(x, b, zp)


################################################################
# BCG priot aux function,
#################################################################
def F_BCG(x, b=0.4, zp=0.5):

    #print "will use zp:", zp
    #print "will use b:", b
    # Recenter at 50% (0.5) or at 68.2% (0.682)
    dx = math.log10(-math.log(zp)) / b
    u = x + dx
    phi = numpy.exp(-10**(b * u))
    return phi


#######################################################################
# Modified Schechter magnitude function from Postman et al (2001)
# uses alpha+2 rather than alpha+1 because of the extra 10^-0.4(m-m*)
# phi = (10^(-0.4(m-m*)))^(alpha+1) * exp[-10^(-0.4(m-m*))]
# PHI = phi*10^(-0.4(m-m*))
#######################################################################
def PHI(m, mstar, alpha):
    exp = numpy.exp
    a = 10**(-0.4 * (m - mstar))
    # Note (alpha+2) normally is just (alpha+1)
    phi = a**(alpha + 2) * exp(-a)
    return phi


###########################################
# Get m_star aparent mangnitude in i-band
###########################################
def mi_star(z, cosmo=(0.3, 0.7, 0.7)):

    # Set the cosmology
    c = cosmopy.set(cosmo)
    dlum = c.dlum(z)[0]

    Mb_star = -19.43 - 1.01 * z
    Mi_star = cosmology.reobs('El_Benitez2003',
                              m=Mb_star,
                              oldfilter="B_Johnson",
                              newfilter="i_MOSAICII")
    return Mi_star + 5.0 * math.log10(dlum) + 25

####################################################
# Fake an ellipse using an N-sided polygon
#####################################################
def PEllipse(xxx_todo_changeme,
             xxx_todo_changeme1,
             resolution=100,
             angle=0.0,
             **kwargs):
    (xo, yo) = xxx_todo_changeme
    (A, B) = xxx_todo_changeme1
    pi = math.pi
    cos = math.cos
    sin = math.sin
    angle = -angle * pi / 180.  # hack to make it work, angle=-angle

    t = 2 * pi / resolution * numpy.arange(resolution)
    xtmp = (A + 5) * numpy.cos(t)
    ytmp = (B + 5) * numpy.sin(t)

    x = xtmp * cos(angle) - ytmp * sin(angle) + xo
    y = xtmp * sin(angle) + ytmp * cos(angle) + yo
    return Polygon(list(zip(x, y)), **kwargs)

##############################
# A circle as a polygon too
###############################
def PCircle(xxx_todo_changeme2, radius, resolution=100, **kwargs):
    (xo, yo) = xxx_todo_changeme2
    pi = math.pi
    t = 2 * pi / resolution * numpy.arange(resolution)
    xtmp = radius * numpy.cos(t)
    ytmp = radius * numpy.sin(t)
    x = xtmp + xo
    y = ytmp + yo
    return Polygon(list(zip(x, y)), **kwargs)

#######################################
# make an array with power law growth
########################################
def mklogarray(x1, x2, n):

    if x1 > 0:
        i = numpy.indices((n + 1, ))[0] * 1.0
        x = x1 * (old_div(x2, x1))**(old_div(i, n))
        #dx = x1*( (x2/x1)**(i/n) - (x2/x1)**((i-1)/n))

    elif x1 == 0:
        i = numpy.indices((n, ))[0] * 1.0 + 1
        x = numpy.zeros((n + 1, ))
        #x[1:] = x2**(i/n)
        dx = (x2 + 1)**(old_div(i, n)) - (x2 + 1)**(old_div((i - 1), n))
        x[1:] = dx.cumsum()
    else:
        print("ERROR, x < 0")
        return

    return x


#################################################
# Make histogram using xbin, gives the same
# results as numpy.histogram
#################################################
def histo(x, xbin, center=None):

    n = len(xbin) - 1

    nbin = numpy.zeros(n).astype(int16)
    for i in range(n):
        if i == 0:
            nbin[i] = len(numpy.where(land(x >= xbin[i], x <= xbin[i + 1]))[
                0])
        else:
            nbin[i] = len(numpy.where(land(x > xbin[i], x <= xbin[i + 1]))[
                0])
    # Center and reduce to n-1
    if center:
        ix = numpy.indices(xbin.shape)[0]
        i1 = ix[:-1]
        i2 = ix[1:]
        dx = xbin[i2] - xbin[i1]
        xbin = xbin[:-1] + old_div(dx, 2.0)

    return nbin, xbin


################################################################
# Bin data in y(n) acoording to x(n) using bin spacing in xbin
###############################################################
def bin_data(x, y, xbin, center=None):

    n = len(xbin) - 1
    ybin = numpy.zeros(n).astype(float64)
    for i in range(n):
        if i == 0:
            idx = numpy.where(land(x >= xbin[i], x <= xbin[i + 1]))
        else:
            idx = numpy.where(land(x > xbin[i], x <= xbin[i + 1]))
        ybin[i] = y[idx].sum()
    # Center and reduce to n-1
    if center:
        ix = numpy.indices(xbin.shape)[0]
        i1 = ix[:-1]
        i2 = ix[1:]
        dx = xbin[i2] - xbin[i1]
        xbin = xbin[:-1] + old_div(dx, 2.0)

    return ybin, xbin

def Mass_calib(N200, L200, LBCG, z, h=1.0):

    # The best fit parameters
    if z < 0.23:
        M0_N = 1.27
        M0_L = 1.81
        alphaN = 1.20
        alphaL = 1.27
        gammaN = 0.71
        gammaL = 0.40
        aN = old_div(1.54, h**2)
        #aL     = 7.77/h**2 # bad value
        aL = old_div(0.61, h**2)
        bN = 0.41
        bL = 0.67

    else:
        M0_N = 1.57
        M0_L = 1.76
        alphaN = 1.12
        alphaL = 1.30
        gammaN = 0.34
        gammaL = 0.26
        aN = old_div(1.64, h**2)
        #aL     = 7.92/h**2 # bad value
        aL = old_div(0.58, h**2)
        bN = 0.43
        bL = 0.66

    L200 = L200 * (h**2) / 1.e10
    LBCG = LBCG * (h**2) / 1.e10

    LBCG_N = aN * N200**bN
    LBCG_L = aL * L200**bL

    M_N200 = (old_div(1.e14, h)) * M0_N * (
        (old_div(N200, 20.0))**alphaN) * (old_div(LBCG, LBCG_N))**gammaN
    M_L200 = (old_div(1.e14, h)) * M0_L * (
        (old_div(L200, 40.0))**alphaL) * (old_div(LBCG, LBCG_L))**gammaL
    M_LBCG = (old_div(1.e14, h)) * 1.07 * (old_div(LBCG, 5.0))**1.10

    return M_N200, M_L200, M_LBCG

def cmdline():

    from optparse import OptionParser

    # Read in the command line options
    USAGE = "usage:\t %prog <tilename> [options] \n"
    USAGE += "i.e.: %prog   --path ./PSZ2_G137.24+53.93/mosaic3/resampled/"
    parser = OptionParser(usage=USAGE)

    parser.add_option("--path", dest="path", default='./', help="Path to data")
    parser.add_option("--radius", dest="radius", default=1000.0,
                        help="Radius in kpc")
    parser.add_option("--zo", dest="zo", default=None, help="zo of cluster")
    parser.add_option("--dz", dest="dz", default=0.08, help="dz of shell")
    parser.add_option("--zuse", dest="zuse", default="ZB", help="use ZB or ML")
    parser.add_option("--dx", dest="dx", default=-1, help="cutout width")
    parser.add_option("--dy", dest="dy", default=-1, help="cutout heigth")
    parser.add_option("--RA", dest="RA", default=None,
                      help="Center Right Ascension - x pixel or SEX RA")
    parser.add_option("--DEC", dest="DEC", default=None,
                      help="Center Declination - y pixel or SEX DEC")

    # add a bit to figure out Mosaic1/mosaic3
    parser.add_option("--pixelscale", dest='pixelscale', default=-1,
        help='pixel scale of the instrument used MOS1:0.2666 or MOS3:0.25')

    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.error("Must supply at least one arguments required")

    return options, args

def main():

    opt, arg = cmdline()

    #print(opt)
    #print(arg)

    ctile = arg[0]
    radius = float(opt.radius)
    if opt.zo:
        zo = float(opt.zo)
    else:
        zo = None

    if float(opt.pixelscale) < 0:
        print('--pixelscale not defined.... assuming 0.25')
        pixelscale = 0.25
    else:
        pixelscale = float(opt.pixelscale)

    if opt.dx < 0 and opt.dy < 0:
        opt.dx = 5.1 * 60 / pixelscale
        opt.dy = 5.1 * 60 / pixelscale

    if not os.path.isfile('{}/{}irg.tiff'.format(opt.path, ctile)):
        print('RGB image does not exist -- probably only kband data')
        sys.exit()

    f = finder(ctile,
                maglim=26.0,
                pixscale=pixelscale,
                zlim=1.8,
                zo=zo,
                dz=float(opt.dz),
                radius=radius,
                cosmo=(0.3, 0.7, 0.7),
                #zuse="ZB", # Use ZB (Bayesian) or ML (Max Like)
                zuse=opt.zuse, # Use ZB (Bayesian) or ML (Max Like)
                outpath='plots',
                path=opt.path,
                evolfile="0_1gyr_hr_m62_salp_Ks.color",
                p_lim=0.4,
                verb='yes')

    f.jpg_read(dx=opt.dx, dy=opt.dy, RA=opt.RA, DEC=opt.DEC)
    f.jpg_display()

    return


main()
