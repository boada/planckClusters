#!/usr/bin/env python

import os, sys

# Decided between numarray and numpy
if os.environ['NUMERIX'] == 'numpy':
    import numpy as numarray
    import numpy.random as nrandom
    float32 = numarray.float32
    float64 = numarray.float64
    int16 = numarray.int16
    nstr = numarray.char
    print("Will use numpy")
elif os.environ['NUMERIX'] == 'numarray':
    import numarray
    import numarray.random_array as nrandom
    float32 = numarray.Float32
    float64 = numarray.Float64
    int16 = numarray.Int16
    import numarray.strings as nstr
    print("Will use numarray :(")
else:
    sys.exit("Must define NUMERIX")
import glob
import math
from pyfits import getheader, getval
#import astrometry_hack as astrometry
import astrometry
import re
import time
import extras
#import bpz_mix
import scipy
#import scipy.misc.pilutil as pilutil
import scipy.misc as sci_misc
import pylab
import cosmopy
import tableio
import cosmology
import aux

land = numarray.logical_and
lge = numarray.greater_equal
lle = numarray.less_equal
lor = numarray.logical_or

sout = sys.stderr


class finder:

    def __init__(self, ctile,maglim=25.0,
                 pixscale = 0.2666,
                 zlim =1.8,
                 zo=None,
                 dz=0.05,
                 radius = 1000.0, # Radius in kpc
                 cosmo=(0.3,0.7,0.7),
                 zuse="ZB", # Use ZB (Bayesian) or ML (Max Like)
                 outpath='plots',
                 path = os.path.join(os.environ['HOME'],"SOAR-data/COMB"),
                 evolfile = "0_1gyr_hr_m62_salp.color",
                 p_lim = 0.4,
                 verb='yes'):

        # Check for environ vars
        self.home = os.environ['HOME']
        if not os.getenv('SOARpipe'):
            os.environ['SOARpipe'] = os.path.join(self.home, 'SOARpipe')
        self.SOARpipe = os.getenv('SOARpipe')

        self.zlim = zlim
        self.cosmo = cosmo
        self.evolfile = os.path.join(self.SOARpipe, "LIB/evol", evolfile)
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

        return

    #########################################
    # Read in the big catalog of photometry
    #########################################
    def read_cat(self):

        cols = (1,
                2,
                23,
                26,
                27,
                28,
                29,
                30,
                3,
                4,  #5,
                6,
                7,  #8,
                9,
                10,  #11,
                #12,13,#14,
                15,
                16,
                17,
                18,
                19,
                20,
                31,
                32,
                33,
                34,
                35,
                36)

        t1 = time.time()
        sout.write("# Reading cols:%s\n# Reading cats from: %s... \n" %
                   (cols, self.catsfile))

        (ra,
         dec,
         z_b,
         t_b,
         odds,
         z_ml,
         t_ml,
         chi,
         g,
         g_err,  #g_sn,
         r,
         r_err,  #r_sn,
         i,
         i_err,  #i_sn, 
         g_bpz,
         g_berr,
         r_bpz,
         r_berr,
         i_bpz,
         i_berr,
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
            t = t_ml
        elif self.zuse == "ZB":
            z_ph = z_b
            t = t_b

        i_lim = self.maglim
        odds_lim = 0.80
        star_lim = 0.80

        # Clean up according to BPZ
        #sout.write( "# Avoiding magnitudes -99 and 99 in BPZ \n") 
        #g_mask   = numarray.where(lor( g_bpz == 99,g_bpz == -99),0,1)
        #r_mask   = numarray.where(lor( r_bpz == 99,r_bpz == -99),0,1)
        #i_mask   = numarray.where(lor( i_bpz == 99,i_bpz == -99),0,1)
        #bpz_mask = g_mask*r_mask*i_mask

        # Clean up to avoid 99 values and very faint i_mag values
        #sout.write( "# Avoiding magnitudes 99 in MAG_AUTO \n") 
        #g_mask = numarray.where( g >= 99,    0 , 1)
        #r_mask = numarray.where( r >= 99,    0 , 1)
        #i_mask = numarray.where( i >= i_lim, 0 , 1)
        #sout.write( "# Avoiding magnitudes i > %s in MAG_AUTO \n" % i_lim)

        # Clean by class_star
        #sout.write( "# Avoiding CLASS_STAR > %s \n" % star_lim) 
        #mask_star = numarray.where( class_star > star_lim, 0 , 1)

        # Clean up by odds
        #sout.write( "# Avoiding ODDS < %s in BPZ \n" % odds_lim) 
        #odds_mask = numarray.where( odds > odds_lim, 1, 0)

        # Avoid z> zlim objects too.
        sout.write("# Avoiding objects with z > %s " % self.zlim)
        zp_mask = numarray.where(z_ph > self.zlim, 0, 1)

        # The final 'good' mask
        #mask_good = g_mask*r_mask*i_mask*zp_mask * odds_mask * mask_star
        mask_good = zp_mask
        idx = numarray.where(mask_good == 1)

        # Make ids a Char String in numarray
        self.id = nstr.array(id)[idx]

        # Only keep the 'good' one, avoid -99 and 99 values in BPZ mags
        self.ra = ra[idx]
        self.dec = dec[idx]
        self.z_b = z_b[idx]
        self.odds = odds[idx]
        self.chi = chi[idx]

        self.z_ml = z_ml[idx]
        self.t_ml = t_ml[idx]
        self.t_b = t_b[idx]
        self.t_ml = t_ml[idx]

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
        self.g_err = g_err[idx]
        self.r_err = r_err[idx]
        self.i_err = i_err[idx]

        self.g_bpz = g_bpz[idx]
        self.r_bpz = r_bpz[idx]
        self.i_bpz = i_bpz[idx]
        self.g_berr = g_berr[idx]
        self.r_berr = r_berr[idx]
        self.i_berr = i_berr[idx]

        self.class_star = class_star[idx]
        self.a_image = a_image[idx]
        self.b_image = b_image[idx]
        self.theta = theta[idx]
        self.x_image = x_image[idx]
        self.y_image = y_image[idx]

        # Color of selected galaxies
        self.gr = self.g_bpz - self.r_bpz
        self.ri = self.r_bpz - self.i_bpz

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
        sout.write("# Reading probs from :%s... " % self.probsfile)

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
                    zx = numarray.arange(z1, z2, dz)
                continue
            ID = fields[0]
            probs.append(numarray.asarray(list(map(float, fields[1:]))))

        # Transform the list into an N array
        p_z = numarray.asarray(probs)

        # select same galaxies as in catalogs we just read
        self.p_z = p_z[self.idx_cat][:]
        self.zx = zx
        sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t0))

        t1 = time.time()
        # Get the 1-sigma z1, z2 limits for each galaxy
        # Cumulatibe P(<z) function for each selected galaxy
        self.Psum = numarray.cumsum(self.p_z, axis=1)
        sout.write("# Getting +/- 1sigma (z1,z2) limits for each galaxy ")
        self.z1 = self.ra * 0.0
        self.z2 = self.ra * 0.0

        # One by one in the list
        for i in range(len(self.ra)):
            i1 = numarray.where(self.Psum[i, :] >= 0.159)[0][0]
            i2 = numarray.where(self.Psum[i, :] > 0.842)[0][0]
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
        self.DM = 25.0 + 5.0 * numarray.log10(self.dlum)

        t0 = time.time()
        # Get the absolute magnitudes, *** not including evolution ***, only Kcorr
        # We use a BPZ E's template for Kcorr
        #sout.write("# Computing absolute magnitudes interpolating Kcorr ")
        #k = Kcorr_fit(sed='El_Benitez2003')

        # Alternatibely we can get both the kcorr and the evol from
        # the *.color file from BC03 *.ised file
        sout.write(
            "# Computing absolute magnitudes interpolating konly from BC03 model \n")
        k, ev = KEfit(self.evolfile)

        self.Mg = self.g - self.DM - k['g'](self.z_ph)
        self.Mr = self.r - self.DM - k['r'](self.z_ph)
        self.Mi = self.i - self.DM - k['i'](self.z_ph)

        sout.write("# Computing evolution ev(z) for each galaxy ")
        self.ev_g = ev['g'](self.z_ph)
        self.ev_r = ev['r'](self.z_ph)
        self.ev_i = ev['i'](self.z_ph)

        # Also get the luminosities in Msun
        # taken from http://www.ucolick.org/~cnaw/sun.html
        self.Msun = {}
        self.Msun['g'] = 5.11
        self.Msun['r'] = 4.65
        self.Msun['i'] = 4.54

        # Mags k-corrected to z=0.25 as done in Reyes el al 2009
        #Mg = self.g - self.DM  -  k['g'](self.z_ph) + k['g'](0.25) 
        #Mr = self.r - self.DM  -  k['r'](self.z_ph) + k['r'](0.25) 
        #Mi = self.i - self.DM  -  k['i'](self.z_ph) + k['i'](0.25)

        # Mags k-corrected
        Mg = self.g - self.DM - k['g'](self.z_ph)
        Mr = self.r - self.DM - k['r'](self.z_ph)
        Mi = self.i - self.DM - k['i'](self.z_ph)

        self.Lg = 10.0**(-0.4 * (Mg - self.Msun['g']))
        self.Lr = 10.0**(-0.4 * (Mr - self.Msun['r']))
        self.Li = 10.0**(-0.4 * (Mi - self.Msun['i']))
        self.Lg_err = self.Lg * self.g_err / 1.0857
        self.Lr_err = self.Lr * self.r_err / 1.0857
        self.Li_err = self.Li * self.i_err / 1.0857

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
                0.1)  #+ self.DM_factor
            Mi_BCG_limit = Mi_limit + self.ev_i - self.evf['i'](
                0.1)  #+ self.DM_factor
            # Evaluate the BCG Probability function, we get the limit for each object
            self.p = p_BCG(self.Mr, Mr_BCG_limit)
            self.BCG_probs = True

            i_lim = 25.0
            star_lim = 0.5
            mask_p = numarray.where(self.p >= p_lim, 1, 0)
            mask_g = numarray.where(self.g < i_lim + 5, 1, 0)
            mask_r = numarray.where(self.r < i_lim + 5, 1, 0)
            mask_i = numarray.where(self.i < i_lim, 1, 0)
            mask_t = numarray.where(self.type < 2.0, 1, 0)

            # Avoid freakishly bright objects, 2.5 mags brighter than the M_BCG_limit
            mask_br = numarray.where(self.Mr > Mr_BCG_limit - 3.5, 1, 0)
            mask_bi = numarray.where(self.Mi > Mi_BCG_limit - 3.5, 1, 0)

            # Put a more strict cut in class_star for bcg candidates
            sout.write("# Avoiding CLASS_STAR > %s in BGCs\n" % star_lim)
            mask_star = numarray.where(self.class_star <= star_lim, 1, 0)

            # Construct the final mask now
            self.mask_BCG = mask_t * mask_g * mask_r * mask_i * mask_br * mask_bi * mask_p
            self.BCG_masked = True
            sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t0))

        # Select the candidates now
        idx = numarray.where(self.mask_BCG == 1)

        # And pass up to to class        
        # The index number
        self.idx_BCG = idx
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
        self.x_BCG = self.x_image[idx]
        self.x_BCG = self.x_image[idx]

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
    def jpg_read(self, dx=1200, dy=1200):

        # The fitsfile with the wcs information
        self.fitsfile = os.path.join(self.datapath, self.ctile,
                                     self.ctile + 'i.fits')
        self.jpgfile = os.path.join(self.datapath, self.ctile,
                                    self.ctile + '.tiff')
        t0 = time.time()
        print("Reading %s" % self.jpgfile, file=sys.stderr)
        #self.jpg_array  = pilutil.imread(self.jpgfile)
        self.jpg_array = sci_misc.imread(self.jpgfile)
        print("Done in %.3f sec." % (time.time() - t0), file=sys.stderr)

        # Get the shape of the array
        print(self.jpg_array.shape)
        (self.ny, self.nx, self.nz) = self.jpg_array.shape

        self.dx = self.nx / 2.0
        self.dy = self.ny / 2.0
        self.xcenter = self.nx / 2.0
        self.ycenter = self.ny / 2.0
        self.xo = self.nx / 2.0
        self.yo = self.ny / 2.0

        # Select the limits of the image to display
        #self.dx = dx
        #self.dy = dy
        #yo = self.yo 
        #xo = self.xo
        #x1 = int(xo - dx)
        #x2 = int(xo + dx)
        #y1 = int(yo - dy)
        #y2 = int(yo + dy)

        # Get the region to use for plotting
        #self.jpg_region = self.jpg_array[x1:x2, y1:y2, :]
        self.jpg_region = self.jpg_array

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

        sx = numarray.arange(s1, s2 + 0.05, ds)

        xtext = []
        xtick = []
        for s in sx:
            x = xo + s / scale  #+ ds/scale
            xtick.append(x)
            xtext.append("%.1f" % s)

        s1 = int((-dy / 2.0) * scale)
        s2 = int((+dy / 2.0) * scale)
        sy = numarray.arange(s1, s2 + 0.05, ds)

        ytext = []
        ytick = []
        for s in sy:
            y = yo + s / scale  #+ ds/scale
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
        t0 = time.time()
        print("Displaying... be patient", file=sys.stderr)
        # Display
        self.ax = pylab.figure(1, figsize=(12, 12))
        #self.ax = pylab.figure(1,figsize=(9,9))
        pylab.imshow(self.jpg_region)
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
        self.draw_zone(n=8)
        # register this function with the event handler
        pylab.connect('button_press_event', self.get_object)
        pylab.show()
        return

    #####################################
    # Draw the zone where to select 
    ######################################
    def draw_zone(self, n=8):
        # Draw the region to select from
        n = float(n)
        (nx, ny, nz) = self.jpg_region.shape
        dx = nx / n
        dy = ny / n
        xx = [dx, (n - 1) * dx, (n - 1) * dx, dx, dx]
        yy = [dy, dy, (n - 1) * dy, (n - 1) * dy, dy]
        pylab.plot(xx, yy, 'w--', linewidth=0.05)
        return

    #######################################
    # This is the loop routine interactive
    #######################################
    def get_object(self, event):

        if event.key == 'q' or event.key == 'Q':
            sys.exit()
            return

        # Remap to right positions
        ximage = event.xdata
        yimage = self.ny - event.ydata

        # Fins the closest one
        self.get_nearest(ximage, yimage)
        iclose = self.iclose

        # Plot BCG candidates
        if event.key == 'b':
            self.ellipse_BCGs()
            return

        # Plot around ID's redshift
        if event.key == 'c':
            radius = self.radius
            if self.zo:
                self.z_ph[iclose] = self.zo
            # Select members using distance and 3sigma clipping
            z_cl, z_err = self.select_members_radius(self.iclose,
                                                     radius=radius,
                                                     zo=self.zo)
            z_cl, z_err = self.select_members_radius(self.iclose,
                                                     radius=radius,
                                                     zo=z_cl)
            self.iBCG = self.iclose
            print("\t Ngal: %d" % self.Ngal)
            print("\t z_cl: %.3f +/- %.3f" % (self.z_cl, self.z_clerr))
            print("\t L   : %.3e [Lsun]" % self.Lsum)
            print("\t Mi  : %6.2f " % self.Mi[iclose])
            print("\t Mr  : %6.2f " % self.Mr[iclose])
            self.ellipse_members()
            return

        # Plot probs
        if event.key == 'p':
            self.plot_probs()
            return

        if event.key == 'v':
            iclose = self.iclose
            text = " z_cl = %.3f\n zBCG = %.3f\n Ngal = %d\n R  = %d[kpc]" % (
                self.z_cl, self.z_ph[iclose], self.Ngal, self.radius)
            xo = 80
            yo = 80
            pylab.text(xo + 2,
                       self.ny - yo + 2,
                       text,
                       color='black',
                       fontsize=18)
            pylab.text(xo, self.ny - yo, text, color='white', fontsize=18)
            pylab.draw()
            pylab.show()
            return

        if event.key == 'w':
            self.write_info()
            self.write_redshift()
            self.write_members()
            pylab.savefig(self.figBname,
                          transparent=True,
                          dpi=100,
                          bbox_inches='tight')
            sys.exit()
            return
        # Print info
        self.print_info()
        self.handle_ellipses(event)

        if event.key == 'p':
            self.print_selection()
            return

        pylab.draw()
        pylab.show()

        return

    ########################################################
    # Modified/updated from find_clusters_ext_auto.py
    # Select galaxies around ID galaxy un redshift range
    ########################################################
    def select_members(self, i, Mi_lim=-20.25):

        # Width of the redshift shell
        dz = self.dz

        t0 = time.time()
        sout.write("# Selecting Cluster members... Ngal, N200, R200 \n")
        # Get the relevant info for ith BCG
        zo = self.z_ph[i]
        ra0 = self.ra[i]
        dec0 = self.dec[i]
        Mi_BCG = self.Mi[i]
        DM = self.DM[i]
        ID_BCG = self.id[i]

        # 1 - Select in position around ra0,dec0
        # Define 1h^-1 Mpc radius in degress @ zo
        R1Mpc = 1000 * 1.0 / self.h  # in kpc
        r1Mpc = astrometry.kpc2arc(zo, R1Mpc,
                                   self.cosmo) / 3600.  # in degrees.
        rcore = r1Mpc / 2.0
        dist = astrometry.circle_distance(ra0,
                                          dec0,
                                          self.ra,
                                          self.dec,
                                          units='deg')
        mask_R1Mpc = numarray.where(dist <= r1Mpc, 1, 0)
        mask_rcore = numarray.where(dist <= rcore, 1, 0)
        arcmin2Mpc = astrometry.arc2kpc(
            zo, 60.0, self.cosmo) / 1000.0  # scale between arcmin and Mpc

        # 2 - Select in redshift
        z1 = zo - dz
        z2 = zo + dz
        mask_z = numarray.where(land(self.z_ph >= z1, self.z_ph <= z2), 1, 0)

        # 3 - Select in brightness
        Mi_lim_zo = Mi_lim + self.evf['i'](zo) - self.evf['i'](0.1)
        mask_L1 = numarray.where(self.Mi <= Mi_lim_zo, 1,
                                 0)  # Faint  cut > 0.4L*
        mask_L2 = numarray.where(self.Mi >= Mi_BCG, 1, 0)  # Bright cut < L_BCG

        # The final selection mask, position x redshift x Luminosity
        idx = numarray.where(mask_R1Mpc * mask_L1 * mask_L2 * mask_z == 1)
        idc = numarray.where(mask_rcore * mask_L1 * mask_L2 * mask_z == 1)

        # Shot versions handles
        gr = self.gr
        ri = self.ri

        # Some simple 3-sigma clipping defined using r< rcore
        Nsigma = 3.0
        loop = 1
        converge = False
        while not converge:
            # The conditions to apply
            c1 = numarray.abs(gr[idc] - gr[idc].mean(
            )) > Nsigma * numarray.std(gr[idc], ddof=1)
            c2 = numarray.abs(ri[idc] - ri[idc].mean(
            )) > Nsigma * numarray.std(ri[idc], ddof=1)
            iclip = numarray.where(lor(
                c1, c2))  # where any of the conditions fails
            if len(iclip[0]) > 0:
                idc = numarray.delete(idc, iclip[0])  # Removed failed ones
                converge = False
            else:
                converge = True
            loop = loop + 1

        sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t0))

        # Or we can make a new mask where the condition's are true
        c1 = numarray.abs(self.gr - gr[idc].mean()) > Nsigma * numarray.std(
            gr[idc], ddof=1)
        c2 = numarray.abs(self.ri - ri[idc].mean()) > Nsigma * numarray.std(
            ri[idc], ddof=1)
        mask_cm = numarray.where(lor(c1, c2), 0, 1)  # where condition fails
        iR1Mpc = numarray.where(
            mask_R1Mpc * mask_L1 * mask_L2 * mask_z * mask_cm == 1)
        Ngal = len(iR1Mpc[0])
        sout.write("# Total: %s objects selected in 1h^-1Mpc around %s\n" %
                   (Ngal, self.ID))

        #############################################################################
        # We'll skip 200 measurement as they depend on the corrected values of Ngal
        # Now let's get R200 and N200
        R200 = 0.156 * (Ngal**0.6) / self.h  # In Mpc
        r200 = astrometry.kpc2arc(zo, R200 * 1000.0,
                                  self.cosmo) / 3600.  # in degrees.
        mask_r200 = numarray.where(dist <= r200, 1, 0)
        i200 = numarray.where(
            mask_r200 * mask_L1 * mask_L2 * mask_z * mask_cm == 1)
        N200 = len(i200[0])
        self.i200 = i200
        self.N200 = N200
        self.R200 = R200
        self.r200 = r200
        self.L200 = self.Lr[i200].sum()
        ############################################################################

        # And the value for all galaxies up NxR1Mpc -- change if required.
        mask_R = numarray.where(dist <= 10 * r1Mpc, 1, 0)
        iR = numarray.where(mask_R * mask_L1 * mask_L2 * mask_z * mask_cm == 1)

        # Pass up
        self.iR = iR
        self.iR1Mpc = iR1Mpc
        self.N1Mpc = Ngal
        self.r1Mpc = r1Mpc  # in degress
        self.dist2BCG = dist
        self.arcmin2Mpc = arcmin2Mpc

        # Sort indices radially for galaxies < N*R1Mpc, will be used later
        i = numarray.argsort(self.dist2BCG[iR])
        self.ix_radial = iR[0][i]

        # We want to keep i200 and iR1Mpc to write out members.
        return Ngal, N200, R200  # iR1Mpc,i200

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
        DM = self.DM[i]
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
        mask_R = numarray.where(dist <= r, 1, 0)
        mask_rcore = numarray.where(dist <= rcore, 1, 0)
        arcmin2Mpc = astrometry.arc2kpc(
            zo, 60.0, self.cosmo) / 1000.0  # scale between arcmin and Mpc

        # 2 - Select in redshift
        z1 = zo - dz
        z2 = zo + dz
        mask_z = numarray.where(land(self.z_ph >= z1, self.z_ph <= z2), 1, 0)

        # 3 - Select in brightness
        Mi_lim_zo = Mi_lim + self.evf['i'](zo) - self.evf['i'](0.1)
        mask_L1 = numarray.where(self.Mi <= Mi_lim_zo, 1,
                                 0)  # Faint  cut > 0.4L*
        mask_L2 = numarray.where(self.Mi >= Mi_BCG, 1, 0)  # Bright cut < L_BCG

        # The final selection mask, position x redshift x Luminosity
        idx = numarray.where(mask_R * mask_L1 * mask_L2 * mask_z == 1)[0]
        idc = numarray.where(mask_rcore * mask_L1 * mask_L2 * mask_z == 1)[0]

        # Shot versions handles
        gr = self.gr
        ri = self.ri

        # Some simple 3-sigma clipping defined using r< rcore
        Nsigma = 3.0
        loop = 1
        converge = False
        while not converge:
            # The conditions to apply
            c1 = numarray.abs(gr[idc] - gr[idc].mean(
            )) > Nsigma * numarray.std(gr[idc], ddof=1)
            c2 = numarray.abs(ri[idc] - ri[idc].mean(
            )) > Nsigma * numarray.std(ri[idc], ddof=1)
            iclip = numarray.where(lor(c1, c2))[
                0]  # where any of the conditions fails
            if len(iclip) > 0:
                idc = numarray.delete(idc, iclip)  # Removed failed ones
                converge = False
            else:
                converge = True
            loop = loop + 1

            # Get the average redshift within the core:
            #z_cl    = self.z_ph[idc].mean()
            #z_clrms = self.z_ph[idc].std()

        print(idc)
        print(self.z_ph[idc])

        # Compute the weighted average and rms
        dz = 0.5 * numarray.abs(self.z2[idc] - self.z1[idc])
        # Fix zeros
        dz[dz == 0] = 1e-5
        z_cl, z_clrms = aux.statsw(self.z_ph[idc], weight=1.0 / dz)
        sout.write(" \t Done: %s\n" % extras.elapsed_time_str(t0))

        # Or we can make a new mask where the condition's are true
        c1 = numarray.abs(self.gr - gr[idc].mean()) > Nsigma * numarray.std(
            gr[idc], ddof=1)
        c2 = numarray.abs(self.ri - ri[idc].mean()) > Nsigma * numarray.std(
            ri[idc], ddof=1)
        mask_cm = numarray.where(lor(c1, c2), 0, 1)  # where condition fails
        iRadius = numarray.where(
            mask_R * mask_L1 * mask_L2 * mask_z * mask_cm == 1)
        Ngal = len(iRadius[0])
        sout.write("# Total: %s objects selected in %s [kpc] around %s\n" %
                   (Ngal, radius, self.ID))

        # Pass up
        self.iRadius = iRadius
        self.arcmin2Mpc = arcmin2Mpc
        self.Lsum = self.Lr[iRadius].sum()
        self.Ngal = Ngal
        self.z_cl = z_cl
        self.z_clerr = z_clrms
        self.rdeg = r  # in degress
        self.idc = idc  # galaxies used for mean redshift
        self.ID_BCG = ID_BCG

        return z_cl, z_clrms

    #####################################
    # Draw the BGCs and cluster members
    #####################################
    def ellipse_members(self, k=0):

        nx = self.nz
        ny = self.ny
        nz = self.nz

        iclose = self.iclose

        ax = pylab.gca()
        # Delete all patches, reset ellipses before redraw
        del ax.patches[:]
        self.ellipses = {}
        zo = self.z_ph[iclose]
        pylab.title("%s" % (self.ctile))

        # construct the ellipses for each members
        for i in self.iRadius[0]:
            ra = self.ra[i]
            dec = self.dec[i]
            a = self.a_image[i]
            b = self.b_image[i]
            theta = self.theta[i]  #*math.pi/180.0
            xo = self.x_image[i]
            yo = self.y_image[i]
            #(xo,yo) = astrometry.rd2xy(ra,dec,self.fitsfile)
            # Change the referece pixel to reflect jpg standards where the
            # origin is at (0,ny), is the upper left corner
            yo = ny - yo
            if i == self.iclose:
                ec = 'yellow'
            else:
                ec = 'red'
            E = PEllipse((xo, yo), (a, b),
                         resolution=80,
                         angle=theta,
                         fill=0,
                         edgecolor=ec,
                         linewidth=1.0)
            self.ellipse[i] = E
            ax.add_patch(E)

        Xo = self.x_image[iclose]
        Yo = ny - self.y_image[iclose]
        ## And a circle of [kpc] in radius
        r_pixels = self.rdeg * 3600.0 / self.pixscale
        C = PCircle((Xo, Yo),
                    r_pixels,
                    resolution=80,
                    fill=0,
                    edgecolor="white",
                    linestyle='dashed',
                    linewidth=0.5)
        ax.add_patch(C)

        # Get the area coverage
        self.area_in_circle(Xo, Yo, r_pixels)

        pylab.draw()
        pylab.show()
        return

    # Get the area inside circle
    def area_in_circle(self, xo, yo, r_pixels):
        (ix, iy) = numarray.indices((self.nx, self.ny))
        d = numarray.sqrt((xo - ix)**2 + (yo - iy)**2)
        pix_in = numarray.where(d <= r_pixels)
        area_in = float(len(pix_in[0]))
        area_tot = math.pi * r_pixels**2
        self.area_fraction = area_in / area_tot
        return

        #####################################
        # Draw the BGCs candidates
        #####################################
    def ellipse_BCGs(self):

        nx = self.nz
        ny = self.ny
        nz = self.nz
        ax = pylab.gca()
        # Delete all patches, reset ellipses before redraw
        del ax.patches[:]
        self.ellipses = {}
        # construct the ellipses for each members
        for i in self.idx_BCG[0]:
            ra = self.ra[i]
            dec = self.dec[i]
            a = self.a_image[i]
            b = self.b_image[i]
            theta = self.theta[i]  #*math.pi/180.0
            xo = self.x_image[i]
            yo = self.y_image[i]
            # Change the referece pixel to reflect jpg standards where the
            # origin is at (0,ny), is the upper left corner
            yo = ny - yo
            ec = 'red'
            E = PEllipse((xo, yo), (a, b),
                         resolution=80,
                         angle=theta,
                         fill=0,
                         edgecolor=ec,
                         linewidth=0.5)
            self.ellipse[i] = E
            ax.add_patch(E)

        pylab.draw()
        pylab.show()
        return

    # Get the nearest object in catalog
    def get_nearest(self, x, y):
        distance = numarray.sqrt((self.x_image - x)**2 + (self.y_image - y)**2)
        self.iclose = numarray.argmin(distance)
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
        print(" RA,DEC:%s %s" % (ra, dec))
        print(" T_B:   %6.3f " % self.t_b[i])
        print(" Z_B:   %6.3f (%.3f)" % (self.z_b[i], self.odds[i]))
        print(" Z_ML:  %6.3f (%.3f)" % (self.z_ml[i], self.chi[i]))
        print(" Mi:    %6.2f " % self.Mi[i])
        print(" Mr:    %6.2f " % self.Mr[i])
        print(" p_BCG: %6.3f " % self.p[i])
        print(" i :    %6.2f (%.3f) " % (self.i[i], self.i_err[i]))
        print(" r :    %6.2f (%.3f) " % (self.r[i], self.r_err[i]))
        print(" g :    %6.2f (%.3f) " % (self.g[i], self.g_err[i]))
        print(" g-r:   %6.2f " % gr)
        print(" r-i:   %6.2f " % ri)
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
        head = "# %-18s %12s %12s %7s %7s %5s %10s %10s %8s %8s %8s %8s %8s %8s %8s\n" % (
            'ID_BCG', 'RA', 'DEC', 'zBCG', 'z_cl', 'Ngal', 'L_i', 'L_iBCG',
            'Mr', 'Mi', 'r', 'i', 'p_BCG', 'R[kpc]', 'area[%]')
        format = "%20s %12s %12s %7.3f %7.3f %5d %10.3e %10.3e %8.2f %8.2f %8.2f %8.2f %8.3f %8.1f %8.2f\n"
        vars = (self.ID_BCG, RA, DEC, self.z_ph[i], self.z_cl, self.Ngal,
                self.Lsum, self.Lr[i], self.Mr[i], self.Mi[i], self.r[i],
                self.i[i], self.p[i], self.radius, self.area_fraction)
        s.write(head)
        s.write(format % vars)
        s.close()
        return

    # Click for testing x,y recovery
    def click(self, event):

        xo = self.xo
        yo = self.yo

        ximage = event.xdata
        yimage = self.ny - event.ydata

        print('you clicked', event.xdata, event.ydata)
        print('xo,yo      ', xo, yo)
        print('ximage,yimage', ximage, yimage)

        ra, dec = astrometry.xy2rd(ximage, yimage, self.fitsfile)
        RA = astrometry.dec2deg(ra / 15)
        DEC = astrometry.dec2deg(dec)
        print("ra,dec,filename", RA, DEC, self.fitsfile)
        return  #event.xdata,event.ydata

    def handle_ellipses(self, event, figure=1):

        i = self.iclose
        # Grab the plot
        pylab.figure(figure)
        ax = pylab.gca()

        # Delete all patches, before redraw
        del ax.patches[:]

        # construct the ellipse for the current display
        a = self.a_image[i]
        b = self.b_image[i]
        theta = self.theta[i]  #*math.pi/180.0
        xo = self.x_image[i]
        yo = self.y_image[i]
        # Change the referece pixel to reflect jpg standards where the
        # origin is at (0,ny), is the upper left corner
        yo = self.ny - yo
        ec = 'white'
        E = PEllipse((xo, yo), (a, b),
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
            except:
                print("Ellipse for ID:%s is not defined" % ID)
            current_ell = None
        else:
            self.ellipse[i] = PEllipse((xo, yo), (a, b),
                                       resolution=100,
                                       angle=theta,
                                       fill=0,
                                       edgecolor="yellow")
            current_ell = PEllipse((xo, yo), (a, b),
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
        pylab.show()
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
        x = numarray.arange(zx[2], zx[-3], ds)
        y = scipy.interpolate.splev(x, tck, der=0)
        y = numarray.where(y < 0, 0, y)

        # Get the 1 sigma error bars (68.2%)
        Psum = numarray.cumsum(y * ds / dz)
        i1 = numarray.where(Psum >= 0.159)[0][0]
        i2 = numarray.where(Psum > 0.841)[0][0]
        z1_68 = x[i1]
        z2_68 = x[i2]
        dz1 = zo - x[i1]
        dz2 = x[i2] - zo

        # And the 2-sigma (95.4%)
        i1 = numarray.where(Psum >= 0.023)[0][0]
        i2 = numarray.where(Psum > 0.977)[0][0]
        z1_95 = x[i1]
        z2_95 = x[i2]
        dz1 = zo - x[i1]
        dz2 = x[i2] - zo

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
        pylab.plot(zx, p_zBCG, 'k-')
        for k in self.idc:
            pylab.plot(zx, self.p_z[k], 'r-', lw=0.2)
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
        pylab.show()
        return

    # Write the members information out
    def write_members(self):

        filename = os.path.join(self.datapath, self.ctile,
                                self.ctile + ".members")
        m = open(filename, "w")
        print("Will write members to %s" % filename)

        head = "# %-23s %15s %15s  %6s %6s %6s %6s %6s  %6s  %6s  %6s  %6s  %6s \n" % (
            "ID", "RA", "DEC", "ZB", "TB", "ZML", "TML", "g_mag", "g_err",
            "r_mag", "r_err", "i_mag", "i_err")
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

    (z, k_g, k_r, k_i, k_z, e_g, e_r, e_i,
     e_z) = tableio.get_data(modelfile,
                             cols=(0, 10, 11, 12, 13, 14, 15, 16, 17))

    # K-only correction at each age SED,
    k['g'] = scipy.interpolate.interp1d(z, k_g)
    k['r'] = scipy.interpolate.interp1d(z, k_r)
    k['i'] = scipy.interpolate.interp1d(z, k_i)

    # Evolution term alone
    e['g'] = scipy.interpolate.interp1d(z, e_g)
    e['r'] = scipy.interpolate.interp1d(z, e_r)
    e['i'] = scipy.interpolate.interp1d(z, e_i)

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
    phi = numarray.exp(-10**(b * u))
    return phi


#######################################################################    
# Modified Schechter magnitude function from Postman et al (2001)
# uses alpha+2 rather than alpha+1 because of the extra 10^-0.4(m-m*)
# phi = (10^(-0.4(m-m*)))^(alpha+1) * exp[-10^(-0.4(m-m*))]
# PHI = phi*10^(-0.4(m-m*))
#######################################################################
def PHI(m, mstar, alpha):
    exp = numarray.exp
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
import matplotlib.patches
import math
Polygon = matplotlib.patches.Polygon


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

    t = 2 * pi / resolution * numarray.arange(resolution)
    xtmp = A * numarray.cos(t)
    ytmp = B * numarray.sin(t)

    x = xtmp * cos(angle) - ytmp * sin(angle) + xo
    y = xtmp * sin(angle) + ytmp * cos(angle) + yo
    return Polygon(list(zip(x, y)), **kwargs)


##############################
# A circle as a polygon too
###############################
def PCircle(xxx_todo_changeme2, radius, resolution=100, **kwargs):
    (xo, yo) = xxx_todo_changeme2
    pi = math.pi
    cos = math.cos
    sin = math.sin
    t = 2 * pi / resolution * numarray.arange(resolution)
    xtmp = radius * numarray.cos(t)
    ytmp = radius * numarray.sin(t)
    x = xtmp + xo
    y = ytmp + yo
    return Polygon(list(zip(x, y)), **kwargs)


def cmdline():

    from optparse import OptionParser

    # Read in the command line options
    USAGE = " usage:\t %prog <tilename> [options] \n i.e.: %prog ACT-J0102-4915  --path ~/SOAR-data/COMB"
    parser = OptionParser(usage=USAGE)

    parser.add_option(
        "--path",
        dest="path",
        default=os.path.join(os.environ['HOME'], "SOAR-data/COMB"),
        help="Path to data")

    parser.add_option("--radius",
                      dest="radius",
                      default=1000.0,
                      help="Radius in kpc")

    parser.add_option("--zo", dest="zo", default=None, help="zo of cluster")

    parser.add_option("--dz", dest="dz", default=0.08, help="dz of shell")

    parser.add_option("--zuse", dest="zuse", default="ZB", help="use ZB or ML")

    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.error("Must supply at least one arguments required")

    return options, args


def main():

    opt, arg = cmdline()

    print(opt)
    print(arg)

    ctile = arg[0]
    radius = float(opt.radius)
    if opt.zo:
        zo = float(opt.zo)
    else:
        zo = None

    f = finder(ctile,maglim=26.0,
               zlim =1.8,
               zo=zo,
               dz=float(opt.dz),
               radius=radius,
               cosmo=(0.3,0.7,0.7),
               #zuse="ZB", # Use ZB (Bayesian) or ML (Max Like)
               zuse=opt.zuse, # Use ZB (Bayesian) or ML (Max Like)
               outpath='plots',
               path = opt.path,
               evolfile = "0_1gyr_hr_m62_salp.color",
               p_lim = 0.4,
               verb='yes')

    f.jpg_read()
    f.jpg_display()

    return


main()
