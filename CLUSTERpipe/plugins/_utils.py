import numpy
import time
import sys
import pylab
import math
from utils import KEfit
try:
    import extras
    import astrometry
except ImportError:
    sys.path.append('/home/boada/Projects/planckClusters/MOSAICpipe/pipe_utils')
    import extras
    import astrometry

sout = sys.stderr


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

# Get the area inside circle
def area_in_circle(self, xo, yo, r_pixels):
    (ix, iy) = numpy.indices((self.nx, self.ny))
    d = numpy.sqrt((xo - ix)**2 + (yo - iy)**2)
    pix_in = numpy.where(d <= r_pixels)
    area_in = float(len(pix_in[0]))
    area_tot = math.pi * r_pixels**2
    self.area_fraction = area_in / area_tot
    return

# Get the nearest object in catalog
def get_nearest(self, x, y):
    distance = numpy.sqrt((self.x_image - x)**2 + (self.y_image - y)**2)
    self.iclose = numpy.argmin(distance)
    self.ID = self.id[self.iclose]
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
