import pylab
import time
import sys
import os

# get the utils from the parent directory
try:
    from cluster_utils import (PCircle, PEllipse)
except ImportError:
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from cluster_utils import (PCircle, PEllipse, color)

# fix large image error
import PIL
PIL.Image.MAX_IMAGE_PIXELS = None


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
        pylab.close('all')
        sys.exit()
        return

    try:
        # handle clicking outside the image
        ximage = event.xdata + (self.xo - self.dx)
        yimage = (self.ny - event.ydata) + (self.yo - self.dy)
    except TypeError:
        return

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
        self.background_map()
        #self.background()
        print("")
        print(color(f"\tNgal: {self.Ngal}", 36, 1))
        print(color(f"\tNgal_c: {self.Ngal_c}", 36, 1))
        print(color(f"\tz_cl: {self.z_cl:.3f} +/- {self.z_clerr:.3f}", 36, 1))
        print(f"\tL   : {self.Lsum:.3e} [Lsun]")
        print(f"\tMi  : {self.Mi[iclose]:6.2f}")
        print(f"\tMr  : {self.Mr[iclose]:6.2f}")

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
                                    fontsize=18,
                                    ha='right')
        self.conf_front = pylab.text(xo,
                                     self.ny - yo,
                                     text,
                                     color='white',
                                     fontsize=18,
                                     ha='right')
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
        self.txt_front = pylab.text(xo,
                                    self.ny - yo,
                                    text,
                                    color='white',
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
        self.jpg_read(dx=self.dx,
                      dy=self.dy,
                      RA=self.RA,
                      DEC=self.DEC,
                      toggle=self.toggled)
        self.jpg_display()

    # Print info
    self.click(event)
    self.print_info()
    self.handle_ellipses(event)

    pylab.draw()
    pylab.show()

    return


####################################
# The loop click-n-display routine
####################################
def jpg_display(self):

    # Measure time
    pylab.close('all')
    t0 = time.time()
    print("# Displaying... be patient", file=sys.stderr)
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
    print("\tDone in %s sec." % (time.time() - t0), file=sys.stderr)
    print("# Wrote Fig A, %s " % (self.figAname), file=sys.stderr)

    # save blank info.
    self.write_info(blank=True)
    # Draw the search zone
    self.mark_PSZ()
    self.draw_zone2()
    pylab.savefig(self.figAname,
                  transparent=True,
                  dpi=100,
                  bbox_inches='tight')
    # register this function with the event handler
    pylab.connect('button_press_event', self.get_object)
    pylab.show()
    return


#####################################
# Draw the zone where to select
######################################
def mark_PSZ(self):
    ax = pylab.gca()
    (nx, ny, nz) = self.jpg_region.shape
    ax.scatter(nx / 2,
               ny / 2,
               s=100,
               marker='*',
               facecolor='None',
               edgecolor='yellow')

    return


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

    r_pixels = 2 * self.rdeg * 3600.0 / self.pixscale
    C = PCircle((Xo, Yo),
                r_pixels,
                resolution=80,
                fill=0,
                edgecolor="yellow",
                linestyle='solid',
                linewidth=0.5)
    ax.add_patch(C)

    # Get the area coverage
    self.area_in_circle(Xo, Yo, r_pixels)

    pylab.draw()
    #pylab.show()
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


# Plot the redshift distribution of the members
def redshift_members(self):
    pylab.figure(3)
    pylab.hist(self.z_ph[self.iRadius])
    pylab.draw()
    pylab.show()
    return
