import os
from astropy.io import fits
import subprocess
import shlex

def read_assoc(self):
    ''' Read the association files for a given tile'''

    self.filters = []
    self.filelist = []
    self.infiles = []

    self.exptime = {}
    self.airmass = {}

    # Make image list per filter
    self.files = {}
    self.files_weight = {}
    # Exptimes per filter
    self.exptimes = {}

    # The the DQmaks files
    self.dqmask = {}

    print("# Will read %s" % self.assocfile)

    with open(self.assocfile, 'r') as assocfile:
        lines = assocfile.readlines()
        subprocs = []
        for line in lines:
            if line[0] == "#":
                continue

            vals = line.split()
            fname = f'tiles/{os.path.basename(vals[0])}'

            # hack to deal with compressed files in the association
            if '.fz' in fname:
                fname = fname.rstrip('.fz')
                if not os.path.isfile(fname) and not self.noSWarp:
                    subprocs.append(subprocess.Popen(
                        shlex.split('funpack -v {}.fz'.format(fname))))

            # Figure out the dqmask
            if 'tu' in fname:
                i = 1
                while True:
                    nid = int(os.path.splitext(fname)[0][8:]) + i
                    ext = os.path.splitext(fname)[1]
                    pre = fname[0:8]
                    dqmask = "%s%s%s" % (pre, nid, ext)
                    if os.path.isfile('{}.fz'.format(dqmask)):
                        break
                    else:
                        print('Looking...')
                        i += 1
                if not os.path.isfile(dqmask) and not self.noSWarp:
                    subprocs.append(subprocess.Popen(
                        shlex.split('funpack -v {}.fz'.format(dqmask))))
            elif 'k4' in fname:
                dqmask = fname.replace('opi', 'opd')
                if not os.path.isfile(dqmask) and not self.noSWarp:
                    subprocs.append(subprocess.Popen(
                        shlex.split('funpack -v {}.fz'.format(dqmask))))

        [p.wait() for p in subprocs]
        [i.kill() for i in subprocs]

    with open(self.assocfile, 'r') as assocfile:
        lines = assocfile.readlines()
        for line in lines:
            if line[0] == "#":
                continue

            vals = line.split()
            fname = f'tiles/{os.path.basename(vals[0])}'

            # hack to deal with compressed files in the association
            if '.fz' in fname:
                fname = fname.rstrip('.fz')
                #if not os.path.isfile(fname) and not self.noSWarp:
                #os.system('funpack -v {}.fz'.format(fname))

            filtername = vals[1]
            self.exptime[fname] = float(vals[2])
            self.airmass[fname] = float(vals[3])

            # Figure out the dqmask
            if 'tu' in fname:
                i = 1
                while True:
                    nid = int(os.path.splitext(fname)[0][8:]) + i
                    ext = os.path.splitext(fname)[1]
                    pre = fname[0:8]
                    dqmask = "%s%s%s" % (pre, nid, ext)
                    if os.path.isfile('{}.fz'.format(dqmask)):
                        break
                    else:
                        print('Looking...')
                        i += 1
                if not os.path.isfile(dqmask) and not self.noSWarp:
                    pass
                    #os.system('funpack -v {}.fz'.format(dqmask))
                    #subprocs.append('funpack -v {}.fz'.format(dqmask))
            elif 'k4' in fname:
                dqmask = fname.replace('opi', 'opd')
                if not os.path.isfile(dqmask) and not self.noSWarp:
                    pass
                    #os.system('funpack -v {}.fz'.format(dqmask))
                    #subprocs.append('funpack -v {}.fz'.format(dqmask))

            # make a list of the file names
            self.filelist.append(fname)
            self.dqmask[fname] = dqmask
            # this is just a list of all files used
            self.infiles.append(fname)
            self.infiles.append(dqmask)

            # update the list of filters
            if filtername not in self.filters:
                self.filters.append(vals[1])
                self.files[filtername] = []
                self.files_weight[filtername] = []
                self.exptimes[filtername] = []

            # append the filename to the right filter
            self.files[filtername].append(fname)
            if not self.noSWarp:
                with fits.open(dqmask, mode='readonly') as inHDU:
                    if len(inHDU) > 2:
                        self.files_weight[filtername].append(
                            "%s.weight.fits" % os.path.splitext(fname)[0])
                    else:
                        self.files_weight[filtername].append(
                            "%s.weight.fits'[0]'" %
                            os.path.splitext(fname)[0])

            self.exptimes[filtername].append(self.exptime[fname])

    return

def get_filenames(self):
    ''' This script will get called if we choose not to swarp the images.
    This makes sure that all of the variables are defined.

    '''
    # this is to convert the hms/dms coordinates back into decimal degrees.
    # I don't really want to use this, but I am not rewriting things.
    from astLib import astCoords

    self.combima = {}
    self.combcat = {}
    self.weightima = {}
    self.outdir = os.path.join(self.outpath, self.tilename)

    self.center_dither(dryrun=True)
    self.xo = astCoords.hms2decimal(self.xo, ':')
    self.yo = astCoords.dms2decimal(self.yo, ':')

    for filter in self.filters:
        # Store the names
        self.combima[filter] = "%s%s" % (self.tilename, filter)
        self.weightima[filter] = "%s%s_weight.fits" % (self.tilename,
                                                       filter)
        # make the combined catalog filenames
        if os.path.isfile('{}{}.cat'.format(self.tilename, filter)):
            self.combcat[filter] = '{}{}.cat'.format(self.tilename, filter)
        elif os.path.isfile('{}{}_cal.cat'.format(self.tilename, filter)):
            self.combcat[filter] = '{}{}_cal.cat'.format(self.tilename,
                                                         filter)
    # The default output names
    self.colorCat = self.tilename + ".color"
    self.columnsFile = self.tilename + ".columns"
    return


def update_header_projection(self, proj='TAN'):
    for fname in self.infiles:
        print('# Updating {} '.format(fname), end='...')
        with fits.open(fname, mode='update') as hdulist:
            for i, hdu in enumerate(hdulist):
                if not i:
                    continue
                elif proj in hdu.header['CTYPE1']:
                    continue
                else:
                    hdu.header['CTYPE1'] = 'RA---{}'.format(proj)
                    hdu.header['CTYPE2'] = 'DEC--{}'.format(proj)
        print(' to {}'.format(proj))
    return

# Generate the mask from the weight
def generate_masks(self, filters, dryrun):

    self.mask = {}
    for filter in filters:

        self.mask[filter] = "%s%s_mask.fits" % (self.tilename, filter)
        weight = self.weightima[filter]
        mask = self.mask[filter]

        if dryrun:
            print("# Skipping Mask Generation...", file=sys.stderr)
            continue

        print("# Generating mask image from %s --> %s" % (
            weight, mask), file=sys.stderr)
        mask_from_weight(weight, mask, value=0)

    return

# Make the detection Image using the mask
def makeDetectionIma(self, filter='i'):

    mask = self.mask[filter]
    image = self.combima[filter] + ".fits"
    DetImage = "%s_detection.fits" % self.tilename

    print("# Generating Detection image for filter:%s" % filter,
          file=sys.stderr)
    print("# \t\t%s x %s --> %s" % (image, mask, DetImage),
          file=sys.stderr)

    # Read in the mask and the data
    m_data, m_hdr = fits.getdata(mask, header=True)
    i_data, i_hdr = fits.getdata(image, header=True)

    # Multiply to get the detection
    d_data = m_data * i_data

    # Make a fits file out of it
    newfits = fits.HDUList()
    hdu = fits.PrimaryHDU()
    hdu.data = d_data
    hdu.header = i_hdr

    # Remove old version of file before
    if os.path.isfile(DetImage):
        os.remove(DetImage)

    newfits.append(hdu)
    newfits.writeto(DetImage)
    newfits.close

    # Make it visible outside
    self.DetImage = DetImage

    return
