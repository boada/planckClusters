import os
import sys
from astropy.io.fits import getheader
from astropy.io import fits
import numpy as np
import subprocess
import shlex

# get the utils from the parent directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import check_exe


def get_astrometry(self, newfirm=False):
    ''' This calls the pp script to astrometrically correct the images. I'm
    breaking it appart from the rest of the script so I can control when
    things happen. It is mostly for testing.

    '''

    # we have to prepare the images first.
    print()
    check_exe('pp_prepare')
    check_exe('pp_register')

    # little patch to keep it from crashing
    try:
        os.mkdir('.diagnostics')
    except FileExistsError:
        pass

    # make the diagnostics file too
    open('diagnostics.html', 'a').close()

    try:
        if ':' in self.xo or ':' in self.yo:
            from astLib import astCoords
            self.xo = astCoords.hms2decimal(self.xo, ':')
            self.yo = astCoords.dms2decimal(self.yo, ':')
    except TypeError:
        pass

    subprocs = []
    for filter in self.filters:
        if not newfirm and filter == 'K':
            continue
        mosaic = '{}.fits'.format(self.combima[filter])
        if self.pixscale == 0.25:
            cmd = 'pp_prepare -ra {} -dec {} -telescope {} {}'.format(
                self.xo, self.yo, 'KPNO4MOS3', mosaic)
        elif self.pixscale == 0.2666:
            cmd = 'pp_prepare -ra {} -dec {} -telescope {} {}'.format(
                self.xo, self.yo, 'KPNO4MOS1', mosaic)
        else:
            cmd = 'pp_prepare -ra {} -dec {} {}'.format(
                self.xo, self.yo, mosaic)

        if not self.dryrun:
            print(cmd)
            subprocs.append(subprocess.Popen(shlex.split(cmd)))

    [p.wait(timeout=600) for p in subprocs]
    [i.kill() for i in subprocs]

    for filter in self.filters:
        if not newfirm and filter == 'K':
            continue
        mosaic = '{}.fits'.format(self.combima[filter])

        cmd = 'pp_register -snr 10 -minarea 12 {}'.format(mosaic)

        if not self.dryrun:
            subprocs.append(subprocess.Popen(shlex.split(cmd)))

    [p.wait(timeout=600) for p in subprocs]
    [i.kill() for i in subprocs]

    for filter in self.filters:
        if not newfirm and filter == 'K':
            continue
        mosaic = '{}.fits'.format(self.combima[filter])

        # correct the header information to make sure floats are floats and
        # not strings
        with fits.open(mosaic, mode='update') as f:
            header = f[0].header
            for key, val in list(header.items()):
                if 'CD1_' in key or 'CD2_' in key or \
                    'CRVAL' in key or 'CRPIX' in key or \
                        'EQUINOX' in key:
                    f[0].header[key] = float(val)
                if 'PV' in key:
                    f[0].header[key] = str(val)

    return

def get_zeropt(self, newfirm=False):
    from astropy.io import fits
    ''' This is going to call the photometrypipline script that lives in
    the projects directory. It should both astrometrically correct the
    mosaics and calculate the overall zeropoint of the mosaic. Currently,
    this has only been tested on the MOSAIC camera, and should NOT work
    with NEWFIRM, yet. It's on my TODO list.

    '''

    print()
    check_exe('pp_photometry')
    check_exe('pp_calibrate')
    check_exe('pp_distill')

    try:
        os.mkdir('.diagnostics')
    except FileExistsError:
        pass

    # make the diagnostics file too
    with open('diagnostics.html', 'a') as f:
        pass

    subprocs = []

    for filter in self.filters:
        if not newfirm and filter == 'K':
            continue
        if os.path.isfile('photometry_control_star_{}.dat'.format(filter)):
            print('# remove old photometry')
            os.remove('photometry_control_star_{}.dat'.format(filter))

        mosaic = '{}.fits'.format(self.combima[filter])
        # check to make sure it has the right header keywords
        with fits.open(mosaic, mode='update') as mos:
            keywords = ['TEL_KEYW', 'TELINSTR', 'MIDTIMJD', 'SECPIXY',
                        'SECPIXX', 'PHOT_K']

            if self.pixscale == 0.25:
                values = ['KPNO4MOS3', 'KPNO4m/MOSAIC', 0.0, 0.25, 0.25,
                              0.5]
            elif self.pixscale == 0.2666:
                values = ['KPNO4MOS1', 'KPNO4m/MOSAIC', 0.0, 0.2666, 0.2666,
                              0.5]
            else:
                values = ['KPNO4NEWF', 'KPNO4m/NEWFIRM', 0.0, 0.4, 0.4,
                              0.5]

            for kw, val in zip(keywords, values):
                # only update if they don't exist
                try:
                    mos[0].header[kw]
                except KeyError:
                    mos[0].header[kw] = val

            #####################################################
            ##### HACK TO MAKE THE ZEROPOINT ERRORS SMALLER #####
            #####################################################
            if filter != 'K':
                mos[0].header['GAIN'] = 1.

        if filter != 'K':
            cmd = 'pp_photometry -snr 10 -aprad 5.7 {}'.format(mosaic)
        elif self.pixscale == 0.25:
            cmd = 'pp_photometry -snr 10 -aprad 16 {}'.format(mosaic)
        elif self.pixscale == 0.2666:
            cmd = 'pp_photometry -snr 10 -aprad 15 {}'.format(mosaic)
        else:
            cmd = 'pp_photometry -snr 10 -aprad 10 {}'.format(mosaic)

        if not self.dryrun:
            subprocs.append(subprocess.Popen(shlex.split(cmd)))

    [p.wait() for p in subprocs]

    for filter in self.filters:
        if not newfirm and filter == 'K':
            continue
        mosaic = '{}.fits'.format(self.combima[filter])
        cmd = 'pp_calibrate {}'.format(mosaic)

        if not self.dryrun:
            if filter != 'K':
                os.system(cmd)
            else:
                cmd = 'pp_calibrate -cat 2MASS {}'.format(mosaic)
                os.system(cmd)
            #subprocs.append(subprocess.Popen(shlex.split(cmd)))

    #[p.wait() for p in subprocs]
    #[i.kill() for i in subprocs]

    for filter in self.filters:
        if not newfirm and filter == 'K':
            continue
        mosaic = '{}.fits'.format(self.combima[filter])
        # the positions file fails -- it's just used to make sure it doesn't
        # query JPL because our stuff isn't an asteroid.
        cmd = 'pp_distill {}'.format(mosaic)
        os.system(cmd) # gotta call it to make it work

        if os.path.isfile('photometry_control_star.dat'):
            os.rename('photometry_control_star.dat',
                        'photometry_control_star_{}.dat'.format(filter))
        else:
            print('Photometric calibration failed')

        # correct the header information to make sure floats are floats and
        # not strings
        with fits.open(mosaic, mode='update') as f:
            header = f[0].header
            for key, val in list(header.items()):
                if 'CD1_' in key or 'CD2_' in key or \
                    'CRVAL' in key or 'CRPIX' in key or \
                        'EQUINOX' in key:
                    f[0].header[key] = float(val)
                if 'PV' in key:
                    f[0].header[key] = str(val)

            # write the zeropt info to the header
            _tmp = np.genfromtxt(
                'photometry_control_star_{}.dat'.format(filter), dtype=None)

            f[0].header['MAGZERO'] = float(_tmp['f11'])
            f[0].header['MAGZSIG'] = float(_tmp['f12'])
            f[0].header['MAGZCAT'] = _tmp['f15'].flatten()[0].decode()

    return
