from astropy.io import fits
from astLib import astImages
from astLib import astCoords
from astLib import astCalc
from astLib.astWCS import WCS
from numpy import genfromtxt
import os

data_dir = '../data/proc2'

def make_RGB(name, conf='stiff-common.conf', suffix=''):

    _conf = os.path.join('/home', 'boada', 'Projects', 'planckClusters',
                         'MOSAICpipe', 'confs', conf)

    # input files
    red = '{}/{}/{}Red_cutout{}.fits'.format(data_dir, name, name, suffix)
    green = '{}/{}/{}Green_cutout{}.fits'.format(data_dir, name, name, suffix)
    blue = '{}/{}/{}Blue_cutout{}.fits'.format(data_dir, name, name, suffix)

    # output file
    output = '{}/{}/{}_cutout{}.tiff'.format(data_dir, name, name, suffix)

    # options
    opts = ['-MIN_LEVEL', '0.01',
            '-MAX_LEVEL', '0.993,0.993,0.993',
            '-MIN_TYPE', 'GREYLEVEL',
            '-MAX_TYPE', 'QUANTILE',
            '-COLOUR_SAT', '2.0',
            '-DESCRIPTION', "'{} RGB'".format(name),
            '-WRITE_XML', 'N',
            '-COPYRIGHT', "'Steven Boada'"]

    # build the command -- a space is always last
    cmd = 'stiff {} {} {} -c {} '.format(red, green, blue, _conf)
    cmd += '-OUTFILE_NAME {} '.format(output)
    # append the options
    for opt in opts:
        cmd += '{} '.format(opt)

    print(cmd)
    os.system(cmd)

    return

def cutouts(name, ra, dec, z=False, suffix=''):

    if isinstance(ra, (str, bytes)):
        ra = astCoords.hms2decimal(ra.decode(), ':')
    if isinstance(dec, (str, bytes)):
        dec = astCoords.dms2decimal(dec.decode(), ':')

    if os.path.isfile('{}/{}/{}Red.fits'.format(data_dir, name, name)):
        filters = ['Red', 'Green', 'Blue']
    else:
        filters = ['i', 'r', 'g']

    for color in filters:
        try:
            wcs = WCS('{}/{}/{}{}.fits'.format(data_dir, name, name, color))
        except FileNotFoundError:
            print('{}/{}/{}{}.fits missing!'.format(data_dir, name, name, color))
            return

        try:
            img_data = fits.getdata('{}/{}/{}{}.fits'.format(data_dir, name, name,
                                                         color))
        except TypeError:
            print('{}/{}/{}{}.fits problem!'.format(data_dir, name, name, color))
            continue

        if z:
            print(z)
            s = 206265. / astCalc.da(z) # arcseconds for 1 Mpc
            s = round(s) * 2 # want 1 Mpc radius
        else:
            s = 12 * 60 # arcseconds

        # clip and write the new image
        img_clipped = astImages.clipImageSectionWCS(img_data, wcs, ra, dec,
                                                   clipSizeDeg=float(s / 3600))
        ra2, dec2 = wcs.wcs2pix(ra, dec)
        d = astImages.clipImageSectionPix(img_data, ra2, dec2,
                                                    float(s / 0.25))

        astImages.saveFITS('{}/{}/{}{}_cutout{}.fits'.format(data_dir, name,
                                                           name,
                                                           color, suffix),
                           #img_clipped['data'],
                           d,
                           img_clipped['wcs'])
    print('Done {}'.format(name))

    return


if __name__ == "__main__":
    # get catalog data
    data = genfromtxt('../catalogs/PSZ2_unconfirmed_catalog - current.csv',
    #data = genfromtxt('../results/preconfirmed.csv',
            delimiter=',', names=True, dtype=None)

    for i, (ra, dec, name) in enumerate(zip(data['RA'], data['DEC'],
                                            data['Name'])):
        #for i, (ra, dec, name, z) in enumerate(zip(data['RA_EX'], data['DEC_EX'],
        #                                        data['Name'], data['z_EX'])):
        print(name, ra, dec, end='...')

        if os.path.isdir('{}/{}'.format(data_dir, name.decode())):
            if ra == b'':
                continue
            cutouts(name.decode(), ra, dec)
            #cutouts(name.decode(), ra, dec, z, suffix='_paper')
            make_RGB(name.decode())
            #make_RGB(name.decode(), suffix='_paper')
        else:
            print('')
