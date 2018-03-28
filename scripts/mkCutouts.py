from astropy.io import fits
from astLib import astImages
from astLib import astCoords
from astLib.astWCS import WCS
from numpy import genfromtxt
import os

data_dir = '../data/proc2'

def make_RGB(name, conf='stiff-common.conf'):

    _conf = os.path.join('/home', 'boada', 'Projects', 'planckClusters',
                         'MOSAICpipe', 'confs', conf)

    # input files
    red = '{}/{}/{}Red_cutout_paper.fits'.format(data_dir, name, name)
    green = '{}/{}/{}Green_cutout_paper.fits'.format(data_dir, name, name)
    blue = '{}/{}/{}Blue_cutout_paper.fits'.format(data_dir, name, name)

    # output file
    output = '{}/{}/{}_cutout_paper.tiff'.format(data_dir, name, name)

    # options
    opts = ['-MIN_LEVEL', '0.01',
            '-MAX_LEVEL', '0.993',
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

def cutouts(name, ra, dec):

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
            return

        try:
            img_data = fits.getdata('{}/{}/{}{}.fits'.format(data_dir, name, name,
                                                         color))
        except TypeError:
            continue
        # clip and write the new image
        img_clipped = astImages.clipImageSectionWCS(img_data, wcs, ra, dec,
                                                    8 / 60)
        ra2, dec2 = wcs.wcs2pix(ra, dec)
        d = astImages.clipImageSectionPix(img_data, ra2, dec2,
                                                    1920)
        astImages.saveFITS('{}/{}/{}{}_cutout_paper.fits'.format(data_dir, name,
                                                           name,
                                                           color),
                           #img_clipped['data'],
                           d,
                           img_clipped['wcs'])
    print('Done {}'.format(name))

    return


if __name__ == "__main__":
    # get catalog data
    data = genfromtxt('../catalogs/PSZ2_unconfirmed_catalog - matches.csv',
            delimiter=',', names=True, dtype=None)

    for i, (ra, dec, name) in enumerate(zip(data['RA_x'], data['DEC_x'],
                                            data['Name'])):
        print(name, ra, dec, end='...')

        if os.path.isdir('{}/{}'.format(data_dir, name.decode())):
            if ra == b'':
                continue
            cutouts(name.decode(), ra, dec)
            make_RGB(name.decode())
        else:
            print('')
