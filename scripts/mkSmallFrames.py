from astropy.io import fits
from astLib import astImages
from astLib import astCoords
from astLib import astCalc
from astLib.astWCS import WCS
from numpy import genfromtxt
import os

data_dir = '../data/proc2_small'


def cutouts(name, ra, dec, z=False, suffix=''):

    if isinstance(ra, (str, bytes)):
        ra = astCoords.hms2decimal(ra.decode(), ':')
    if isinstance(dec, (str, bytes)):
        dec = astCoords.dms2decimal(dec.decode(), ':')

    filters = ['g', 'r', 'i', 'z']
    for color in filters:
        try:
            wcs = WCS('{}/{}/{}{}.fits'.format(data_dir, name, name, color))
        except FileNotFoundError:
            print('{}/{}/{}{}.fits missing!'.format(data_dir, name, name,
                                                    color))
            return

        try:
            img_data = fits.getdata('{}/{}/{}{}.fits'.format(
                data_dir, name, name, color))
            wht_data = fits.getdata('{}/{}/{}{}_weight.fits'.format(
                data_dir, name, name, color))
        except TypeError:
            print('{}/{}/{}{}.fits problem!'.format(data_dir, name, name,
                                                    color))
            continue

        if z:
            print(z)
            s = 206265. / astCalc.da(z)  # arcseconds for 1 Mpc
            s = round(s) * 2  # want 1 Mpc radius
        else:
            s = 20 * 60  # arcseconds

        # clip and write the new image
        img_clipped = astImages.clipImageSectionWCS(
            img_data, wcs, ra, dec, clipSizeDeg=float(s / 3600))
        wht_clipped = astImages.clipImageSectionWCS(
            wht_data, wcs, ra, dec, clipSizeDeg=float(s / 3600))
        ra2, dec2 = wcs.wcs2pix(ra, dec)
        d_img = astImages.clipImageSectionPix(img_data, ra2, dec2,
                                              float(s / 0.25))
        d_wht = astImages.clipImageSectionPix(wht_data, ra2, dec2,
                                              float(s / 0.25))

        astImages.saveFITS('{}/{}/{}{}{}.fits'.format(
            data_dir, name, name, color, suffix), d_img, img_clipped['wcs'])
        astImages.saveFITS('{}/{}/{}{}{}.fits'.format(
            data_dir, name, name, color, '_weight'), d_wht, wht_clipped['wcs'])
    print('Done {}'.format(name))

    return


if __name__ == "__main__":
    # get catalog data
    data = genfromtxt(
        '../catalogs/PSZ2_unconfirmed_catalog - current.csv',
        #data = genfromtxt('../results/preconfirmed.csv',
        delimiter=',',
        names=True,
        dtype=None)

    for i, (ra, dec, name) in enumerate(
            zip(data['RA'], data['DEC'], data['Name'])):
        #for i, (ra, dec, name, z) in enumerate(zip(data['RA_EX'], data['DEC_EX'],
        #                                        data['Name'], data['z_EX'])):
        print(name, ra, dec, end='...')

        if os.path.isdir('{}/{}'.format(data_dir, name.decode())):
            if ra == b'':
                continue
            cutouts(name.decode(), ra, dec)
            #cutouts(name.decode(), ra, dec, z, suffix='_paper')
            #make_RGB(name.decode())
            #make_RGB(name.decode(), suffix='_paper')
        else:
            print('')
