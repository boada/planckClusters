from glob import glob
import pylab as pyl
from astropy.io import ascii
import numpy as np
import scipy
import os

filters = ['g', 'r', 'i', 'z', 'K']
#filters = ['g']

data = {}
data_sdss = {}
data_sdss_type = {}

dirs = [dirs for _, dirs, _ in os.walk('./')][0] # only want top level
cwd = os.getcwd()
for d in dirs[:50]:
    if 'PSZ' not in d:
        continue
    os.chdir(d)
    catalogs = glob('*_cal.cat')
    if len(catalogs) == 0:
        os.chdir(cwd)
        continue
    else:
        print(d, catalogs)

    # make the combined catalog Dictionary and filter information

    tilename = catalogs[0].split('_cal.cat')[0][:-1]
    complete_cat = ascii.read('{}_complete.catalog'.format(tilename))

    for c in catalogs:
        filter = c.split('_cal.cat')[0][-1]
        print(c, filter)
        tmp = ascii.read(c)

        try:
            photo_sdss = tmp['sdss_{}'.format(filter)]
            type_sdss = tmp['sdss_type']
        except KeyError:
            try:
                photo_sdss = tmp['2MASS_{}'.format(filter)]
                #type_sdss = tmp['sdss_type']
            except KeyError:
                continue
        try:
            photo = complete_cat['{}_MOSAICII_MAG_ISO'.format(filter)]
        except KeyError:
            try:
                photo = complete_cat['{}_KittPeak_MAG_ISO'.format(filter)]
            except KeyError:
                continue

        try:
            data_sdss['{}'.format(filter)] = np.append(data_sdss[
                                            '{}'.format(filter)], photo_sdss)
        except KeyError:
            data_sdss['{}'.format(filter)] = photo_sdss

        try:
            data_sdss_type['{}'.format(filter)] = np.append(data_sdss_type[
                                            '{}'.format(filter)], type_sdss)
        except KeyError:
            data_sdss_type['{}'.format(filter)] = type_sdss

        try:
            data['{}'.format(filter)] = np.append(data['{}'.format(filter)],
                                                  photo)
        except KeyError:
            data['{}'.format(filter)] = photo

    os.chdir(cwd)

fig, axes = pyl.subplots(2, 3)

#histogram definition
xyrange = [[14, 30], [14, 30]]  # data range
bins = [25, 25]  # number of bins
thresh = 5  # density threshold

for filt, ax in zip(filters, axes.ravel()):

    mask1 = data[filt] != 99
    mask2 = data[filt] != -99
    mask3 = data_sdss[filt] != 99
    mask = mask1 & mask2 & mask3

    xdat = data_sdss[filt][mask]
    ydat = data[filt][mask]

    ax.scatter(xdat, ydat, alpha=0.7, c='0.7')

    '''

    # histogram the data
    hh, locx, locy = scipy.histogram2d(xdat, ydat, range=xyrange, bins=bins)
    posx = np.digitize(xdat, locx)
    posy = np.digitize(ydat, locy)

    #select points within the histogram
    ind = (posx > 0) & (posx <= bins[0]) & (posy > 0) & (posy <= bins[1])
    hhsub = hh[posx[ind] - 1,
            posy[ind] - 1]  # values of the histogram where the points are
    xdat1 = xdat[ind][hhsub < thresh]  # low density points
    ydat1 = ydat[ind][hhsub < thresh]
    hh[hh < thresh] = np.nan  # fill the areas with low density by NaNs

    ax.scatter(xdat1, ydat1, alpha=0.7, c='0.7')
    ax.imshow(hh.T,
            cmap='binary',
            extent=np.array(xyrange).flatten(),
            interpolation='none')
    '''
