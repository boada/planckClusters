import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import sys
from astropy.io import ascii
from sklearn.metrics import median_absolute_error, mean_squared_error

def main(tilename, filter=None):

    filters = ['g', 'r', 'i', 'z']

    if filter:
        cat = ascii.read('{}{}_cal.cat'.format(tilename, filter))
        mask_obs = (cat['sdss_{}'.format(filter)] != 99.0) & \
                                                    (cat['MAG_BEST'] != 99)
        mask_fwhm = (cat['FWHM_IMAGE'] != 0.0)
        mask_type = (cat['sdss_type'] == 3) # 3 = galaxies, 6 = stars
        mask = mask_obs & mask_fwhm & mask_type
        dm = cat['sdss_{}'.format(filter)][mask] - cat['MAG_BEST'][mask]
        bad_idx = np.where(abs(dm) > 1)[0]
        ids = cat['NUMBER'][mask][bad_idx]

    # create a figure and the subplots
    f, axes = plt.subplots(4, 1)

    for ax, f in zip(axes, filters):
        # read the data
        cat = ascii.read('{}{}_cal.cat'.format(tilename, f))
        # mask out bad data
        mask_obs = (cat['sdss_{}'.format(f)] != 99.0) & (cat['MAG_BEST'] != 99)
        # clean the data which lie outside the detection image
        mask_fwhm = (cat['FWHM_IMAGE'] != 0.0)

        mask_type = (cat['sdss_type'] == 3) # 3 = galaxies, 6 = stars

        mask = mask_obs & mask_fwhm & mask_type

        # make the magnitude difference
        dm = cat['sdss_{}'.format(f)][mask] - cat['MAG_BEST'][mask]

        ax.scatter(cat['sdss_{}'.format(f)][mask], dm, alpha=0.5)
        ax.axhline(0.0)

        dm_full = cat['sdss_{}'.format(f)] - cat['MAG_BEST']
        # only use the ids we want to use
        dm_ids = dm_full[ids - 1]
        cat_ids = cat[ids - 1]

        # now you have to mask out the crazy values
        mask_obs = (cat_ids['sdss_{}'.format(f)] != 99.0) & \
                                        (cat_ids['MAG_BEST'] != 99)
        # clean the data which lie outside the detection image
        mask_fwhm = (cat_ids['FWHM_IMAGE'] != 0.0)
        if filter == f:
            ax.scatter(cat_ids['sdss_{}'.format(f)][mask_obs & mask_fwhm],
                    dm_ids[mask_obs & mask_fwhm], alpha=0.5, c='#CF4451')
        else:
            ax.scatter(cat_ids['sdss_{}'.format(f)][mask_obs & mask_fwhm],
                    dm_ids[mask_obs & mask_fwhm], alpha=0.5, c='#467821')
        mask_spec = (cat_ids['sdss_specz'] != 99) & (cat_ids['sdss_specz'] > 0)

        ax.scatter(cat_ids['sdss_{}'.format(f)][mask_obs & mask_fwhm & mask_spec],
                    dm_ids[mask_obs & mask_fwhm & mask_spec], alpha=0.5,
                    marker='*', c='k')
        for i in range(len(cat_ids['sdss_objid'][mask_obs & mask_fwhm])):
            print(cat_ids['sdss_objid'][mask_obs & mask_fwhm][i],
              dm_ids[mask_obs & mask_fwhm][i])

        fit = stats.linregress(cat['sdss_{}'.format(f)][mask], dm)
        x = np.linspace(15, 30)
        ax.plot(x, fit.slope * x + fit.intercept, c='#e24a33')

        ax.set_ylabel('sdss$_{}$-me'.format(f))
        print(median_absolute_error(cat['sdss_{}'.format(f)][mask],
                                    cat['MAG_BEST'][mask]), end=' ')
        print(np.sqrt(mean_squared_error(cat['sdss_{}'.format(f)][mask],
                                    cat['MAG_BEST'][mask])))

        ax.set_xlim(15.5, 27)

    ax.set_xlabel('SDSS MAG')

    plt.suptitle(tilename, fontweight='medium')


if __name__ == "__main__":
    # call with the name of the tile that we are working with
    main(sys.argv[1], sys.argv[2])
