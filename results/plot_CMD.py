import matplotlib.pyplot as plt
from get_results import loadClusters, loadMembers
from get_catalogs import loadCatalogs
from astropy.coordinates import SkyCoord
import astropy.units as u
from math import sqrt

results = loadClusters()

users = ['boada', 'felipe', 'doze']

for user in users:
    r_user = results[~results['Confidence_{}'.format(user)].isnull()]

    # now we loop over all of the results
    for i, (indx, c) in enumerate(r_user.Cluster.iteritems()):
        print(user, '--', c)
        bcg = r_user.iloc[i]['BCG_{}'.format(user)]
        bcg = int(bcg)
        # load the member data
        try:
            members = loadMembers(user, c)
        except FileNotFoundError:
            print('missing --', c)
            continue
        # load the catalog data
        cat = loadCatalogs(user, c)

        # coords of BCG_
        bcg_coord = SkyCoord(ra=cat.iloc[bcg - 1]['RA'] * u.degree,
                         dec=cat.iloc[bcg - 1]['DEC'] * u.degree)
        cat_coord = SkyCoord(ra=cat.RA.values * u.degree,
                             dec=cat.DEC.values * u.degree)

        # find all of the close objects
        mask = bcg_coord.separation(cat_coord) <= 5 * u.arcmin
        close = cat[mask]

        # now it's time to make a figure
        f, ax = plt.subplots(3, 1, figsize=(7 * (sqrt(5.) - 1.0) / 2.0, 7))
        filters = ['g_MOSAICII_MAG_ISO', 'r_MOSAICII_MAG_ISO',
                   'i_MOSAICII_MAG_ISO', 'z_MOSAICII_MAG_ISO']

        # loop over the filters
        for j in range(3):
            # make sure the data exists.
            if not filters[j] in close.columns.values:
                continue
            if not filters[j + 1] in close.columns.values:
                continue

            # plot the background galaxies
            ax[j].scatter(close[filters[j + 1]], close[filters[j]] -
                       close[filters[j + 1]], marker='.', label='BG Galaxies',
                         c='k')
            # plot the cluster galaxies
            ax[j].scatter(cat.iloc[members.ID.astype('int') - 1][filters[j + 1]],
                    cat.iloc[members.ID.astype('int') - 1][filters[j]] -
                    cat.iloc[members.ID.astype('int') - 1][filters[j + 1]],
                    c='#a60628', label='Cluster Galaxies')

            ax[j].set_xlabel(filters[j + 1][0]) # only the filter letter
            ax[j].set_ylabel('{} - {}'.format(filters[j][0], filters[j +
                                                                     1][0]))

        # add extra labels
        ax[0].set_title('{} - {}'.format(c,
                            r_user.iloc[i]['Confidence_{}'.format(user)]))

        plt.tight_layout()
        plt.savefig('{}/{}/{}/{}_CMD.png'.format(user, c, c, c), bbox='tight')
        plt.close('all')

