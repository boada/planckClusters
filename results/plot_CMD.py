import matplotlib.pyplot as plt
from matplotlib import cm
from get_results import loadClusters, loadMembers
from get_catalogs import loadCatalogs
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import ascii
from numpy import abs, linspace
import ezgal # BC03 model maker

def setup_models(zf):
    model = ezgal.model('bc03_exp_0.1_z_0.02_salp.model')
    # set cosmology
    model.set_cosmology(Om=0.3, Ol=0.7, h=0.7, w=-1)
    # set the model normalization to Dai et al 2009 (ApJ, 697, 506)
    model.set_normalization('ch1', 0.24, -25.06, vega=True)
    model.add_filter('g_MOSAICII.res', name='g_MOSAICII_MAG_ISO')
    model.add_filter('r_MOSAICII.res', name='r_MOSAICII_MAG_ISO')
    model.add_filter('SLOAN-SDSS.i.res', name='i_MOSAICII_MAG_ISO')
    model.add_filter('SLOAN-SDSS.z.res', name='z_MOSAICII_MAG_ISO')
    model.add_filter('K_KittPeak.res', name='K_KittPeak_MAG_ISO')
    model.set_zfs(zf)

    return model

def find_nearest(array, value):
    ''' returns the index and value of the nearest value in an array. '''
    idx = (abs(array - value)).argmin()
    return idx, array[idx]


# get the models for BC03
zf = [1, 1.5, 2, 2.5, 3]
models = setup_models(zf)

# make some colored lines for the plots
cm_subset = linspace(0.2, 0.8, len(zf))
colors = [cm.bone(x) for x in cm_subset]

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

        # what is the redshift of the cluster?
        z_cl = r_user.iloc[i]['zBCG_{}'.format(user)]

        # load the catalog data
        cat = loadCatalogs(user, c)

        # coords of BCG_
        bcg_coord = SkyCoord(ra=cat.iloc[bcg - 1]['RA'] * u.degree,
                         dec=cat.iloc[bcg - 1]['DEC'] * u.degree)
        cat_coord = SkyCoord(ra=cat.RA.values * u.degree,
                             dec=cat.DEC.values * u.degree)

        # find all of the close objects
        mask = bcg_coord.separation(cat_coord) <= 1.5 * u.arcmin
        close = cat[mask]

        # now it's time to make a figure
        #f, ax = plt.subplots(2, 2, figsize=(7 * (sqrt(5.) - 1.0) / 2.0, 7))
        f, ax = plt.subplots(2, 2, figsize=(7, 7))
        ax = ax.flatten()
        filters = ['g_MOSAICII_MAG_ISO', 'r_MOSAICII_MAG_ISO',
                   'i_MOSAICII_MAG_ISO', 'z_MOSAICII_MAG_ISO',
                   'K_KittPeak_MAG_ISO']

        # loop over the filters
        for j in range(len(filters) - 1):
            # make sure the data exists.
            if not filters[j] in close.columns.values:
                continue
            if not filters[j + 1] in close.columns.values:
                continue

            # plot the background galaxies
            ax[j].scatter(close[filters[j + 1]], close[filters[j]] -
                       close[filters[j + 1]], marker='.', # label='BG Galaxies',
                         c='k')

            # plot the cluster galaxies
            ax[j].scatter(cat.iloc[members.ID.astype('int') - 1][filters[j + 1]],
                    cat.iloc[members.ID.astype('int') - 1][filters[j]] -
                    cat.iloc[members.ID.astype('int') - 1][filters[j + 1]],
                    c='#a60628', )# label='Cluster Galaxies')

            # add the model colors
            for ii, color in zip(zf, colors):
                mag1 = models.get_observed_absolute_mags(ii,
                                                        filters=filters[j],
                                                        zs=z_cl)
                mag2 = models.get_observed_absolute_mags(ii,
                                                        filters=filters[j + 1],
                                                        zs=z_cl)
                ax[j].axhline(mag1 - mag2, color=color,
                                  label='zf={}'.format(ii))

            ax[j].set_xlabel(filters[j + 1][0], fontsize=16) # only the filter letter
            ax[j].set_ylabel('{} - {}'.format(filters[j][0], filters[j +
                                                                     1][0]),
                             fontsize=16)

            ax[j].set_ylim(-2, 3)

        # add extra labels
        ax[0].set_title('{} - {} - {}'.format(c,
                            r_user.iloc[i]['Confidence_{}'.format(user)], z_cl))

        ax[0].legend(loc='lower left', ncol=2)

        plt.tight_layout()
        plt.savefig('{}/{}/{}/{}_CMD.png'.format(user, c, c, c), bbox='tight')
        plt.close('all')

