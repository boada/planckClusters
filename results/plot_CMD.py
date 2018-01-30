import matplotlib.pyplot as plt
from get_results import loadClusters, loadMembers
from get_catalogs import loadCatalogs
from astropy.coordinates import SkyCoord
import astropy.units as u


results = loadClusters()

users = ['boada', 'felipe', 'doze']

for user in users:
    r_user = results[~results['Confidence_{}'.format(user)].isnull()]

    # now we loop over all of the results
    for i, c in r_user.Cluster.iteritems():
        bcg = r_user.iloc[i]['BCG_{}'.format(user)]
        bcg = int(bcg)
        # load the member data
        members = loadMembers(user, c)
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



