import pandas as pd
from astropy.io import ascii
from numpy import nan

''' specify which round of cluster finding we want use. This can be changed by
calling

>>> loadClusters(round=2)

for example.

'''

def loadClusters(confirmed=False, round=1):

    users = ['boada', 'felipe', 'doze', 'jph']

    results = ['round{}/{}/{}_results.csv'.format(round, u, u) for u in users]

    # read in the files
    tables = [pd.read_csv(r) for r in results]
    #tables = [pd.read_csv(r) for r in results if os.path.isfile(r)]

    # clean off the extra columns
    cols = ['RA', 'DEC', 'Ngal', 'L_i', 'L_iBCG', 'Mr', 'Mi',
            'r', 'i', 'p_BCG', 'R[kpc]', 'area[%]']

    tables = [t.drop(cols, axis=1) for t in tables]

    # add the top level cluster names
    for t in tables:
        clusters = [c.rsplit('_', maxsplit=1)[:-1][0] for c in t['ID_BCG']]
        bcgs = [c.rsplit('_', maxsplit=1)[-1] for c in t['ID_BCG']]
        t['Cluster'] = clusters
        t['ID_BCG'] = bcgs

    # rename the columns
    for i, u in enumerate(users):
        tables[i] = tables[i].rename(columns={'ID_BCG': 'BCG_{}'.format(u),
                                             'Confidence': 'Conf_{}'.format(u),
                                             'zBCG': 'zBCG_{}'.format(u),
                                             'z_cl': 'z_cl_{}'.format(u),
                                             'z_clerr': 'z_clerr_{}'.format(u),
                                             'Ngal_c': 'Ngal_c_{}'.format(u)})

    # merge the first two together
    for i in range(len(users) - 1):
        if not i:
            df = pd.merge(tables[i], tables[i + 1], how='outer', on=['Cluster'])
        else:
            df = pd.merge(df, tables[i + 1], how='outer', on=['Cluster'])

    if confirmed:
        high_conf = ['PSZ1_G206.45+13.89',
                    'PSZ1_G224.82+13.62',
                    'PSZ2_G029.66-47.63',
                    'PSZ2_G043.44-41.27',
                    'PSZ2_G096.43-20.89',
                    'PSZ2_G120.76+44.14',
                    'PSZ2_G125.55+32.72',
                    'PSZ2_G137.24+53.93',
                    'PSZ2_G305.76+44.79',
                    'PSZ2_G107.83-45.45',
                    'PSZ2_G098.38+77.22',
                    # added later
                    'PSZ1_G084.62-15.86',
                    'PSZ2_G106.11+24.11',
                    'PSZ2_G173.76+22.92',
                    'PSZ2_G191.82-26.64']

        df = df.loc[df.Cluster.isin(high_conf)]

    return df

def loadMembers(user, cluster, round=1):

    data_dir = 'round{}/{}/{}/{}/'.format(round, user, cluster, cluster)

    table = ascii.read('{}/{}.members'.format(data_dir, cluster))

    # convert to pandas
    table = table.to_pandas()

    # clean things up a bit
    ids = [c.rsplit('_', maxsplit=1)[-1] for c in table['ID']]
    table['ID'] = ids
    table.replace(99, nan, inplace=True)

    return table


if __name__ == "__main__":
    print("Don't call this script. Import it into other scripts!")
