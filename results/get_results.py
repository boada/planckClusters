import pandas as pd
from astropy.io import ascii
from numpy import nan
import os

''' specify which round of cluster finding we want use. This can be changed by
calling

>>> loadClusters(round=2)

for example.

'''

def loadClusters(confirmed=False, round=1):

    users = ['boada', 'felipe', 'doze']

    results = ['round{}/boada/boada_results.csv'.format(round),
            'round{}/felipe/felipe_results.csv'.format(round),
            'round{}/doze/doze_results.csv'.format(round)]
    # read in the files
    tables = [pd.read_csv(r) for r in results if os.path.isfile(r)]

    # clean off the extra columns
    cols = ['RA', 'DEC', 'z_cl', 'Ngal', 'L_i', 'L_iBCG', 'Mr', 'Mi',
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
                                             'Confidence':
                                             'Confidence_{}'.format(u),
                                             'zBCG':
                                             'zBCG_{}'.format(u)})

    # merge the first two together
    df = pd.merge(tables[0], tables[1], how='outer', on=['Cluster'])
    df = pd.merge(df, tables[2], how='outer', on=['Cluster'])

    if confirmed:
        high_conf = ['PSZ2_G145.25+50.84',
             'PSZ2_G120.76+44.14',
             'PSZ2_G305.76+44.79',
             'PSZ2_G029.66-47.63',
             'PSZ2_G173.76+22.92',
             'PSZ1_G224.82+13.62',
             'PSZ2_G048.47+34.86',
             'PSZ2_G106.11+24.11',
             'PSZ1_G084.62-15.86',
             'PSZ2_G125.55+32.72',
             'PSZ2_G043.44-41.27',
             'PSZ2_G096.43-20.89']

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
