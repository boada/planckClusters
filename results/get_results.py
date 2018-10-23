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

    users = ['boada', 'felipe', 'doze', 'jph']
    data_dir = os.path.dirname(os.path.abspath(__file__))

    tables = []
    i = -1
    for u in users:
        if not os.path.isfile('{}/round{}/{}/{}_results.csv'.format(
                                                    data_dir, round, u, u)):
            continue

        i += 1
        results = '{}/round{}/{}/{}_results.csv'.format(data_dir, round, u, u)

        table = pd.read_csv(results)

        # clean off the extra columns
        cols = ['RA', 'DEC', 'Ngal', 'L_i', 'L_iBCG', 'Mr', 'Mi',
                'r', 'i', 'p_BCG', 'R[kpc]', 'area[%]']

        table = table.drop(cols, axis=1)

        clusters = [c.rsplit('_', maxsplit=1)[:-1][0] for c in table['ID_BCG']]
        bcgs = [c.rsplit('_', maxsplit=1)[-1] for c in table['ID_BCG']]
        table['Cluster'] = clusters
        table['ID_BCG'] = bcgs

        tables.append(table)

        tables[i] = tables[i].rename(columns={'ID_BCG': 'BCG_{}'.format(u),
                                             'Confidence': 'Conf_{}'.format(u),
                                             'zBCG': 'zBCG_{}'.format(u),
                                             'z_cl': 'z_cl_{}'.format(u),
                                             'z_clerr': 'z_clerr_{}'.format(u),
                                             'Ngal_c': 'Ngal_c_{}'.format(u)})

    # start counting from 1
    i += 1

    # merge the tables together -- if there are at least two tables
    if i >= 2:
        for j in range(i):
            if not j:
                df = pd.merge(tables[i], tables[j + 1], how='outer',
                              on=['Cluster'])
            elif j + 2 <= i:
                df = pd.merge(df, tables[j + 2], how='outer', on=['Cluster'])
            else:
                break

    else:
        df = pd.DataFrame(tables[0])

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
