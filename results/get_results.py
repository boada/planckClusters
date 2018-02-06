import pandas as pd
from astropy.io import ascii
from numpy import nan

def loadClusters():

    users = ['boada', 'felipe', 'doze']

    results = ['boada/boada_results.txt',
            'felipe/felipe_results.txt',
            'doze/doze_results.txt']
    # read in the files
    tables = [pd.read_table(r, delim_whitespace=True) for r in results]

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

    return df

def loadMembers(user, cluster):

    data_dir = '{}/{}/{}/'.format(user, cluster, cluster)

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
