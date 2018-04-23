from numpy import nan
from astropy.io import ascii

''' specify which round of cluster finding you want to use by calling the function
with

>>> loadCatalogs(round=1)

for example.

'''

def loadCatalogs(user, cluster, round=1):
    data_dir = 'round{}/{}/{}/'.format(round, user, cluster)
    table = ascii.read('{}/{}.color'.format(data_dir, cluster))

    # convert to pandas
    table = table.to_pandas()
    # clean up things
    table.replace(-99., nan, inplace=True)
    table.replace(99., nan, inplace=True)

    # clean off the extra columns
    cols = ['X_IMAGE', 'Y_IMAGE']
    table.drop(cols, axis=1, inplace=True)

    table.rename(columns={'NUMBER': 'ID'}, inplace=True)

    return table


if __name__ == "__main__":
    print("Don't call this script. Import it into other scripts!")
