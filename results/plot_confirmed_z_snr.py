import pandas as pd
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from get_results import loadClusters


def bin_by(x, y, nbins=30, bins=None):
    """
    Divide the x axis into sections and return groups of y based on its x value
    Taken from : github.com/tommlogan
    """
    if bins is None:
        bins = np.linspace(x.min(), x.max(), nbins)

    bin_space = (bins[-1] - bins[0]) / (len(bins) - 1) / 2

    indicies = np.digitize(x, bins + bin_space)

    output = []
    for i in range(0, len(bins)):
        output.append(y[indicies == i])
    #
    # prepare a dataframe with cols: median; mean; 1up, 1dn, 2up, 2dn, 3up, 3dn
    df_names = [
        'mean', 'median', '5th', '95th', '10th', '90th', '25th', '75th'
    ]
    df = pd.DataFrame(columns=df_names)
    to_delete = []
    # for each bin, determine the std ranges
    for y_set in output:
        if y_set.size > 0:
            av = y_set.mean()
            intervals = np.percentile(y_set, q=[50, 5, 95, 10, 90, 25, 75])
            res = [av] + list(intervals)
            df = df.append(pd.DataFrame([res], columns=df_names))
        else:
            # just in case there are no elements in the bin
            to_delete.append(len(df) + 1 + len(to_delete))

    # add x values
    bins = np.delete(bins, to_delete)
    df['x'] = bins

    return df


# read all of the data
ps1 = Table.read('../catalogs/PSZ1v2.1.fits')
ps2 = Table.read('../catalogs/PSZ2v1.fits')
Bpaper = Table.read('../papers/1803.05764/Barrena_tbl3.csv')
# the confirmed = True gets the 12 confirmed clusters
results = loadClusters(round=2, confirmed=True)

# convert to pandas
df1 = ps1.to_pandas()
df2 = ps2.to_pandas()
df_paper = Bpaper.to_pandas()

# clean up strings -- not required
df1 = df1.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)
df2 = df2.applymap(lambda x: x.decode() if isinstance(x, bytes) else x)

# clean up underscores in bpaper and results
df_paper['Planck Name'] = df_paper['Planck Name'].str.replace('_', ' ')
results['Cluster'] = results['Cluster'].str.replace('_', ' ')
# drop empty values
df_paper.dropna(subset=['Planck Name', 'z_cl'], inplace=True)

# merge catalogs
df_m = df1.merge(
    df2,
    how='outer',
    left_on='INDEX',
    right_on='PSZ',
    suffixes=('_PSZ1', '_PSZ2'))

# merge in the Barrena paper
df_mp = df_m.merge(
    df_paper[['Planck Name', 'z_cl']],
    how='outer',
    left_on='NAME_PSZ1',
    right_on='Planck Name')

df_mp['redshift'] = df_mp['REDSHIFT_PSZ2'].fillna(df_mp['REDSHIFT_PSZ1'])
df_mp['redshift'] = df_mp['redshift'].fillna(df_mp['z_cl'])
df_mp['SNR'] = df_mp['SNR_PSZ2'].fillna(df_mp['SNR_PSZ1'])

# messy bit to add the SNR for our results. The problem is that the PSZ1
# sources are actually in the PSZ2 catalog. We've not updated the names. I
# might do that and I might not.
results_m = results.merge(df_mp, how='left', left_on='Cluster',
                        right_on='NAME_PSZ1')
results_m = results_m.merge(df_mp, how='left', left_on='Cluster',
                        right_on='NAME_PSZ2')
results_m['SNR'] = results_m['SNR_y'].fillna(results_m['SNR_x'])

mask = df_mp['redshift'] > 0
df = bin_by(df_mp.loc[mask, 'redshift'], df_mp.loc[mask, 'SNR'])
x = df_mp.loc[mask, 'redshift']
y = df_mp.loc[mask, 'SNR']

col1 = list(
    plt.style.library['fivethirtyeight']['axes.prop_cycle'])[0]['color']
col2 = list(
    plt.style.library['fivethirtyeight']['axes.prop_cycle'])[1]['color']
fig_width = 30  # 15.24#6.5 # if columns==1 else 6.9 # width in cm
font_size = 20

golden_mean = (np.sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
fig_height = fig_width * golden_mean  # height in cm,
params = {
    'axes.labelsize': font_size,  # fontsize for x and y labels (was 10)
    'font.size': font_size,  # was 10
    'legend.fontsize': font_size,  # was 10
    'xtick.labelsize': font_size,
    'lines.linewidth': 2,
    'figure.autolayout': True,
    'figure.figsize': [fig_width / 2.54, fig_height / 2.54]
}

with plt.style.context('fivethirtyeight'):
    mpl.rcParams.update(params)
    # plot the points
    plt.scatter(x, y, facecolors=col1, label='Known Clusters')
    # plot the 1st band
    plt.fill_between(df.x, df['5th'], df['95th'], alpha=0.1, color=col2,
                     label='1st band')
    # plot the 2nd band
    plt.fill_between(df.x, df['10th'], df['90th'], alpha=0.5, color=col2,
                     label='2nd band')
    # plot the 3rd band
    plt.fill_between(df.x, df['25th'], df['75th'], alpha=0.7, color=col2,
                     label='3rd band')
    # plt the line
    plt.plot(df.x, df['median'], color='1', alpha=0.7, linewidth=1)

plt.scatter(results_m['z_cl_boada'], results_m['SNR'], marker='*', s=150,
            facecolor='#188487', edgecolor='k', linewidths=0.5,
            label='This Work', zorder=10)

plt.xlabel('redshift')
plt.ylabel('SNR')
plt.show()
