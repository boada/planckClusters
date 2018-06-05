import aplpy
import matplotlib.pyplot as plt
from get_results import loadClusters, loadMembers
from astLib import astCalc
from astropy.table import Table


data_dir = '/home/boada/Projects/planckClusters/data/proc2_small/'

# the confirmed = True gets the 15 confirmed clusters
results = loadClusters(round=2, confirmed=True)

# load the master spreadsheet
t_ss = Table.read('../catalogs/PSZ2_unconfirmed_catalog - current.csv')
df_ss = t_ss.to_pandas()

observed = df_ss.loc[df_ss['MOSAIC Imaging'].notnull()]

for i, row in results.iterrows():
    mems = loadMembers('boada', row['Cluster'], round=2)
    try:
        ra = mems.loc[mems['ID'] == row['BCG_boada'], 'RA'].values[0]
    except IndexError:
        continue
    dec = mems.loc[mems['ID'] == row['BCG_boada'], 'DEC'].values[0]
    # convert to sexigesimal -- don't need
#    ra_sex = astCoords.decimal2hms(ra, ':')
#    dec_sex = astCoords.decimal2dms(dec, ':')

    # write it back into the main frame
    results.loc[i, 'RA BCG'] = ra
    results.loc[i, 'DEC BCG'] = dec

results = results.sort_values('Cluster')
# tack on all of the original SS info.
confirmed = results.merge(observed, left_on='Cluster', right_on='Name',
                            how='left')

# number of panels per figure
panels = 4
no_figures = len(results) // panels + 1

#figures = [plt.figure(i, figsize=(20, 20)) for i in range(no_figures)]


for i, cluster in enumerate(confirmed['Cluster']):

    fits = '{}/{}/{}i.fits'.format(data_dir, cluster, cluster)
    png = '{}/{}/{}irg.tiff'.format(data_dir, cluster, cluster)

    gc = aplpy.FITSFigure(fits)
    try:
        gc.show_rgb(png)
    except FileNotFoundError:
        gc.show_grayscale(stretch='arcsinh', pmin=1, pmax=98)
        gc.set_theme('publication')

    gc.set_tick_labels_format(xformat='hh:mm:ss', yformat='dd:mm')
    #gc.set_tick_labels_size('small')

    ###
    # move things around and draw the labels
    ###

    # recenter
    window = 206265. / astCalc.da(results.iloc[i]['zBCG_boada'])
    gc.recenter(results.iloc[i]['RA BCG'], results.iloc[i]['DEC BCG'],
                window / 3600)

    # add the circles
    gc.show_circles(confirmed.iloc[i]['RA'], confirmed.iloc[i]['DEC'], 2 / 60,
                    linestyle='--', edgecolor='#188487', facecolor='none')
    gc.show_circles(confirmed.iloc[i]['RA'], confirmed.iloc[i]['DEC'], 5 / 60,
                    linestyle='-', edgecolor='#188487', facecolor='none')

    # add the source position
    gc.show_markers(confirmed.iloc[i]['RA'], confirmed.iloc[i]['DEC'],
                    marker='*', s=150, layer='psz', color='#188487')

    # add extra info
    gc.add_scalebar(1 / 60, color='w', label="$1'$")
    text = ("%s\n"
            "z$_{phot}$ = %.3f\n" % (cluster, confirmed.iloc[i]['z_cl_boada']))

    ax = plt.gca()
    txt_front = plt.text(0.1, 0.97, text, ha='left', va='top',
                     transform=ax.transAxes, color='white', fontsize=20)

    gc.axis_labels.set_xtext('Right Ascension (J2000)')
    gc.axis_labels.set_ytext('Declination (J2000)')

    # now for the axis lables
    if i % 4 < 2:
        gc.axis_labels.hide_x()
    if i % 4 == 1 or i % 4 == 3:
        gc.axis_labels.hide_y()

    plt.tight_layout()
    plt.savefig(r'{}.pdf'.format(cluster), bbox='tight')
    plt.savefig(r'{}.png'.format(cluster), bbox='tight')


