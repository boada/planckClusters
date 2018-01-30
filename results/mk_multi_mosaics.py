import aplpy
from glob import glob
import pylab as pyl

files = glob('*.fits')
fig = pyl.figure(figsize=(13.2, 10.5))

for i, f in enumerate(files):
    cluster = f.split('_')[0]
    if cluster == 'c203p83+41p0':
        continue
    gc = aplpy.FITSFigure(f, figure=fig, north=True, subplot=(3, 3, i + 1))
    gc.show_grayscale(stretch='arcsinh', pmin=1, pmax=99.9)
    gc.set_tick_labels_format(xformat='hh:mm:ss', yformat='dd:mm:ss')
    gc.set_theme('publication')
    gc.set_tick_labels_size('small')

    # now for the axis lables
    if not i % 3 == 0:
        gc.axis_labels.hide_y()
    if i < 6:
        gc.axis_labels.hide_x()

    ax = fig.axes[-1]
    ax.set_title(cluster)

    data = pyl.genfromtxt('./../analysis_all/redshifts/'
            cluster.split('_')[0] + '_redshifts.csv', delimiter=',',
                          names=True, dtype=None)

    try:
        # filter out the specz's
        x = pyl.isnan(data['Specz'])
        # draw the specz's
        gc.show_markers(data['ra'][~x], data['dec'][~x], edgecolor='#ffbf00',
                        facecolor='none', marker='D', s=50)

    except ValueError:
        print 'no Speczs found'

    # draw observed but not redshifted
    x = data['Q'] == 2
    gc.show_markers(data['ra'][x], data['dec'][x], edgecolor='#a60628',
                    facecolor='none', marker='s', s=30)

    # draw redshifted
    x = (data['Q'] == 0) | (data['Q'] == 1)
    gc.show_markers(data['ra'][x], data['dec'][x], edgecolor='#188487',
                    facecolor='none', marker='o', s=30)

pyl.tight_layout()
