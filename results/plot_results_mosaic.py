import aplpy
import matplotlib.pyplot as plt
from numpy import sort, sqrt
from astropy.io import ascii

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

data_dir = '/home/boada/Projects/planckClusters/data/proc2'
results_dir = '/home/boada/Projects/planckClusters/results/boada'

fig = plt.figure(1, figsize=(10 * (sqrt(5.) - 1.0) / 2.0, 10))

for i, c in enumerate(sort(high_conf)):

    fits = '{}/{}/{}Blue_cutout.fits'.format(data_dir, c, c)
    png = '{}/pngs_cutout/{}_cutout.tiff'.format(data_dir, c)
    mems = '{}/{}/{}/{}.members'.format(results_dir, c, c, c)

    mems = ascii.read(mems)
    gc = aplpy.FITSFigure(fits, figure=fig, subplot=(4, 3, i + 1))
    try:
        gc.show_rgb(png)
    except FileNotFoundError:
        gc.show_grayscale(stretch='arcsinh', pmin=1, pmax=99.9)
        gc.set_theme('publication')

    gc.set_tick_labels_format(xformat='hh:mm:ss', yformat='dd:mm:ss')
    gc.set_tick_labels_size('small')

    # now for the axis lables
    if not i % 3 == 0:
        gc.axis_labels.hide_y()
    if i < 9:
        gc.axis_labels.hide_x()

    ax = fig.axes[-1]
    ax.set_title(c, fontsize='small')
