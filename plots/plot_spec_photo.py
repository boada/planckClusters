from glob import glob
import pylab as pyl
from astropy.io import ascii
import numpy as np
import scipy

files = glob('./**/*.bpz', recursive=True)
#files = glob('./*.bpz')
files.sort()
ax1 = pyl.subplot(121)
ax2 = pyl.subplot(122)

gals = 0

for i, f in enumerate(files):
    tmp = ascii.read(f)
    try:
        zspec = tmp['Z_S']
    except KeyError:
        continue
    zphot = tmp['Z_B']
    zml = tmp['Z_ML']

    mask1 = tmp['Z_S'] != -99.
    mask2 = tmp['T_B'] < 1.5
    mask3 = tmp['ODDS'] > 0.7
    lowz = tmp['Z_B'] <= 0.03

    mask = mask1 & mask2 & mask3 & ~lowz

    #    ax1.scatter(zspec[mask], zphot[mask], c='#a60628', alpha=0.6)
    #    ax2.scatter(zspec[mask], zml[mask], c='#188487', alpha=0.6)

    print(f, zspec[mask].size)

    gals += zspec[mask].size

    lowz = tmp['Z_B'][mask] <= 0.03

    print(tmp['ID'][mask][lowz][:5])

    # make some big arrays so we can calculate some statistics
    try:
        specs_all = np.append(specs_all, zspec[mask])
    except NameError:
        specs_all = zspec[mask]
    try:
        phots_all = np.append(phots_all, zphot[mask])
    except NameError:
        phots_all = zphot[mask]
    try:
        ml_all = np.append(ml_all, zml[mask])
    except NameError:
        ml_all = zml[mask]

#histogram definition
xyrange = [[0, 1], [0, 1]]  # data range
bins = [25, 25]  # number of bins
thresh = 5  # density threshold

xdat = specs_all

for ax, ydat, c in zip([ax1, ax2], [phots_all, ml_all],
                       ['#a60628', '#188487']):

    # histogram the data
    hh, locx, locy = scipy.histogram2d(xdat, ydat, range=xyrange, bins=bins)
    posx = np.digitize(xdat, locx)
    posy = np.digitize(ydat, locy)

    #select points within the histogram
    ind = (posx > 0) & (posx <= bins[0]) & (posy > 0) & (posy <= bins[1])
    hhsub = hh[posx[ind] - 1,
               posy[ind] - 1]  # values of the histogram where the points are
    xdat1 = xdat[ind][hhsub < thresh]  # low density points
    ydat1 = ydat[ind][hhsub < thresh]
    hh[hh < thresh] = np.nan  # fill the areas with low density by NaNs

    #    ax.scatter(xdat1, ydat1, alpha=0.7, c=c)
    ax.scatter(xdat, ydat, alpha=0.7, c=c, edgecolors='k')
#    ax.imshow(hh.T,
#            cmap='binary',
#            extent=np.array(xyrange).flatten(),
#            interpolation='none')

zml = pyl.Line2D(
    (0, 1), (0, 0), color='#188487', marker='o', linestyle='', label='Z_ML')
zb = pyl.Line2D(
    (0, 1), (0, 0), color='#a60628', marker='o', linestyle='', label='Z_B')

ax1.legend([zb], ['z_b'], frameon=True, loc='upper left')
ax2.legend([zml], ['z_ml'], frameon=True, loc='upper left')

ax1.set_xlim(0, 1)
ax1.set_ylim(0, 1)
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 1)

ax1.plot([0, 1], [0, 1], c='k')
ax2.plot([0, 1], [0, 1], c='k')

ax1.set_xlabel('z-spec')
ax2.set_xlabel('z-spec')
ax1.set_ylabel('z-bpz')

ax2.set_yticklabels([])
ax2.set_xticks([0.25, 0.5, 0.75, 1.0])

dz = abs(specs_all - phots_all) / (1 + specs_all)
sf = np.sqrt(np.mean(dz**2))
nmad = 1.48 * np.median(dz)
olf = np.where(abs(dz) > 5 * nmad)[0].size / dz.size

print(sf, nmad, olf)

dz = abs(specs_all - ml_all) / (1 + specs_all)
sf = np.sqrt(np.mean(dz**2))
nmad = 1.48 * np.median(dz)
olf = np.where(abs(dz) > 5 * nmad)[0].size / dz.size

print(sf, nmad, olf)

pyl.show()
