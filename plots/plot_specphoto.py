from glob import glob
import pylab as pyl
from astropy.io import ascii
import numpy as np

# find all of the fields we have hunted
imgs = glob('./../cluster_search/round2/PSZ*/**/*A.png', recursive=True)
fields = [i.split('/')[-2] for i in imgs]

f = pyl.figure(figsize=(7 * (pyl.sqrt(5.) - 1.0) / 2.0, 7))
ax2 = pyl.subplot2grid((3, 1), (2, 0))
ax1 = pyl.subplot2grid((3, 1), (0, 0), rowspan=2)

gals = 0

for i, f in enumerate(fields):
    f = './../data/proc2/{}/{}.bpz'.format(f, f)

    cat = ascii.read(f)

    cat = cat.to_pandas()
    try:
        cat['Z_S'].replace(-99.0, np.nan, inplace=True)
    except KeyError:
        # no specz's?
        continue

    # remove non-specz's
    cat = cat.loc[cat['Z_S'].notnull()]

    # only ellipticals
    cat = cat.loc[cat['T_B'] < 2]

    # create a new column to store the "final" Redshift
    cat['redshift'] = np.nan

    # put the Z_MLs where T_ML and T_B are both type 1.000. Where they are not
    # both class 1.000 then we should use the Z_B redshift
    cat['redshift'] = np.where((cat.T_B == 1) & (cat.T_ML == 1), cat.Z_ML,
                               cat.Z_B)

    # originally just Z_B
    #cat['redshift'] = cat.Z_B

    #zphot = tmp['Z_B']
    #mask1 = tmp['Z_S'] != -99.
    #mask2 = tmp['T_B'] < 2
    #mask3 = tmp['ODDS'] > 0.7
    #lowz = tmp['Z_B'] <= 0.03

    zspec = cat.loc[cat['Z_S'].notnull(), 'Z_S']
    zphot = cat.loc[cat['Z_S'].notnull(), 'redshift']

    #mask = mask1 & mask2 & mask3 & ~lowz

    ax1.scatter(zspec, zphot, c='#e24a33', alpha=0.6,
                edgecolor='none')
    ax2.scatter(zspec, zspec - zphot, c='#e24a33', alpha=0.6,
               edgecolor='none')

    print(f, zspec.size)

    gals += zspec.size

    #lowz = tmp['Z_B'][mask] <= 0.03

    #badz = zspec[mask] - zphot[mask] > 0.2

    #print(tmp['ID'][mask][badz][:5])

    # make some big arrays so we can calculate some statistics
    try:
        specs_all = np.append(specs_all, zspec)
    except NameError:
        specs_all = zspec
    try:
        phots_all = np.append(phots_all, zphot)
    except NameError:
        phots_all = zphot

zb = pyl.Line2D((0, 1), (0, 0), color='#e24a33', marker='o', linestyle='',
                 label='Z_B')

ax1.legend([zb], ['z_b'], frameon=True, loc='upper left')

ax1.set_xlim(0, 1.0)
ax1.set_ylim(0, 1.0)
ax2.set_xlim(0, 1.0)
ax2.set_ylim(-0.5, 0.5)

ax1.plot([0, 1], [0, 1], c='k')
ax2.axhline(0)
ax2.set_xlabel('Redshift $z_{SPEC}$')
ax1.set_ylabel('Redshift $z_{BPZ}$')
ax2.set_ylabel('$z_{SPEC} - z_{BPZ}$')

ax1.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
ax1.set_yticks([0.25, 0.5, 0.75, 1.0])
ax1.set_xticklabels([])
ax2.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])

dz = abs(specs_all - phots_all) / (1 + specs_all)
sf = np.sqrt(np.mean(dz**2))
nmad = 1.48 * np.median(dz)
olf = np.where(abs(dz) > 5 * nmad)[0].size / dz.size

print(sf, nmad, olf)

pyl.show()
