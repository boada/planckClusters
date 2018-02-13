from glob import glob
import pylab as pyl
from astropy.io import ascii
import numpy as np

files = glob('./**/*.bpz', recursive=True)
files.sort()

f = pyl.figure(1, figsize=(7 * (pyl.sqrt(5.) - 1.0) / 2.0, 7))
ax2 = pyl.subplot2grid((3, 1), (2, 0))
ax1 = pyl.subplot2grid((3, 1), (0, 0), rowspan=2)

gals = 0

for i, f in enumerate(files):
    tmp = ascii.read(f)
    try:
        zspec = tmp['Z_S']
    except KeyError:
        continue
    zphot = tmp['Z_B']

    mask1 = tmp['Z_S'] != -99.
    mask2 = tmp['T_B'] < 2
    mask3 = tmp['ODDS'] > 0.7
    lowz = tmp['Z_B'] <= 0.03

    mask = mask1 & mask2 & mask3 & ~lowz

    ax1.scatter(zspec[mask], zphot[mask], c='#a60628', alpha=0.6,
                edgecolor='none')
    ax2.scatter(zspec[mask], zspec[mask] - zphot[mask], c='#a60628', alpha=0.6,
               edgecolor='none')

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

zb = pyl.Line2D((0, 1), (0, 0), color='#a60628', marker='o', linestyle='',
                 label='Z_B')

ax1.legend([zb], ['z_b'], frameon=True, loc='upper left')

ax1.set_xlim(0, 1.0)
ax1.set_ylim(0, 1.0)
ax2.set_xlim(0, 1.0)
ax2.set_ylim(-1, 1.)

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
