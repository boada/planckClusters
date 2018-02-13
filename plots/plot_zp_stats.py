from glob import glob
import pylab as pyl

zpts = glob('./**/photometry*.dat', recursive=True)
zpts.sort()
ax = pyl.subplot(111)

for i, zpt in enumerate(zpts):
    if not 'PS' in zpt:
        continue
    if '_g' in zpt:
        c = 'g'
    if '_r' in zpt:
        c = 'r'
    if '_i' in zpt:
        c = 'k'
    if '_z' in zpt:
        c = 'b'
    if '_K' in zpt:
        c = 'orange'
    tmp = pyl.genfromtxt(zpt, names=True, dtype=None)
    zp = tmp['ZP']
    zp_err = tmp['ZP_sig']
    if zp_err > 0.3:
        print(zpt)
    if not zp == 0.0 or zp_err == 0.0:
        cat = tmp['6'] # the data file reference by column number
        if 'sdss' in cat.tostring().decode().lower():
            marker = 'o'
        if 'panstarrs' in cat.tostring().decode().lower():
            marker = 'v'
        if 'apass' in cat.tostring().decode().lower():
            marker = 'P'
        if '2mass' in cat.tostring().decode().lower():
            marker = 'D'

        ax.scatter(zp_err, zp, c=c, marker=marker)
pyl.semilogx()
pyl.xlim(0.02, 4.5)
pyl.ylim(22, 33)
pyl.xlabel('zp error')
pyl.ylabel('zp')

# make a legend

g = pyl.Line2D((0, 1), (0, 0), color='g', marker='o', linestyle='', label='g')
r = pyl.Line2D((0, 1), (0, 0), color='r', marker='o', linestyle='', label='r')
i = pyl.Line2D((0, 1), (0, 0), color='k', marker='o', linestyle='', label='i')
z = pyl.Line2D((0, 1), (0, 0), color='b', marker='o', linestyle='', label='z')
K = pyl.Line2D((0, 1), (0, 0), color='orange', marker='o', linestyle='',
               label='k')
# surveys
sloan = pyl.Line2D((0, 1), (0, 0), color='0.6', marker='o', linestyle='', label='z')
ps = pyl.Line2D((0, 1), (0, 0), color='0.6', marker='v', linestyle='', label='z')
apa = pyl.Line2D((0, 1), (0, 0), color='0.6', marker='P', linestyle='', label='z')
twomass = pyl.Line2D((0, 1), (0, 0), color='0.6', marker='D', linestyle='', label='z')

pyl.legend([g, r, i, z, K, sloan, ps, apa, twomass], ['g', 'r', 'i', 'z', 'K', 'sdss', 'PS',
                                          'apass', '2mass'], ncol=2, frameon=True)

pyl.legend()

pyl.show()
