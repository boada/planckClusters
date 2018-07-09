#!/usr/bin/env python3
from astropy.table import Table
import numpy
import matplotlib.pyplot as plt
from pyTinker import tinker

def generate_Tinker(mass, z1, z2, dz):

    h = 0.7

    # Initialize and generate the dn/dm
    t = tinker.counts(
        Om=0.3,
        OL=0.7,
        Ob=0.0456,
        H0=100.0,
        h0=0.70,
        ITRANS=5,
        sigma8=0.82,
        spectral_index=0.963)

    t.get_dndM(z1=z1, z2=z2, dz=dz, delta=500)

    nm = len(mass)
    nz = len(t.zx)

    dndz = numpy.zeros((nm, nz))
    Ns = numpy.zeros((nm, nz))

    for i, m in enumerate(mass):
        dn, N = t.get_dndz(m / h)
        dndz[i, :] = dn
        Ns[i, :] = N
    return t.zx, dndz, Ns


if __name__ == "__main__":

    h = 0.7
    z1 = 0
    z2 = 1.5
    dz = 0.025

    # build the mass array
    zarr = numpy.arange(z1, z2 + dz, dz)
    mass = numpy.ones_like(zarr) * 1e14

    ps2 = Table.read('../catalogs/PSZ2v1.fits')
    df2 = ps2.to_pandas()
    data = df2[['REDSHIFT', 'MSZ']]
    data['REDSHIFT'].replace(-1, numpy.nan, inplace=True)

    # redshift bins
    zbins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 3]

    big_mass = []
    for j in range(100):
        mass = numpy.ones_like(zarr) * 1e14
        for i in range(len(zbins) - 1):
            mask = (zbins[i] <= zarr) & (zarr < zbins[i + 1])

            mass[mask] *= float(data.loc[(zbins[i] <= data['REDSHIFT']) &
                            (data['REDSHIFT'] < zbins[i + 1]),
                                'MSZ'].sample()) * h
        big_mass.append(mass)

    mass = numpy.vstack(big_mass)

    #mass = numpy.array([3e14, 6e14, 9e14]) * h

    z, dN, Ns = generate_Tinker(mass, z1, z2, dz)

    # Plot the cumulative distributions
    f, ax = plt.subplots(1)

    #style = ['k-','k--','k:']
    lws = [0.5, 1.0, 2.0, 2.5]

    labels = [r'$3\times10^{14}M_{\odot}$', r'$6\times10^{14}M_{\odot}$',
              r'$9\times10^{14}M_{\odot}$']

    for i in range(100):
        # Normalize to unity
        Ns[i] = Ns[i] / Ns[i][-1]
        #Ns[i] = Ns[i] * 41000
        ax.plot(z, Ns[i], c='0.7', lw=lws[0], alpha=0.7)
        #ax.plot(z, Ns[i], 'k-', lw=lws[i], label=label)

    ax.plot(z, numpy.median(Ns, axis=0), 'r-', lw=lws[1])
    ax.plot(z, numpy.percentile(Ns, 50 + 34.1, axis=0), 'r--', lw=lws[1])
    ax.plot(z, numpy.percentile(Ns, 50 - 34.1, axis=0), 'r--', lw=lws[1])

    bins = numpy.linspace(0, 1.5, 50)

    n2, bins2, patches = ax.hist(df2.loc[df2['REDSHIFT'].notnull(),
                                         'REDSHIFT'], bins, histtype='step',
                                 cumulative=True, normed=True, zorder=2)

    ax.set_xlim(0, z2)
    ax.set_ylim(0, 1.01)
    ax.set_xlabel('Redshift')
    ax.set_ylabel('N(<z) ')
    # The 50% and 90% levels
    ax.axhline(0.9, lw=0.8)
    ax.axhline(0.5, lw=0.8)
    ax.text(0.10, 0.92, "90%")
    ax.text(0.10, 0.52, "50%")
    ax.grid(alpha=0.25)
    # Make the legend
    ax.legend()
    plt.show()


