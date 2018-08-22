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
    z2 = 1.2
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

    nMasses = 100
    big_mass = []
    for j in range(nMasses):
        mass = numpy.ones_like(zarr) * 1e14
        for i in range(len(zbins) - 1):
            mask = (zbins[i] <= zarr) & (zarr < zbins[i + 1])

            mass[mask] *= float(data.loc[(zbins[i] <= data['REDSHIFT']) &
                            (data['REDSHIFT'] < zbins[i + 1]),
                                'MSZ'].sample()) * h
        big_mass.append(mass)

    mass = numpy.vstack(big_mass)

    z, dN, Ns = generate_Tinker(mass, z1, z2, dz)

    # Plot the cumulative distributions
    #f, ax = plt.subplots(1)

    f = plt.figure(figsize=(7 * (numpy.sqrt(5.) - 1.0) / 2.0, 7))
    ax = plt.subplot2grid((3, 1), (0, 0), rowspan=2)
    axs = plt.subplot2grid((3, 1), (2, 0))

    lws = [0.5, 1.5, 2.5]

    # plot the masses at the bottom
    for i in range(nMasses):
        axs.plot(z, big_mass[i] / 1e14, c='0.7', lw=lws[0], alpha=0.7)

    for i in range(nMasses):
        # Normalize to unity
        Ns[i] = Ns[i] / Ns[i][-1]
        #axs.plot(z, Ns[i], c='0.7', lw=lws[0], alpha=0.7, zorder=1)

    # add median and 68% lines
    ax.plot(z, numpy.median(Ns, axis=0), lw=lws[1], zorder=2, color='#e24a33')
    ax.plot(z, numpy.percentile(Ns, 50 + 34.1, axis=0), ls='--', lw=lws[1],
            zorder=2, color='#e24a33')
    ax.plot(z, numpy.percentile(Ns, 50 - 34.1, axis=0), ls='--', lw=lws[1],
            zorder=2, color='#e24a33')

    # bottom plot
    # make all the curves in one shot. The 68% are harder to make, so that's
    # why we are doing it all in one go, versus one at a time.
    # get unique masses in each set of bins
    umasses = [numpy.unique(numpy.unique(mass, axis=0)[:, i])
                            for i in range(zarr.size)]

    bounds = [numpy.percentile(umasses[i] / 1e14, [50 - 34.1, 50, 50 + 34.1])
              for i in range(zarr.size)]

    # convert to array
    bounds = numpy.array(bounds)

    axs.plot(z, bounds[:, 0], ls='--', lw=lws[1], color='#e24a33')
    axs.plot(z, bounds[:, 1], lw=lws[1], color='#e24a33')
    axs.plot(z, bounds[:, 2], ls='--', lw=lws[1], color='#e24a33')

    bins = numpy.linspace(0, 1.5, 50)

    n2, bins2, patches = ax.hist(df2.loc[df2['REDSHIFT'].notnull(),
                                         'REDSHIFT'], bins, histtype='step',
                                 cumulative=True, density=True, zorder=3,
                                 lw=lws[-1])

    # finish the plot
    ax.set_xlim(0, z2)
    ax.set_xticklabels([])
    axs.set_xlim(0, z2)
    ax.set_ylim(0, 1.01)
    axs.set_xlabel('Redshift (z)')
    ax.set_ylabel('N (<z)')
    axs.set_ylabel('M$_{500}$ (10$^{14} M_\odot$)')
    # The 50% and 90% levels
    ax.axhline(0.9, lw=0.8)
    ax.axhline(0.5, lw=0.8)
    ax.text(0.10, 0.92, "90%")
    ax.text(0.10, 0.52, "50%")
    ax.grid(alpha=0.25)
    # Make a custom legend
    indv_runs = plt.Line2D((0, 0), (0, 1), color='0.7', alpha=0.7, lw=lws[0])
    med = plt.Line2D((0, 0), (0, 1), color='#e24a33', lw=lws[1])
    quartile = plt.Line2D((0, 0), (0, 1), color='#e24a33', ls='--', lw=lws[1])
    psz_hist = plt.Line2D((0, 0), (0, 1), color='#348abd', lw=lws[-1])

    ax.legend([med, quartile, psz_hist], ['Median', '68%',
                                          'PSZ Confirmed (<z)'])
    axs.legend([indv_runs], ['M(z)'], loc='lower right')

    plt.show()


