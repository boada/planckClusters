#!/usr/bin/env python3

import numpy
import matplotlib.pyplot as plt
from pyTinker import tinker

zmax = 1.5

def main():

    h = 0.7

    mass = numpy.array([3e14, 6e14, 9e14]) * h

    z, dN, Ns = generate_Tinker(mass)

    # Plot the cumulative distributions
    f, ax = plt.subplots(1)

    #style = ['k-','k--','k:']
    lws = [0.5, 1.0, 2.0]

    labels = [r'$3\times10^{14}M_{\odot}$', r'$6\times10^{14}M_{\odot}$',
              r'$9\times10^{14}M_{\odot}$']

    for i, label in enumerate(labels):
        # Normalize to unity
        Ns[i] = Ns[i] / Ns[i][-1]
        ax.plot(z, Ns[i], 'k-', lw=lws[i], label=label)

    ax.set_xlim(0, zmax)
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

    return

def generate_Tinker(mass):

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

    t.get_dndM(z1=0.05, z2=1.50, dz=0.025, delta=200)

    nm = len(mass)
    nz = len(t.zx)

    dndz = numpy.zeros((nm, nz))
    Ns = numpy.zeros((nm, nz))

    for i in range(nm):
        dn, N = t.get_dndz(mass[i] / h)
        dndz[i, :] = dn
        Ns[i, :] = N
    return t.zx, dndz, Ns


if __name__ == "__main__":
    main()
