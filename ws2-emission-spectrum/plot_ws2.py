#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys

def main():
    filename = sys.argv[1]
    x, y, _, _, _, _ = np.loadtxt(filename, unpack=True) # get columns instead of rows
    plt.scatter(x, y)
    plt.xlabel("Photon Energy (eV)")
    plt.ylabel("Intensity (Counts)")
    plt.savefig("WS2_monolayer_spectrum.png", dpi=300)
    plt.savefig("WS2_monolayer_spectrum.eps")
    plt.clf()

    hc = 1239.8 # eV nm
    x_nm = hc/x # nm
    plt.scatter(x_nm, y, color="blue")
    plt.xlabel("Photon wavelength (nm)")
    plt.ylabel("Intensity (Counts)")
    plt.savefig("WS2_monolayer_spectrum_nm.png", dpi=300)
    plt.savefig("WS2_monolayer_spectrum_nm.eps")


if __name__ == '__main__':
    main()
