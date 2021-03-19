#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys

def main():
    filename = sys.argv[1]
    x, y, _, _, _, _ = np.loadtxt(filename, unpack=True) # get columns instead of rows
    plt.scatter(x, y)
    plt.xlabel("Photon Energy (nm)")
    plt.ylabel("Emission intensity (a.u.)")
    plt.savefig("NV_center_spectrum.png")
    plt.clf()

    hc = 1239.8 # eV nm
    x_eV = hc/x # eV
    plt.scatter(x_eV, y)
    plt.xlabel("Photon wavelength (eV)")
    plt.ylabel("Emission intensity (a.u.)")
    plt.savefig("NV_center_spectrum_eV.png")


if __name__ == '__main__':
    main()
