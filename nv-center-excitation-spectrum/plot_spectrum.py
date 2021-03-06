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


if __name__ == '__main__':
    main()
