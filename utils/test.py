"""
This file is primarily used to compute the MSE between two rasters.
It expects a .asc format and CAN NOT deal with raw tiff files.
You can either pass two files as parameters to the script or you can define their paths in the paths variable.

Usage: 
`test.py path/to/raster_1.asc path/to/raster_2.asc`
"""

import sys
import numpy as np

def main(argv):
    a, b, = [np.loadtxt(x, skiprows=6) for x in argv]

    print(f"MSE: {((a - b) ** 2).mean()}")

if __name__ == '__main__':
    paths = []

    if len(sys.argv) > 1:
        main(sys.argv[1:])
    else:
        main(paths)