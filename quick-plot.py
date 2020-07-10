#!/usr/bin/env python
""" A script that uses xarray, numpy, matplotlib, and cmocean to draw a contour plot
"""

import os
import xarray as xr
import numpy as np

# --------------------------------------------- #

def _floor_1sigfig(x):
    """ Round a number down to it's first significant digit

        Examples
        ----
        _floor_1sigfig(45.2) = 40
        _floor_1sigfig(-45.2) = -50
        _floor_1sigfig(0.0223) = 0.02
    """
    if x == 0:
        return 0
    sgnx = np.sign(x)
    absx = sgnx*x
    pwr_of_10 = 10**np.floor(np.log10(absx))
    return pwr_of_10 * np.floor(x/pwr_of_10)

# --------------------------------------------- #

def _ceil_1sigfig(x):
    """ Round a number up to it's first significant digit

        Examples
        ----
        _ceil_1sigfig(45.2) = 50
        _ceil_1sigfig(-45.2) = -40
        _ceil_1sigfig(0.0223) = 0.03
    """
    if x == 0:
        return 0
    sgnx = np.sign(x)
    absx = sgnx*x
    pwr_of_10 = 10**np.floor(np.log10(absx))
    return pwr_of_10 * np.ceil(x/pwr_of_10)

# --------------------------------------------- #

def _make_plot(da, var, contour_levels):
    """ Generate the requested plot
    """
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import cmocean

    # valid strings are ['GTK3Agg', 'GTK3Cairo', 'MacOSX', 'nbAgg', 'Qt4Agg', 'Qt4Cairo', 'Qt5Agg', 'Qt5Cairo', 'TkAgg', 'TkCairo', 'WebAgg', 'WX', 'WXAgg', 'WXCairo', 'agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template']
    matplotlib.use('TkAgg')
    long_name = da.attrs['long_name']
    minval = da.min().data
    maxval = da.max().data
    print(f'Plotting {long_name}...')
    print(f'Min: {minval}')
    print(f'Max: {maxval}')
    if contour_levels is not None:
        print(f'User specified contour levels: {contour_levels}')
        levels = np.array(contour_levels)
    else:
        # using levels from cesm2-marbl repo
        levels_dict = dict(
                        NO3=[0, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 4, 8, 12., 16, 20, 24, 28, 32,],
                        PO4=[0, 0.02, 0.04, 0.08, 0.12, 0.16, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.2, 1.4, 1.6, 1.8, 2.0],
                        SiO3=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90,],
                        )
        if var in levels_dict:
            levels = levels_dict[var]
        else:
            levels = np.linspace(_floor_1sigfig(minval), _ceil_1sigfig(maxval), 16)

        print(f'Contour levels run from {levels[0]} to {levels[-1]}')

    da.plot(levels=levels, cmap=cmocean.cm.dense, norm=colors.BoundaryNorm(levels, ncolors=len(levels)))
    plt.title(f'{long_name}\nmin: {minval:.2f}, max: {maxval:.2f}');
    plt.show()

# --------------------------------------------- #

def _parse_args():
    """ Parse command line arguments
    """

    import argparse

    parser = argparse.ArgumentParser(description='Generate a quick contour plot',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Which run to look at
    parser.add_argument('-r', '--run', action='store', dest='run', default='004', choices=['004', '005'],
                        help='Which of the two runs to look at')

    # Which stream to look at
    parser.add_argument('-s', '--stream', action='store', dest='stream', default='monthly', choices=['monthly', 'daily', '5day'],
                        help='Which stream to look at')

    # Which variable to plot
    parser.add_argument('-v', '--var', action='store', dest='var', required=True,
                        help='Which of the variables to plot')

    # Which level to plot
    parser.add_argument('-l', '--level', action='store', dest='level', default=0,
                        help='Which level to plot (index, not physical depth)')

    # What contour levels to plot
    parser.add_argument('-c', '--contour-levels', nargs='+', type=float, action='store', dest='contour_levels', default=None,
                        help='What contour levels to use')

    return parser.parse_args()

# --------------------------------------------- #

if __name__ == '__main__':
    args = _parse_args()
    case = f'g.e22a06b.G1850ECOIAF_JRA_HR.TL319_t13.first_month.{args.run}'
    suffixes = {'monthly' : '0001-01.nc',
                'daily' : 'nday1.0001-01-01.nc',
                '5day' : 'nday5.0001-01-05.nc'
               }
    infile = os.path.join(case, f'{case}.pop.h.{suffixes[args.stream]}')
    ds = xr.open_dataset(infile).isel(z_t=args.level, z_t_150m=args.level, z_w=args.level, z_w_top=args.level, z_w_bot=args.level)
    _make_plot(ds[args.var], args.var, args.contour_levels)
