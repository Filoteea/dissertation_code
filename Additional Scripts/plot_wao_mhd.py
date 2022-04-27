# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 12:27:48 2022

@author: filot
"""

# Standard Library imports
import argparse
import gzip
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import os
import pandas as pd
import sys
import xarray as xr

# Third party imports
from collections import OrderedDict
from datetime import datetime

# Semi-local imports
import name_qch4_couple.io
import name_qch4_couple.name
import name_qch4_couple.plot_h2


# Local imports
import chem_co


date = "2018-01"

# Dates

dates_tHour = pd.date_range(
    pd.to_datetime(date),
    pd.to_datetime(date) + pd.DateOffset(months=12),
    closed='left',
    freq='1H'
    )


# observations
obs_H2, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "WAO", 0)   # read WAO dataset

obs_H2_mhd, sigma_obs_H2_mhd = chem_co.read_obs(dates_tHour, "MHD_10magl", 0)   # read MHD dataset



figs = {}
axs = {}
pobjs = {}


colours = {

    'p_nw_all': '#A6DAF3',
    'background': '#000000',
    'final': '#0000FF'
    }
zorder = {
    'f_local': 4,
    'f_europe': 3,
    'f_south': 2,
    'p_nw_all': 1,
    'background': 5,
    'final': 6
    }
fig_param = {
    'mw': 6, 'mh': 5,
    'mpw': 5.3, 'mph': 4.5,
    'mgap': 0.05,
    'mlmargin': 0.7, 'mbmargin': 0.5,
    'ylblx': 0.05, 'ylbly': 1.5,  # left, centre aligned
    'fontsize': 8,
    }
plt.close('all')

figs['main'] = plt.figure(figsize=(fig_param['mw'], fig_param['mh']), dpi=300)
axs['main'] = {}
pobjs['main'] = {}
figs['main'].clf()

year = 2018

# =============================================================================
# Plot observed H2 concentrations
# =============================================================================

ylabel = u'$\chi$ H$_{2}$ (nmol mol$^{-1}$)'
ylim = [460., 570.]
yticks = np.arange(460., 570., 20.)

name_qch4_couple.plot_h2.generic3(
    fig=figs['main'],
    axs=axs['main'],
    pobjs=pobjs['main'],
    new_axs={
        'date1': [
            dict(
                rect=[
                    (1 * fig_param['mlmargin']
                        + fig_param['mgap'])
                    / fig_param['mw'],
                    (fig_param['mh']
                        - fig_param['mph']
                        + fig_param['mgap'])
                    / fig_param['mh'],
                    (fig_param['mpw']
                        - 2*fig_param['mgap'])
                    / fig_param['mw'],
                    (fig_param['mph']
                        - 2*fig_param['mgap'])
                    / fig_param['mh']
                    ],
                label='date1',
                projection=None
                ),
            {
                "set_yticks": [[yticks], {}],
                "set_xlim": [[
                    dates_tHour[dates_tHour.year == year][0] - pd.to_timedelta('1 day'),
                    dates_tHour[dates_tHour.year == year][-1] + pd.to_timedelta('1 day')
                    ], {}],
                "set_ylim": [[ylim], {}],
                "tick_params": [[], dict(
                    axis='both', which='major', direction='in',
                    labelsize=fig_param['fontsize'],
                    left=True, bottom=True,
                    right=False, top=False,
                    labelleft=True, labelbottom=False,
                    labelright=False, labeltop=False,
                    )],
                "xaxis.set_major_locator": [
                    [mdates.DayLocator(bymonthday=1)], {}
                    ],
                "xaxis.set_major_formatter": [
                    [mdates.DateFormatter('%Y-%m-%d')], {}
                    ],
                },
            dict(
                patch_alpha=0.0
                )
            ]
        },
    new_pobjs={

        'mhd': [
            'date1', 'plot', [obs_H2_mhd.index, obs_H2_mhd.values, 'o'],
            {'c': "#000000", 'ms': 1., 'mew': 0.,
             'zorder': zorder['final'],
             'label': 'MHD'}
            ],
        'wao': [
            'date1', 'plot', [obs_H2.index, obs_H2.values, 'o'],
            {'c': "#0012FF", 'ms': 1., 'mew': 0.,
             'zorder': zorder['final'],
             'label': 'WAO'}
            ],
    
        'legend':[
            'date1', 'legend', [],
            dict(
                loc='upper right',
                numpoints=2, fontsize=fig_param['fontsize'], ncol=3,
                markerscale=5.0/1.5, handletextpad=0.2, columnspacing=1.0,
                borderpad=0.2, borderaxespad=0.2
                )
            ]
        },
    texts=[
        {
            'x': fig_param['mgap'] / fig_param['mw'],
            'y': (fig_param['mh']
                    - 1/2*fig_param['mph'])
                / fig_param['mh'],
            's': ylabel,
            'ha': 'left', 'va': 'center',
            'size': fig_param['fontsize'], 'rotation': 90
            }
        ],
    legend_params=[
        [],
        [],
        {}
        ]
    )
for l in axs['main']['date1'].get_xticklabels():
    l.set_ha("right")
    l.set_rotation(30)
#figs['main'].savefig(os.path.join('outputs/observations/redisual.png'))


# =============================================================================
# Plot difference in observed H2 mixing ratios
# =============================================================================

ylabel = u'MHD - WAO (nmol mol$^{-1}$)'
ylim = [-50., 90.]
yticks = np.arange(-50., 90., 20.)


year = 2018

name_qch4_couple.plot_h2.generic3(
    fig=figs['main'],
    axs=axs['main'],
    pobjs=pobjs['main'],
    new_axs={
        'date1': [
            dict(
                rect=[
                    (1 * fig_param['mlmargin']
                        + fig_param['mgap'])
                    / fig_param['mw'],
                    (fig_param['mh']
                        - fig_param['mph']
                        + fig_param['mgap'])
                    / fig_param['mh'],
                    (fig_param['mpw']
                        - 2*fig_param['mgap'])
                    / fig_param['mw'],
                    (fig_param['mph']
                        - 2*fig_param['mgap'])
                    / fig_param['mh']
                    ],
                label='date1',
                projection=None
                ),
            {
                "set_yticks": [[yticks], {}],
                "set_xlim": [[
                    dates_tHour[dates_tHour.year == year][0] - pd.to_timedelta('1 day'),
                    dates_tHour[dates_tHour.year == year][-1] + pd.to_timedelta('1 day')
                    ], {}],
                "set_ylim": [[ylim], {}],
                "tick_params": [[], dict(
                    axis='both', which='major', direction='in',
                    labelsize=fig_param['fontsize'],
                    left=True, bottom=True,
                    right=False, top=False,
                    labelleft=True, labelbottom=False,
                    labelright=False, labeltop=False,
                    )],
                "xaxis.set_major_locator": [
                    [mdates.DayLocator(bymonthday=1)], {}
                    ],
                "xaxis.set_major_formatter": [
                    [mdates.DateFormatter('%Y-%m-%d')], {}
                    ],
                },
            dict(
                patch_alpha=0.0
                )
            ]
        },
    new_pobjs={

        'mhd': [
            'date1', 'plot', [obs_H2_mhd.index, obs_H2_mhd.values - obs_H2.values, 'o'],
            {'c': "#000000", 'ms': 1., 'mew': 0.,
             'zorder': zorder['final'],
             }
            ],

        'legend':[
            'date1', 'legend', [],
            dict(
                loc='upper right',
                numpoints=2, fontsize=fig_param['fontsize'], ncol=3,
                markerscale=5.0/1.5, handletextpad=0.2, columnspacing=1.0,
                borderpad=0.2, borderaxespad=0.2
                )
            ]
        },
    texts=[
        {
            'x': fig_param['mgap'] / fig_param['mw'],
            'y': (fig_param['mh']
                    - 1/2*fig_param['mph'])
                / fig_param['mh'],
            's': ylabel,
            'ha': 'left', 'va': 'center',
            'size': fig_param['fontsize'], 'rotation': 90
            }
        ],
    legend_params=[
        [],
        [],
        {}
        ]
    )
for l in axs['main']['date1'].get_xticklabels():
    l.set_ha("right")
    l.set_rotation(30)
#figs['main'].savefig(os.path.join('outputs/observations/redisual.png'))