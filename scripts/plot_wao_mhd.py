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
import routines
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
obs_H2, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "WAO")   # could add st dev

obs_H2_mhd, sigma_obs_H2_mhd = chem_co.read_obs(dates_tHour, "MHD_10magl")   # could add st dev

def read_baseline(timestamps):
    year = timestamps[0].strftime('%Y')
    date = timestamps[0].strftime('%Y-%m')
    chi0file = (
                'outputs/baseline/baseline-MHD_10magl-h2-2018.nc'
                )
    with xr.open_dataset(chi0file) as ds_read:     #put as
        with ds_read.load() as ds:
                chi0 = ds.chi_H2.sel(time=year).to_series()
                
    return chi0

chi0 = read_baseline(dates_tHour)


# residual
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


ylabel = u'MHD - WAO (nmol mol$^{-1}$)'
ylim = [-50., 100.]
yticks = np.arange(-50., 100., 20.)

#ylabel = u'H$_{2}$ (nmol mol$^{-1}$)'
#ylim = [470., 580.]
#yticks = np.arange(470., 580., 20.)

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
                    labelleft=True, labelbottom=True,
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
       # 'mhd': [
       #     'date1', 'plot', [obs_H2_mhd.index, obs_H2_mhd.values, 'o'],
       #     {'c': "#a50026", 'ms': 1., 'mew': 0.,
       #      'zorder': zorder['final'],
       #      'label': 'MHD'}
       #     ],
       # 'wao': [
       #     'date1', 'plot', [obs_H2.index, obs_H2.values, 'o'],
       #     {'c': "#fdae61", 'ms': 1., 'mew': 0.,
       #      'zorder': zorder['final'],
       #      'label': 'WAO'}
       #     ],
       # 'chi0': [
       #     'date1', 'plot', [chi0.index, chi0.values, '-'],
       #     {'c': "#d73027", 'ms': 1., 'mew': 0.,
       #      'zorder': zorder['final'],
       #      'label': 'baseline'}
       #     ],
    
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
figs['main'].savefig(os.path.join('outputs/observations/redisual.png'))


'''
# Plot
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

ylabel = u'$\chi$ H$_{2}$ (nmol mol$^{-1}$)'
ylim = [460., 600.]
yticks = np.arange(460., 600., 30.)
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
                    labelleft=True, labelbottom=True,
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
            'date1', 'plot', [obs_H2_mhd.index, obs_H2_mhd.values, '-'],
            {'c': "#540303", 'ms': 1., 'mew': 0.,
             'zorder': zorder['final'],
             'label': 'MHD observations'}
            ],

        'wao': [
            'date1', 'plot', [obs_H2.index, obs_H2.values, '-'],
            {'c': "#FC8E02", 'ms': 1., 'mew': 0.,
             'zorder': zorder['final'],
             'label': 'WAO observations'}
            ],
        'processsed': [
            'date1', 'plot', [chi0.index, chi0.values, '-'],
            {'c': "#FC0202", 'ms': 1., 'mew': 0.,
             'zorder': zorder['final'],
             'label': 'Porocessed MHD'}
            ],

        'legend':[
            'date1', 'legend', [],
            dict(
                loc='upper right',
                numpoints=2, fontsize=fig_param['fontsize'], ncol=3,
                markerscale=5.0/3.5, handletextpad=0.2, columnspacing=1.0,
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
figs['main'].savefig(os.path.join('outputs/baseline/baseline-mhd-wao-processed-2018-01.png'))
'''



# Plots
'''
fig = {}
ax = {}
colours = {
    'obs': '#000000',
    'mo1': '#8888FF',
    'bas': '#886600',
    'mo2': '#FF8800',
    }
fig_param = {
    'w': 6, 'h': 3,
    'px0': 0.80, 'py0': 0.50,
    'pw': 5.15, 'ph': 2.45,
    'ylblx': 0.05, 'ylbly': 1.5,  # left, centre aligned
    'fontsize': 8,
    }
plt.close('all')

# Concentration
fig['main'] = plt.figure(figsize=(fig_param['w'], fig_param['h']), dpi=300)
for i in ['H2']:
    fig['main'].clf()
    ax['main'] =  name_qch4_couple.plot_h2.generic(
        fig=fig['main'],
        idata={
            'obs': [
                'line',
                [obs_H2.index, np.array(obs_H2), '-'],
                {'c': colours['obs'], 'lw': 0.5, 'label': 'Measured WAO'}
                ],
            'obs2': [
                'line',
                [obs_H2_mhd.index, np.array(obs_H2_mhd), '-'],
                {'c': colours['bas'], 'lw': 0.5, 'label': 'Measured MHD'}
                ],
            'bas': [
                'line',
                [chi0.index, chi0.values, '-'],
                {'c': colours['bas'], 'lw': 0.5, 'label': 'Processed Baseline'}
                ],                
              
              },

        texts=[
            {
                'x': fig_param['ylblx'] / fig_param['w'],
                'y': fig_param['ylbly'] / fig_param['h'],
                's': (
                    u'$\chi$ H$_{2}$ (nmol mol$^{-1}$)'
                    ),
                'ha': 'left', 'va': 'center',
                'size': fig_param['fontsize'], 'rotation': 90
                }
            ],
        
        xlim=[
            pd.to_datetime(date), 
            pd.to_datetime(date) + pd.DateOffset(months=12), 
            ],
        ylim=(
            [450., 680.]
            ),
        yticks=(
            np.arange(450., 680., 20.)
            ),
        tick_fontsize=fig_param['fontsize'],
        loc_plot=[
            fig_param['px0'] / fig_param['w'],
            fig_param['py0'] / fig_param['h'],
            fig_param['pw'] / fig_param['w'],
            fig_param['ph'] / fig_param['h']
            ],
        xtick_params=[
            True,
            mdates.DateFormatter('%Y-%m-%d'),
            np.linspace(0,365,13),
            ]
        )
    for l in ax['main'].get_xticklabels():
        l.set_ha("right")
        l.set_rotation(30)
    ax['main'].legend(
        loc='upper right', ncol=4, fontsize=fig_param['fontsize']
        )
    fig['main'].savefig(os.path.join('outputs/simplified/baseline-mhd-wao-processed-2018-01.png'))

'''