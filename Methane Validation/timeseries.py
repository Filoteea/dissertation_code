# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 00:44:36 2022

@author: filot
Create timeseries
"""
import pandas as pd
import glob
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
from sklearn import linear_model
from scipy.stats import spearmanr
from scipy.stats import pearsonr

# Semi-local imports
import name_qch4_couple.io
import name_qch4_couple.name
import name_qch4_couple.plot_h2

# Local imports
import routines
import chem_ch4_validation
import chem_co


#os.chdir('C:/Users/filot/Desktop/YEAR_4/Dissertation/Ed_new_script/scripts')

date = "2018-01"
# Dates
dates_tHour = pd.date_range(
    pd.to_datetime(date),
    pd.to_datetime(date) + pd.DateOffset(months=12),
    closed='left',
    freq='1H'
    )

# create timeseries!!
'''
files = ["hfd_ch4_2018-01.csv", "hfd_ch4_2018-02.csv", "hfd_ch4_2018-03.csv", "hfd_ch4_2018-04.csv",
         "hfd_ch4_2018-05.csv", "hfd_ch4_2018-06.csv", "hfd_ch4_2018-07.csv", "hfd_ch4_2018-08.csv",
         "hfd_ch4_2018-09.csv", "hfd_ch4_2018-10.csv", "hfd_ch4_2018-11.csv", "hfd_ch4_2018-12.csv"]

# Combine all three CSV files using the concat method
a_b_c = pd.concat([pd.read_csv(f) for f in files])
# Export to csv
a_b_c.to_csv( "ABC.csv", index=False, encoding='utf-8-sig')
'''

# correlation
obs, sigma_obs_H2_mhd = chem_ch4_validation.read_obs(dates_tHour, "hfd")   # could add st dev

mod_ch4, sigma_obs_H2_mhd = chem_ch4_validation.read_obs(dates_tHour, "mod")   # could add st dev


        
    #imterpolate missing data with average
for i in range(0,len(obs)):
    if np.isnan(obs[i]):
        obs[i]=523.024

corr, p_value = pearsonr(obs.values, mod_ch4.values)

obs_ch4, sigma_obs_H2_mhd = chem_ch4_validation.read_obs(dates_tHour, "hfd")   # could add st dev


# obs mhd
#obs_ch4, sigma_obs_H2_mhd = chem_ch4_validation.read_obs(dates_tHour, "MHD_10magl")   # could add st dev
#mod_ch4, sigma_obs_H2_mhd = chem_ch4_validation.read_obs(dates_tHour, "mod")   # could add st dev
#obs_ch4, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "WAO")   # could add st dev


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


ylabel = u'$\chi$ CH$_{4}$ (nmol mol$^{-1}$)'
ylim = [1850., 2320.]
yticks = np.arange(1850., 2320., 100.)
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

        'obs': [
            'date1', 'plot', [obs_ch4.index, obs_ch4.values, 'o'],
            {'c': "#000000", 'ms': 1., 'mew': 0.,
             'zorder': zorder['final'],
             'label': 'Observed HFD'}
            ],
        'mod': [
            'date1', 'plot', [mod_ch4.index, mod_ch4.values, 'o'],
            {'c': "#0000FF", 'ms': 1., 'mew': 0.,
             'zorder': zorder['final'],
             'label': 'Modelled HFD'}
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
    l.set_rotation(35)
figs['main'].savefig(os.path.join('outputs/validation_ch4/timeseries_2018_chi.png'))

'''
fig = {}
ax = {}
colours = {
    'obs_mhd': '#000000',
    'mod_ch4': '#0000FF',
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
for i in ['CH4']:
    fig['main'].clf()
    ax['main'] =  name_qch4_couple.plot_h2.generic(
        fig=fig['main'],
        idata={
            'hfd': [
                'line',
                [obs_ch4.index, np.array(obs_ch4), 'o'],
                {'c': colours['obs_mhd'], 'lw': 0.5, 'label': 'Measured HFD'}
                ],

            'mod': [
                'line',
                [mod_ch4.index, np.array(mod_ch4), 'o'],
                {'c': colours['mod_ch4'], 'lw': 0.5, 'label': 'Modelled HFD'}
                ],
            },

        texts=[
            {
                'x': fig_param['ylblx'] / fig_param['w'],
                'y': fig_param['ylbly'] / fig_param['h'],
                's': (
                    u'CH$_{4}$ (nmol mol$^{-1}$)'
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
            [1800., 2320.]
            ),
        yticks=(
            np.arange(1850., 2400., 100.)
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
            mdates.MonthLocator(bymonthday=1),
            ]
        )
    for l in ax['main'].get_xticklabels():
        l.set_ha("right")
        l.set_rotation(30)
    ax['main'].legend(
        loc='upper right', ncol=4, fontsize=fig_param['fontsize']
        )
    fig['main'].savefig('outputs/validation_ch4/timeseries_2018_{i}.png')
'''