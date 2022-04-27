# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 18:15:35 2022

Used for plottinf future H2 scenarios for Section 3.4

@author: Filoteea Moldovan
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
import csv
import random
import matplotlib.cm as cm
import scipy.linalg
import scipy.stats
from scipy.stats import pearsonr
from numpy import genfromtxt


# Third party imports
from collections import OrderedDict
from datetime import datetime

# Semi-local imports
import name_qch4_couple.io
import name_qch4_couple.name
import name_qch4_couple.plot_h2


# Local imports
import chem_co



# Plots

date = '2018-04'
# Dates
dates_tHour = pd.date_range(
    pd.to_datetime(date),
    pd.to_datetime(date) + pd.DateOffset(months=12),
    closed='left',
    freq='1H'
    )

# import scenarios        

mod_0, sigma_obs_H2 = chem_co.read_obs(dates_tHour, 'mod', 0)  
 
mod_4, sigma_obs_H2 = chem_co.read_obs(dates_tHour, 'mod', 0.035)  

mod_8, sigma_obs_H2 = chem_co.read_obs(dates_tHour, 'mod', 0.165)  

mod_12, sigma_obs_H2 = chem_co.read_obs(dates_tHour, 'mod', 0.241)  

# import modelled 'baselines'
 
bas_mhd, sigma_obs_H2 = chem_co.read_obs(dates_tHour, 'bas_mhd', 0)  

bas_wao, sigma_obs_H2 = chem_co.read_obs(dates_tHour, 'bas_wao', 0)  


# plot H2 concentraion scenarios
fig = {}
ax = {}

fig_param = {
    'w': 10, 'h': 3,
    'px0': 0.80, 'py0': 0.50,
    'pw': 9.15, 'ph': 2.45,
    'ylblx': 0.05, 'ylbly': 1.5,  # left, centre aligned
    'fontsize': 6,
    }
plt.close('all')


# Concentration
fig['main'] = plt.figure(figsize=(fig_param['w'], fig_param['h']), dpi=300)
for i in ['H2']:
    fig['main'].clf()
    ax['main'] =  name_qch4_couple.plot_h2.generic(
        fig=fig['main'],
        idata={

            'bas_mhd': [
                'fill',
                [dates_tHour, np.array(bas_mhd), np.array(bas_wao)],
                { 'facecolor': '#9D9D9D', 'lw': 0.5, 'label': '', 'ls':'-'}
                ],

            'mod13': [
                'line',
                [dates_tHour, np.array(mod_12), '-'],
                {'c': '#d73027', 'lw': 0.5, 'label': ''}
                ],
            'mod9': [
                'line',
                [dates_tHour, np.array(mod_8), '-'],
                {'c': '#fc8d59', 'lw': 0.5, 'label': ''}
                ],
            'mod5': [
                'line',
                [dates_tHour, np.array(mod_4), '-'],
                {'c': '#91bfdb', 'lw': 0.5, 'label': ''}
                ],
           # 'mod1': [
           #     'line',
           #     [dates_tHour, np.array(mod_0), '--'],
           #     {'c': '#fee090', 'lw': 0.5, 'label': ''}
           #     ],
           
            
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
            pd.to_datetime(date) + pd.DateOffset(months=3), 
            ],
        ylim=(
            [470., 590.]
            ),
        yticks=(
            np.arange(470., 590., 20.)
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
            mdates.DateFormatter('%m-%d'),
            mdates.WeekdayLocator(byweekday=6),
            ]
        )
        
    for l in ax['main'].get_xticklabels():
        l.set_ha("right")
        l.set_rotation(30)

    ax['main'].legend(
        loc='upper right', ncol=7, fontsize=fig_param['fontsize']
        )

   # fig['main'].savefig(f'outputs/scenarios/new_figures/scenario_{date}_hfd.png')

