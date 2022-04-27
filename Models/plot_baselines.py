# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 03:47:03 2022

@author: Filoteea Moldovan
Adapted from E. Chung

Code for plotting modelled baselines with tweaked emissions and sinks

"""

# plots models for 2018

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
from scipy.stats import spearmanr
from scipy.stats import pearsonr



# Third party imports
from collections import OrderedDict
from datetime import datetime
from numpy import genfromtxt


# Semi-local imports
import name_qch4_couple.io
import name_qch4_couple.name
import name_qch4_couple.plot_h2


# Local imports
import chem_co


'''
# Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument("-date", required=True)    # yyyy-mm
args = parser.parse_args()

date = args.date
'''
date = '2018-01'

'''
# import models
if date == '2018-01':    
    bas_mhd = genfromtxt('outputs/models/baselines/higher_emm/merged_01_mhd.csv', delimiter=',')
    bas_wao = genfromtxt('outputs/models/baselines/higher_emm/merged_01_wao.csv', delimiter=',')
else:
    bas_mhd = genfromtxt('outputs/models/baselines/higher_emm/merged_07_mhd.csv', delimiter=',')
    bas_wao = genfromtxt('outputs/models/baselines/higher_emm/merged_07_wao.csv', delimiter=',')

    
bas_mhd = np.transpose(bas_mhd)
bas_wao = np.transpose(bas_wao)
'''

# Dates
dates_tHour = pd.date_range(
    pd.to_datetime(date),
    pd.to_datetime(date) + pd.DateOffset(months=12),
    closed='left',
    freq='1H'
    )

# observations
obs_mhd, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "MHD_10magl", 0)  #read MHD observations
obs_wao, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "WAO", 0)   #read WAO observations 




# Plot settings

# generic 3 
figs = {}
axs = {}
pobjs = {}

zorder = {
    'background': 1,
    'final': 2
    }
fig_param = {
    'mw': 10.5, 'mh': 7,
    'mpw': 8.3, 'mph': 5.5,
    'mgap': 0.05,
    'mlmargin': 1.2, 'mbmargin': 0.5,
    'ylblx': 0.05, 'ylbly': 1.5,  # left, centre aligned
    'fontsize': 15,
    }
plt.close('all')

ylabel = u'H$_{2}$ (nmol mol$^{-1}$)'
ylabel2 = 'Residuals'
ylim = [400., 550.]
ylim2 = [-50., 400.]

yticks = np.arange(400., 550., 20.)
yticks2 = np.arange(-50., 400., 100.)

var_long_name = 'mole_fraction_of_hydrogen'
var_units = 'nmol mol-1'


# =============================================================================
# Plot observed concentrations and their difference and SD
# =============================================================================        


figs['main'] = plt.figure(figsize=(fig_param['mw'], fig_param['mh']), dpi=300)
axs['main'] = {}
pobjs['main'] = {}
figs['main'].clf()
i = 0

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
                pd.to_datetime(date), 
                pd.to_datetime(date) + pd.DateOffset(months=12), 
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

        'wao': [
            'date1', 'plot', [dates_tHour, np.array(obs_wao[i]), '-'],
            {'c': '#000000', 'ms': 1., 'mew': 0.,
             'zorder': zorder['background'],
             'label': 'Modelled WAO baseline'}
            ],
        'mhd': [
            'date1', 'plot', [dates_tHour, np.array(obs_mhd[i]), '-'],
            {'c': '#0012FF', 'ms': 1., 'mew': 0.,
             'zorder': zorder['final'],
             'label': 'Modelled MHD baseline'}
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
                "set_yticks": [[yticks2], {}],
                "set_xlim": [[
                pd.to_datetime(date), 
                pd.to_datetime(date) + pd.DateOffset(months=12), 
                    ], {}],
                "set_ylim": [[ylim2], {}],
                "tick_params": [[], dict(
                    axis='both', which='major', direction='in',
                    labelsize=fig_param['fontsize'],
                    left=False, bottom=True,
                    right=True, top=False,
                    labelleft=False, labelbottom=True,
                    labelright=True, labeltop=False,
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

        'residual': [
            'date1', 'plot', [dates_tHour, np.array(obs_mhd[i]) - np.array(obs_wao[i]), '-'],
            {'c': '#FF0000', 'ms': 1., 'mew': 0.,
             'zorder': zorder['final'],
             'label': 'MHD-WAO'}
            ],

        'legend':[
            'date1', 'legend', [],
            dict(
                loc='upper left',
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
            's': ylabel2,
            'ha': 'right', 'va': 'center',
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

#figs['main'].savefig(f'outputs/models/baselines/higher_emm/twiny/trial3ndy_{i}_{date}.png')



        
        
# =============================================================================
# Plot all modelled 'baselines' and their difference with the observed concentrations and SD
# =============================================================================        

'''
mhd = obs_mhd
       
    #remove nans
for i in range(0,len(mhd)):
    if np.isnan(mhd[i]):
        mhd[i]=523.024
        
wao = obs_wao
       
    #remove nans
for i in range(0,len(wao)):
    if np.isnan(wao[i]):
        wao[i]=510.400
        
obs_wao, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "WAO")   # read WAO observations

obs_mhd, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "MHD_10magl")   # read MHD observations


    # create correlation list
corr_mhd = np.ndarray(shape = (25))
i = 0
while i < 25:
    corr_mhd[i], _ = pearsonr(mhd, bas_mhd[i])
    i += 1

corr_wao = np.ndarray(shape = (25))
i = 0
while i < 25:
    corr_wao[i], _ = pearsonr(wao, bas_wao[i])
    i += 1

#create the 25 plots

for i in range(0, 24):
    fig['main'] = plt.figure(figsize=(fig_param['w'], fig_param['h']), dpi=300)
    for j in ['H2']:
        fig['main'].clf()
        ax['main'] =  name_qch4_couple.plot_h2.generic(
            fig=fig['main'],
            idata={
                'mhd': [
                    'line',
                    [dates_tHour, np.array(obs_mhd), '-'],
                    {'c': '#000000', 'lw': 0.5, 'label': 'Observed MHD'}
                    ],
                'model': [
                    'line',
                    [dates_tHour, np.array(bas_mhd[i]), '-'],
                    {'c': '#0012FF', 'lw': 0.5, 'label': 'Modelled MHD baseline: {:.2f}'.format(corr_mhd[i])}
                    ],
                
                  },
    
            texts=[
                {
                    'x': fig_param['ylblx'] / fig_param['w'],
                    'y': fig_param['ylbly'] / fig_param['h'],
                    's': (
                        u'H$_{2}$ (nmol mol$^{-1}$)'
                        ),
                    'ha': 'left', 'va': 'center',
                    'size': fig_param['fontsize'], 'rotation': 90
                    }
                ],
            xlim=[
                pd.to_datetime(date), 
                pd.to_datetime(date) + pd.DateOffset(months=6), 
                ],
            ylim=(
                [400., 660.]
                ),
            yticks=(
                np.arange(400., 660., 40.)
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
            loc='upper right', ncol=7, fontsize=fig_param['fontsize1']
            )
        fig['main'].savefig(f'outputs/models/baselines/mhd/mhd_{i}_{date}.png')
        
        
        
        
    fig['main'] = plt.figure(figsize=(fig_param['w'], fig_param['h']), dpi=300)
    for j in ['H2']:
        fig['main'].clf()
        ax['main'] =  name_qch4_couple.plot_h2.generic(
            fig=fig['main'],
            idata={
                'wao': [
                    'line',
                    [dates_tHour, np.array(obs_wao), '-'],
                    {'c': '#000000', 'lw': 2, 'label': 'Observed WAO'}
                    ],
    
                'model8': [
                    'line',
                    [dates_tHour, np.array(bas_wao[i]), '-'],
                    {'c': '#0012FF', 'lw': 0.5, 'label': 'Modelled WAO baseline: {:.2f}'.format(corr_wao[i])}
                    ],
                
                  },
    
            texts=[
                {
                    'x': fig_param['ylblx'] / fig_param['w'],
                    'y': fig_param['ylbly'] / fig_param['h'],
                    's': (
                        u'H$_{2}$ (nmol mol$^{-1}$)'
                        ),
                    'ha': 'left', 'va': 'center',
                    'size': fig_param['fontsize'], 'rotation': 90
                    }
                ],
            xlim=[
                pd.to_datetime(date), 
                pd.to_datetime(date) + pd.DateOffset(months=6), 
                ],
            ylim=(
                [400., 660.]
                ),
            yticks=(
                np.arange(400., 660., 40.)
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
            loc='upper right', ncol=7, fontsize=fig_param['fontsize1']
            )
        fig['main'].savefig(f'outputs/models/baselines/wao/wao_{i}_{date}.png')


'''