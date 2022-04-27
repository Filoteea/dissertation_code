# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 00:45:50 2022

@author: Filoteea Moldovan
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
import routines
import chem_co
import chem_ch4_validation
import merge_csv
import create_model
import corr_try

# Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument("-date", required=True)    # yyyy-mm
args = parser.parse_args()

date = args.date

# Dates
dates_tHour = pd.date_range(
    pd.to_datetime(date),
    pd.to_datetime(date) + pd.DateOffset(months=12),
    closed='left',
    freq='1H'
    )


# import modelled 'baselines' for either MHD or WAO
#bas_mhd = genfromtxt('outputs/models/baselines/lower_emm/2018/2018_mhd.csv', delimiter=',')
bas_wao = genfromtxt('outputs/models/baselines/lower_emm/2018/2018_wao.csv', delimiter=',')
#bas_mhd = np.transpose(bas_mhd)
bas_wao = np.transpose(bas_wao)


# import obs
#obs_mhd, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "MHD_10magl")  
obs_wao, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "WAO") 


# calculate difference between modelled 'baseline' and observations
dif = []
z = 0
count = 0
for i, j in zip(obs_mhd, bas_mhd):
    if np.isnan(i) or np.isnan(j):
        count += 1
    else:
        dif.append(i - j)
        z += 1
#calculate SD
dev = np.std(dif)


# =============================================================================
# Plot modelled 'baseline' along observations and their difference and the difference SD
# =============================================================================

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
    'mpw': 8.5, 'mph': 5.7,
    'mgap': 0.05,
    'mlmargin': 1.2, 'mbmargin': 1.5,
    'ylblx': 0.05, 'ylbly': 1.5,  # left, centre aligned
    'fontsize': 15,
    'fontsize2': 12,

    }
plt.close('all')

ylabel = u'H$_{2}$ (nmol mol$^{-1}$)'
#ylabel = u'$\chi$ H$_{2}$ (nmol mol$^{-1}$)'
ylabel2 = 'Residuals (nmol mol$^{-1}$)'
ylim = [300., 600.]
ylim2 = [-70., 600.]

yticks = np.arange(300., 600., 50.)
yticks2 = np.arange(-70., 600., 100.)

var_long_name = 'mole_fraction_of_hydrogen'
var_units = 'nmol mol-1'



# plot modelled against observed
for i in range(12, 13):
    figs['main'] = plt.figure(figsize=(fig_param['mw'], fig_param['mh']), dpi=300)
    axs['main'] = {}
    pobjs['main'] = {}
    figs['main'].clf()
 #   dev = np.std(np.array(bas_mhd[i]) - np.array(bas_wao[i]))
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
                        left=True, bottom=False,
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
  
            'mhd observed': [
                'date1', 'plot', [dates_tHour, np.array(obs_wao), '-'],
                {'c': '#525255', 'ms': 1., 'mew': 0.,
                 'zorder': zorder['final'],
                 'label': 'Observed WAO '}
                ],
          
            'mhd': [
                'date1', 'plot', [dates_tHour, np.array(bas_wao[i]), 'o'],
                {'c': '#0012FF', 'ms': 1., 'mew': 0.,
                 'zorder': zorder['final'],
                 'label': f'Modelled WAO baseline M{i+1}'}
                ],
            

    
            'legend':[
                'date1', 'legend', [],
                dict(
                    loc='upper left',
                    numpoints=2, fontsize=fig_param['fontsize2'], ncol=3,
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
                'date1', 'plot', [dates_tHour, np.array(obs_wao) - np.array(bas_wao[i]), '--'],
                {'c': '#767676', 'ms': 1., 'mew': 0.,
                 'zorder': zorder['final'],
                 'label': 'Observed - Modelled'}
                ],
            'legend':[
                'date1', 'legend', [],
                dict(
                    loc='upper right',
                    numpoints=2, fontsize=fig_param['fontsize2'], ncol=3,
                    markerscale=5.0/3.5, handletextpad=0.2, columnspacing=1.0,
                    borderpad=0.2, borderaxespad=0.2
                    )
                ]
            },
        
        texts=[
            {
                'x': 1,
                'y': (fig_param['mh']
                        - 1/2*fig_param['mph'])
                    / fig_param['mh'],
                's': ylabel2,
                'ha': 'right', 'va': 'center',
                'size': fig_param['fontsize'], 'rotation': 270
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

    figs['main'].savefig(f'outputs/models/baselines/lower_emm/2018/wao_{i+1}.png')

'''
# plot baseline scenarios
for i in range(0, 25):
    figs['main'] = plt.figure(figsize=(fig_param['mw'], fig_param['mh']), dpi=300)
    axs['main'] = {}
    pobjs['main'] = {}
    figs['main'].clf()
    dev = np.std(np.array(bas_mhd[i]) - np.array(bas_wao[i]))
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
                        left=True, bottom=False,
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
    
            'wao': [
                'date1', 'plot', [dates_tHour, np.array(bas_wao[i]), '-'],
                {'c': '#000000', 'ms': 1., 'mew': 0.,
                 'zorder': zorder['background'],
                 'label': 'Modelled WAO baseline'}
                ],
            'mhd': [
                'date1', 'plot', [dates_tHour, np.array(bas_mhd[i]), '-'],
                {'c': '#0012FF', 'ms': 1., 'mew': 0.,
                 'zorder': zorder['final'],
                 'label': 'Modelled MHD baseline'}
                ],
    
            'legend':[
                'date1', 'legend', [],
                dict(
                    loc='upper left',
                    numpoints=2, fontsize=fig_param['fontsize2'], ncol=3,
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
                'date1', 'plot', [dates_tHour, np.array(bas_mhd[i]) - np.array(bas_wao[i]), '--'],
                {'c': '#767676', 'ms': 1., 'mew': 0.,
                 'zorder': zorder['final'],
                 'label': 'MHD-WAO  SD:{:.2f}'.format(dev)}
                ],
            'legend':[
                'date1', 'legend', [],
                dict(
                    loc='upper right',
                    numpoints=2, fontsize=fig_param['fontsize2'], ncol=3,
                    markerscale=5.0/3.5, handletextpad=0.2, columnspacing=1.0,
                    borderpad=0.2, borderaxespad=0.2
                    )
                ]
            },
        
        texts=[
            {
                'x': 1,
                'y': (fig_param['mh']
                        - 1/2*fig_param['mph'])
                    / fig_param['mh'],
                's': ylabel2,
                'ha': 'right', 'va': 'center',
                'size': fig_param['fontsize'], 'rotation': 270
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

    figs['main'].savefig(f'outputs/models/baselines/lower_emm/2018/bas_{i+1}.png')
'''