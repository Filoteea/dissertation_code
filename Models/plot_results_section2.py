# -*- coding: utf-8 -*-
"""
Created on March 17 04:50:04 2022

@author: Filoteea Moldovan
Adapted from: Edward Chung

Script used for creating the plots in section 3.3.
"""
# Standard Library imports

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy import genfromtxt

# Semi-local imports

import name_qch4_couple.plot_h2

# Local imports

import chem_co


# Set the date-time variable
date = '2018-01'

dates_tHour = pd.date_range(
    pd.to_datetime(date),
    pd.to_datetime(date) + pd.DateOffset(months=12),
    closed='left',
    freq='1H'
    )

# =============================================================================
#   Plotting the MHD and WAO observations and their difference and
#  calculating the SD of this difference
# =============================================================================

'''
Inputs:
    - H2 observations at WAO for 2018
    - H2 observations at MHD for 2018
'''

# read observations
obs_mhd, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "MHD_10magl", 0)  
obs_wao, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "WAO", 0) 

'''
# calculate SD
dif = []
z = 0
count = 0
for i, j in zip(obs_mhd, obs_wao):
    if np.isnan(i) or np.isnan(j):
        count += 1
    else:
        dif.append(i - j)
        z += 1
dev = np.std(dif)


# Plot

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

ylabel = u'$\chi$ H$_{2}$ (nmol mol$^{-1}$)'
ylabel2 = 'Differences (nmol mol$^{-1}$)'
ylim = [350., 600.]
ylim2 = [-70., 600.]

yticks = np.arange(300., 600., 50.)
yticks2 = np.arange(-70., 600., 100.)

var_long_name = 'mole_fraction_of_hydrogen'
var_units = 'nmol mol-1'

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
  
            'mhd': [
                'date1', 'plot', [dates_tHour, np.array(obs_mhd), 'o'],
                {'c': '#000000', 'ms': 1., 'mew': 1.,
                 'zorder': zorder['final'],
                 'label': 'Observed MHD'}
                ],
          
            'wao': [
                'date1', 'plot', [dates_tHour, np.array(obs_wao), 'o'],
                {'c': '#0012FF', 'ms': 1., 'mew': 1.,
                 'zorder': zorder['final'],
                 'label': 'Observed WAO'}
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
                'date1', 'plot', [dates_tHour, np.array(obs_mhd)-np.array(obs_wao), '--'],
                {'c': '#767676', 'ms': 1., 'mew': 0.,
                 'zorder': zorder['final'],
                 'label': 'MHD-WAO - SD: {:.2f}'.format(dev)}
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

  #  figs['main'].savefig(f'outputs/obs_dif_sd.png')
''' 

# =============================================================================
#   Plotting the MHD and WAO modelled 'baselines' and their difference and
#  calculating the SD of this difference
# =============================================================================

'''
Inputs:
    - MHD modelled 'baseline' calculated in create_baseline.py
    - WAO modelled 'baseline' calculated in create_baseline.py
    S13 (Appendix 7.2.1.) was used for the final result
'''

# import modelled 'baselines'
   
bas_mhd = genfromtxt('outputs/models/baselines/lower_emm/2018/2018_mhd.csv', delimiter=',')
bas_wao = genfromtxt('outputs/models/baselines/lower_emm/2018/2018_wao.csv', delimiter=',')
bas_mhd = np.transpose(bas_mhd)
bas_wao = np.transpose(bas_wao)

# calculate SD
dev = np.std(bas_mhd[12] - bas_wao[12])

# remove data point where observations are missing
dif = []
z = 0
count = 0
for i in range(0, len(obs_mhd)):
    if np.isnan(obs_mhd[i]) or np.isnan(obs_wao[i]):
        dif.append(np.nan)
    else:
        dif.append(bas_mhd[12][i] - bas_wao[12][i])



mhd = bas_mhd[12]
wao = bas_wao[12]
for i in range(0, len(obs_mhd)):
    if np.isnan(obs_mhd[i]):
        mhd[i] = 0
    if np.isnan(obs_wao[i]):
        wao[i] = 0 


# Plot

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

ylabel = u'$\chi$ H$_{2}$ (nmol mol$^{-1}$)'
ylabel2 = 'Differences (nmol mol$^{-1}$)'
ylim = [350., 600.]
ylim2 = [-70., 600.]

yticks = np.arange(300., 600., 50.)
yticks2 = np.arange(-70., 600., 100.)

var_long_name = 'mole_fraction_of_hydrogen'
var_units = 'nmol mol-1'

for i in range(12, 13):         # allows to create the plots for all the scenarios in one run
                                # in this case only plotting the modelled 'baseline' with the lowest SD
    figs['main'] = plt.figure(figsize=(fig_param['mw'], fig_param['mh']), dpi=300)
    axs['main'] = {}
    pobjs['main'] = {}
    figs['main'].clf()
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
  
            'mhd': [
                'date1', 'plot', [dates_tHour, np.array(mhd), 'o'],
                {'c': '#000000', 'ms': 1., 'mew': 1.,
                 'zorder': zorder['final'],
                 'label': 'Mobelled MHD baseline'}
                ],
          
            'wao': [
                'date1', 'plot', [dates_tHour, np.array(wao), 'o'],
                {'c': '#0012FF', 'ms': 1., 'mew': 1.,
                 'zorder': zorder['final'],
                 'label': 'Modelled WAO baseline'}
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
                'date1', 'plot', [dates_tHour, np.array(dif), '--'],
                {'c': '#767676', 'ms': 1., 'mew': 0.,
                 'zorder': zorder['final'],
                 'label': 'MHD-WAO - SD: {:.2f}'.format(dev)}
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

  #  figs['main'].savefig(f'outputs/obs_dif_sd.png')
