# -*- coding: utf-8 -*-
"""
Created on Sun March 20 13:23:31 2022

Script used to plot the final model used for future scenarios
@author: Filoteea Moldovan
Adapted from Edward Chung
"""

# Standard Library imports
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from numpy import genfromtxt


# Semi-local imports
import name_qch4_couple.plot_h2


# Local imports
import chem_co



'''
Inputs:
    - WAO observations
    - WAO modelled mole fractions calculated in create_models.py

    S13 (Appendix 7.2.1.) was used     
'''

# Plots

date = '2018-01'
# Dates
dates_tHour = pd.date_range(
    pd.to_datetime(date),
    pd.to_datetime(date) + pd.DateOffset(months=12),
    closed='left',
    freq='1H'
    )

        
# import model output
mod_0, sigma_obs_H2 = chem_co.read_obs(dates_tHour, 'mod', 0)  
 
# import WAO observations
obs_wao, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "WAO", 0)   # could add st dev

# calculate correlation

    #remove nans
for i in range(0,len(obs_wao)):
    if np.isnan(obs_wao[i]):
        obs_wao[i]=523.024
        
 
corr, p_value = pearsonr(obs_wao, mod_0)


obs_wao, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "WAO", 0)   # could add st dev

for i in range(0, len(obs_wao)):
    if np.isnan(obs_wao[i]):
        mod_0[i] = 0

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
ylim = [450., 600.]
yticks = np.arange(450., 600., 20.)

var_long_name = 'mole_fraction_of_hydrogen'
var_units = 'nmol mol-1'

for i in range(12, 13):         # allows to create the plots for all the scenarios in one run
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
                'date1', 'plot', [dates_tHour, np.array(obs_wao), 'o'],
                {'c': '#000000', 'ms': 1., 'mew': 1.,
                 'zorder': zorder['final'],
                 'label': 'Observed WAO'}
                ],
          
            'mod': [
                'date1', 'plot', [dates_tHour, np.array(mod_0), 'o'],
                {'c': '#0012FF', 'ms': 1., 'mew': 1.,
                 'zorder': zorder['final'],
                 'label': 'Modelled WAO'}
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

  #  figs['main'].savefig(f'outputs/obs_dif_sd.png')


