# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 13:47:47 2022

@author: filot
"""

r"""

CH4 Forward Model verification

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


# Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument("-date", required=True)    # yyyy-mm
parser.add_argument("-odir", required=True)
args = parser.parse_args()


date = args.date
odir = args.odir

# Dates
dates_tHour = pd.date_range(
    pd.to_datetime(date),
    pd.to_datetime(date) + pd.DateOffset(months=1),
    closed='left',
    freq='1H'
    )

# Grid
grid_info = routines.define_grid()
inv_reg_map0 = grid_info['inv_reg_map']
nlat = grid_info['nlat']
nlon = grid_info['nlon']
area = grid_info['area']
grid_centre = grid_info['grid_centre']
grid_vertex = grid_info['grid_vertex']
inv_reg_uniq = grid_info['inv_reg_uniq']

M_CH4 = 16.043  # g mol-1 - IUPAC
M_H2 = 2.016


# using chem_co
#Q = chem_co.read_Q(dates_tHour)


# Prior emission Q 
Qfiles_H2 = OrderedDict([
    (0, [
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_AWB.nc', 'CO_emissions', '1M'],
        ]),
    (1, [
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_ENE.nc', 'CO_emissions', '1M'],   
        ]),
    (2, [
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_REF.nc', 'CO_emissions', '1M'],        
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_IND.nc', 'CO_emissions', '1M'],    
        ]),
    (3, [
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_CDS.nc', 'CO_emissions', '1M'],        
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_CRS.nc', 'CO_emissions', '1M'],        
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_TRO.nc', 'CO_emissions', '1M'],        
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_LTO.nc', 'CO_emissions', '1M'], 
        ]),
    (4, [
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_RCO.nc', 'CO_emissions', '1M'],        
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_PRO.nc', 'CO_emissions', '1M'], 
        ]),
    (5, [
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_FFF.nc', 'CO_emissions', '1M'], 
        ]),
    (6, [
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_SWD.nc', 'CO_emissions', '1M'],        
        ]),
    (7, [
        ['inputs/emissions/biomass/gfed_2015.nc', 'H2_emissions', '1M'],        
        ]),    
    ])
    

Q_factor = {
        0: [0.0357],      
        1: [0.0143],
        2: [0.0143],
        3: [0.0357],
        4: [0.0217], 
        5: [0.0143],
        6: [0.005],
        7: [1],        
        }

Q = np.zeros((nlat, nlon))    # stil anthropogenic
for s, vs in Qfiles_H2.items():
    for v in vs:
        with xr.open_dataset(v[0]) as ds_read:
            with ds_read.load() as Q_in:
                t_Q = Q_in['time']
                if v[2] == '1Y':
                    t = np.datetime64(
                        dates_tHour[0].floor('d').replace(month=1, day=1)
                        )
                    t_in = min(t_Q, key=lambda x: abs(x - t))
                else:
                    t = np.datetime64(
                        dates_tHour[0].floor('d').replace(day=1)
                        )
                    t_in = min(
                        t_Q[t_Q.dt.month==dates_tHour[0].month],
                        key=lambda x: abs(x - t)
                        )
                Q += Q_in[v[1]].sel(time=t_in).values * Q_factor[s] * 1.e3  # kg -> g
                
# Scenario emission Qs

Qfiles_H2 = OrderedDict([
    (0, [
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_AWB.nc', 'CO_emissions', '1M'],
        ]),
    (1, [
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_ENE.nc', 'CO_emissions', '1M'],   
        ]),
    (2, [
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_REF.nc', 'CO_emissions', '1M'],        
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_IND.nc', 'CO_emissions', '1M'],    
        ]),
    (3, [
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_CDS.nc', 'CO_emissions', '1M'],        
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_CRS.nc', 'CO_emissions', '1M'],        
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_TRO.nc', 'CO_emissions', '1M'],        
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_LTO.nc', 'CO_emissions', '1M'], 
        ]),
    (4, [
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_RCO.nc', 'CO_emissions', '1M'],        
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_PRO.nc', 'CO_emissions', '1M'], 
        ]),
    (5, [
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_FFF.nc', 'CO_emissions', '1M'], 
        ]),
    (6, [
        ['inputs/emissions/CO_emissions_processed/prior_edgar_v6_0_co_SWD.nc', 'CO_emissions', '1M'],        
        ]),
    (7, [
        ['inputs/emissions/biomass/gfed_2015.nc', 'H2_emissions', '1M'],        
        ]),    
    (8, [
        ['inputs/emissions/prior_edgar_v6_0_PRO_GAS.nc', 'CH4_emissions', '1M'],        
        ]), 
    ])
    

Q_factor = {
        0: [0.0357],      
        1: [0.0143],
        2: [0.0143],
        3: [0.0357],
        4: [0.0217], 
        5: [0.0143],
        6: [0.005],
        7: [1],        
        8: [0.368],        
        }

Qs = np.zeros((nlat, nlon))    # stil anthropogenic
for s, vs in Qfiles_H2.items():
    for v in vs:
        with xr.open_dataset(v[0]) as ds_read:
            with ds_read.load() as Q_in:
                t_Q = Q_in['time']
                if v[2] == '1Y':
                    t = np.datetime64(
                        dates_tHour[0].floor('d').replace(month=1, day=1)
                        )
                    t_in = min(t_Q, key=lambda x: abs(x - t))
                else:
                    t = np.datetime64(
                        dates_tHour[0].floor('d').replace(day=1)
                        )
                    t_in = min(
                        t_Q[t_Q.dt.month==dates_tHour[0].month],
                        key=lambda x: abs(x - t)
                        )
                Qs += Q_in[v[1]].sel(time=t_in).values * Q_factor[s] * 1.e3  # kg -> g

# Dilution matrix - CH4
Dfile_CH4 = (
    'inputs/footprints_wao/'
    f'WAO-20magl_UKV_EUROPE_{date}.nc'
    )
with xr.open_dataset(Dfile_CH4) as ds_read:
    with ds_read.load() as Din:
        D = Din.fp.transpose('time', 'lat', 'lon').values

        
# baseline
chi0 = chem_co.read_baseline(dates_tHour, "MHD-10magl")


# concentration H2

'''
conversion documentation:
    1 kg CH4 to 1 kg H2: 0.368
    H2 demand (by 2018 natural gas consumption in 2018):
        2035 -> 0.147
        2050 -> 0.687
        ideal -> 1
    leakage: 
        1% -> 0.01
        10% -> 0.1
        50% -> 0.5
'''
# prior
mod_H2 = pd.Series(
    523.024 + (D * Q).sum((1, 2)) / M_H2 * 1e9,
    index=dates_tHour
    )
'''
(D * Q.sum(0)).sum((1,  2))

mod_H2 = pd.Series(
    523.024 + (D * Q.sum(0)).sum((1,  2)) / M_H2 * 1e9,
    index=dates_tHour
    )
'''

#scenario
mod_H2_s = pd.Series(
    523.024 + (D * Qs).sum((1, 2)) / M_H2 * 1e9,
    index=dates_tHour
    )

# observations
obs_H2, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "WAO")   # could add st dev

# to save 

pd.concat([
    pd.Series(mod_H2, index=dates_tHour, name='chi_H2'),

    ], axis=1).to_csv(os.path.join(odir, f'{date}H2_emission_prior_withprogas.csv'))

# Plots

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
                {'c': colours['obs'], 'lw': 0.5, 'label': 'Measured'}
                ],
            'mod': [
                'line',
                [mod_H2.index, np.array(mod_H2), '-'],
                {'c': colours['mo1'], 'lw': 0.5, 'label': 'Modelled'}
                ],
            'sce': [
                'line',
                [mod_H2_s.index, np.array(mod_H2_s), '-'],
                {'c': colours['mo2'], 'lw': 0.5, 'label': 'Modelled Progas Scenario'}
                ],        
            'bas': [
                'line',
                [chi0.index, np.array(chi0), '-'],
                {'c': colours['bas'], 'lw': 0.5, 'label': 'Baseline'}
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
            pd.to_datetime(date) + pd.DateOffset(months=1), 
            ],
        ylim=(
            [450., 650.]
            ),
        yticks=(
            np.arange(450., 600., 50.)
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
            mdates.WeekdayLocator(byweekday=6),
            ]
        )
    for l in ax['main'].get_xticklabels():
        l.set_ha("right")
        l.set_rotation(30)
    ax['main'].legend(
        loc='upper right', ncol=4, fontsize=fig_param['fontsize']
        )
    fig['main'].savefig(os.path.join(odir, f'emissions_co_biomass{i}_progas.png'))



