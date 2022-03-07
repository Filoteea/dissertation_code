
"""
Created on Wed Jan 26 13:47:47 2022

@author: filoteea
Adapted from Edward Chung
"""

r"""

H2 Forward Model verification

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


# Prior emission
Q = chem_co.read_Q(dates_tHour)

# Soil sink
Q_sink = chem_co.read_Qsink(dates_tHour)
              
# Scenario emission Qs

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
factor = 0.368

Qs = chem_co.read_Qs(dates_tHour, factor)


# Dilution matrix - CH4
Dfile_CH4 = (
    'inputs/footprints_wao/'
    f'WAO-20magl_UKV_EUROPE_{date}.nc'
    )
with xr.open_dataset(Dfile_CH4) as ds_read:
    with ds_read.load() as Din:
        D = Din.fp.transpose('time', 'lat', 'lon').values

        
# baseline

def read_baseline(timestamps):
    date = timestamps[0].strftime('%Y-%m')
    year = timestamps[0].strftime('%Y')
    chi0file = (
                'outputs/baseline/baseline-MHD_10magl-h2-2018.nc'
                )
    with xr.open_dataset(chi0file) as ds_read:     #put as
        with ds_read.load() as ds:
                chi0 = ds.chi_H2.sel(time=date).to_series()
                
    return chi0

chi0 = read_baseline(dates_tHour)


# prior
'''
mod_H2 = pd.Series(
    523.024 + (D * Q).sum((1, 2)) / M_H2 * 1e9,
    index=dates_tHour
    )
'''

mod_H2_bas = pd.Series(
    chi0 + (D * Q).sum((1, 2)) / M_H2 * 1e9,
    index=dates_tHour
    )

mod_sink = pd.Series(
    chi0 + (D * Q_sink).sum((1, 2)) / M_H2 * 1e9,
    index=dates_tHour
    )
'''
#scenario
mod_H2_s = pd.Series(
    chi0 + (D * Qs).sum((1, 2)) / M_H2 * 1e9,
    index=dates_tHour
    )
'''

# observations
obs_H2, sigma_obs_H2 = chem_co.read_obs(dates_tHour, "WAO")   # could add st dev

obs_H2_mhd, sigma_obs_H2_mhd = chem_co.read_obs(dates_tHour, "MHD_10magl")   # could add st dev

'''
# to save 
pd.concat([
    pd.Series(mod_H2_s, index=dates_tHour, name='chi_H2'),
    ], axis=1).to_csv(os.path.join(odir, f'{date}simplified.csv'))
'''

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
            'wao': [
                'line',
                [obs_H2.index, np.array(obs_H2), '-'],
                {'c': colours['obs'], 'lw': 0.5, 'label': 'Measured WAO'}
                ],
            'mhd': [
                'line',
                [obs_H2_mhd.index, np.array(obs_H2_mhd), '-'],
                {'c': '#0139FF', 'lw': 0.5, 'label': 'Measured MHD'}
                ],
            'obs2': [
                'line',
                [mod_H2_bas.index, np.array(mod_H2_bas), '-'],
                {'c': colours['bas'], 'lw': 0.5, 'label': 'Modelled with baseline'}
                ],
          #  'mod': [
          #      'line',
          #      [mod_H2.index, np.array(mod_H2), '-'],
          #      {'c': '#FF0000', 'lw': 0.5, 'label': 'Modelled'}
          #      ],
         #   'sce': [
         #       'line',
         #       [mod_H2_s.index, np.array(mod_H2_s), '-'],
         #       {'c': colours['mo2'], 'lw': 0.5, 'label': 'Modelled Progas Scenario'}
         #       ],        
            'bas': [
                'line',
                [chi0.index, chi0.values, '-'],
                {'c': colours['bas'], 'lw': 0.5, 'label': 'Baseline'}
                ], 
            'sink': [
                'line',
                [mod_sink.index, np.array(mod_sink), '-'],
                {'c': '#D201FF', 'lw': 0.5, 'label': 'Soil sink'}
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
            [460., 650.]
            ),
        yticks=(
            np.arange(460., 650., 20.)
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
    fig['main'].savefig(os.path.join(odir, f'soil sink{i}_{date}.png'))


