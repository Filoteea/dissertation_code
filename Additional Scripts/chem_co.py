# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 00:10:23 2022

@author: Filoteea Moldovan
Adapted from Edward Chung

"""

r"""

Forward Model

"""
# Standard Library imports
import gzip
import numpy as np
import pandas as pd
import xarray as xr

# Third party imports
from collections import OrderedDict

# Semi-local imports
import name_qch4_couple.io

# Local imports
import routines



# Function used to create the emissions map
# For Section 2.9.: factor_q = a (Table B.1); factor_s = b (Table B.1); factor_sce = 0
# For Section 2.11.: factor_q = 1; factos_s = 1; factor_sce - Table 2.4

def read_Qsink(dates_tHour, factor_q, factor_s, factor_sce):
    grid_info = routines.define_grid()
    nlat = grid_info['nlat']
    nlon = grid_info['nlon']

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
            ['inputs/emissions/biomass/gfed_2012.nc', 'H2_emissions', '1M'],        
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
            8: [factor_sce],        
            }
    Q = np.zeros((nlat, nlon))    
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
                    Q += Q_in[v[1]].sel(time=t_in).values * Q_factor[s] * 1.e3 * factor_q # kg -> g
    lwfile = 'inputs/sink/land_mask.nc'
    with xr.open_dataset(lwfile) as ds_read:
        lwin = ds_read.load()
        soil_sink = np.array(lwin.lo_land)/100 * -0.000000005 * factor_s
    Q += soil_sink 
    return Q



def r_decc(fpath):
    odata = pd.read_csv(
                fpath,
                usecols=lambda x: x.lower() in ['time', 'h2_ppb'],
                index_col=['time'],
                skipinitialspace=True,
                parse_dates=['time']
                ).dropna()
    odata.columns = odata.columns.str.lower()
    return odata


def read_obs(timestamps, site, factor, resample='1H'):
    date = timestamps[0].strftime('%Y-%m')
    t0 = timestamps[0].strftime('%Y-%m-%d %H')
    t1 = timestamps[-1].strftime('%Y-%m-%d %H')
    if site == 'WAO':
        ifile = 'inputs/obs/WAO_H2_oct2021.csv'
        col_or_no = 'h2_ppb'
        sigma_col_or_no = 0.2
    elif site == 'MHD_10magl':
        ifile = 'inputs/baseline/MHD_2018.csv'
        col_or_no = 'h2_ppb'
        sigma_col_or_no = 0.2
    elif site == 'bas':
        ifile = 'outputs/models/mhd_bas/chi0_proc.csv'
        col_or_no = 'chi0p_H2'
        sigma_col_or_no = 0.2
    elif site == 'mod':
        ifile = f'outputs/scenarios/new_merged/merged_wao_scenario_{factor}.csv'
        col_or_no = 'h2_ppb'
        sigma_col_or_no = 0.2
    elif site == 'bas_mhd':
        ifile = 'outputs/scenarios/merged/merged_bas_mhd.csv'
        col_or_no = 'h2_ppb'
        sigma_col_or_no = 0.2
    elif site == 'bas_wao':
        ifile = 'outputs/scenarios/merged/merged_bas_wao.csv'
        col_or_no = 'h2_ppb'
        sigma_col_or_no = 0.2
    else:
        ifile = False
        col_or_no = np.nan
        sigma_col_or_no = np.nan
    if ifile:
        all_obs_raw = r_decc(ifile).sort_index().loc[t0:t1]
        obs_raw = all_obs_raw[col_or_no]
        sigma_obs_raw = (all_obs_raw[sigma_col_or_no]
                         if isinstance(sigma_col_or_no, str) else
                         pd.Series(sigma_col_or_no, index=all_obs_raw.index))
    if isinstance(col_or_no, str):
        obs = (obs_raw
               if resample is False else
               obs_raw.resample('1H').mean().reindex(timestamps))
    else:
        obs = pd.Series(col_or_no, index=timestamps)
    if isinstance(sigma_col_or_no, str) or isinstance(col_or_no, str):
        sigma_obs = (
            sigma_obs_raw
            if resample is False else
            sigma_obs_raw.resample('1H').apply(
                lambda x: np.sum(x**2)).reindex(timestamps))
    else:
        sigma_obs = pd.Series(sigma_col_or_no, index=timestamps)
    return obs, sigma_obs 


def read_baseline(timestamps, site, btype="default"):
    date = timestamps[0].strftime('%Y-%m')
    year = timestamps[0].strftime('%Y')
    if site == 'MHD_10magl':
        if btype == 'default':
            chi0file = (
                'outputs/baseline/baseline-MHD_10magl-h2-2018.nc'
                )
            with xr.open_dataset(chi0file) as ds_read:     #put as
                with ds_read.load() as ds:
                    chi0 = ds.chi_H2.sel(time=date).to_series()
                    var_chi0 = ds.var_chi_H2.sel(time=date).to_series()
        elif btype == 'intem':
            if timestamps[0] < pd.to_datetime('2020-07'):
                bmonth = '2020-07'
                bflag = 1
            elif timestamps[0] > pd.to_datetime('2020-12'):
                bmonth = '2020-12'
                bflag = 2
            else:
                bmonth = date
                bflag = 0
            m_sta = (pd.to_datetime(bmonth)).date().strftime('%Y%m%d')
            m_end = (
                pd.to_datetime(bmonth)+pd.tseries.offsets.MonthEnd(0)
                ).date().strftime('%Y%m%d')
            chi0file = (
                '/home/ec5/hpc-work/data_archive/decc/'
                'EUROPE_UKV_HFD_100magl/Pos_CH4/'
                'H1_C_MHT1T2R1ANH1B2CBWB_ch4_OBUSEXL_4h_Fnc10_'
                f'{m_sta}-{m_end}_average_f.gz'
                )
            with gzip.open(chi0file, mode='r') as chi0in:
                chi0_0all = pd.read_csv(
                    chi0in, sep=' ', skipinitialspace=True,
                    skiprows=[5], header=4,
                    parse_dates={'datetime': ['YYYY', 'MM', 'DD', 'HH', 'MI']},
                    #converters={
                    #    'datetime': lambda Y, m, d, H, M:
                    #        pd.to_datetime(f'{Y} {m} {d} {H} {M}',
                    #                       format='%Y %m %d %H %M'),
                    #    },
                    #index_col='datetime'
                    )
                chi0_0all.index = pd.to_datetime(
                        chi0_0all['datetime'], format='%Y %m %d %H %M')
                chi0_0 = chi0_0all['BasePos']
                var_chi0_0 = chi0_0all['BasePos']
            if bflag == 1:
                chi0 = pd.Series(chi0_0.iloc[-1], index=timestamps)
                var_chi0 = pd.Series(var_chi0_0.iloc[-1], index=timestamps)
            elif bflag == 2:
                chi0 = pd.Series(chi0_0.iloc[0], index=timestamps)
                var_chi0 = pd.Series(var_chi0_0.iloc[0], index=timestamps)
            else:
                chi0 = chi0_0.reindex(timestamps).interpolate(
                    method='time', limit_direction='both'
                    )
                var_chi0 = var_chi0_0.reindex(timestamps).interpolate(
                    method='time', limit_direction='both'
                    )
    else:
        chi0 = pd.Series(1900., index=timestamps)
        var_chi0 = pd.Series(0.1, index=timestamps)
    return chi0

