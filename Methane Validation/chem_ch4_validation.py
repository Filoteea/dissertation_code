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


def read_Q(dates_tHour):
    grid_info = routines.define_grid()
    nlat = grid_info['nlat']
    nlon = grid_info['nlon']
    Qfiles_CH4 = OrderedDict([
        (0, [
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_ENE.nc', 'CH4_emissions', '1M'],
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_FFF.nc', 'CH4_emissions', '1M'],
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_REF_TRF.nc', 'CH4_emissions', '1M'],
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_RCO.nc', 'CH4_emissions', '1M'],
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_IND.nc', 'CH4_emissions', '1M'],
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_CHE.nc', 'CH4_emissions', '1M'],
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_IRO.nc', 'CH4_emissions', '1M'],
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_PRO.nc', 'CH4_emissions', '1M'],
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_TRO_noRES.nc', 'CH4_emissions', '1M'],
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_TNR_Aviation_CDS.nc', 'CH4_emissions', '1M'],
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_TNR_Aviation_CRS.nc', 'CH4_emissions', '1M'],
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_TNR_Other.nc', 'CH4_emissions', '1M'],
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_TNR_Ship.nc', 'CH4_emissions', '1M']
            ]),
        (1, [
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_SWD_LDF.nc', 'CH4_emissions', '1M'],
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_SWD_INC.nc', 'CH4_emissions', '1M'],
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_WWT.nc', 'CH4_emissions', '1M']
            ]),
        (2, [
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_ENF.nc', 'CH4_emissions', '1M'],
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_MNM.nc', 'CH4_emissions', '1M'],
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_AWB.nc', 'CH4_emissions', '1M'],
            ['inputs/prior_edgar_v6_0_ch4/prior_edgar_v6_0_AGS.nc', 'CH4_emissions', '1M']
            ]),
        (3, [
            ['inputs/prior_edgar_v6_0_ch4/prior_gfed.nc', 'CH4_emissions', '1M'],
            ]),
        (4, [
            ['inputs/prior_edgar_v6_0_ch4/prior_wetcharts.nc', 'CH4_emissions', '1M']
            ]),
        ])
    
    Q = np.zeros((nlat, nlon))    
    for s, vs in Qfiles_CH4.items():
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
                    Q += Q_in[v[1]].sel(time=t_in).values * 1.e3  # kg -> g
    return Q


def r_decc(fpath):
    odata = pd.read_csv(
                fpath,
                usecols=lambda x: x.lower() in ['time', 'ch4_c'],
                index_col=['time'],
                skipinitialspace=True,
                parse_dates=['time']
                ).dropna()
    odata.columns = odata.columns.str.lower()
    return odata


def read_obs(timestamps, site, resample='1H'):
    date = timestamps[0].strftime('%Y-%m')
    t0 = timestamps[0].strftime('%Y-%m-%d %H')
    t1 = timestamps[-1].strftime('%Y-%m-%d %H')
    if site == 'mod':
        ifile = 'outputs/validation_ch4/ABC.csv'
        col_or_no = 'ch4_c'
        sigma_col_or_no = 0.2
    elif site == 'hfd':
        ifile = 'inputs/obs/hfd-g2401_10.csv'
        col_or_no = 'ch4_c'
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
