import numpy as np
import xarray as xr
import pandas as pd
import polars as pl
import json
import matplotlib.pyplot as plt
import glob

#LOAD FILES

test_folder = 'test'
instance_id = 'lockdown_escenarios'

simulation_fname = f"../{test_folder}/{instance_id}/output/compartments_full.nc"
config_fname =f"../{test_folder}/{instance_id}/config_prueba.json"
sim_xa = xr.load_dataarray(simulation_fname)
#sim_xa = sim_xa.transpose()
with open(config_fname) as fh:
    config_dict = json.load(fh)
    
#DEFINE FUNCTION

def post_processing(sim_xa, config_dict):
    #We create a dictionary
    G_dic = sim_xa.coords['G'].values
    M_dic = sim_xa.coords['M'].values

    age_labels = config_dict['population_params']['age_labels']

    epi_params = config_dict['epidemic_params']
    alphas = {}
    h = np.array(epi_params['μᵍ'])*(1-np.array(epi_params['θᵍ']))*np.array(epi_params['γᵍ'])
    hosp_rates = {}
    for i,g in enumerate(age_labels):
        alphas[g] = epi_params['αᵍ'][i]
        hosp_rates[g] = h[i]

    alpha_mean = np.mean(epi_params['αᵍ'])
    hosp_rate_mean = h.mean()

    observable_labels = ['A', 'I', 'D']
    sim_observables_xa = sim_xa.loc[:,:,:,observable_labels]

    alphas_vec = np.array([alphas[g] for g in sim_xa.coords['G'].values])
    hosp_rates_vec = np.array([hosp_rates[g] for g in sim_xa.coords['G'].values]) 

    sim_observables_xa = sim_observables_xa.transpose('epi_states', 'M', 'T', 'G')

    sim_observables_xa.loc['A', :, :, :] = np.multiply(sim_observables_xa.loc['A', :, :, :], alphas_vec)
    sim_observables_xa.loc['I', :, :, :] = np.multiply(sim_observables_xa.loc['I', :, :, :], hosp_rates_vec)

    sim_observables_xa = sim_observables_xa.transpose('epi_states', 'G', 'M', 'T')
    dims = sim_observables_xa.dims
    data = sim_observables_xa.values
    coords = sim_observables_xa.coords

    coords['epi_states'] = ['I', 'H', 'D']
    sim_observables_xa = xr.DataArray(data, coords=coords, dims=dims)
    sim_observables_xa.loc['D',:] = np.diff(sim_observables_xa.loc['D',:], prepend =0)
    return sim_observables_xa



#AGGREGATE

agg_sim_xa = sim_xa.sum(['V'])

#CONVERT THE XARRAY TO A DATAFRAME

df_sim = agg_sim_xa.to_dataframe()
df_sim = df_sim.reset_index(["T", "G", "M", "epi_states"])
df_sim = df_sim.set_index(["M"])

#LOAD THE LIST OF MOBILITY AREAS WITH THE PROVINCE TO WHICH THEY BELONG

path     = 'data/patch_agg.csv'
patch_agg = pd.read_csv(path, dtype = {'id_map': 'str'},keep_default_na = False)
patch_agg = patch_agg.set_index(['id'])

#GROUP THE SIMULATION DATA BY PROVINCES

df_sim = df_sim.join(patch_agg, how = 'outer')
df_sim = df_sim.reset_index()
df_sim = df_sim.drop(['M'], axis=1)
df_sim = df_sim.rename(columns = {'id_map': 'M'})
df_sim = df_sim.groupby(['T','G','epi_states','M']).sum()
df_sim = df_sim.unstack(level='epi_states')
df_sim = df_sim.droplevel(level = 0, axis=1)


#CONVERT THE DATAFRAME IN A XARRAY

sim_xa = df_sim.to_xarray().to_array()
sim_xa = sim_xa.rename({'variable': 'epi_states'})
sim_xa = sim_xa.transpose("G","M","T","epi_states")

xa = post_processing(sim_xa, config_dict)

#SAVE THE XARRAY IN A NETCDF FILE

xa.to_netcdf('test/no_vaccination/output/IHD_agg.nc')






