import h5py
import pandas as pd
import xarray as xa

with h5py.File('./test/output/compartments_full.h5','r') as file:
    nda = np.array(file.get('compartments'))
    nda = nda.T

df_patches = pd.read_csv("data/patch_ids.csv")

patch_coords = df_patches['id'].values
age_coords = ["Y", "M", "O"]
compartments_coords = ["S", "E", "A", "I", "PH", "PD", "HR", "HD", "R", "D"]
vaccine_coords = ["NV", "V"]
time_coords = range(nda.shape[2])

infected_compartments = ['A', 'I', 'PH', 'PD', 'HR', 'HD', 'D']

xr = xa.DataArray(nda, coords={'age':age_coords,'patch': patch_coords, 'date': time_coords, 'vacc': vaccine_coords, 'compartment': compartments_coords})
incidence_xr = xr.loc[:, :, :, :, infected_compartments].sum('compartment')
tot_incidence = incidence_xr.sum('age').sum('patch').sum('vacc')
tot_incidence.loc[1:]
tot_incidence.loc[1:] - tot_incidence[:-1]
tot_incidence.loc[1:] - tot_incidence[0:-1]
tot_incidence
tot_incidence.loc[1:] - tot_incidence.loc[0:-1]
tot_incidence.values[1:] - tot_incidence.values[0:-1]
