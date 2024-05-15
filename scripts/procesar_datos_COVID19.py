import pandas as pd
import polars as pl
import xarray as xr
import pyarrow
import netCDF4

#READ THE DATAFILE

fname = "scripts/data/casos_hosp_uci_def_sexo_edad_provres.csv"
data = pl.read_csv(fname, separator = ",")
df_pl = data

#DATA PROCESSING

"""
We drop the sex column, because in the simulations we don't 
have information about it. 
We regroup the columns
"""

df_pl = df_pl.drop('sexo')
df_pl = df_pl.group_by(['provincia_iso','grupo_edad','fecha']).agg(pl.col(['num_casos','num_hosp','num_uci','num_def']).sum())

"""
We group by age strata:
    - Young people:  (0-9, 10-19)
    - Adults:        (20-29, 30-39, 40-49, 50-59)
    - Elderly people:(60-69, 70-79, 80+)
"""

df_pl = df_pl.with_columns(
                      pl.when(pl.col('grupo_edad')=='0-9').then(pl.lit('Y'))
                     .when(pl.col('grupo_edad')=='10-19').then(pl.lit('Y'))
                     .when(pl.col('grupo_edad')=='20-29').then(pl.lit('M'))
                     .when(pl.col('grupo_edad')=='30-39').then(pl.lit('M'))
                     .when(pl.col('grupo_edad')=='40-49').then(pl.lit('M'))
                     .when(pl.col('grupo_edad')=='50-59').then(pl.lit('M'))
                     .when(pl.col('grupo_edad')=='60-69').then(pl.lit('O'))
                     .when(pl.col('grupo_edad')=='70-79').then(pl.lit('O'))
                     .when(pl.col('grupo_edad')=='80+').then(pl.lit('O'))
                     .otherwise(pl.col('grupo_edad')).alias('grupo_edad')
                     )
df_pl = df_pl.group_by(['provincia_iso','grupo_edad','fecha']).agg(pl.col(['num_casos','num_hosp','num_uci','num_def']).sum())

"""
The porcentage of NC data is very low
We drop all the rows with NC data
"""

df_pl = df_pl.filter((pl.col('grupo_edad') != 'NC')&(pl.col('provincia_iso') != 'NC'))
df_pl = df_pl.drop('total_casos')

#DATA STANDARIZATION

df = df_pl.to_pandas()

"""
In the simulations the provinces are identified by the 
INE code
We change from ISO to INE code
"""

#We set as the index the provinces codes
df = df.set_index(['provincia_iso'])

#We read the data that links the codes with the provinces
#and we use it as a dictionary
path1     = "scripts/data/id_provincias_INE-ISO.csv"
cod_provincias = pd.read_csv(path1, dtype = {'codigo_INE': 'str','codigo_ISO': 'str', 'distrito_mitma': 'str'},keep_default_na = False)
cod_provincias = cod_provincias.drop(['Provincia','distrito_mitma'],axis = 1)
cod_provincias = cod_provincias.set_index(['codigo_ISO'])


#Change the province code from ISO to INE

df = df.join(cod_provincias, how = 'outer')



df = df.reset_index(drop=True)

"""
We rename the columns so they coincide with the results of the simulations
"""
df['H'] = df.loc[:,['num_hosp','num_uci']].sum(axis=1)
df = df.drop(['num_hosp', 'num_uci'], axis = 1)
df = df.rename(columns={"fecha": "T", "codigo_INE": "M", "grupo_edad":"G", "num_casos":"I", "num_def":"D"})

#CREATE xarray

df = df.set_index(['M', 'G', 'T'])
xa = df.to_xarray().to_array()
xa = xa.rename({'variable': 'epi_states'})
xa = xa.transpose("M", "G", "T", "epi_states")

#SAVE THE XARRAY

xa.to_netcdf('data/casos_hosp_def_edad_provres.nc')





