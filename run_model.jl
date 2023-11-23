#using DelimitedFiles


using Printf
using ArgParse
using Logging

using Dates
using CSV
using NPZ
using JSON
using HDF5
using DataStructures
using DelimitedFiles
using DataFrames



base_folder = ""

include(joinpath(base_folder, "MMCAcovid19_vac/markov_vac_aux.jl"))
include(joinpath(base_folder, "MMCAcovid19_vac/markov_vac_io.jl"))
include(joinpath(base_folder, "MMCAcovid19_vac/markov_vac.jl"))


####################################################
##############   FUNCTIONS     #####################
####################################################

function set_compartments!(epi_params, initial_compartments)
    @assert size(initial_compartments) == (size(epi_params.ρˢᵍᵥ)[1], size(epi_params.ρˢᵍᵥ)[2], size(epi_params.ρˢᵍᵥ)[4], 10)
    total_population = sum(initial_compartments, dims=(3))[:,:,1]
    
    # Index of the initial condition
    T0 = 1
    
    epi_params.ρˢᵍᵥ[:,:,T0,:]  .= initial_compartments[:, :, T0, :, 1] ./ total_population
    epi_params.ρᴱᵍᵥ[:,:,T0,:]  .= initial_compartments[:, :, T0, :, 2] ./ total_population
    epi_params.ρᴬᵍᵥ[:,:,T0,:]  .= initial_compartments[:, :, T0, :, 3] ./ total_population
    epi_params.ρᴵᵍᵥ[:,:,T0,:]  .= initial_compartments[:, :, T0, :, 4] ./ total_population
    epi_params.ρᴾᴴᵍᵥ[:,:,T0,:] .= initial_compartments[:, :, T0, :, 5] ./ total_population
    epi_params.ρᴾᴰᵍᵥ[:,:,T0,:] .= initial_compartments[:, :, T0, :, 6] ./ total_population
    epi_params.ρᴴᴿᵍᵥ[:,:,T0,:] .= initial_compartments[:, :, T0, :, 7] ./ total_population
    epi_params.ρᴴᴰᵍᵥ[:,:,T0,:] .= initial_compartments[:, :, T0, :, 8] ./ total_population
    epi_params.ρᴿᵍᵥ[:,:,T0,:]  .= initial_compartments[:, :, T0, :, 9] ./ total_population
    epi_params.ρᴰᵍᵥ[:,:,T0,:]  .= initial_compartments[:, :, T0, :, 10] ./ total_population

    epi_params.ρˢᵍᵥ[isnan.(epi_params.ρˢᵍᵥ)] .= 0
    epi_params.ρᴱᵍᵥ[isnan.(epi_params.ρᴱᵍᵥ)] .= 0
    epi_params.ρᴬᵍᵥ[isnan.(epi_params.ρᴬᵍᵥ)] .= 0
    epi_params.ρᴵᵍᵥ[isnan.(epi_params.ρᴵᵍᵥ)] .= 0
    epi_params.ρᴾᴴᵍᵥ[isnan.(epi_params.ρᴾᴴᵍᵥ)] .= 0
    epi_params.ρᴾᴰᵍᵥ[isnan.(epi_params.ρᴾᴰᵍᵥ)] .= 0
    epi_params.ρᴴᴿᵍᵥ[isnan.(epi_params.ρᴴᴿᵍᵥ)] .= 0
    epi_params.ρᴴᴰᵍᵥ[isnan.(epi_params.ρᴴᴰᵍᵥ)] .= 0
    epi_params.ρᴿᵍᵥ[isnan.(epi_params.ρᴿᵍᵥ)] .= 0
    epi_params.ρᴰᵍᵥ[isnan.(epi_params.ρᴰᵍᵥ)] .= 0
end

###########################################
############# FILE READING ################
###########################################

args = parse_commandline()

config_fname  = args["config"]
data_path     = args["data-folder"]
instance_path = args["instance-folder"]

config = JSON.parsefile(config_fname);
update_config!(config, args)

# Output simulation
output_path = joinpath(instance_path, "output")
if !isdir(output_path)
    println("Creating output folder: $output_path")
    mkpath(output_path)
end

# Reading simulation start and end dates
first_day = Date(config["simulation"]["first_day_simulation"])
last_day  = Date(config["simulation"]["last_day_simulation"])
# Converting dates to time steps
T = (last_day - first_day).value + 1

T_coords = string.(collect(first_day:last_day))

A0_instance_filename = get(config["simulation"], "A0_filename", nothing)
A0_instance_filename = joinpath(instance_path, A0_instance_filename)
initial_compartments_path = get(config["simulation"], "initial_compartments", nothing)

if A0_instance_filename !== nothing && initial_compartments_path !== nothing
    println("ERROR!!!")
end


#########################
# Simulation output 
#########################

export_compartments_full = get(config["simulation"], "export_compartments_full", false)
export_compartments_time_t = get(config["simulation"], "export_compartments_time_t", nothing)

println("first_day_simulation = ", first_day)
println("last_day_simulation = ", last_day)
println("export_compartments_full = ", export_compartments_full)
println("export_compartments_time_t = ", export_compartments_time_t)
println("initial_compartments = ", initial_compartments_path)

########################################
####### VARIABLES INITIALIZATION #######
########################################

data_dict        = config["data"]
epi_params_dict  = config["epidemic_params"]
pop_params_dict  = config["population_params"]
vac_params_dict  = config["vaccination"]
npi_params_dict  = config["NPI"]

# Population info

# Loading metapopulation patches info (surface, label, population by age)
metapop_data_filename    = joinpath(data_path, data_dict["metapopulation_data_filename"])
metapop_df = CSV.read(metapop_data_filename, DataFrame)

# Loading mobility network
mobility_matrix_filename = joinpath(data_path, data_dict["mobility_matrix_filename"])
network_df  = CSV.read(mobility_matrix_filename, DataFrame)

# Metapopulations patches coordinates (labels)
M_coords = map(String,metapop_df[:, "id"])
M = length(M_coords)

# Coordinates for each age strata (labels)
G_coords = map(String, pop_params_dict["age_labels"])
G = length(G_coords)

# Num. of vaccination statuses Vaccinated/Non-vaccinated
V = length(epi_params_dict["kᵥ"])

##################################################
####### INITIALIZATION OF THE EPIDEMICS ##########
##################################################

## POPULATION PARAMETERS
population = init_pop_param_struct(G, M, G_coords, 
                                   pop_params_dict, 
                                   metapop_df, network_df)
## EPIDEMIC PARAMETERS 
epi_params = init_epi_parameters_struct(G, M, G_coords, epi_params_dict)


total_population = sum(population.nᵢᵍ)

#########################################################
# Vaccination parameters
#########################################################

# vaccionation dates
start_vacc = vac_params_dict["start_vacc"]
dur_vacc   = vac_params_dict["dur_vacc"]
end_vacc   = start_vacc + dur_vacc

# total vaccinations per age strata
ϵᵍ = vac_params_dict["ϵᵍ"] * round( total_population * vac_params_dict["percentage_of_vacc_per_day"] )

tᵛs = [start_vacc, end_vacc, T]
ϵᵍs = ϵᵍ .* [0  Int(vac_params_dict["are_there_vaccines"])  0] 

#########################################################
# Containment measures
#########################################################

# Mobility reduction
κ₀_df = CSV.read(joinpath(data_path, config["data"]["kappa0_filename"]), DataFrame);

# syncronize containment measures with simulation
κ₀_df.time = map(x -> (x .- first_day).value + 1, κ₀_df.date)

# Timesteps when the containment measures will be applied
tᶜs = Int64.(npi_params_dict["tᶜs"])

# Array of level of confinement
# κ₀s = κ₀_df.reduction[:]
κ₀s = Float64.(npi_params_dict["κ₀s"])
# Array of premeabilities of confined households
ϕs = Float64.(npi_params_dict["ϕs"])
# Array of social distancing measures
δs = Float64.(npi_params_dict["δs"])




println("M = ", M)
println("G = ", G)
println("T = ", T)
println("V = ", V)
println("N. of epi compartments = ", epi_params.NumComps)

# vac_parms = Vaccination_Params(tᵛs, ϵᵍs)
# npi_params = NPI_Params(tᶜs, κ₀s, ϕs, δs)
# run_epidemic_spreading_mmca!(epi_params, population, npi_params, vac_parms; verbose = true )

# Initial seeds (intial condition at the begining of the pandemic)
# Load initial full conditions
if initial_compartments_path !== nothing
    # use initial compartments matrix to initialize simulations
    initial_compartments = h5open(initial_compartments_path, "r") do file
        read(file, "compartments")
    end
    # set the full initial condition o a user defined
    set_compartments!(epi_params, initial_compartments)
else
    Sᵛ₀ = zeros(Float64, G, M)
    E₀  = zeros(Float64, G, M)
    A₀  = zeros(Float64, G, M)
    I₀  = zeros(Float64, G, M)
    H₀  = zeros(Float64, G, M)
    R₀  = zeros(Float64, G, M)
    if A0_instance_filename !== nothing
        # Initial number of infectious asymptomatic individuals
        # use seeds to initialize simulations
        conditions₀ = CSV.read(A0_instance_filename, DataFrame)        
        A₀[1, Int.(conditions₀[:,"idx"])] .= 0.12 .* conditions₀[:,"seed"]
        A₀[2, Int.(conditions₀[:,"idx"])] .= 0.16 .* conditions₀[:,"seed"]
        A₀[3, Int.(conditions₀[:,"idx"])] .= 0.72 .* conditions₀[:,"seed"]    
    else
        # Initial set custom number of infected
        E₀ = nᵢᵍ / total_population * 1000
        A₀ = nᵢᵍ / total_population * 1000
        I₀ = nᵢᵍ / total_population * 1000    
    end
    set_initial_conditions!(epi_params, population, Sᵛ₀, E₀, A₀, I₀, H₀, R₀)
end

########################################################
################ RUN THE SIMULATION ####################
########################################################

run_epidemic_spreading_mmca!(epi_params, population, tᶜs, tᵛs, κ₀s, ϕs, δs, ϵᵍs; verbose = true )

##############################################################
################## STORING THE RESULTS #######################
##############################################################

if export_compartments_full
    filename = joinpath(output_path, "compartments_full.h5")
    println("Storing full simulation output")
    println("\t- filename: $(filename)")
    save_simulation_hdf5(epi_params, population, filename)
    filename = joinpath(output_path, "compartments_full.nc")
    save_simulation_netCDF(epi_params, population, filename;G_coords=G_coords, M_coords=M_coords, T_coords=T_coords)
end

if export_compartments_time_t !== nothing
    export_compartments_date = first_day + Day(export_compartments_time_t - 1)
    filename = joinpath(output_path, "compartments_t_$(export_compartments_date).h5")
    println("Storing compartments at single date $(export_compartments_date):")
    println("\t- Simulation step: $(export_compartments_time_t)")
    println("\t- filename: $(filename)")
    save_simulation_hdf5(epi_params, population, filename; 
                         export_time_t = export_compartments_time_t)
end
