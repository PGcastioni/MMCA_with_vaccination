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

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--config", "-c"
            help = "config file (json file)"
            required = true
        "--data-folder", "-d"
            help = "data folder"
            required = true
        "--instance-folder", "-i"
            help = "instance folder (experiment folder)"
            required = true 
        "--export-compartments-full"
            help = "export compartments of simulations"
            action = :store_true
        "--export-compartments-time-t"
            help = "export compartments of simulations at a given time"
            default = nothing
            arg_type = Int
        "--initial-compartments"
            help = "compartments to initialize simulation. If missing, use the seeds to initialize the simulations"
            default = nothing
        "--start-date"
            help = "starting date of simulation. Overwrites the one provided in config.json"
            default = nothing
        "--end-date"
            help = "end date of simulation. Overwrites the one provided in config.json"
            default = nothing
    end

    return parse_args(s)
end

function set_compartments!(epi_params, population, 
                            initial_compartments; normalize=true)
    
    # Index of the initial condition
    t₀ = 1
    if normalize
        for i in 1:2
            epi_params.ρˢᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 1] ./ population.nᵢᵍ
            epi_params.ρᴱᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 2] ./ population.nᵢᵍ
            epi_params.ρᴬᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 3] ./ population.nᵢᵍ
            epi_params.ρᴵᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 4] ./ population.nᵢᵍ
            epi_params.ρᴾᴴᵍᵥ[:,:,t₀,i] .= initial_compartments[:, :, i, 5] ./ population.nᵢᵍ
            epi_params.ρᴾᴰᵍᵥ[:,:,t₀,i] .= initial_compartments[:, :, i, 6] ./ population.nᵢᵍ
            epi_params.ρᴴᴿᵍᵥ[:,:,t₀,i] .= initial_compartments[:, :, i, 7] ./ population.nᵢᵍ
            epi_params.ρᴴᴰᵍᵥ[:,:,t₀,i] .= initial_compartments[:, :, i, 8] ./ population.nᵢᵍ
            epi_params.ρᴿᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 9] ./ population.nᵢᵍ
            epi_params.ρᴰᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 10] ./ population.nᵢᵍ
        end
    else
        for i in 1:2
            epi_params.ρˢᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 1] 
            epi_params.ρᴱᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 2] 
            epi_params.ρᴬᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 3] 
            epi_params.ρᴵᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 4] 
            epi_params.ρᴾᴴᵍᵥ[:,:,t₀,i] .= initial_compartments[:, :, i, 5] 
            epi_params.ρᴾᴰᵍᵥ[:,:,t₀,i] .= initial_compartments[:, :, i, 6] 
            epi_params.ρᴴᴿᵍᵥ[:,:,t₀,i] .= initial_compartments[:, :, i, 7] 
            epi_params.ρᴴᴰᵍᵥ[:,:,t₀,i] .= initial_compartments[:, :, i, 8] 
            epi_params.ρᴿᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 9] 
            epi_params.ρᴰᵍᵥ[:,:,t₀,i]  .= initial_compartments[:, :, i, 10]
        end
    end

    epi_params.ρˢᵍᵥ[isnan.(epi_params.ρˢᵍᵥ)]   .= 0
    epi_params.ρᴱᵍᵥ[isnan.(epi_params.ρᴱᵍᵥ)]   .= 0
    epi_params.ρᴬᵍᵥ[isnan.(epi_params.ρᴬᵍᵥ)]   .= 0
    epi_params.ρᴵᵍᵥ[isnan.(epi_params.ρᴵᵍᵥ)]   .= 0
    epi_params.ρᴾᴴᵍᵥ[isnan.(epi_params.ρᴾᴴᵍᵥ)] .= 0
    epi_params.ρᴾᴰᵍᵥ[isnan.(epi_params.ρᴾᴰᵍᵥ)] .= 0
    epi_params.ρᴴᴿᵍᵥ[isnan.(epi_params.ρᴴᴿᵍᵥ)] .= 0
    epi_params.ρᴴᴰᵍᵥ[isnan.(epi_params.ρᴴᴰᵍᵥ)] .= 0
    epi_params.ρᴿᵍᵥ[isnan.(epi_params.ρᴿᵍᵥ)]   .= 0
    epi_params.ρᴰᵍᵥ[isnan.(epi_params.ρᴰᵍᵥ)]   .= 0
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

simulation_dict = config["simulation"]
data_dict       = config["data"]
epi_params_dict = config["epidemic_params"]
pop_params_dict = config["population_params"]
vac_params_dict = config["vaccination"]
npi_params_dict = config["NPI"]

#########################
# Simulation output 
#########################
output_path = joinpath(instance_path, "output")
if !isdir(output_path)
    println("Creating output folder: $output_path")
    mkpath(output_path)
end

output_format    = simulation_dict["output_format"]
save_full_output = get(simulation_dict, "save_full_output", false)
save_time_step   = get(simulation_dict, "save_time_step", nothing)
init_format      = get(simulation_dict, "init_format", "netcdf")
initial_compartments_path = get(data_dict, "initial_condition_filename", nothing)

#########################
# Initial Condition
#########################


if isnothing(initial_compartments_path)
    println("ERROR. Missing initial condition file")
end

########################################
####### VARIABLES INITIALIZATION #######
########################################

# Reading simulation start and end dates
first_day = Date(simulation_dict["first_day_simulation"])
last_day  = Date(simulation_dict["last_day_simulation"])
# Converting dates to time steps
T = (last_day - first_day).value + 1
# Array with time coordinates (dates)
T_coords  = string.(collect(first_day:last_day))

# Loading metapopulation patches info (surface, label, population by age)
metapop_data_filename = joinpath(data_path, data_dict["metapopulation_data_filename"])
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

####################################################
#####   INITIALIZATION OF DATA Structures   ########
####################################################

## POPULATION PARAMETERS
population       = init_pop_param_struct(G, M, G_coords, pop_params_dict, metapop_df, network_df)
total_population = sum(population.nᵢᵍ)

## EPIDEMIC PARAMETERS 
epi_params       = init_epi_parameters_struct(G, M, T, G_coords, epi_params_dict)

##################################################

println("- Initializing MMCA epidemic simulations")
println("\t- first_day_simulation = ", first_day)
println("\t- last_day_simulation = ", last_day)
println("\t- G (agent class) = ", G)
println("\t- M (n. of metapopulations) = ", M)
println("\t- T (simulation steps) = ", T)
println("\t- V (vaccination states) = ", V)
println("\t- N. of epi compartments = ", epi_params.NumComps)

# println("\t- Initial file = ", initial_compartments_path)
println("\t- Save full output = ", save_full_output)
if save_time_step !== nothing
    println("\t- Save time step at t=", save_time_step)
end

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

# Daily Mobility reduction
kappa0_filename = get(data_dict, "kappa0_filename", nothing)

if !isnothing(kappa0_filename)
    kappa0_filename = joinpath(data_path, kappa0_filename)
    @info "- Loading κ₀ time series from $(kappa0_filename)"
    κ₀_df = CSV.read(kappa0_filename, DataFrame);
    # syncronize containment measures with simulation
    @info "- Synchronizing to dates"
    κ₀_df.time = map(x -> (x .- first_day).value + 1, κ₀_df.date)
    # Timesteps when the containment measures will be applied
    tᶜs = κ₀_df.time[:]
    # Array of level of confinement
    κ₀s = κ₀_df.reduction[:]
    # Array of premeabilities of confined households
    ϕs = ones(Float64, length(tᶜs))
    # Array of social distancing measures
    δs = ones(Float64, length(tᶜs))
else
    # Timesteps when the containment measures will be applied
    tᶜs = Int64.(npi_params_dict["tᶜs"])
    # Array of level of confinement
    κ₀s = Float64.(npi_params_dict["κ₀s"])
    # Array of premeabilities of confined households
    ϕs = Float64.(npi_params_dict["ϕs"])
    # Array of social distancing measures
    δs = Float64.(npi_params_dict["δs"])
end


# vac_parms = Vaccination_Params(tᵛs, ϵᵍs)
# npi_params = NPI_Params(tᶜs, κ₀s, ϕs, δs)
# run_epidemic_spreading_mmca!(epi_params, population, npi_params, vac_parms; verbose = true )

# Initial seeds (intial condition at the begining of the pandemic)
# Load initial full conditions
if initial_compartments_path !== nothing
    # use initial compartments matrix to initialize simulations
    if init_format == "netcdf"
        initial_compartments = ncread(initial_compartments_path, "data")
    elseif init_format == "hdf5"
        initial_compartments = h5open(initial_compartments_path, "r") do file
            read(file, "data")
        end
    else
        @error "init_format must be one of : netcdf/hdf5"
    end
end


τ_inc    = 5.2
scale_ea = 0.5775
τᵢ       = 3.9131

# epi_params.ηᵍ .= 1.0/(τ_inc * (1.0 - scale_ea))
# epi_params.αᵍ .= [1.0/(τᵢ - 1 + τ_inc * scale_ea),
#                   1.0/(τ_inc * scale_ea),
#                   1.0/(τ_inc * scale_ea)]
# epi_params.μᵍ .= [1.0, 1.0/τᵢ, 1.0/τᵢ]

scale₀   = 0.8
initial_compartments[:, :, :, 2] .= initial_compartments[:, :, :, 2] * scale₀

@assert size(initial_compartments) == (G, M, V, epi_params.NumComps)
set_compartments!(epi_params, population, initial_compartments)

########################################################
################ RUN THE SIMULATION ####################
########################################################

run_epidemic_spreading_mmca!(epi_params, population, tᶜs, tᵛs, κ₀s, ϕs, δs, ϵᵍs; verbose = true )

##############################################################
################## STORING THE RESULTS #######################
##############################################################


if save_full_output
    println("Storing full simulation output in $(output_format)")
    if output_format == "netcdf"
        filename = joinpath(output_path, "compartments_full.nc")
        println("\t- Output filename: $(filename)")
        save_simulation_netCDF(epi_params, population, filename;G_coords=G_coords, M_coords=M_coords, T_coords=T_coords)
    elseif output_format == "hdf5"
        filename = joinpath(output_path, "compartments_full.h5")
        println("\t- Output filename: $(filename)")
        save_simulation_hdf5(epi_params, population, filename)
    end
end

if save_time_step !== nothing
    export_compartments_date = first_day + Day(export_compartments_time_t - 1)
    filename = joinpath(output_path, "compartments_t_$(export_compartments_date).h5")
    println("Storing compartments at single date $(export_compartments_date):")
    println("\t- Simulation step: $(export_compartments_time_t)")
    println("\t- filename: $(filename)")
    save_simulation_hdf5(epi_params, population, filename; 
                         export_time_t = export_compartments_time_t)
end
