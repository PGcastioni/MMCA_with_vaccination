using DelimitedFiles
using Statistics
using DataFrames
using CSV
using Printf
using DataStructures
using Base.Threads
using ProgressMeter
using NPZ
using ArgParse
using JSON
using Dates
using ArgParse
using HDF5

base_folder = ""

include(joinpath(base_folder, "MMCAcovid19_vac/markov_vac_aux.jl"))
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

function run_simu_params!(epi_params::Epidemic_Params,
    population::Population_Params,
    E₀::Array{Float64, 2},
    A₀::Array{Float64, 2},
    I₀::Array{Float64, 2},
    H₀::Array{Float64, 2},
    R₀::Array{Float64, 2},
    Sᵛ₀::Array{Float64, 2},
    compartments::Union{Nothing, Array{Float64, 5}})


# Reset compartments
reset_params!(epi_params, population)

if initial_compartments != nothing
set_compartments!(epi_params, initial_compartments)
else
set_initial_conditions!(epi_params, population, Sᵛ₀, E₀, A₀, I₀, H₀, R₀)
end


## RUN EPIDEMIC SPREADING
run_epidemic_spreading_mmca!(epi_params, population, tᶜs, tᵛs, κ₀s, ϕs, δs, ϵᵍs; verbose = true )

# Store compartments to later export (can't write to disk here, hdf5 is not thread safe)
if export_compartments
compartments[:, :, :, :, 1]  .= epi_params.ρˢᵍᵥ .* population.nᵢᵍ
compartments[:, :, :, :, 2]  .= epi_params.ρᴱᵍᵥ .* population.nᵢᵍ
compartments[:, :, :, :, 3]  .= epi_params.ρᴬᵍᵥ .* population.nᵢᵍ
compartments[:, :, :, :, 4]  .= epi_params.ρᴵᵍᵥ .* population.nᵢᵍ
compartments[:, :, :, :, 5]  .= epi_params.ρᴾᴴᵍᵥ .* population.nᵢᵍ
compartments[:, :, :, :, 6]  .= epi_params.ρᴾᴰᵍᵥ .* population.nᵢᵍ
compartments[:, :, :, :, 7]  .= epi_params.ρᴴᴿᵍᵥ .* population.nᵢᵍ
compartments[:, :, :, :, 8]  .= epi_params.ρᴴᴰᵍᵥ .* population.nᵢᵍ
compartments[:, :, :, :, 9]  .= epi_params.ρᴿᵍᵥ .* population.nᵢᵍ
compartments[:, :, :, :, 10] .= epi_params.ρᴰᵍᵥ .* population.nᵢᵍ
end
# PIER: if we want to save just the age, location and time infos we have to use
# compartments[:, :, :, :, 1] .= sum( epi_params.ρˢᵍᵥ .* population.nᵢᵍ, dims=(4) )[:,:,:,1]
end

###########################################
############# FILE READING ################
###########################################

data_path = joinpath(base_folder,"data/")
instance_path = joinpath(base_folder,"test")

# Output simulation
output_path = joinpath(instance_path, "output")


config_fname = joinpath(instance_path, "config.json")
config = JSON.parsefile(config_fname);

fitting_params_fname = joinpath(instance_path, "fitting_parameters.json")
fitting_params = JSON.parsefile(fitting_params_fname)

vacparams_dict = config["vaccination"]
npiparams_dict = config["NPI"]
epiparams_dict = config["model"]

first_day = Date(config["simulation"]["first_day_simulation"])
last_day = Date(config["simulation"]["last_day_simulation"])


#T: time steps
T = (last_day - first_day).value + 1

initial_compartments_path = get(config["simulation"], "initial_compartments", nothing)
export_compartments_full = get(config["simulation"], "export_compartments_full", false)
export_compartments_time_t = get(config["simulation"], "export_compartments_time_t", nothing)

if export_compartments_time_t !== nothing
    export_compartments_date = first_day + Day(export_compartments_time_t - 1)
end

if export_compartments_full || export_compartments_time_t != nothing
    export_compartments = true
else
    export_compartments = false
end

println("first_day_simulation = ", first_day)
println("last_day_simulation = ", last_day)
println("export_compartments = ", export_compartments)
println("export_compartments_full = ", export_compartments_full)
println("export_compartments_time_t = ", export_compartments_time_t)
println("initial_compartments = ", initial_compartments_path)

num_compartments = 10

########################################
####### VARIABLES INITIALIZATION #######
########################################

# Patch surface
sᵢ = CSV.read(joinpath(data_path, config["data"]["surface_filename"]), DataFrame)[:,"area"]
# Population info
nᵢ_ages = CSV.read(joinpath(data_path, config["data"]["population_age_filename"], ), DataFrame);
# Patch population by age
nᵢᵍ = copy(transpose(Array{Float64,2}(nᵢ_ages[:,epiparams_dict["age_labels"]])))
# Total patch population
nᵢ = Array{Float64,1}(nᵢ_ages[:,"Total"])
# Total population
total_population = sum(nᵢ)

# Age Contact Matrix
C = readdlm(joinpath(data_path, config["data"]["contact_matrix_filename"]), ',', Float64)
# Num. of patches
M = length(nᵢ)
# Num of stratas
G = size(C)[1]
# Num. of vaccination statuses Vaccinated/Non-vaccinated
V = length(epiparams_dict["kᵥ"])


# Loading mobility network
network = CSV.read(joinpath(data_path, config["data"]["mobility_matrix_filename"]), DataFrame)
edgelist = Array{Int64, 2}(network[:, 1:2])
Rᵢⱼ = copy(network[:, 3])
# Correcting Self Loops
edgelist, Rᵢⱼ = correct_self_loops(edgelist, Rᵢⱼ, M)


## EPIDEMIC PARAMETERS HUMAN BEHAVIOUR
# Average number of contacts per strata
kᵍ = Float64.(epiparams_dict["kᵍ"])
# Average number of contacts at home per strata
kᵍ_h = Float64.(epiparams_dict["kᵍ_h"])
# Average number of contacts at work per strata
kᵍ_w = Float64.(epiparams_dict["kᵍ_w"])
# Degree of mobility per strata
pᵍ = Float64.(epiparams_dict["pᵍ"])
# Density factor
ξ = epiparams_dict["σ"]
# Average household size
σ = epiparams_dict["σ"]
# Check network structure and self-loop correction


## EPIDEMIC PARAMETERS TRANSITION RATES
# Scaling of the asymptomatic infectivity
scale_β = epiparams_dict["scale_β"]
# Infectivity of Symptomatic
βᴵ = epiparams_dict["βᴵ"]
# Infectivity of Asymptomatic
βᴬ = scale_β * βᴵ
# Exposed rate
ηᵍ = Float64.(epiparams_dict["ηᵍ"])
# Asymptomatic rate
αᵍ = Float64.(epiparams_dict["αᵍ"])
# Infectious rate
μᵍ = Float64.(epiparams_dict["μᵍ"])


## EPIDEMIC PARAMETERS TRANSITION RATES VACCINATION
# Direct death probability
θᵍ = Float64.(reduce(hcat, [epiparams_dict["θᵍ"], epiparams_dict["θᵍ"] * epiparams_dict["risk_reduction_dd"]]) )
# Hospitalization probability
γᵍ = Float64.(reduce(hcat, [epiparams_dict["γᵍ"], epiparams_dict["γᵍ"] * epiparams_dict["risk_reduction_h"]]) )
# Fatality probability in ICU
ωᵍ = Float64.(reduce(hcat, [epiparams_dict["ωᵍ"], epiparams_dict["ωᵍ"] * epiparams_dict["risk_reduction_d"]]) )
# Pre-deceased rate
ζᵍ = Float64.(epiparams_dict["ζᵍ"])
# Pre-hospitalized in ICU rate
λᵍ = Float64.(epiparams_dict["λᵍ"])
# Death rate in ICU
ψᵍ = Float64.(epiparams_dict["ψᵍ"])
# ICU discharge rate
χᵍ = Float64.(epiparams_dict["χᵍ"])

# Waning immunity rate 
Λ = epiparams_dict["Λ"] 
# Reinfection rate
Γ = epiparams_dict["Γ"] 
# Relative risk reduction of the probability of infection
rᵥ = Float64.(epiparams_dict["rᵥ"])
# Relative risk reduction of the probability of transmission
kᵥ = Float64.(epiparams_dict["kᵥ"])


println("M = ", M)
println("G = ", G)
println("T = ", T)
println("V = ", V)
println("num_compartments : ", num_compartments)

#########################################################
# Vaccination parameters
#########################################################

# vaccionation dates
start_vacc = vacparams_dict["start_vacc"]
dur_vacc   = vacparams_dict["dur_vacc"]
end_vacc   = start_vacc + dur_vacc

# total vaccinations per age strata
ϵᵍ = vacparams_dict["ϵᵍ"] * round( total_population * vacparams_dict["percentage_of_vacc_per_day"] )

tᵛs = [start_vacc, end_vacc, T]
ϵᵍs = ϵᵍ .* [0  Int(vacparams_dict["are_there_vaccines"])  0] 


#########################################################
# Containement measures
#########################################################

# Mobility reduction
κ₀_df = CSV.read(joinpath(data_path, config["data"]["kappa0_filename"]), DataFrame);

# syncronize containment measures with simulation
κ₀_df.time = map(x -> (x .- first_day).value + 1, κ₀_df.date)

# Timesteps when the containment measures will be applied
tᶜs = Int64.(npiparams_dict["tᶜs"])

# Array of level of confinement
# κ₀s = κ₀_df.reduction[:]
κ₀s = Float64.(npiparams_dict["κ₀s"])
# Array of premeabilities of confined households
ϕs = Float64.(npiparams_dict["ϕs"])
# Array of social distancing measures
δs = Float64.(npiparams_dict["δs"])

##################################################
####### INITIALIZATION OF THE EPIDEMICS ##########
##################################################

# Load initial conditions
if initial_compartments_path != nothing
    # use initial compartments matrix to initialize simulations
    initial_compartments = h5open(initial_compartments_path, "r") do file
        read(file, "compartments")
    end
else
    # use seeds to initialize simulations
    initial_compartments = nothing
end

# Initial seeds (intial condition at the begining of the pandemic)
A0_instance_filename = joinpath(instance_path, config["simulation"]["A0_filename"])
conditions₀ = CSV.read(A0_instance_filename, DataFrame)


## INITIALIZATION OF THE EPIDEMICS

# Initial number of exposed individuals
E₀ = zeros(G, M)
# Initial number of infectious asymptomatic individuals
A₀ = zeros(Float64, G, M)

println("Total population = ", total_population)
# Distribution of the intial infected individuals per strata
# WARN: ni idea de por que es necesario el .+ 1 aqui (las semillas ya vienen indexadas partiendo de 1) 
# pero si no lo pongo la simulacion simplemente no funciona

#descomentado
A₀[1, Int.(conditions₀[:,"idx"])] .= 0.12 .* conditions₀[:,"seed"]
A₀[2, Int.(conditions₀[:,"idx"])] .= 0.16 .* conditions₀[:,"seed"]
A₀[3, Int.(conditions₀[:,"idx"])] .= 0.72 .* conditions₀[:,"seed"]

#A₀ = A₀ / total_population

# Initial number of infectious symptomatic individuals
I₀ = zeros(Float64, G, M)

E₀ = nᵢᵍ / total_population * 1000
A₀ = nᵢᵍ / total_population * 1000
I₀ = nᵢᵍ / total_population * 1000

H₀ = nᵢᵍ * 0
# R₀ = population.nᵢᵍ / total_population * 23e5
R₀ = nᵢᵍ * 0
#Sᵛ₀ = (population.nᵢᵍ .- E₀ .- A₀ .- I₀ .- H₀ .- R₀) .* 0.5
Sᵛ₀ = nᵢᵍ * 0

#########################################################
######### SETTING UP SIMULATION VARIABLES ###############
#########################################################

# Parameters to simulate


compartments = zeros(Float64, G, M, T, V, num_compartments);
size(compartments)


population = Population_Params(G, M, nᵢᵍ, kᵍ, kᵍ_h, kᵍ_w, C, pᵍ, edgelist, Rᵢⱼ, sᵢ, ξ, σ)
epi_param = Epidemic_Params(βᴵ,  βᴬ, ηᵍ, αᵍ, μᵍ, θᵍ, γᵍ, ζᵍ, λᵍ, ωᵍ, ψᵍ, χᵍ,  Λ, Γ, rᵥ, kᵥ, G, M, T, V)



# Patch to CCAA mapping matrix
PatchToCCAA = npzread(joinpath(data_path, fitting_params["patch_to_ccaa_filename"]))
n_ccaa = size(PatchToCCAA)[1]
n_patches = size(PatchToCCAA)[2]

########################################################
################ RUN THE SIMULATION ####################
########################################################

# Run the simulation for all the parameters
    
# Run simu
run_simu_params!(epi_param,
                 population,
                 E₀,
                 A₀,
                 I₀,
                 H₀,
                 R₀,
                 Sᵛ₀,
                compartments)   


##############################################################
################## STORING THE RESULTS #######################
##############################################################

println("Storing results to: $output_path")

if !isdir(output_path)
  mkpath(output_path)
end



# store compartments
if export_compartments
    println("Storing compartments")

    if export_compartments_time_t != nothing
        filename = joinpath(output_path, "compartments_$(export_compartments_date).h5")
        # filename = joinpath(output_path, "compartments_$(export_compartments_date):$(export_compartments_time_t)_$id.h5")
        h5open(filename, "w") do file
            write(file, "compartments", compartments[:,:,export_compartments_time_t,:,:])
        end
    end

    if export_compartments_full
        filename = joinpath(output_path, "compartments_full.h5")
        h5open(filename, "w") do file
            write(file, "compartments", compartments[:,:,:,:,:])
        end
    end
end