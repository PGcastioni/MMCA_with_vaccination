using DelimitedFiles
# using Statistics
using DataFrames
using CSV
using Printf
using DataStructures
using Base.Threads
using ProgressMeter
using NPZ
# using MMCAcovid19
using ArgParse
using JSON
using Dates
using ArgParse
using HDF5


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
        "--params", "-p"
            help = "parameters file (params.csv)"
            default = nothing
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


include("MMCAcovid19_vac/markov_vac_aux.jl")
include("MMCAcovid19_vac/markov_vac.jl")

args = parse_commandline()
data_path = args["data-folder"]
instance_path = args["instance-folder"]

# Config file
config = JSON.parsefile(args["config"])
if !haskey(config, "simulation")
    config["simulation"] = Dict()
end

# overwrite config with command line
if args["start-date"] != nothing
    config["simulation"]["first_day_simulation"] = args["start-date"]
end
if args["end-date"] != nothing
    config["simulation"]["last_day_simulation"] = args["end-date"]
end
if args["initial-compartments"] != nothing
    config["simulation"]["initial_compartments"] = args["initial-compartments"]
end
if args["export-compartments-time-t"] != nothing
    config["simulation"]["export_compartments_time_t"] = args["export-compartments-time-t"]
end
if args["export-compartments-full"] == true
    config["simulation"]["export_compartments_full"] = true
end


first_day = Date(get(config["simulation"], "first_day_simulation", "2020-02-09"))
last_day = Date(get(config["simulation"], "last_day_simulation", "2020-03-09"))
initial_compartments_path = get(config["simulation"], "initial_compartments", nothing)
export_compartments_full = get(config["simulation"], "export_compartments_full", false)
export_compartments_time_t = get(config["simulation"], "export_compartments_time_t", nothing)

if export_compartments_time_t != nothing
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

## -----------------------------------------------------------------------------
## LOADING DATA
## -----------------------------------------------------------------------------


# Parameters to simulate
if args["params"] != nothing
    paramsDF = CSV.read(args["params"], DataFrame)
else
    paramsDF = CSV.read(joinpath(instance_path, "params.csv"), DataFrame)
end

# Output simulation
output_path = joinpath(instance_path, "output")

# Loading mobility network
network = CSV.read(joinpath(data_path, "R_mobility_matrix.csv"), DataFrame)
edgelist = Array{Int64, 2}(network[:, 1:2])
Rᵢⱼ = copy(network[:, 3])

# Patch surface
sᵢ = CSV.read(joinpath(data_path, "S_patch_surface_area.csv"), DataFrame)[:,"area"]

# Population info
nᵢ_ages = CSV.read(@sprintf("%s/n_population_patch_age.csv", data_path), DataFrame);

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
A0_instance_filename = @sprintf("%s/%s", instance_path, config["data"]["A0_filename"])
A0_data_filename = @sprintf("%s/%s", data_path, config["data"]["A0_filename"])
if isfile(A0_instance_filename)
  println("Loading A0 from instance/ folder")
  conditions₀ = CSV.read(A0_instance_filename, DataFrame)
else
  println("Loading A0 from data/ folder")
  conditions₀ = CSV.read(A0_data_filename, DataFrame)
end

C = readdlm(@sprintf("%s/C_age_contact_matrix.csv", data_path),
            ',', Float64)

# Containement measures
κ₀_df = CSV.read(@sprintf("%s/%s", data_path, config["data"]["kappa0_filename"]), DataFrame)

# Patch to CCAA mapping matrix
PatchToCCAA = npzread(@sprintf("%s/patch_to_ccaa.npy", data_path))
n_ccaa = size(PatchToCCAA)[1]
n_patches = size(PatchToCCAA)[2]


## -----------------------------------------------------------------------------
## SETTING PARAMETERS (some of them are just placeholders for the values on DF)
## -----------------------------------------------------------------------------


## POPULATION

# Patch population for each strata
nᵢᵍ = copy(transpose(Array{Float64,2}(nᵢ_ages[:,["Y", "M", "O"]])))

# Total patch population
nᵢ = Array{Float64,1}(nᵢ_ages[:,"Total"])

# Num. of patches
M = length(nᵢ)

# Num of stratas
G = size(C)[1]

# Num. of vaccination statuses Vaccinated/Non-vaccinated
V = 2

# Average number of contacts per strata
kᵍ = [11.8, 13.3, 6.76]

# Average number of contacts at home per strata
kᵍ_h = [3.15, 3.17, 3.28]

# Average number of contacts at work per strata
kᵍ_w = [1.72, 5.18, 0.0]

# Degree of mobility per strata
pᵍ = [0.0, 1.0, 0.00]

# Density factor
ξ = 0.01

# Average household size
σ = 2.5

# Check network structure and self-loop correction
edgelist, Rᵢⱼ = correct_self_loops(edgelist, Rᵢⱼ, M)


## EPIDEMIC PARAMETERS

# Scaling of the asymptomatic infectivity
scale_β = 0.51
# Infectivity of Symptomatic
βᴵ = 0.0903
# Infectivity of Asymptomatic
βᴬ = scale_β * βᴵ
# Exposed rate
ηᵍ = [1/3.64, 1/3.64, 1/3.64]
# Asymptomatic rate
αᵍ = [1/3.76, 1/1.56, 1/1.56]
# Infectious rate
μᵍ = [1/1.0, 1/3.2, 1/3.2]


############################################
# CHANGE CODE TO RUN WITH VACCINATION MODEL
############################################

# Direct death probability
# θᵍ = [0.00, 0.008, 0.047]
θᵍ = [0.     0.0*(1-0.95);
      0.0  0.0*(1-0.95); 
      0.0  0.0*(1-0.95)]

# Direct death probability
# γᵍ = [0.0003, 0.003, 0.026]
γᵍ = [0.003  0.003*(1-0.9);
      0.01   0.01*(1-0.9); 
      0.08   0.08*(1-0.9)] 

# Fatality probability in ICU
# ωᵍ = [0.30, 0.30, 0.30]
ωᵍ = [0.  0.0*(1-0.95); 
      0.04  0.04*(1-0.95); 
      0.3  0.3*(1-0.95)]

############################################

# Pre-deceased rate
ζᵍ = [1/7.8, 1/7.8, 1/7.8]
# Pre-hospitalized in ICU rate
λᵍ = ones(Float64, 3)
# Death rate in ICU
ψᵍ = [1/7.0, 1/7.0, 1/7.0]
# ICU discharge rate
χᵍ = [1/21.0, 1/21.0, 1/21.0]


# Number of timesteps
# dia inicial: 9 Feb
# dia final: 14 Abril (66 dias)
# T = 123
# T = 66
T = (last_day - first_day).value + 1


############################################
# CHANGE CODE TO RUN WITH VACCINATION MODEL
############################################

# Epidemic parameters
# Epidemic_Params(βᴵ,  βᴬ, ηᵍ, αᵍ, μᵍ, θᵍ, γᵍ, ζᵍ, λᵍ, ωᵍ, ψᵍ, χᵍ,  Λ, Γ, rᵥ, kᵥ, G, M, T, V)


Λ = 1 / 50 # Waning immunity rate 
Γ = 1 / 100 # Reinfection rate
rᵥ = [0., 0.6 ] # Relative risk reduction of the probability of infection
kᵥ = [0., 0.4 ] # Relative risk reduction of the probability of transmission

# We will need to define tᵛs for the vaccination times
# in order to have times for confinementmeasures (tᶜs) and 
# vaccinations (tᵛs) --> also update run_epidemic_spreading_mmca!
# so it accept both parameters (PG now howto)

total_population = sum(nᵢ)
percentage_of_vacc_per_day = 0.05
# total vaccinations per age strata
ϵᵍ = [0.1 , 0.4 , 0.5] * round( total_population * percentage_of_vacc_per_day )

start_vacc = 10
dur_vacc   = 20
end_vacc   = start_vacc + dur_vacc
tᶜs = [start_vacc, end_vacc, T] # rename to tᵛs
ϵᵍs = ϵᵍ .* [0 1 0] # reshape(ϵᵍ, (3, 1))

# syncronize containment measures with simulation
κ₀_df.time = map(x -> (x .- first_day).value + 1, κ₀_df.date)
# Timesteps when the containment measures will be applied
# tᶜs = κ₀_df.time[:]

# Array of level of confinement
# κ₀s = κ₀_df.reduction[:]
κ₀s = zeros(length(tᶜs))
# Array of premeabilities of confined households
ϕs = ones(Float64, length(tᶜs))
# Array of social distancing measures
#δs = ones(Float64, length(tᶜs))
δs = zeros(Float64, length(tᶜs))
############################################

## INITIALIZATION OF THE EPIDEMICS

# Initial number of exposed individuals
E₀ = zeros(G, M)
# Initial number of infectious asymptomatic individuals
A₀ = zeros(Float64, G, M)

println("Total population = ", total_population)
# Distribution of the intial infected individuals per strata
# WARN: ni idea de por que es necesario el .+ 1 aqui (las semillas ya vienen indexadas partiendo de 1) 
# pero si no lo pongo la simulacion simplemente no funciona
"""
A₀[1, Int.(conditions₀[:,"idx"])] .= 0.12 .* conditions₀[:,"seed"]
A₀[2, Int.(conditions₀[:,"idx"])] .= 0.16 .* conditions₀[:,"seed"]
A₀[3, Int.(conditions₀[:,"idx"])] .= 0.72 .* conditions₀[:,"seed"]
"""
#A₀ = A₀ / total_population

# Initial number of infectious symptomatic individuals
I₀ = zeros(Float64, G, M)


E₀ = nᵢᵍ / total_population * 1000
A₀ = nᵢᵍ / total_population * 1000
I₀ = nᵢᵍ / total_population * 1000
H₀ = nᵢᵍ * 0
# R₀ = population.nᵢᵍ / total_population * 23e5 ### NUMERO REALISTICO?
R₀ = nᵢᵍ * 0
#S₁ = (population.nᵢᵍ .- E₀ .- A₀ .- I₀ .- H₀ .- R₀) .* 0.5
S₁ = nᵢᵍ * 0

## -----------------------------------------------------------------------------
## SETTING UP SIMULATION VARIABLES
## -----------------------------------------------------------------------------

# Num. of set of epidemic parameters
total_simulations = length(paramsDF.id)
# total_simulations = 3

println("Allocating output matrices (size = n_ccaa x T x total_simulations)")
println("n_patches: ", n_patches)
println("n_ccaa: ", n_ccaa)
println("T: ", T)
println("total_simulations: ", total_simulations)

# Arrays where the output variables will be stored
incidence  = zeros(Float64, n_ccaa, T - 1, total_simulations)
prevalence = zeros(Float64, n_ccaa, T, total_simulations)
deaths = zeros(Float64, n_ccaa, T, total_simulations)
deaths_new = zeros(Float64, n_ccaa, T - 1, total_simulations)

if export_compartments
  num_compartments = 10
  compartments = zeros(Float64, G, M, T, num_compartments, total_simulations)
else
  compartments = nothing
end


## SETTING UP THE THREADING VARIABLES

# Number of threads used to parallelize the execution
# nThreads = 48;
nThreads = nthreads();
# nThreads = 1;
println("nThreads: ", nThreads)
flush(stdout)

# Circular queues to optimize and reuse the population and epidemic structures
populations = CircularDeque{Population_Params}(nThreads)
epi_params = CircularDeque{Epidemic_Params}(nThreads)

# Populate the circular deque
for t in 1:nThreads
    push!(populations, Population_Params(G, M, nᵢᵍ, kᵍ, kᵍ_h, kᵍ_w, C, pᵍ, edgelist, Rᵢⱼ, sᵢ, ξ, σ))
    push!(epi_params, Epidemic_Params(βᴵ,  βᴬ, ηᵍ, αᵍ, μᵍ, θᵍ, γᵍ, ζᵍ, λᵍ, ωᵍ, ψᵍ, χᵍ,  Λ, Γ, rᵥ, kᵥ, G, M, T, V))
end


## -----------------------------------------------------------------------------
## FUNCTIONS
## -----------------------------------------------------------------------------


function set_compartments!(epi_params, compartments)
    @assert size(compartments) == (size(epi_params.ρˢᵍ)[1], size(epi_params.ρˢᵍ)[2], 10)
    total_population = sum(compartments, dims=(3))[:,:,1]
    epi_params.ρˢᵍ[:,:,1] .= compartments[:, :, 1] ./ total_population
    epi_params.ρᴱᵍ[:,:,1] .= compartments[:, :, 2] ./ total_population
    epi_params.ρᴬᵍ[:,:,1] .= compartments[:, :, 3] ./ total_population
    epi_params.ρᴵᵍ[:,:,1] .= compartments[:, :, 4] ./ total_population
    epi_params.ρᴾᴴᵍ[:,:,1] .= compartments[:, :, 5] ./ total_population
    epi_params.ρᴾᴰᵍ[:,:,1] .= compartments[:, :, 6] ./ total_population
    epi_params.ρᴴᴿᵍ[:,:,1] .= compartments[:, :, 7] ./ total_population
    epi_params.ρᴴᴰᵍ[:,:,1] .= compartments[:, :, 8] ./ total_population
    epi_params.ρᴰᵍ[:,:,1] .= compartments[:, :, 9] ./ total_population
    epi_params.ρᴿᵍ[:,:,1] .= compartments[:, :, 10] ./ total_population

    epi_params.ρˢᵍ[isnan.(epi_params.ρˢᵍ)] .= 0
    epi_params.ρᴱᵍ[isnan.(epi_params.ρᴱᵍ)] .= 0
    epi_params.ρᴬᵍ[isnan.(epi_params.ρᴬᵍ)] .= 0
    epi_params.ρᴵᵍ[isnan.(epi_params.ρᴵᵍ)] .= 0
    epi_params.ρᴾᴴᵍ[isnan.(epi_params.ρᴾᴴᵍ)] .= 0
    epi_params.ρᴾᴰᵍ[isnan.(epi_params.ρᴾᴰᵍ)] .= 0
    epi_params.ρᴴᴿᵍ[isnan.(epi_params.ρᴴᴿᵍ)] .= 0
    epi_params.ρᴴᴰᵍ[isnan.(epi_params.ρᴴᴰᵍ)] .= 0
    epi_params.ρᴰᵍ[isnan.(epi_params.ρᴰᵍ)] .= 0
    epi_params.ρᴿᵍ[isnan.(epi_params.ρᴿᵍ)] .= 0
end


"""
    run_simu_params!(epi_params::Epidemic_Params,
                     population::Population_Params,
                     paramsdf::dataframe,
                     indx_id::int64,
                     a₀::array{float64, 2},
                     i₀::array{float64, 2},
                     incidence::array{float64, 2},
                     prevalence::array{float64, 2},
                     deaths_new::array{float64, 2},
                     deaths::array{float64, 2})

Runs a SEAIRHD simulation on a specific set of parameters stored in the DF
updateing the tables of prevalence, incidence, num. of deaths and daily deaths.
"""
function run_simu_params!(epi_params::Epidemic_Params,
                          population::Population_Params,
                          paramsDF::DataFrame,
                          indx_id::Int64,
                          A₀::Array{Float64, 2},
                          I₀::Array{Float64, 2},
                          incidence::Array{Float64, 3},
                          prevalence::Array{Float64, 3},
                          deaths_new::Array{Float64, 3},
                          deaths::Array{Float64, 3},
                          compartments::Union{Nothing, Array{Float64, 5}})

    # Loading parameters from DF
    id, P, β, scale_β, τ_inc, scale_ea, τᵢ, δ, ϕ, scale₀ = paramsDF[indx_id, :]

    # Set epidemic params to the ones speficied on the DF
    epi_params.βᴬ .= scale_β * β
    epi_params.βᴵ .= β
    epi_params.ηᵍ .= 1.0/(τ_inc * (1.0 - scale_ea))
    epi_params.αᵍ .= [1.0/(τᵢ - 1 + τ_inc * scale_ea),
                      1.0/(τ_inc * scale_ea),
                      1.0/(τ_inc * scale_ea)]
    epi_params.μᵍ .= [1.0, 1.0/τᵢ, 1.0/τᵢ]

    # Reset compartments
    reset_params!(epi_params, population)
    """
    if initial_compartments != nothing
        set_compartments!(epi_params, initial_compartments)
    else
        set_initial_infected!(epi_params, population, E₀, scale₀ .* A₀, I₀)
    end
    """
    set_initial_conditions!(epi_params, population, S₁, E₀, A₀, I₀, H₀, R₀)


    # Set containment parameters
    ϕs .= ϕ
    δs .= δ

    ## RUN EPIDEMIC SPREADING
    run_epidemic_spreading_mmca!(epi_params, population, tᶜs, κ₀s, ϕs, δs, ϵᵍs; verbose=false)

    # Compute the prevalence
    prevalence[:, :, indx_id] = PatchToCCAA * sum((epi_params.ρᴵᵍ[:, :, 1:epi_params.T] .+
                                   epi_params.ρᴴᴰᵍ[:, :, 1:epi_params.T] .+
                                   epi_params.ρᴴᴿᵍ[:, :, 1:epi_params.T] .+
                                   epi_params.ρᴾᴴᵍ[:, :, 1:epi_params.T] .+
                                   epi_params.ρᴾᴰᵍ[:, :, 1:epi_params.T] .+
                                   epi_params.ρᴰᵍ[:, :, 1:epi_params.T] .+
                                   epi_params.ρᴿᵍ[:, :, 1:epi_params.T] .+
                                   epi_params.ρᴬᵍ[:, :, 1:epi_params.T] .+
                                   epi_params.ρᴱᵍ[:, :, 1:epi_params.T]) .*
                                  population.nᵢᵍ, dims = (1))[1,:,:]  # sum by CCAA

    # Compute the incidence
    incidence[:, :, indx_id]  = diff(prevalence[:, :, indx_id], dims=2)

    # Compute total number of deaths
    deaths[:, :, indx_id] = PatchToCCAA * sum(epi_params.ρᴰᵍ[:, :, 1:epi_params.T] .*
                                          population.nᵢᵍ, dims = (1))[1,:,:]  # sum by CCAA

    # Compute daily new deaths
    deaths_new[:, :, indx_id] = diff(deaths[:, :, indx_id], dims=2)

    # Store compartments to later export (can't write to disk here, hdf5 is not thread safe)
    if export_compartments
      compartments[:, :, :, 1, indx_id] = epi_params.ρˢᵍ .* population.nᵢᵍ
      compartments[:, :, :, 2, indx_id] = epi_params.ρᴱᵍ .* population.nᵢᵍ
      compartments[:, :, :, 3, indx_id] = epi_params.ρᴬᵍ .* population.nᵢᵍ
      compartments[:, :, :, 4, indx_id] = epi_params.ρᴵᵍ .* population.nᵢᵍ
      compartments[:, :, :, 5, indx_id] = epi_params.ρᴾᴴᵍ .* population.nᵢᵍ
      compartments[:, :, :, 6, indx_id] = epi_params.ρᴾᴰᵍ .* population.nᵢᵍ
      compartments[:, :, :, 7, indx_id] = epi_params.ρᴴᴿᵍ .* population.nᵢᵍ
      compartments[:, :, :, 8, indx_id] = epi_params.ρᴴᴰᵍ .* population.nᵢᵍ
      compartments[:, :, :, 9, indx_id] = epi_params.ρᴰᵍ .* population.nᵢᵍ
      compartments[:, :, :, 10, indx_id] = epi_params.ρᴿᵍ .* population.nᵢᵍ
    end
end


## -----------------------------------------------------------------------------
## RUN THE SIMULATION
## -----------------------------------------------------------------------------

# Setup progress bar
p = Progress(total_simulations)

# Circular Deque Lock
lockS = SpinLock()

# Run the simulation for all the parameters
@threads for indx_id in 1:total_simulations
    
    # println("Thread: ", threadid(), " row: ", indx_id)

    # Lock the data structure
    lock(lockS)

    # Recover population and epidemic parameters
    population = pop!(populations)
    epi_param = pop!(epi_params)

    # Unlock data structure
    unlock(lockS)

    # Run simu
    run_simu_params!(epi_param,
                     population,
                     paramsDF,
                     indx_id,
                     A₀,
                     I₀,
                     incidence,
                     prevalence,
                     deaths_new,
                     deaths,
                     compartments)

    # Lock the data structure
    lock(lockS)

    # Free population and epidemic parameters
    push!(populations, population)
    push!(epi_params, epi_param)

    # Unlock data structure
    unlock(lockS)

    # Update progress bar
    next!(p)
    
end


## -----------------------------------------------------------------------------
## STORING THE RESULTS
## -----------------------------------------------------------------------------

println("Storing results to: $output_path")

if !isdir(output_path)
  mkpath(output_path)
end

# npzwrite(joinpath(output_path, "incidence.npy"), incidence)
# npzwrite(joinpath(output_path, "prevalence.npy"), prevalence)
# npzwrite(joinpath(output_path, "deaths.npy"), deaths)
# npzwrite(joinpath(output_path, "deaths_new.npy"), deaths_new)

# store incidence and deaths
for i in 1:total_simulations
    id = paramsDF[i, "id"]
    npzwrite(joinpath(output_path, "incidence_$id.npy"), incidence[:,:,i])
    npzwrite(joinpath(output_path, "deaths_new_$id.npy"), deaths_new[:,:,i])
    # npzwrite(joinpath(output_path, "deaths_$id.npy"), deaths[:,:,i])
    # npzwrite(joinpath(output_path, "prevalence_$id.npy"), prevalence[:,:,i])
end

# store compartments
if export_compartments
    println("Storing compartments")
    for i in 1:total_simulations

        id = paramsDF[i, "id"]

        if export_compartments_time_t != nothing
            filename = joinpath(output_path, "compartments_$(export_compartments_date)_$id.h5")
            # filename = joinpath(output_path, "compartments_$(export_compartments_date):$(export_compartments_time_t)_$id.h5")
            h5open(filename, "w") do file
                write(file, "compartments", compartments[:,:,export_compartments_time_t,:,i])
            end
        end

        if export_compartments_full
            filename = joinpath(output_path, "compartments_full_$id.h5")
            h5open(filename, "w") do file
                write(file, "compartments", compartments[:,:,:,:,i])
            end
        end
    end
end

# store configuration
open(joinpath(output_path, "config.json"), "w") do f
    JSON.print(f, config, 4)
end
