using Printf
using ArgParse
using Logging

using CSV
using NPZ
using JSON
using HDF5
using DataStructures
using DelimitedFiles
using DataFrames
using NetCDF

base_folder = ".."

include(joinpath(base_folder, "MMCAcovid19_vac/markov_vac_aux.jl"))
include(joinpath(base_folder, "MMCAcovid19_vac/markov_vac_io.jl"))
include(joinpath(base_folder, "MMCAcovid19_vac/markov_vac.jl"))

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--config", "-c"
            help = "config file (json file)"
            required = true
        "--data-folder", "-d"
            help = "data folder"
            required = true
        "--output-folder", "-o"
            help = "output folder (to store intial condition)"
            default = ""
        "--seed-fname", "-s"
            help = "file name with initiall seeds)"
            default=""
        "--init-type", "-t"
            help = "Initialization option (A0 seeds / uniform)"
            default = "PG"
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()

    config_fname  = args["config"]
    data_folder   = args["data-folder"]
    output_folder = args["output-folder"]
    seeds_fname   = args["seed-fname"]

    output_format = "netcdf"

    config          = JSON.parsefile(config_fname);
    data_dict       = config["data"]
    pop_params_dict = config["population_params"]
    epi_params_dict = config["epidemic_params"]

    # Reading metapopulation Dataframe
    metapop_data_filename = joinpath(data_folder, data_dict["metapopulation_data_filename"])
    metapop_df            = CSV.read(metapop_data_filename, DataFrame)

    # Loading metapopulation patches info (surface, label, population by age)
    metapop_data_filename = joinpath(data_folder, data_dict["metapopulation_data_filename"])
    metapop_df = CSV.read(metapop_data_filename, DataFrame)

    # Loading mobility network
    mobility_matrix_filename = joinpath(data_folder, data_dict["mobility_matrix_filename"])
    network_df  = CSV.read(mobility_matrix_filename, DataFrame)

    # Metapopulations patches coordinates (labels)
    M_coords = map(String,metapop_df[:, "id"])
    M = length(M_coords)

    # Coordinates for each age strata (labels)
    G_coords = map(String, pop_params_dict["age_labels"])
    G = length(G_coords)

    T = 1

    population = init_pop_param_struct(G, M, G_coords, pop_params_dict, metapop_df, network_df)
    epi_params = init_epi_parameters_struct(G, M, T, G_coords, epi_params_dict)


    V = epi_params.V
    N = epi_params.NumComps
    comp_coords = epi_params.CompLabels
    V_coords = ["NV", "V"]

    nᵢ = Array{Float64,1}(metapop_df[:,"Total"])
    total_population = sum(nᵢ)

    E₀ = population.nᵢᵍ / total_population * 10
    A₀ = population.nᵢᵍ / total_population * 10
    I₀ = population.nᵢᵍ / total_population * 10


    S₀_idx = 1
    E₀_idx = 2
    A₀_idx = 3
    I₀_idx = 4
    V_idx  = 1
    
    S₀ = population.nᵢᵍ .- (E₀ + A₀ + I₀)

    compartments = zeros(Float64, G, M, V, N);
    
    compartments[:, :, V_idx, S₀_idx] .= S₀
    compartments[:, :, V_idx, E₀_idx] .= E₀
    compartments[:, :, V_idx, A₀_idx] .= A₀
    compartments[:, :, V_idx, I₀_idx] .= I₀

    if output_format == "netcdf"
        output_fname = joinpath(output_folder, "generic_initial_condition.nc")
        isfile(output_fname) && rm(output_fname)
        nccreate(output_fname, "data", "G", G_coords, "M", M_coords, "V", V_coords, "epi_states", comp_coords)
        ncwrite(compartments, output_fname, "data")
    elseif output_format == "hdf5"
        output_fname = joinpath(output_folder, "generic_initial_condition.h5")
        isfile(output_fname) && rm(output_fname)
        h5open(output_fname, "w") do file
            write(file, "data", compartments)
        end
    end
end

main()