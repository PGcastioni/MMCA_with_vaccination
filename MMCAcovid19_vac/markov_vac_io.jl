using Logging


function create_default_epiparameters()
    epiparams_dict = Dict()
    epiparams_dict["kᵍ"] = [11.8, 13.3, 6.76]
    epiparams_dict["kᵍ_h"] = [3.15, 3.17, 3.28]
    epiparams_dict["kᵍ_w"] = [1.72, 5.18, 0.0]
    epiparams_dict["pᵍ"] = [0.0, 1.0, 0.00]
    epiparams_dict["ξ"] = 0.01
    epiparams_dict["σ"] = 2.5
    epiparams_dict["scale_β"] = 0.51
    epiparams_dict["βᴵ"] = 0.0903
    epiparams_dict["ηᵍ"] = [0.2747252747252747, 0.2747252747252747, 0.2747252747252747]
    epiparams_dict["αᵍ"] = [0.26595744680851063, 0.641025641025641, 0.641025641025641]
    epiparams_dict["μᵍ"] = [1.0, 0.3125, 0.3125]
    epiparams_dict["θᵍ"] = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
    epiparams_dict["γᵍ"] = [0.003, 0.01, 0.08]
    epiparams_dict["risk_reduction_h"] = 0.1
    epiparams_dict["ωᵍ"] = [0.0, 0.04, 0.3]
    epiparams_dict["risk_reduction_d"] = 0.05
    epiparams_dict["ζᵍ"] = [0.12820512820512822, 0.12820512820512822, 0.12820512820512822]
    epiparams_dict["λᵍ"] = [1.0, 1.0, 1.0]
    epiparams_dict["ψᵍ"] = [0.14285714285714285, 0.14285714285714285, 0.14285714285714285]
    epiparams_dict["χᵍ"] = [0.047619047619047616, 0.047619047619047616, 0.047619047619047616]
    epiparams_dict["Λ"] = 0.02
    epiparams_dict["Γ"] = 0.01
    epiparams_dict["rᵥ"] = [0.0, 0.6]
    epiparams_dict["kᵥ"] = [0.0, 0.4]
    epiparams_dict["age_labels"] = ["Y", "M", "O"]

    return epiparams_dict
end

function create_default_vacparameters()
    vacparams_dict = Dict()
    vacparams_dict["ϵᵍ"] = [0.1 , 0.4 , 0.5]
    vacparams_dict["percentage_of_vacc_per_day"] = 0.005
    vacparams_dict["start_vacc"] = 2
    vacparams_dict["dur_vacc"] = 8
    vacparams_dict["are_there_vaccines"] = false

    return vacparams_dict
end

function create_default_npiparameters()
    npiparams_dict = Dict()
    ## It's important that the default parameters are those of the absence of 
    ## lockdowns, because they are the one the code refers to if the key "are_there_npi" = false
    npiparams_dict["κ₀s"] = [0.0]
    npiparams_dict["ϕs"] = [1.0]
    npiparams_dict["δs"] = [0.0]
    npiparams_dict["tᶜs"] =  [1]

    return npiparams_dict
end

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


function update_config!(config, cmd_line_args)
    # Define dictionary containing epidemic parameters
    if !haskey(config, "model")
        config["model"] = create_default_epiparameters()
    end

    # Define dictionary containing vaccination parameters
    if !haskey(config, "vaccination")
        config["vaccination"] = create_default_vacparameters()
    end

    # Define dictionary containing npi parameters
    if (!haskey(config, "NPI") | !config["NPI"]["are_there_npi"] )
        config["NPI"] = create_default_npiparameters()
    end

    # overwrite config with command line
    if cmd_line_args["start-date"] != nothing
        config["simulation"]["first_day_simulation"] = cmd_line_args["start-date"]
    end
    if cmd_line_args["end-date"] != nothing
        config["simulation"]["last_day_simulation"] = cmd_line_args["end-date"]
    end
    if cmd_line_args["initial-compartments"] != nothing
        config["simulation"]["initial_compartments"] = cmd_line_args["initial-compartments"]
    end
    if cmd_line_args["export-compartments-time-t"] != nothing
        config["simulation"]["export_compartments_time_t"] = cmd_line_args["export-compartments-time-t"]
    end
    if cmd_line_args["export-compartments-full"] == true
        config["simulation"]["export_compartments_full"] = true
    end

    nothing
end