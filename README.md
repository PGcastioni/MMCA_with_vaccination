
# MMCA_with_vaccination
Further development of the Microscopic Markov Chain Approach to epidemics. More specifically this code is the extension of the MMCAcovid19.jl code, created exactly to deal with COVID19.
Each of the previous compartment has being doubled to account for the fact that there are vaccinated and unvaccinated individuals. Having received the vaccine changes some of the transition probabilities related with the most negative aspects of the disease, such as transmission, hospitalization and death.

The file markov_vac.jl contains all the function needed to numerically solve the difference equation of the model, while the markov_vac_aux.jl contains additional functions related to the definition of structures of data, initial condition, print functions etc.

Also a list Jupyter notebooks is present. The file Working_model.ipynb shows how to initialize the parameters and how to run simulations, together with a example of a plot with multiple scenarios with different vaccination strategies. In Scenarios.ipynb we run different scenarios of the vaccination processes, compare them against the baseline simulation with no vaccination and plot the results. In Pareto.ipynb we run a big number of simulations that will then just be plotted in a scatter plot to identify the Pareto front. In Calibration_model.ipynb we show how the model is calibrated using the Turing package of Julia.


## run_simulations.jl

Main script to run model simulations, with multiple configuration options.

Examples:

	julia --project=scripts/ scripts/run_simulations.jl

	# start and end date can be provided in the config.json
	julia --project=scripts/ scripts/run_simulations.jl -i experiments/test20 -d data -c data/config.json

	# or can be provided as command line arguments
	julia --project=scripts/ scripts/run_simulations.jl -d data -c data/config.json -i experiments/test20 --start-date 2020-02-09 --end-date 2020-05-01

	# to export the compartments at a given time t use the option --export-compartments-time-t
	julia --project=scripts/ scripts/run_simulations.jl -d data -c data/config.json -i experiments/test20 --start-date 2020-02-09 --end-date 2020-05-01 --export-compartments-time-t 5

	# to export the full compartments time series (for all dates), use the option --export-compartments
	julia --project=scripts/ scripts/run_simulations.jl -d data -c data/config.json -i experiments/test20 --start-date 2020-02-09 --end-date 2020-05-01 --export-compartments

	# to use the compartments from other simulation as initial conditions, use the option --initial-compartments
	julia --project=scripts/ scripts/run_simulations.jl -d data -c data/config.json -i experiments/test20 --start-date 2020-02-13 --end-date 2020-05-01 --initial-compartments experiments/test20/output/compartments_2020-02-13_10.h5


Usage:
	
	usage: run_simulations.jl -c CONFIG -d DATA-FOLDER -i INSTANCE-FOLDER
	                        [--export-compartments-full]
	                        [--export-compartments-time-t EXPORT-COMPARTMENTS-TIME-T]
	                        [--initial-compartments INITIAL-COMPARTMENTS]
	                        [--start-date START-DATE]
	                        [--end-date END-DATE]


## run_parallel_simulation.jl

Basic script that runs multiple simulation in parallel. 


Note: the <instance_folder> must contain a `params.csv` file, with one set of parameters per row.


Usage:
	
	julia --project=scripts/. scripts/run_parallel_simulation.jl <data_folder> <instance_folder>



## run_parallel_simulation_daily_mobility.jl

Script that runs multiple simulation in parallel. It uses a modified version of the original package, that uses a different mobility matrix (Rij) for each day.

The modified package is in `$ROOT_FOLDER/julia/MMCACovid19custom`

The script looks at `data/config.json` for the keyword `daily_mobility_matrices`, which defines a folder that contains the mobility matrices for each day.

If the mobility matrix for one day is missing, it can be configured whether to use a default matrix, or previous day matrix. 


Usage:
	
	julia --project=scripts/. scripts/run_parallel_simulation_daily_mobility.jl <data_folder> <instance_folder>


Note: the <instance_folder> must contain a `params.csv` file, with one set of parameters per row.



## run_parallel_simulation_seeds.jl

This script is a modified version of `run_parallel_simulation.jl` which uses a different set of seeds for each simulation.


It requires a `seeds.csv` file to be present in the instance/ folder (appart from `params.csv`), with one row per simulation, defining the set of seeds to be used in each simulation.


Usage:
	
	julia --project=scripts/. scripts/run_parallel_simulation_seeds.jl <data_folder> <instance_folder>



## sample_params.py

Sample parameters from a `priors.json` file:

Usage:

	python python/sample_params.py --priors data/priors.json --size 1000 --output experiments/test1000/instance_1/params.csv



## evaluate.py

Evaluates one instance, and outputs a `scores.csv` file, wich is the same as the `params.csv` with an additional `score` column.


Example:
	
	python python/evaluate.py -i <instance_folder> -d <data_folder>

	python python/evaluate.py -i experiments/test1000/instance_1/ -d data/

	# fit with incidence
	python python/evaluate.py -i experiments/test1000/instance_1/ -d data/ --fit incidence

	# fit with deaths
	python python/evaluate.py -i experiments/test1000/instance_1/ -d data/ --fit deaths


	
Usage:

	usage: evaluate.py [-h] --instance-folder INSTANCE_PATH --data-folder DATA_PATH [--config CONFIG_PATH] [--first-day-train FIRST_DAY_TRAIN]
                   [--last-day-train LAST_DAY_TRAIN] [--metric METRIC] [--weights WEIGHTS] [--fit FIT]


## summarize.py

Script that summarizes an instance, and outputs several plots.

Examples:

	python python/summarize.py -i <instance_folder> -d <data_folder>

	# draw main plots
	python python/summarize.py -i experiments/test1000/instance_1/ -d data/
	
	# draw all plots (takes time)
	python python/summarize.py --all -i experiments/test1000/instance_1/ -d data/

	# fit with incidence
	python python/summarize.py -i experiments/test1000/instance_1/ -d data/ --fit incidence

	# fit with deaths
	python python/summarize.py -i experiments/test1000/instance_1/ -d data/ --fit deaths

Usage:

	usage: summarize.py [-h] --instance-folder INSTANCE_PATH --data-folder DATA_PATH [--config CONFIG_PATH] [--output-folder OUTPUT_PATH]
	                    [--simulation SIMULATION] [--first-day-train FIRST_DAY_TRAIN] [--last-day-train LAST_DAY_TRAIN] [--metric METRIC] [--weights WEIGHTS]
	                    [--fit FIT] [--all] [--run RUN]

