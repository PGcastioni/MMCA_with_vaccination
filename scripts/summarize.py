#!/usr/bin/env python
import argparse
import sys
import os
import json
import math
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pylab as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import seaborn as sns
from datetime import datetime, timedelta

scripts_directory = os.path.dirname(os.path.realpath(__file__))
python_directory = os.path.join(scripts_directory, '../python')
sys.path.insert(0, python_directory)

import metrics
import utils
import evaluate


def dump_json(s, break_each=100):
    class NpEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            else:
                return super(NpEncoder, self).default(obj)
    s = json.dumps(s, cls=NpEncoder)
    return "".join(s[i:i+break_each] + "\n" for i in range(0,len(s),break_each))


def save_xarray(a, filename):
    with open(filename, 'w') as f:
        json.dump(a.to_dict(), f)


def plot_grid(R, S_train, S_test, S_full, sim, title, filename, weights=None, seeds=None):
    rows, cols = 5, 4
    fig, axs = plt.subplots(rows, cols, figsize=(20,20))
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    sim_params = params_indexed.loc[sim].to_dict()
    fig.suptitle(f"{title}. sim={sim}, metric={config['evaluate']['metric']}, timestamp={now}, \nparams={dump_json(sim_params)}")

    # plot for spain
    ax = axs[0, 0]
    R.sum(dim='patch').plot(ax=ax, label='Real data')
    S_full.sel(simulation=sim).sum(dim='patch').plot.line(x='date', ax=ax, label='Simulation (full)')
    S_train.sel(simulation=sim).sum(dim='patch').plot.line(x='date', ax=ax, label='Simulation (fitted date range)', linewidth=2, color='darkred')
    ax.set_title(f'Spain, sim={sim}')

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.2, 0.95))

    # plot for each CCAA
    for i, patch in enumerate(ccaa_names, 1):
        row, col = int(i/cols), i%cols
        ax = axs[row, col]
        R.sel(patch=patch).plot(ax=ax)
        S_full.sel(simulation=sim, patch=patch).plot.line(x='date', ax=ax)
        S_train.sel(simulation=sim, patch=patch).plot.line(x='date', ax=ax, linewidth=2, color='darkred')
        title = f'{patch[:22]} sim={sim}'
        if weights:
            title += f' w={weights.get(patch, 0)}'
        if seeds:
            title += f' seed={seeds.get(patch, None)}'
        ax.set_title(title)

    for ax in axs.flat:
        ax.label_outer()

    # save figure
    fig_filename = os.path.join(output_path, filename)
    plt.savefig(fig_filename)
    plt.close()


def plot_top_n_trajectories_grid(config, top_simulations, R, S_full, S_train, S_test, filename, weights=None):
    print(f'- Plot: {filename}')

    # prepare figure
    rows, cols = 5, 4
    fig, axs = plt.subplots(rows, cols, figsize=(20,20))
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    fig.suptitle(f"Top {len(top_simulations)} trajectories, metric={config['evaluate']['metric']}, timestamp={now}, fit_range={config['evaluate']['first_day_train'], config['evaluate']['last_day_train']}")

    # plot for spain
    ax = axs[0, 0]
    ax.axvspan(S_train.date.min().values, S_train.date.max().values, alpha=0.2, color='grey')
    S_full.sel(simulation=top_simulations).sum(dim='patch').plot.line(x='date', ax=ax, linewidth=1, alpha=0.5)
    S_train.sel(simulation=top_simulations).sum(dim='patch').plot.line(x='date', ax=ax, linewidth=1, alpha=0.5, color='darkred')
    R.sum(dim='patch').plot(ax=ax, linewidth=2, color='darkblue')
    if weights:
        ax.set_title(f'Spain, weight={weights.get("Spain", 0)}')
    else:
        ax.set_title(f'Spain')
    ax.get_legend().remove()

    # plot for each CCAA
    for i, patch in enumerate(ccaa_names, 1):
        row, col = int(i/cols), i%cols
        ax = axs[row, col]
        ax.axvspan(S_train.date.min().values, S_train.date.max().values, alpha=0.2, color='grey')
        S_full.sel(simulation=top_simulations, patch=patch).plot.line(x='date', ax=ax, linewidth=1, alpha=0.5)
        S_train.sel(simulation=top_simulations, patch=patch).plot.line(x='date', ax=ax, linewidth=1, alpha=0.5, color='darkred')
        R.sel(patch=patch).plot(ax=ax, linewidth=2, color='darkblue')
        if weights:
            ax.set_title(f'{patch}, weight={weights.get(patch, 0)}')
        else:
            ax.set_title(f'{patch}')
        ax.get_legend().remove()

    for ax in axs.flat:
        ax.label_outer()

    custom_lines = [Line2D([0], [0], color='darkblue', lw=4, label='Real data'),
                    Line2D([0], [0], color='darkred', lw=4, label='Top fitting simulations'),
                    Patch(facecolor='grey', edgecolor='grey', label='Fitting date range')]
    fig.legend(handles=custom_lines, loc='upper center', bbox_to_anchor=(0.2, 0.95))

    # save figure
    fig_filename = os.path.join(output_path, filename)
    plt.savefig(fig_filename)
    plt.close()


def plot_confidence_intervals_boxplot(ax, X):
    """Calculate confidence intervals using median and quartiles (boxplot-like)"""
    low = X.quantile(0.25, dim='simulation')
    mid = X.quantile(0.5, dim='simulation')
    high = X.quantile(0.75, dim='simulation')
    low.plot(ax=ax, color='#ff9999', alpha=0.1)
    mid.plot(ax=ax, color='darkred')
    high.plot(ax=ax, color='#ff9999', alpha=0.1)
    ax.fill_between(low.date, low, high, color='#ff9999', alpha=0.5)


def plot_confidence_intervals_sns(ax, X):
    data = X.to_dataframe('value').reset_index()
    sns.lineplot(data=data, x='date', y='value', ax=ax, color='lightcoral', ci=95)


def plot_confidence_intervals_mean_std(ax, X):
    """faster version of sns.lineplot, look at: https://stackoverflow.com/questions/46125182/is-seaborn-confidence-interval-computed-correctly"""
    mean = X.mean(dim='simulation')
    s = 1.96 * X.std(dim='simulation') / math.sqrt(len(X.simulation))
    high = mean + s
    low = mean - s
    low.plot(ax=ax, color='#ff9999', alpha=0.1)
    mean.plot(ax=ax, color='darkred')
    high.plot(ax=ax, color='#ff9999', alpha=0.1)
    ax.fill_between(low.date, low, high, color='#ff9999', alpha=0.5)


def plot_top_n_prediction_intervals(config, top_simulations, R, S_full, S_train, S_test, filename, plot_conf_interv_fn=None, weights=None):
    print(f'- Plot: {filename}')

    # prepare figure
    rows, cols = 5, 4
    fig, axs = plt.subplots(rows, cols, figsize=(20,20))
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    fig.suptitle(f"Prediction with confidence intervals using top {len(top_simulations)} trajectories, metric={config['evaluate']['metric']}, timestamp={now}, fit_range={config['evaluate']['first_day_train'], config['evaluate']['last_day_train']}")

    # plot for spain
    ax = axs[0, 0]
    ax.axvspan(S_train.date.min().values, S_train.date.max().values, alpha=0.2, color='grey')
    plot_conf_interv_fn(ax, S_full.sel(simulation=top_simulations).sum(dim='patch'))
    R.sum(dim='patch').plot(ax=ax, linewidth=2, color='darkblue')
    if weights:
        ax.set_title(f'Spain, weight={weights.get("Spain", 0)}')
    else:
        ax.set_title(f'Spain')

    # plot for each CCAA
    for i, patch in enumerate(ccaa_names, 1):
        row, col = int(i/cols), i%cols
        ax = axs[row, col]
        ax.axvspan(S_train.date.min().values, S_train.date.max().values, alpha=0.2, color='grey')
        plot_conf_interv_fn(ax, S_full.sel(simulation=top_simulations, patch=patch))
        R.sel(patch=patch).plot(ax=ax, linewidth=2, color='darkblue')
        if weights:
            ax.set_title(f'{patch}, weight={weights.get(patch, 0)}')
        else:
            ax.set_title(f'{patch}')

    for ax in axs.flat:
        ax.label_outer()

    custom_lines = [Line2D([0], [0], color='darkblue', lw=4, label='Real data'),
                    Line2D([0], [0], color='darkred', lw=4, label='Prediction'),
                    Patch(facecolor='grey', edgecolor='grey', label='Fitting date range')]
    fig.legend(handles=custom_lines, loc='upper center', bbox_to_anchor=(0.2, 0.95))

    # save figure
    fig_filename = os.path.join(output_path, filename)
    plt.savefig(fig_filename)
    plt.close()



def plot_ccaa_score_pairplot():

    print('- Plot: ccaa_score_pairplot.png')

    scores = metric(R, S_train)

    df = scores.to_dataframe('score').reset_index().pivot(columns='patch', index='simulation', values='score')

    sns.pairplot(df, vars=ccaa_subset)
    filename = os.path.join(output_path, 'ccaa_score_pairplot.png')
    plt.savefig(filename)



def evaluate_scores():
    print('- Evaluating scores')
    scores_df, seeds_df = evaluate.evaluate_instance(config, instance_path, data_path)
    if seeds_df is not None:
        seeds_df.to_csv(os.path.join(output_path, 'scores_seeds.csv'), index=False)
    if scores_df is not None:
        scores_df = scores_df.sort_values('score')
        scores_df.to_csv(os.path.join(output_path, 'scores.csv'), index=False)
        scores_df.nsmallest(10, 'score').to_csv(os.path.join(output_path, 'scores_top10.csv'), index=False)
        scores_df.nsmallest(100, 'score').to_csv(os.path.join(output_path, 'scores_top100.csv'), index=False)
        scores_df.nsmallest(500, 'score').to_csv(os.path.join(output_path, 'scores_top500.csv'), index=False)
        scores_df.nsmallest(1000, 'score').to_csv(os.path.join(output_path, 'scores_top1000.csv'), index=False)


def plot_accepted_parameters_distrib_pairplot():

    scores_df = pd.read_csv(os.path.join(output_path, 'scores.csv'))

    top_n_params = config.get('top_n_trajectories_params_pairplot', 1000)
    accepted = scores_df.sort_values('score').nsmallest(top_n_params, 'score')

    # plot parameters distributions
    print('- Plot: top_n_params_distributions.png')
    accepted.hist(figsize=(12, 10))
    # plt.tight_layout()
    fig_filename = os.path.join(output_path, 'top_n_params_distributions.png')
    plt.savefig(fig_filename)
    plt.close()

    # plot parameters pairplot
    print('- Plot: top_n_params_pairplot.png')
    plt.figure(figsize=(20,20))
    sns.pairplot(accepted[parameters_columns])
    fig_filename = os.path.join(output_path, 'top_n_params_pairplot.png')
    plt.savefig(fig_filename)
    plt.close()


def plot_seeds_distrib_pairplot():
    try:
        seeds_df = pd.read_csv(os.path.join(output_path, 'scores_seeds.csv'))
    except FileNotFoundError:
        print('- Skipping plot: top_n_seeds_distributions.png, scores_seeds.csv missing from instance')
        return

    top_n_seeds = config.get('top_n_trajectories_params_pairplot', 1000)
    accepted = seeds_df.sort_values('score').nsmallest(top_n_seeds, 'score')

    names = {str(d['patch_idx']): d['patch_name'] for d in config['seeds']}

    # plot seeds distributions
    print('- Plot: top_n_seeds_distributions.png')
    accepted.rename(columns=names).hist(figsize=(15, 10))
    plt.tight_layout()
    fig_filename = os.path.join(output_path, 'top_n_seeds_distributions.png')
    plt.savefig(fig_filename)
    plt.close()

    # plot seeds pairplot
    print('- Plot: top_n_seeds_pairplot.png')
    plt.figure(figsize=(20,20))
    sns.pairplot(accepted[ccaa_subset].rename(columns=names))
    fig_filename = os.path.join(output_path, 'top_n_seeds_pairplot.png')
    plt.savefig(fig_filename)
    plt.close()


def plot_preditions_confidence_intervals():
    top_n_plot = config.get('top_n_trajectories_plot', 100)
    top_n_ci = config.get('top_n_trajectories_ci', 500)

    # top 10% fitting for all CCAA using equal weights
    scores = metric(R, S_train)
    weights = dict(zip(ccaa_names, [1]*len(ccaa_names)))
    w = utils.build_xarry_from_dict(weights, 'patch', normalize=True)
    scores = xr.dot(scores, w)
    scores_df = scores.to_dataframe("score").sort_values('score')
    top_sims_plot = scores_df.nsmallest(top_n_plot, 'score').index.tolist()
    top_sims_ci = scores_df.nsmallest(top_n_ci, 'score').index.tolist()

    plot_top_n_prediction_intervals(config, top_sims_ci, R, S_full, S_train, S_test, weights=weights, filename='top_n_fitted_CCAA_equal_weights_prediction_intervals_mean_std.png', plot_conf_interv_fn=plot_confidence_intervals_mean_std)
    plot_top_n_prediction_intervals(config, top_sims_ci, R, S_full, S_train, S_test, weights=weights, filename='top_n_fitted_CCAA_equal_weights_prediction_intervals_boxplot.png', plot_conf_interv_fn=plot_confidence_intervals_boxplot)
    plot_top_n_trajectories_grid(config, top_sims_plot, R, S_full, S_train, S_test, weights=weights, filename='top_n_fitted_CCAA_equal_weights_trajectories.png')


    # top 10% fitting for all CCAA using population weights
    scores = metric(R, S_train)
    weights = utils.get_population_weights(data_path)
    w = utils.build_xarry_from_dict(weights, 'patch', normalize=True)
    scores = xr.dot(scores, w)
    scores_df = scores.to_dataframe("score").sort_values('score')
    top_sims_plot = scores_df.nsmallest(top_n_plot, 'score').index.tolist()
    top_sims_ci = scores_df.nsmallest(top_n_ci, 'score').index.tolist()

    plot_top_n_prediction_intervals(config, top_sims_ci, R, S_full, S_train, S_test, weights=weights, filename='top_n_fitted_CCAA_population_weights_prediction_intervals_mean_std.png', plot_conf_interv_fn=plot_confidence_intervals_mean_std)
    plot_top_n_prediction_intervals(config, top_sims_ci, R, S_full, S_train, S_test, weights=weights, filename='top_n_fitted_CCAA_population_weights_prediction_intervals_boxplot.png', plot_conf_interv_fn=plot_confidence_intervals_boxplot)
    plot_top_n_trajectories_grid(config, top_sims_plot, R, S_full, S_train, S_test, weights=weights, filename='top_n_fitted_CCAA_population_weights_trajectories.png')


    # top 10% fitting for Spain
    scores = metric(R, S_train)
    scores_df = scores.sum('patch').to_dataframe("score").sort_values('score')
    top_sims_plot = scores_df.nsmallest(top_n_plot, 'score').index.tolist()
    top_sims_ci = scores_df.nsmallest(top_n_ci, 'score').index.tolist()

    plot_top_n_prediction_intervals(config, top_sims_ci, R, S_full, S_train, S_test, weights={'Spain': 1}, filename='top_n_fitted_for_Spain_prediction_intervals_mean_std.png', plot_conf_interv_fn=plot_confidence_intervals_mean_std)
    plot_top_n_prediction_intervals(config, top_sims_ci, R, S_full, S_train, S_test, weights={'Spain': 1}, filename='top_n_fitted_for_Spain_prediction_intervals_boxplot.png', plot_conf_interv_fn=plot_confidence_intervals_boxplot)
    plot_top_n_trajectories_grid(config, top_sims_plot, R, S_full, S_train, S_test, weights={'Spain': 1}, filename='top_n_fitted_for_Spain_trajectories.png')


def plot_best_weighed_trajectory():
    print('- Plot: weighted_best_sim.png')

    # calculate score for each ccaa, for each simulation
    # 'scores' has two dimmensions: (ccaa, sim)
    scores = metric(R, S_train)

    # apply weights
    weights = utils.get_weights(config, data_path)
    w = utils.build_xarry_from_dict(weights, 'patch', normalize=True)
    scores = xr.dot(scores, w)

    # best sim
    best_sim = int(min(scores).simulation)
    best_sim_score = float(scores.sel(simulation=best_sim))

    # seeds
    seeds_df = utils.load_seeds(config, instance_path, data_path, rename=True)
    if seeds_df is not None:
        seeds = seeds_df.set_index('id').loc[best_sim][ccaa_names].to_dict()
        seeds = {k: round(v, 2) for k,v in seeds.items()}
    else:
        seeds = None

    # plot
    plot_grid(R, S_train, S_test, S_full, best_sim, filename='weighted_best_sim.png', title=f'One simulation that fits better all CCAA at the same time, using weights. \nweights={weights}\nSCORE={best_sim_score}', weights=weights, seeds=seeds)

    # save weighted scores
    filename = os.path.join(output_path, 'weighted_best_sim.params.json')
    save_xarray(scores, filename)


def plot_best_trajectory_for_each_CCAA_single_plot():

    print('- plotting for each CCAA the best trajectory in one single grid')

    # find best simulation for each CCAA

    best_sim_by_ccaa = {}

    for i, patch in enumerate(ccaa_names, 1):
        x = R.sel(patch=patch)
        y = S_train.sel(patch=patch)

        scores = metric(x, y) # calculate score for each simulation
        best_sim = int(min(scores).simulation)

        best_sim_by_ccaa[patch] = best_sim


    # plot
    rows, cols = 5, 4

    fig, axs = plt.subplots(rows, cols, figsize=(20,20))
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    fig.suptitle(f"For each CCAA, plot the best simulation. metric={config['evaluate']['metric']}, timestamp={now}")

    # plot for Spain (select best overall simulation)
    ax = axs[0, 0]

    x = R.sum(dim='patch') # sum to spain level
    y = S_train.sum(dim='patch')

    scores = metric(x, y) # calculate score for each simulation

    best_sim_spain = int(min(scores).simulation)

    R.sum(dim='patch').plot(ax=ax, label="Real data")
    S_full.sel(simulation=best_sim_spain).sum(dim='patch').plot(ax=ax, label="Simulation (full)")
    S_train.sel(simulation=best_sim_spain).sum(dim='patch').plot(ax=ax, label="Simulation (fitted date range)", linewidth=2, color='darkred')

    ax.set_title(f'Spain, sim={best_sim_spain}')

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.2, 0.95))

    # plot best simulation for each CCAA
    for i, (patch, sim) in enumerate(best_sim_by_ccaa.items(), 1):
        row, col = int(i/cols), i%cols
        ax = axs[row, col]

        R.sel(patch=patch).plot(ax=ax)
        S_full.sel(simulation=sim, patch=patch).plot(ax=ax)
        S_train.sel(simulation=sim, patch=patch).plot(ax=ax, linewidth=2, color='darkred')

        ax.set_title(f'{patch}, sim={sim}')

    for ax in axs.flat:
        ax.label_outer()

    # save figure
    fig_filename = os.path.join(output_path, 'best_sim_for_each_ccaa.png')
    plt.savefig(fig_filename)
    plt.close()


def plot_best_trajectory_for_each_CCAA_multiple_plots():
    print('- plotting one grid for each CCAA best trajectory')

    # find best simulation for each CCAA

    best_sim_by_ccaa = {}

    for i, patch in enumerate(ccaa_names, 1):
        x = R.sel(patch=patch)
        y = S_train.sel(patch=patch)

        scores = metric(x, y) # calculate score for each simulation
        best_sim = int(min(scores).simulation)

        best_sim_by_ccaa[patch] = best_sim


    for patch, sim in best_sim_by_ccaa.items():
        print('   -', patch)
        plot_grid(R, S_train, S_test, S_full, sim=sim, title=f'Best simulation for patch: {patch}', filename=f'best_sim_for_patch_{patch}.png')


def plot_best_trajectory_for_Spain():
    print('- Plot: best_sim_for_Spain.png')

    x = R.sum(dim='patch') # sum to spain level
    y = S_train.sum(dim='patch')
    scores = metric(x, y) # calculate score for each simulation
    best_sim_spain = int(min(scores).simulation)

    plot_grid(R, S_train, S_test, S_full, sim=best_sim_spain, title=f'Best simulation for Spain', filename='best_sim_for_Spain.png')


def plot_mean_trajectory():
    print('- Plot: mean_trajectory.png')

    # plot
    rows, cols = 5, 4

    fig, axs = plt.subplots(rows, cols, figsize=(20,20))
    fig.suptitle(f"Plot mean trajectory. metric={config['evaluate']['metric']}")

    # plot for Spain (select best overall simulation)
    ax = axs[0, 0]

    R.sum(dim='patch').plot(ax=ax, label="Real data")
    S_full.sum(dim='patch').mean(dim='simulation').plot(ax=ax, label="Mean trajectory Simulation (full)")
    S_train.sum(dim='patch').mean(dim='simulation').plot(ax=ax, label="Mean trajectory Simulation (fitted date range)", linewidth=2, color='darkred')

    ax.set_title(f'Spain')

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.2, 0.95))

    # plot best simulation for each CCAA
    for i, patch in enumerate(ccaa_names, 1):
        row, col = int(i/cols), i%cols
        ax = axs[row, col]

        R.sel(patch=patch).plot(ax=ax)
        S_full.sel(patch=patch).mean(dim='simulation').plot(ax=ax)
        S_train.sel(patch=patch).mean(dim='simulation').plot(ax=ax, linewidth=2, color='darkred')

        ax.set_title(f'{patch}')

    for ax in axs.flat:
        ax.label_outer()

    # save figure
    fig_filename = os.path.join(output_path, 'mean_trajectory.png')
    plt.savefig(fig_filename)
    plt.close()


def quick_plots():

    ######################################
    # plot trajectory that fits better all weighted CCAA
    ######################################

    plot_best_weighed_trajectory()


    ######################################
    # plot trajectory that fits better Spain
    ######################################

    plot_best_trajectory_for_Spain()


    ######################################
    # score each simulation
    ######################################

    # evaluate_scores()


    ######################################
    # plot parameters distributions (ABC)
    ######################################

    # plot_accepted_parameters_distrib_pairplot()


def all_plots():

    ######################################
    # plot trajectory that fits better all weighted CCAA
    ######################################

    plot_best_weighed_trajectory()


    ######################################
    # plot trajectory that fits better Spain
    ######################################

    plot_best_trajectory_for_Spain()


    ######################################
    # score each simulation
    ######################################

    evaluate_scores()

    
    ######################################
    # plot parameters distributions (ABC)
    ######################################

    plot_accepted_parameters_distrib_pairplot()


    ######################################
    # plot seeds distributions (ABC)
    ######################################

    plot_seeds_distrib_pairplot()


    ######################################
    # grid plot: mean of all trajectories (check lockdowns effect)
    ######################################

    plot_mean_trajectory()


    ######################################
    # plot the best trajectory for each CCAA in a single grid
    ######################################

    plot_best_trajectory_for_each_CCAA_single_plot()


    ######################################
    # now plot one grid for each CCAA, showing the best simulation 
    # for that CCAA and how it performs for the rest of CCAA
    ######################################

    plot_best_trajectory_for_each_CCAA_multiple_plots()


    ######################################
    # show predictions and confidence intervals using 
    # best simulations for fitting range
    ######################################

    plot_preditions_confidence_intervals()


    ######################################
    # compare performance for different CCAA
    # (score distribution pairplot for several CCAA)
    ######################################

    plot_ccaa_score_pairplot()



def plots_for_simulation(simulation_id):
    plot_grid(R, S_train, S_test, S_full, sim=simulation_id, title=f'Simulation: {simulation_id}', filename=f'simulation_{simulation_id}.png')



if __name__ == '__main__':

    ap = argparse.ArgumentParser(description="")

    ap.add_argument('--instance-folder', '-i',dest="instance_path",required=True,help='Path to the instance folder (experiment folder)')
    ap.add_argument('--data-folder', '-d', dest="data_path",required=True,help='Path to the data folder')
    ap.add_argument('--config', '-c',dest="config_path",required=False,help='Path to the config file (JSON)')
    ap.add_argument('--output-folder', '-o',dest="output_path",required=False,help='Output folder')
    ap.add_argument('--simulation',dest="simulation",default=None,type=int,required=False,help='Draw plots for this simulation id.')
    ap.add_argument('--first-day-train',dest="first_day_train",help='')
    ap.add_argument('--last-day-train',dest="last_day_train",help='')
    ap.add_argument('--metric',dest="metric",default="RMSE",help='Options: RMSE, MSE, MAE, MSRE, RMSRE')
    ap.add_argument('--weights',dest="weights",default="ones",help='Options: ones, population')
    ap.add_argument('--fit',dest="fit",default="deaths",help='Options: deaths, incidence')
    ap.add_argument('--all',dest="all_plots",action="store_true",required=False,help='Draw all plots. If ommited, only the main plots will be plotted.')
    ap.add_argument('--run',dest="run",help="Run specific plot. Provide the function name")
    args = ap.parse_args()


    # folder containing simulation outputs
    instance_path = args.instance_path

    # folder containing ccaa_ids.csv, ccaa_covid19_fallecidos_long.csv
    data_path = args.data_path

    # load config.json
    config_path = args.config_path or os.path.join(args.instance_path, 'output/config.json')
    with open(config_path) as f:
        config = json.load(f)

    config['evaluate'] = {}
    config['evaluate']['metric'] = args.metric
    config['evaluate']['weights'] = args.weights
    config['evaluate']['real_data_filename'] = config['data']['fit_alternatives'][args.fit]['real_data_filename']
    config['evaluate']['sim_data_pattern'] = config['data']['fit_alternatives'][args.fit]['sim_data_pattern']
    config['evaluate']['skip_first_day_simulation'] = config['data']['fit_alternatives'][args.fit]['skip_first_day_simulation']
    config['evaluate']['first_day_train'] = args.first_day_train or config['simulation']['first_day_simulation']
    config['evaluate']['last_day_train'] = args.last_day_train or config['simulation']['last_day_simulation']

    # metric to be used, one of MSRE, RMSE, MSE, RMSRE, etc
    metric = metrics.methods[config['evaluate']['metric']]

    # columns in the params.csv to use as parameters and for which show their distribution
    default_parameters_columns = ['beta', 'scale_beta', 't_inc', 'scale_ea', 't_i', 'delta', 'phi', 'scale_s']
    parameters_columns = config.get('summarize', {}).get('parameters_columns', default_parameters_columns)
    default_seeds_columns = ["1573", "360", "2141", "2309", "2628", "2678", "2493"]
    seeds_columns = config.get('summarize', {}).get('seeds_columns', default_seeds_columns)

    default_ccaa_subset = ["Comunidad de Madrid", "Cataluña", "Canarias", "Andalucía", "Comunitat Valenciana", "Castilla y León", "Castilla-La Mancha"]
    ccaa_subset = config.get('summarize', {}).get('pairplot_ccaa_subset', default_ccaa_subset)

    # folder to output the files generated by this script
    output_path = args.output_path or os.path.join(instance_path, 'output/summarize/')
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    print(f'saving output to: {output_path}')

    # load ccaa data
    ccaa_df, ccaa_names, ccaa_ids, ccaa_pop, ccaa_area = utils.load_ccaa_data(data_path)

    # load parameters csv
    print(instance_path)
    params_df = utils.load_parameters(instance_path)
    params_indexed = params_df.set_index('id')
    params_indexed['id'] = params_indexed.index

    print("Parameters: (n_sim, n_params) =", params_df.shape)

    # load simulation results and real covid cases data
    S_full, S_train, S_test, R = utils.load_simulations_and_real(config, instance_path, data_path)


    if args.simulation is not None:
        plots_for_simulation(args.simulation)
    elif args.run is not None:
        locals()[args.run]()
    elif args.all_plots:
        all_plots()
    else:
        quick_plots()
