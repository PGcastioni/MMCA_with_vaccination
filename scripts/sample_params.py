import sys, os, argparse
import json
import pandas as pd
import numpy as np
from scipy.stats import uniform
from scipy.stats import multivariate_normal


def sample_uniform(loc, scale, size):
	return uniform.rvs(loc=loc, scale=scale, size=size)


def get_sample_df(priors, size):
	parameters = {
		'id': list(range(1, size+1)),
	}
	columns = ['id']

	for param_def in priors:

		if param_def['type'] == 'constant':
			parameters[param_def['name']] = [param_def['args']['value']] * size
			columns.append(param_def['name'])

		elif param_def['type'] == 'uniform':
			parameters[param_def['name']] = sample_uniform(**param_def['args'], size=size)
			columns.append(param_def['name'])

		else:
			print('WARN: unrecognized parameter type, ', param_def['type'])

	return pd.DataFrame(parameters, columns=columns)


def sample_parameters(priors_path, output_path, size):

	with open(priors_path) as f:
		priors = json.load(f)

	df = get_sample_df(priors, size)

	# print(f'Writing {df.shape[0]} lines to: {output_path}')
	df.to_csv(output_path, index=False)


def calculate_mean_cov(scores_file, parameters_columns, accept_n=100):
	# read csv file with the parameters and a score column
    scores_df = pd.read_csv(scores_file)

    # accept best simulations
    accepted_df = scores_df.sort_values('score').nsmallest(accept_n, 'score')

    # generate multivariate normal distribution
    # https://numpy.org/doc/stable/reference/random/generated/numpy.random.multivariate_normal.html
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.multivariate_normal.html
    accepted_matrix = accepted_df[parameters_columns].values.T    # one row per parameter
    mean = accepted_matrix.mean(axis=1)
    cov = np.cov(accepted_matrix)
    return mean, cov


def sample_multivariate_ABC(output_file, parameters_columns, mean, cov, size):
    params_sample = multivariate_normal.rvs(mean=mean, cov=cov, size=size)
    df = pd.DataFrame(params_sample, columns=parameters_columns)
    df['id'] = range(1, size+1)
    df = df[['id', *parameters_columns]]
    df.to_csv(output_file, index=False)


if __name__ == '__main__':

	ap = argparse.ArgumentParser(description="")

	ap.add_argument('--priors',dest="priors_path",required=True,help='Config file with parameter priors')
	ap.add_argument('--size',dest="size",required=True,type=int,help='Number of samples to generate')
	ap.add_argument('--out',dest="output_path",required=True,type=str,help='Output file path')

	args = ap.parse_args()

	sample_parameters(args.priors_path, args.output_path, args.size)
