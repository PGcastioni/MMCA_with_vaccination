import argparse
import json
import sys
import os

scripts_directory = os.path.dirname(os.path.realpath(__file__))
python_directory = os.path.join(scripts_directory, '../python')
sys.path.insert(0, python_directory)

from evaluate import evaluate_instance


def main():

    ap = argparse.ArgumentParser(description="")

    ap.add_argument('--instance-folder', '-i',dest="instance_path",required=True,help='Path to the instance folder (experiment folder)')
    ap.add_argument('--data-folder', '-d', dest="data_path",required=True,help='Path to the data folder')
    ap.add_argument('--config', '-c',dest="config_path",required=False,help='Path to the config file (JSON)')
    ap.add_argument('--first-day-train',dest="first_day_train",help='')
    ap.add_argument('--last-day-train',dest="last_day_train",help='')
    ap.add_argument('--metric',dest="metric",default="RMSE",help='')
    ap.add_argument('--weights',dest="weights",default="ones",help='')
    ap.add_argument('--fit',dest="fit",default="deaths",help='')

    args = ap.parse_args()

    # folder containing simulation outputs
    instance_path = args.instance_path

    # '../data', folder containing ccaa_ids.csv, ccaa_covid19_fallecidos_long.csv
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

    scores_df = evaluate_instance(config, instance_path, data_path)

    if scores_df is not None:
        output_file = os.path.join(instance_path, 'scores.csv')
        scores_df.to_csv(output_file, index=False)


if __name__ == '__main__':
    main()
