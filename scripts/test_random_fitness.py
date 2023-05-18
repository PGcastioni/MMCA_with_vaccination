import sys
import pandas as pd
import numpy as np


def main():
    instance_dir = sys.argv[1]
    csv_fname = instance_dir + "params.csv"
    df = pd.read_csv(csv_fname)
    n = df.shape[0]
    for i in range(20000):
        fitness = np.random.rand(n)
    result = ";".join([str(i) for i in fitness])

    csv_fname = instance_dir + "out.txt"
    with open(csv_fname, 'w') as fh:
        fh.write(result)

main()