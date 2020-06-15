import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

op = sys.argv[1]
data = pd.read_csv(f'out/{op}_output.csv')

plt.plot(data.iloc[:, 0].values, data.iloc[:, 1].values, color='r')
plt.title(f"{' '.join(op.split('_'))} operation")
plt.xlabel("size")
plt.ylabel("time(ms)")
plt.savefig(f"plots/timeplot_{op}.png")
plt.clf()

try:
    plt.plot(data.iloc[:, 0].values, data.iloc[:, 2].values, color='g')
    plt.title(f"{' '.join(op.split('_'))} operation")
    plt.xlabel("size")
    plt.ylabel("GFLOPs/s")
    plt.savefig(f"plots/GFLOPplot_{op}.png")
except:
    print("GFLOPs benchmarks failed! Skipping...")
