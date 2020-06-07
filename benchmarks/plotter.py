import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

op = sys.argv[1]
data = pd.read_csv(f'out/{op}_output.csv')

plt.plot(data.iloc[:, 0].values, data.iloc[:, 1].values)
plt.title(f"{op} operation")
plt.xlabel("size")
plt.ylabel("time(ms)")
plt.savefig(f"plots/timeplot_{op}.png")
