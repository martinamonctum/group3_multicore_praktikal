from matplotlib import pyplot as plt
import numpy as np
import sys
import pandas as pd

seq_2000 = 0.252269
seq_6000 = 2.511786

fileblock = sys.argv[1]
filenonblock = sys.argv[2]
fileothernb = sys.argv[3]

df_block = pd.read_csv(fileblock,delimiter=",",dtype={"config": "string","res": np.int32,"rt":np.float64})
df_nonblock = pd.read_csv(filenonblock,delimiter=",",dtype={"config": "string","res": np.int32,"rt":np.float64})
df_othernb = pd.read_csv(fileothernb,delimiter=",",dtype={"config": "string","res": np.int32,"rt":np.float64})

rt_block = df_block.loc[:,"rt"].to_list()
rt_nonblock = df_nonblock.loc[:,"rt"].to_list()
rt_othernb = df_othernb.loc[:,"rt"].to_list()

speedupblock = [seq_6000/rt_block[i] for i in range(len(rt_block))]
speedupnonblock = [seq_6000/rt_nonblock[i] for i in range(len(rt_nonblock))]
speedupothernb = [seq_6000/rt_othernb[i] for i in range(len(rt_othernb))]

plt.plot(df_block.loc[:,"config"].to_list(),speedupblock, label="blocking")
plt.plot(df_nonblock.loc[:,"config"].to_list(),speedupnonblock, label="non-blocking")
plt.plot(df_othernb.loc[:,"config"].to_list(),speedupothernb, label="other non-blocking")
plt.ylabel("xSpeedup")
plt.xlabel("Process configuration")
plt.grid(True)
plt.title("Blocking vs. Non-Blocking")
plt.legend()
plt.show()
