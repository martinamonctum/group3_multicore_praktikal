import json
import numpy as np
import matplotlib.pyplot as plt 
from numpy import loadtxt
from matplotlib import rc
import pandas as pd
from tabulate import tabulate

def load_file(file):
    O2_data = json.load(file)
    O2_data = O2_data["threads"]["regions"]["computation"]

    O2_l2_percent = float(O2_data["PAPI_L2_TCM"]) / float(O2_data["PAPI_L2_TCA"]) 
    O2_l3_percent = float(O2_data['PAPI_L3_TCM']) / float(O2_data["PAPI_L3_TCA"])

    O2_percents = [O2_l2_percent,O2_l3_percent]
    return O2_percents

# Opening JSON file
# f1 = open('papi_500/papi_hl_output_O2/o2-500')
# f2 = open('papi_1500/o2-1500')
# f3 = open('papi_2500/o2-2500')
# f4 = open('papi_3500/o2-3500')
# f5 = open('papi_4500/o2-4500')
f6 = open('papi_500/papi_hl_output_o3/O3-500')
f7 = open('papi_1500/O3-1500')
f8 = open('papi_2500/O3-2500')
f9 = open('papi_3500/O3-3500')
f10 = open('papi_4500/O3-4500')

data = loadtxt("FlopData-O3.txt")
resolution = data[:,0]
flops = data[:,1]
papi_flops = data[:,2]

g = globals()
l2_list =[]
l3_list =[]
for i in range(6,11,1):
    file_var = 'f{}'.format(i)
    data = load_file(g[file_var])
    l2_list.append(data[0])
    l3_list.append(data[1])

# df = pd.DataFrame([l2_list,l3_list],columns=resolution,index=['L2','L3'])
# df = df.map(lambda x: f'{int(x):,d}' if isinstance(x, float) else x)
# print(tabulate(df, headers='keys', tablefmt='psql'))
fig,ax1= plt.subplots()
ax1.set_xlabel("Resolutions")
ax1.set_ylabel("Percentage of Cache Misses", color='blue')
ax1.plot(resolution, l2_list, color ='blue', label = "L2 miss rate")
ax1.plot(resolution, l3_list, color ='orange', label = "L3 miss rate")
ax1.tick_params(axis='y', labelcolor='blue')

ax2 = ax1.twinx()
ax2.set_ylabel('MFlops', color='red')
ax2.plot(resolution, papi_flops, color ='red', label = "MFlops PAPI")
ax2.plot(resolution, flops, color ='black', label = "MFlops Hardcoded")
ax2.tick_params(axis='y', labelcolor='red')
# fig.tight_layout()
fig.legend()
plt.title("Group 3" + '\n'+ "Best Flag Combination MFlops and Cache Miss Rates with Respect to Different Resolutions", y=1)
plt.show()
