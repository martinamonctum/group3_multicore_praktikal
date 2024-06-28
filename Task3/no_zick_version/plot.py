import json
import numpy as np
import matplotlib.pyplot as plt 
from numpy import loadtxt
from matplotlib import rc
import pandas as pd

optimized_flops = loadtxt("Task3/data/FlopData.txt")
old_flop_data = loadtxt("2nd-task/FlopData-O3.txt")
resolution = optimized_flops[:,0]
op_flops = optimized_flops[:,1]
old_flops = old_flop_data[:,1]
flop_percents = []
for i in range(len(old_flops)):
    flop_percents.append((op_flops[i] / old_flops[i])*100)


def load_file():
    optimized = json.load(optimizedf)
    optimized = optimized["threads"]["regions"]["computation"]
    best = json.load(bestf)
    best = best["threads"]["regions"]["computation"]

    l2_percents.append((float(optimized["PAPI_L2_TCM"]) / float(best["PAPI_L2_TCM"]))*100)
    l3_percents.append((float(optimized['PAPI_L3_TCM']) / float(best["PAPI_L3_TCM"]))*100)

l2_percents = []
l3_percents = []
for i in resolution:
    bestf = open('Task3/data/O3-{}'.format(int(i)))
    optimizedf = open("Task3/data/optimized_{}".format(int(i)))
    load_file()

# set width of bar 
barWidth = 0.25
fig = plt.subplots(figsize =(12, 8)) 

# Set position of bar on X axis 
br1 = np.arange(len(resolution)) 
br2 = [x + barWidth for x in br1] 
br3 = [x + barWidth for x in br2] 
# Make the plot
plt.bar(br1, flop_percents, color ='r', width = barWidth, 
        edgecolor ='grey', label ='MFlops/s Percentages') 
plt.bar(br2, l2_percents, color ='g', width = barWidth, 
        edgecolor ='grey', label ='L2 Miss Percentages') 
plt.bar(br3, l3_percents, color ='b', width = barWidth, 
        edgecolor ='grey', label ='L3 Miss Percentages') 
plt.axhline(y=100, color='grey', linestyle='-')
# Adding Xticks 
plt.xlabel('Resolutions', fontweight ='bold', fontsize = 15) 
plt.ylabel('Percentage of Old Performance', fontweight ='bold', fontsize = 15) 
plt.xticks([r + barWidth for r in range(len(resolution))], 
        ['500', '1500', '2500', '3500', '4500'])
plt.title("Group 3 \n" + "Percentages of L2, L3 Misses and MFlops/s Compared to the Unoptimized Version",fontweight = 'bold', fontsize = 17)
plt.legend()
plt.show() 