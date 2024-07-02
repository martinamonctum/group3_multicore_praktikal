import json
import numpy as np
import matplotlib.pyplot as plt 
from numpy import loadtxt

data = loadtxt("data.txt")

grouped_data = {}
for depth, value in data:
    if depth not in grouped_data:
        grouped_data[depth] = []
    grouped_data[depth].append(value)

sorted_depths = sorted(grouped_data.keys())
values = [grouped_data[depth] for depth in sorted_depths]

fig, ax = plt.subplots()

num_groups = len(sorted_depths)

num_bars_per_group = max(len(v) for v in values)

bar_width = 0.2

indices = np.arange(num_groups)
bar_positions = [indices + i * bar_width for i in range(num_bars_per_group)]
labels = ['start position', 'mid game 1', 'mid game 2', 'end game']
colors = ['blue', 'red', 'orange', 'brown']
for i in range(num_bars_per_group):
    heights = [values[j][i] if i < len(values[j]) else 0 for j in range(num_groups)]
    ax.bar(bar_positions[i], heights, bar_width,color=colors[i], label=labels[i])

ax.set_xticks(indices + bar_width * (num_bars_per_group - 1) / 2)
ax.set_xticklabels(sorted_depths)

# Add labels, legend, and title
ax.set_xlabel('Depth')
ax.set_ylabel('Evaluations per Second')
ax.set_title('Evaluations per Second by Depth and different positions\n Group 3')
ax.legend()

# Show the plot
plt.show()