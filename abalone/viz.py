from matplotlib import pyplot as plt
import numpy as np
import sys

threads = [1,2,4,6,8,12,24]
evals = [587364,857262,1223482,1472316,1934179,1743209,1539588]

speedup = [evals[i]/evals[0] for i in range(len(evals))]

plt.plot(threads,speedup)
plt.ylabel("xSpeedup")
plt.xlabel("Number of threads")
plt.grid(True)
plt.xticks([0,4,8,12,16,20,24])
plt.title("Minimax (Loop) Speedup")
plt.legend()
plt.show()
