from matplotlib import pyplot as plt
import numpy as np
import sys

seq_2000 = 0.252269
seq_6000 = 2.511786

fileblock = sys.argv[1]
filenonblock = sys.argv[2]

procblock,resolutionblock,rtblock = np.loadtxt(fileblock,comments=["#"],delimiter=",",unpack=True)
procnonblock,resolutionnonblock,rtnonblock = np.loadtxt(filenonblock,comments=["#"],delimiter=",",unpack=True)

speedupblock = [seq6000/rtblock[i] for i in range(len(rt1D))]
speedupnonblock = [seq_6000/rtnonblock[i] for i in range(len(rt2D))]

plt.plot(procblock,speedupblock, label="blocking")
plt.plot(procnonblock,speedupnonblock, label="non-blocking")
plt.ylabel("xSpeedup")
plt.xlabel("Process configuration")
plt.grid(True)
plt.xticks([2,24,48,72,96,120,144,168,192])
plt.title("Blocking vs. Non-Blocking")
plt.legend()
plt.show()
