from matplotlib import pyplot as plt
import numpy as np
import sys

file1 = sys.argv[1]

cores,rt_2000,rt_6000 = np.loadtxt(file1,comments=["#"],unpack=True)

plt.plot(cores,rt_2000[0]/rt_2000, label="Resolution=2000")
plt.plot(cores,rt_6000[0]/rt_6000, label="Resolution=6000")
plt.ylabel("xSpeedup")
plt.xlabel("Number of Cores")
plt.axis([1,48,1,plt.ylim()[1]])
plt.grid(True)
plt.xticks([1,2,4,8,16,24,32,48])
plt.legend()
plt.show()
