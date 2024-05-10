from matplotlib import pyplot as plt
import numpy as np
import sys

file1 = sys.argv[1]
file2 = sys.argv[2]

cores,rt_2000,rt_6000 = np.loadtxt(file1,comments=["#"],unpack=True)
cores_2,rt_2000_2,rt_6000_2 = np.loadtxt(file2,comments=["#"],unpack=True)

plt.plot(cores,rt_2000[0]/rt_2000, label="Resolution=2000")
plt.plot(cores,rt_6000[0]/rt_6000, label="Resolution=6000")
plt.plot(cores_2,rt_2000_2[0]/rt_2000_2, label="Resolution=2000 with ft")
plt.plot(cores_2,rt_6000_2[0]/rt_6000_2, label="Resolution=6000 with ft")
plt.ylabel("xSpeedup")
plt.xlabel("Number of Cores")
plt.axis([1,48,1,plt.ylim()[1]])
plt.grid(True)
plt.xticks([1,2,4,8,16,24,32,48])
plt.title("Parallelization of relax_jacobi with OMP")
plt.legend()
plt.show()
