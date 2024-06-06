from matplotlib import pyplot as plt
import numpy as np
import sys

seq_2000 = 0.252269
seq_6000 = 2.511786

file1D = sys.argv[1]
file2D = sys.argv[2]

proc1D,resolution1D,rt1D = np.loadtxt(file1D,comments=["#"],delimiter=",",unpack=True)
proc2D,resolution2D,rt2D = np.loadtxt(file2D,comments=["#"],delimiter=",",unpack=True)

res2000_1D = [seq_2000/rt1D[i] for i in range(len(rt1D)) if resolution1D[i]==2000]
res2000_1D_proc = [proc1D[i] for i in range(len(rt1D)) if resolution1D[i]==2000]
res6000_1D = [seq_6000/rt1D[i] for i in range(len(rt1D)) if resolution1D[i]==6000]
res6000_1D_proc = [proc1D[i] for i in range(len(rt1D)) if resolution1D[i]==6000]
res2000_2D = [seq_2000/rt2D[i] for i in range(len(rt2D)) if resolution2D[i]==2000]
res2000_2D_proc = [proc2D[i] for i in range(len(rt2D)) if resolution2D[i]==2000]
res6000_2D = [seq_6000/rt2D[i] for i in range(len(rt2D)) if resolution2D[i]==6000]
res6000_2D_proc = [proc2D[i] for i in range(len(rt2D)) if resolution2D[i]==6000]

plt.plot(res2000_1D_proc,res2000_1D, label="1D res=2000")
plt.plot(res6000_1D_proc,res6000_1D, label="1D res=6000")
plt.plot(res2000_2D_proc,res2000_2D, label="2D res=2000")
plt.plot(res6000_2D_proc,res6000_2D, label="2D res=6000")
plt.ylabel("xSpeedup")
plt.xlabel("Number of Processes")
plt.grid(True)
plt.xticks([2,24,48,72,96,120,144,168,192])
plt.title("MPI speedup")
plt.legend()
plt.show()
