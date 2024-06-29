from matplotlib import pyplot as plt
import numpy as np
import sys

seq_2000 = 0.252269
seq_6000 = 2.511786

file_hyb = sys.argv[1]

proc,resolution,rt = np.loadtxt(file_hyb,comments=["#"],delimiter=",",unpack=True)

res2000 = [seq_2000/rt[i] for i in range(len(rt)) if resolution[i]==2000]
res2000_proc = [proc[i] for i in range(len(rt)) if resolution[i]==2000]
res6000 = [seq_6000/rt[i] for i in range(len(rt)) if resolution[i]==6000]
res6000_proc = [proc[i] for i in range(len(rt)) if resolution[i]==6000]

plt.plot(res2000_proc,res2000, label="res=2000")
plt.plot(res6000_proc,res6000, label="res=6000")
plt.ylabel("xSpeedup")
plt.xlabel("Number of Processes")
plt.grid(True)
plt.xticks([4,8,16,32,64])
plt.title("MPI hybrid speedup")
plt.legend()
plt.show()
