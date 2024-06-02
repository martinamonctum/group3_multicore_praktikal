from matplotlib import pyplot as plt
import numpy as np
import sys

file1 = sys.argv[1]

byte,core,socket,node = np.loadtxt(file1,comments=["#"],unpack=True)

Mbytes = [2**i/10**6 for i in byte]

core_mb_s = [Mbytes[i]/core[i] for i in range(len(core))]
socket_mb_s = [Mbytes[i]/socket[i] for i in range(len(core))]
node_mb_s = [Mbytes[i]/node[i] for i in range(len(core))]

plt.plot(byte,core_mb_s, label="Different Core")
plt.plot(byte,socket_mb_s, label="Different Socket")
plt.plot(byte,node_mb_s, label="Different Node")
plt.ylabel("Mbyte/s")
plt.xlabel("Byte (Powers of two)")
plt.grid(True)
plt.xticks([0,1,2,4,8,16,24])
plt.title("Ping-Pong for different process assignments")
plt.legend()
plt.show()
