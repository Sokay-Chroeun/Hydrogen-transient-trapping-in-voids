import numpy as np
import matplotlib.pyplot as plt

def readData(FileName):
    txt = [];X1 = [];X2 = [];X3 = [];X4 = [];
    for Line in open(FileName,'r'):
        t = Line.split()
        txt.append(t)
    for i in range(len(txt)-1):
        X1.append(float(txt[i+1][0]))
        X2.append(float(txt[i+1][1]))
        X3.append(float(txt[i+1][2]))
        X4.append(float(txt[i+1][3]))
    return (X1,X2,X3,X4)

time, cL, ct, cH = readData("Job-1-Qtot.txt")

plt.figure(1,figsize=(6,5))
plt.plot(time, cL,label='H in lattice sites')
plt.plot(time, ct,label='H in voids')
plt.plot(time, cH,label='Total H')
plt.legend(fontsize=12)
plt.xlabel('Time (s)',fontsize=14)
plt.ylabel("H retention (H atoms)",fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid()
plt.savefig("mass_conservation.png", bbox_inches='tight',dpi=500)