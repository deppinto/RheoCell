import sys
import numpy as np
from math import *
import os.path
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage

from matplotlib import cm

if len(sys.argv)!=3:
        print(sys.argv[0]," [topology file] [conf file]")
        sys.exit(1)


tfile=open(sys.argv[1],"r")
N=int(tfile.readline().split()[0])
species=[]
line=tfile.readline()
lline=line.split()
for l in lline:
        species.append(int(l))
tfile.close()

numspecies=len(set(species))


cfile=open(sys.argv[2],"r")
header=cfile.readline().split()
t=int(header[2])

header=cfile.readline().split()
lx=int(float(header[2]))
ly=int(float(header[3]))

sim_grid=np.zeros(N*lx*ly)
pt_num=0
x=np.arange(0,lx,1)
y=np.arange(0,ly,1)
Z=[[0 for q in range(lx)] for k in range(ly)]
fig = plt.figure(figsize=(6,6))
start_value = 9
for line in cfile:
    area=0
    out_area=0
    words=line.split()
    Z=[[0 for q in range(lx)] for k in range(ly)]
    track_problem = 0
    for i in range(start_value,len(words),2):
        site=int(float(words[i]))
        value=float(words[i+1])
        sim_grid[pt_num*site]=value
        yy=int(site/lx)
        xx=site-int(yy*lx)
        if value>0.5:
            Z[yy][xx]=value
        else:
            out_area=value*value
        area+=value*value
        if value>1.5 or value<-0.5:
            print("phase field is not in [0,1]!: ",xx,yy,value,pt_num)
            track_problem+=1

    if abs(1-area/(pi*8*8)>0.2):
        print("area is not conserved!")

    if out_area/area > 0.9:
        print("cell is leaking!")
    

    X, Y = np.meshgrid(x, y)
    #norm = cm.colors.Normalize(vmax=abs(Z).max(), vmin=-abs(Z).max())
    if track_problem==0:
        cmap = cm.binary
    else:
        cmap = cm.cool

    #axs = _axs.flatten()

    step = 0.01
    m = np.amax(Z)
    if m<0.000001:
        continue
    levels = np.arange(0.0, m, step) + step
    cset1 = plt.contourf(X, Y, Z, levels, cmap=cmap, alpha=0.5)
    #cset1 = plt.imshow(Z,interpolation='gaussian', cmap=cmap)
    #fig.colorbar(cset1)

    print(pt_num, area)
    #if pt_num==3:
        #print(float(words[2]),float(words[3]))
    pt_num+=1


ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
fig.tight_layout()
plt.show()
