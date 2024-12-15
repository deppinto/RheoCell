import sys
import numpy as np
from math import *
import os.path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgba, to_hex, Normalize
import numpy as np
import scipy.ndimage

from matplotlib import cm

if len(sys.argv)!=4:
    print(sys.argv[0]," [topology file] [conf file] [1:save conf; 2:make plot]")
    sys.exit(1)

variable=int(float(sys.argv[3])) 

def set_walls(lx,ly, walls):
    for y in range(ly):
        for x in range(lx):
            k = x + y * lx
            walls[k] = exp(-float(y)/5) + exp(-float(ly-y-1)/5)


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

walls = [0. for i in range(lx*ly)]
set_walls(lx,ly,walls)

x=np.arange(0,lx,1)
y=np.arange(0,ly,1)
Z_x=[[0 for q in range(lx)] for k in range(ly)]
Z_y=[[0 for q in range(lx)] for k in range(ly)]
Z=[[0 for q in range(lx)] for k in range(ly)]
fig = plt.figure(figsize=(6,6))
start_value=0
for line in cfile:
    words=line.split()
    Z_x=[[0 for q in range(lx)] for k in range(ly)]
    Z_y=[[0 for q in range(lx)] for k in range(ly)]
    for i in range(start_value,len(words),3):
        site=int(float(words[i]))
        value_x=float(words[i+1])
        value_y=float(words[i+2])
        yy=int(site/lx)
        xx=site-int(yy*lx)

        Z_x[yy][xx]=value_x
        Z_y[yy][xx]=value_y
        Z[yy][xx]=sqrt(value_x*value_x+value_y*value_y)
        cset1 = plt.arrow(xx, yy, value_x, value_y, width=0.075, color='k')

    #z_min, z_max = -np.abs(Z).max(), np.abs(Z).max()
    z_min, z_max = 0., np.abs(Z).max()
    X, Y = np.meshgrid(x, y)
    #cset1 = plt.imshow(Z, cmap='hot', interpolation='nearest')
    #cset1 = plt.pcolormesh(X, Y, Z, cmap='RdBu', vmin=z_min, vmax=z_max)



ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
ax.set_xlim([0, lx])
ax.set_ylim([0, ly])
#fig.tight_layout()
if variable==1:
    plt.savefig('frame.png')
if variable==2:
    plt.show()
