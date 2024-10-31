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

sim_grid=np.zeros(N*lx*ly)
pt_num=0
x=np.arange(0,lx,1)
y=np.arange(0,ly,1)
Z=[[0 for q in range(lx)] for k in range(ly)]
fig = plt.figure(figsize=(6,6))
start_value = 9
totphi=[0. for i in range(lx*ly)]
for line in cfile:
    area=0
    out_area=0
    words=line.split()
    Z=[[0 for q in range(lx)] for k in range(ly)]
    track_problem = 0
    for i in range(start_value,len(words),2):
        site=int(float(words[i]))
        value=float(words[i+1])
        sim_grid[site+pt_num*lx*ly]=value
        yy=int(site/lx)
        xx=site-int(yy*lx)

        if value>0.5:
            Z[yy][xx]=value
        else:
            out_area=value*value

        area+=value*value
        if value>1.5 or value<-0.5:
            print("phase field is not in [0,1]!: ", pt_num, xx, yy, value)
            track_problem+=1

    if track_problem==0:
        cmap = cm.binary
    else:
        cmap = cm.cool

    if abs(1-area/(pi*8*8)>0.2):
        print("area is not conserved: ", pt_num, area)
        cmap=cm.winter

    if out_area/area > 0.9:
        print("cell is leaking: ", pt_num, area, out_area)
        cmap=cm.autumn
    

    X, Y = np.meshgrid(x, y)
    #axs = _axs.flatten()

    step = 0.01
    m = np.amax(Z)
    if m<0.000001:
        continue

    levels = np.arange(0.0, m, step) + step
    if pt_num<3:
        cset1 = plt.contour(X, Y, Z, levels, cmap=cmap, alpha=0.5)
    else:
        cset1 = plt.contour(X, Y, Z, levels, levels=[0.5], cmap=cmap, alpha=0.5)

    #print(pt_num, area)
    pt_num+=1


ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
fig.tight_layout()
if variable==1:
    plt.savefig('frame.png')
if variable==2:
    plt.show()


'''
fig = plt.figure(figsize=(6,6))
totphi = [[0. for j in range(lx)] for i in range(ly)]
for k in range(lx*ly):
    yy=int(k/lx)
    xx=k-int(yy*lx)
    for i in range(pt_num):
        totphi[yy][xx] += sim_grid[k+i*lx*ly] * walls[k]
        for j in range(i+1, pt_num):
            totphi[yy][xx] += sim_grid[k+i*lx*ly] * sim_grid[k+j*lx*ly]


cmap = LinearSegmentedColormap.from_list('mycmap', ['grey', 'white'])
plt.imshow(totphi, interpolation='lanczos', cmap=cmap, origin='lower')

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
fig.tight_layout()
plt.savefig('frame.png')
plt.show()
'''
