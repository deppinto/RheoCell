import sys
import numpy as np
from math import *
import os.path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgba, to_hex, Normalize
import numpy as np
import scipy.ndimage

from matplotlib import cm

if len(sys.argv)!=5:
    print(sys.argv[0]," [topology file] [conf file] [1:save conf; 2:make plot] [choose field: 0 - all; 1 - passive; 2 - active]")
    sys.exit(1)

variable=int(float(sys.argv[3])) 
f_variable=int(float(sys.argv[4])) 

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
Z_field=[[0 for q in range(lx)] for k in range(ly)]
Z_pass=[[0 for q in range(lx)] for k in range(ly)]
Z_act=[[0 for q in range(lx)] for k in range(ly)]
fig = plt.figure(figsize=(6,6))
start_value=0
for line in cfile:
    words=line.split()
    Z_field=[[0 for q in range(lx)] for k in range(ly)]
    Z_pass=[[0 for q in range(lx)] for k in range(ly)]
    Z_act=[[0 for q in range(lx)] for k in range(ly)]
    for i in range(start_value,len(words),8):
        #site=int(float(words[i]))
        xx=float(words[i])
        yy=float(words[i+1])

        value_field_x=float(words[i+2])
        value_field_y=float(words[i+3])
        value_pass_x=float(words[i+4])
        value_pass_y=float(words[i+5])
        value_act_x=float(words[i+6])
        value_act_y=float(words[i+7])

        Z_field[int(yy)][int(xx)]=sqrt(value_field_x * value_field_x + value_field_y * value_field_y)
        Z_pass[int(yy)][int(xx)]=sqrt(value_pass_x * value_pass_x + value_pass_y * value_pass_y)
        Z_act[int(yy)][int(xx)]=sqrt(value_act_x * value_act_x + value_act_y * value_act_y)

        #if abs(value_field_x) > 0 or abs(value_field_y) > 0 :
            #cset1 = plt.arrow(xx, yy, value_field_x , value_field_y , width=0.5, color='k')
        #else:
            #cset1 = plt.arrow(xx, yy, value_field_x , value_field_y , width=0.5, color='w')

    #z_min, z_max = -np.abs(Z).max(), np.abs(Z).max()
    z_min, z_max = 0., np.abs(Z_field).max()
    X, Y = np.meshgrid(x, y)
    cset1 = plt.imshow(Z_field, cmap='hot', interpolation='nearest')
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
    plt.clf()


    #z_min, z_max = -np.abs(Z).max(), np.abs(Z).max()
    z_min, z_max = 0., np.abs(Z_pass).max()
    X, Y = np.meshgrid(x, y)
    cset1 = plt.imshow(Z_pass, cmap='hot', interpolation='nearest')
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
    plt.clf()


    #z_min, z_max = -np.abs(Z).max(), np.abs(Z).max()
    z_min, z_max = 0., np.abs(Z_act).max()
    X, Y = np.meshgrid(x, y)
    cset1 = plt.imshow(Z_act, cmap='hot', interpolation='nearest')
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
    plt.clf()
