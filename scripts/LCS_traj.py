import sys
import numpy as np
from math import *
import os.path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgba, to_hex, Normalize
import numpy as np
import scipy.ndimage
from itertools import combinations

from matplotlib import cm
import matplotlib
#matplotlib.use('Agg')


if len(sys.argv)!=7:
    print(sys.argv[0]," [topology file] [frames location] [start frame] [end frame] [1:save conf; 2:make plot] [1:forward; -1:backward]")
    sys.exit(1)

variable=int(float(sys.argv[5])) 
direction=sys.argv[6] 
start_frame=int(float(sys.argv[3])) 
end_frame=int(float(sys.argv[4])) 


tfile=open(sys.argv[1],"r")
N=int(tfile.readline().split()[0])
species=[]
line=tfile.readline()
lline=line.split()
for l in lline:
        species.append(int(l))
tfile.close()

numspecies=len(set(species))
time = int(end_frame - start_frame)

cfile=open(sys.argv[2],"r")
header=cfile.readline().split()
t=int(header[2])

header=cfile.readline().split()
lx=int(float(header[2]))
ly=int(float(header[3]))


start_line=(N+2)*int(float(sys.argv[3])) 
for i in range(start_line):
    cfile.readline()


tv = np.zeros(time)
tv[0] = t

xP = np.zeros((N , time))
yP = np.zeros((N , time))
pt_num = 0
time_line = 0

for line in cfile:
    words=line.split()

    if words[0]=='t':
        t=int(float(words[2]))
        time_line += 1
        if end_frame == start_frame + time_line:
            break
        tv[time_line] = t
        pt_num = 0
    elif words[0]=='b':
        lx=int(float(words[2]))
        ly=int(float(words[3]))
    else:
        CoMX=float(words[2])
        CoMY=float(words[3])

        xP[pt_num, time_line] = CoMX
        yP[pt_num, time_line] = CoMY

        pt_num += 1


# Number of particles
NP = xP.shape[0]
tf, t0 = tv[-1], tv[0]
T = (tf-t0) * 0.1
lam=1e-10

# Extract initial and final time snapshots
if direction == 'forward':
    xP0, yP0 = xP[:,0], yP[:,0]
    xPf, yPf = xP[:,-1], yP[:,-1]
elif direction == 'backward':
    xP0, yP0 = xP[:,-1], yP[:,-1]
    xPf, yPf = xP[:,0], yP[:,0]
else:
    print('Error: direction argument must be \'forward\' or \'backward\'')


eps = 24
# Calculate FTLE values at each particle
FTLE = np.zeros(NP)
for i in range(NP):
    # Compute initial distances
    dxP = xP0 - xP0[i]
    dyP = yP0 - yP0[i]

    for j in range(len(dxP)):
        if dxP[j] >= lx/2:
            dxP[j] -= lx
        if dxP[j] <= -lx/2:
            dxP[j] += lx

        if dyP[j] >= ly/2:
            dyP[j] -= ly
        if dyP[j] <= -ly/2:
            dyP[j] += ly
        
    # Find pairwise combinations of neighbor indices
    neighbors = np.flatnonzero((dxP**2+dyP**2)<eps**2)
    combs = list(combinations(range(len(neighbors)),2))
    ind1 = [comb[0] for comb in combs]
    ind2 = [comb[1] for comb in combs]
        
    # Form X and Y data matrices
    X = np.zeros((2,len(combs)))
    Y = np.zeros((2,len(combs)))
    X[0,:] = xP0[neighbors[ind1]]-xP0[neighbors[ind2]]
    X[1,:] = yP0[neighbors[ind1]]-yP0[neighbors[ind2]]
    Y[0,:] = xPf[neighbors[ind1]]-xPf[neighbors[ind2]]
    Y[1,:] = yPf[neighbors[ind1]]-yPf[neighbors[ind2]]

    for j in range(len(combs)):
        if X[0,j] >= lx/2:
            X[0,j] -= lx/2
        if X[0,j] <= -lx/2:
            X[0,j] += lx/2

        if X[1,j] >= ly/2:
            X[1,j] -= ly/2
        if X[1,j] <= -ly/2:
            X[1,j] += ly/2
        
        if Y[0,j] >= lx/2:
            Y[0,j] -= lx/2
        if Y[0,j] <= -lx/2:
            Y[0,j] += lx/2
            
        if Y[1,j] >= ly/2:
            Y[1,j] -= ly/2
        if Y[1,j] <= -ly/2:
            Y[1,j] += ly/2


    # Least square fit of flow map gradient
    A = Y@X.T + lam*max(1,len(neighbors))*np.eye(2)
    B = X@X.T + lam*max(1,len(neighbors))*np.eye(2)
    DF = A@np.linalg.inv(B)
        
    # Calculate FTLE as the largest singular value of DF
    FTLE[i] = np.log(np.linalg.norm(DF,2))/T


'''
FTLE_ridges = [[0 for q in range(LLx)] for k in range(LLy)]
FTLE_ridges_x = []
FTLE_ridges_y = []
lowest_lambda = 0
for i in range(0, LLx*LLy):
    y = int(i/LLx)
    x = i - y * LLx
    yy = (y + 1) % LLy
    yyy = (y - 1 + LLy) % LLy
    xx = (x + 1) % LLx
    xxx = (x - 1 + LLx) % LLx

    dsdx = 0.5 * (FTLE[y][xx] - FTLE[y][xxx])
    dsdy = 0.5 * (FTLE[yy][x] - FTLE[yyy][x])

    dsdxx = (FTLE[y][xx] + FTLE[y][xxx] - 2 * FTLE[y][x])
    dsdyy = (FTLE[yy][x] + FTLE[yyy][x] - 2 * FTLE[y][x])
    dsdxy = 0.25 * (FTLE[yy][xx] - FTLE[yyy][xx] - FTLE[yy][xxx] + FTLE[yyy][xxx])

    H11 = dsdxx
    H12 = dsdxy
    H22 = dsdyy

    lambda_1 = 0.5 * (H11 + H22) + sqrt((H11 - H22) * (H11 - H22) + 4 * H12 * H12)
    lambda_2 = 0.5 * (H11 + H22) - sqrt((H11 - H22) * (H11 - H22) + 4 * H12 * H12)
    
    vx = 0
    vy = 0
    lambda_f = 0
    if lambda_1 < lambda_2:
        normalization = 1/sqrt(H12 * H12 + (H11 - lambda_1) * (H11 - lambda_1))
        vx = - H12 * normalization
        vy = (H11 - lambda_1) * normalization
        lambda_f = lambda_1
    else:
        normalization = 1/sqrt(H12 * H12 + (H11 - lambda_2) * (H11 - lambda_2))
        vx = - H12 * normalization
        vy = (H11 - lambda_2) * normalization
        lambda_f = lambda_2


    if lambda_f < -0.003: #and abs(dsdx * vx + dsdy * vy) < 1e-7:
        if lambda_f < lowest_lambda:
            lowest_lambda = lambda_f
        FTLE_ridges[y][x] = lambda_f
        FTLE_ridges_x.append(x)
        FTLE_ridges_y.append(y)

#print(lowest_lambda)
#print(FTLE_ridges_x, FTLE_ridges_y)
'''


z_min, z_max = -np.abs(FTLE).max(), np.abs(FTLE).max()
'''
for i in range(0, NP):
    if FTLE[i] < 0.5 * z_max:
        FTLE[i]=0
'''

if variable==2:

    plt.figure(figsize=(6,6))
    x=np.arange(0,lx,1)
    y=np.arange(0,ly,1)
    X, Y = np.meshgrid(x, y)

    #z_min, z_max = -np.abs(FTLE).max(), np.abs(FTLE).max()
    #cset1 = plt.imshow(Z, cmap='hot', interpolation='nearest')
    #cset1 = plt.pcolormesh(X, Y, vorticity, cmap='RdBu', vmin=z_min, vmax=z_max)
    #cset1 = plt.imshow(FTLE, cmap='RdBu', interpolation='nearest', vmin=-z_max, vmax=z_max)
    cset1 = plt.scatter(xP[:,-1], yP[:,-1], s=20, c=FTLE, cmap='RdBu', vmin=-z_max, vmax=z_max)
    #cset1 = plt.plot(FTLE_ridges_x, FTLE_ridges_y, 'o', color='k', markersize=6)
    #cset1 = plt.imshow(vorticity, cmap='RdBu', interpolation='nearest', vmin=-1, vmax=1)

    #plt.ylabel(r'$C_v$', fontsize=18)
    #plt.xlabel('R', fontsize=18)
    #plt.xticks(fontsize=18)
    #plt.yticks(fontsize=18)
    #plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
    #plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/mean_velocity_1.png", transparent=True)
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim([0, lx])
    ax.set_ylim([0, ly])
    plt.show()
