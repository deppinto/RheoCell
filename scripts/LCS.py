import sys
import numpy as np
from math import *
import os.path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgba, to_hex, Normalize
import numpy as np
import scipy.ndimage

from matplotlib import cm
import matplotlib
#matplotlib.use('Agg')


if len(sys.argv)!=7:
    print(sys.argv[0]," [topology file] [frames location] [start frame] [end frame] [1:save conf; 2:make plot] [1:forward; -1:backwards]")
    sys.exit(1)

variable=int(float(sys.argv[5])) 
FTLE_variable=int(float(sys.argv[6])) 
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


file = sys.argv[2] + "velocity_field_000.txt"
lfile=open(file,"r")
header=lfile.readline().split()
t=int(header[2])
header=lfile.readline().split()
lx=int(float(header[2]))
ly=int(float(header[3]))

delta_x = 0.25
delta_y = 0.25
LLx = int(lx / delta_x)
LLy = int(ly / delta_y)
flow_map_x = [[q * delta_x for q in range(LLx)] for k in range(LLy)]
flow_map_y = [[k * delta_y for q in range(LLx)] for k in range(LLy)]
tstep = 500 * 0.1
time = (end_frame - start_frame) * tstep


step = FTLE_variable
for frame in range(start_frame, end_frame, 1):

    if frame<10:
        file = sys.argv[2] + "velocity_field_00" + str(frame) + ".txt"
    elif frame<100:
        file = sys.argv[2] + "velocity_field_0" + str(frame) + ".txt"
    else:
        file = sys.argv[2] + "velocity_field_" + str(frame) + ".txt"
    #print(file)

    cfile=open(file,"r")
    header=cfile.readline().split()
    t=int(header[2])

    header=cfile.readline().split()
    lx=int(float(header[2]))
    ly=int(float(header[3]))

    Z_x=[[0 for q in range(lx)] for k in range(ly)]
    Z_y=[[0 for q in range(lx)] for k in range(ly)]
    start_value=0
    read_line = 0
    for line in cfile:

        if read_line==1:
            read_line+=1
            continue

        words=line.split()
        Z_x=[[0 for q in range(lx)] for k in range(ly)]
        Z_y=[[0 for q in range(lx)] for k in range(ly)]
        for i in range(start_value,len(words),4):
            xx=float(words[i])
            yy=float(words[i+1])
            value_x=float(words[i+2])
            value_y=float(words[i+3])

            Z_x[int(yy)][int(xx)]=value_x
            Z_y[int(yy)][int(xx)]=value_y
        read_line += 1

    for i in range(0, LLx*LLy):
        y = int(i/LLx)
        x = i - y * LLx
        yy = int(flow_map_y[y][x])
        xx = int(flow_map_x[y][x])

        xxx = (xx + 1) % lx
        distx = (flow_map_x[y][x] - xx)
        if distx < -lx/2:
            distx+=lx
        distxx = (xxx - xx)
        if distxx < -lx/2:
            distxx+=lx

        yyy = (yy + 1) % ly
        disty = (flow_map_y[y][x] - yy)
        if disty < -ly/2:
            disty+=ly
        distyy = (yyy - yy)
        if distyy < -ly/2:
            distyy+=ly

        flow_map_x[y][x] += step * tstep * (Z_x[yy][xx] + (Z_x[yy][xxx] - Z_x[yy][xx]) * (flow_map_x[y][x] - xx) / distxx)
        flow_map_y[y][x] += step * tstep * (Z_y[yy][xx] + (Z_y[yyy][xx] - Z_y[yy][xx]) * (flow_map_y[y][x] - yy) / distyy)

        if flow_map_x[y][x] >= lx:
            flow_map_x[y][x] -= lx
        if flow_map_x[y][x] < 0:
            flow_map_x[y][x] += lx

        if flow_map_y[y][x] >= ly:
            flow_map_y[y][x] -= ly
        if flow_map_y[y][x] < 0:
            flow_map_y[y][x] += ly


FTLE = [[0 for q in range(LLx)] for k in range(LLy)]
for i in range(0, LLx*LLy):
    y = int(i/LLx)
    x = i - y * LLx
    yy = (y + 1) % LLy
    yyy = (y - 1 + LLy) % LLy
    xx = (x + 1) % LLx
    xxx = (x - 1 + LLx) % LLx

    dxfdx0 = (flow_map_x[y][xx] - flow_map_x[y][xxx])
    if dxfdx0 >= lx/2:
        dxfdx0 -= lx
    if dxfdx0 <= -lx/2:
        dxfdx0 += lx
    dxfdx0 = dxfdx0 / 2

    dxfdy0 = (flow_map_x[yy][x] - flow_map_x[yyy][x])
    if dxfdy0 >= lx/2:
        dxfdy0 -= lx
    if dxfdy0 <= -lx/2:
        dxfdy0 += lx
    dxfdy0 = dxfdy0 / 2

    dyfdx0 = (flow_map_y[y][xx] - flow_map_y[y][xxx])
    if dyfdx0 >= ly/2:
        dyfdx0 -= ly
    if dyfdx0 <= -ly/2:
        dyfdx0 += ly
    dyfdx0 = dyfdx0 / 2

    dyfdy0 = (flow_map_y[yy][x] - flow_map_y[yyy][x])
    if dyfdy0 >= ly/2:
        dyfdy0 -= ly
    if dyfdy0 <= -ly/2:
        dyfdy0 += ly
    dyfdy0 = dyfdy0 / 2

    C11 = dxfdx0 * dxfdx0 + dyfdx0 * dyfdx0
    C12 = dxfdx0 * dxfdy0 + dyfdx0 * dyfdy0
    C22 = dxfdy0 * dxfdy0 + dyfdy0 * dyfdy0

    #lambda_1 = 0.5 * (C11 + C22) - sqrt((C11 - C22) * (C11 - C22) + 4 * C12 * C12)
    lambda_2 = 0.5 * (C11 + C22) + sqrt((C11 - C22) * (C11 - C22) + 4 * C12 * C12)

    FTLE[y][x] = log(lambda_2) / (2 * abs(time))


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

z_min, z_max = -np.abs(FTLE).max(), np.abs(FTLE).max()
'''
for i in range(0, lx*ly):
    y = int(i/lx)
    x = i - y * lx
    if FTLE[y][x] < z_max - 0.6 * z_max:
        FTLE[y][x]=0
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
    array1 = np.array(flow_map_x)
    array2 = np.array(flow_map_y)
    cset1 = plt.scatter(array1.flatten(), array2.flatten(), s=6, c=FTLE, cmap='RdBu', vmin=-z_max, vmax=z_max)
    cset1 = plt.plot(FTLE_ridges_x, FTLE_ridges_y, 'o', color='k', markersize=6)
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
