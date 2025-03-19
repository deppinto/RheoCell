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


if len(sys.argv)!=6:
    print(sys.argv[0]," [topology file] [frames location] [start frame] [end frame] [1:save conf; 2:make plot]")
    sys.exit(1)

variable=int(float(sys.argv[5])) 
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

total_size = 80
delta_R = 1
size_R = int(total_size / delta_R)
corr_dist = [i for i in range(size_R)]
correlations = [0. for i in range(size_R)]
corr_num = [0 for i in range(size_R)]
corr_dem = 0.


for i in range(start_frame, end_frame, 1):

    if i<10:
        file = sys.argv[2] + "nematic_field_00" + str(i) + ".txt"
    elif i<100:
        file = sys.argv[2] + "nematic_field_0" + str(i) + ".txt"
    else:
        file = sys.argv[2] + "nematic_field_" + str(i) + ".txt"

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

    for i in range(0, lx*ly):
        y1 = int(i/lx)
        x1 = i - y1 * lx
        corr_dem += 2 * (Z_x[y1][x1] * Z_x[y1][x1] + Z_y[y1][x1] * Z_y[y1][x1]) / (lx * ly)

        for j in range(i, lx*ly):
            y2 = int(j/lx)
            x2 = j - y2 * lx

            distx = x2 - x1
            if distx <= -lx/2:
                distx+=lx
            if distx >= lx/2:
                distx-=lx

            disty = y2 - y1
            if disty <= -ly/2:
                disty+=ly
            if disty >= ly/2:
                disty-=ly

            dist = sqrt(distx * distx + disty * disty)
            if int(dist)<lx/2:
                correlations[int(dist)] += 2 * (Z_x[y1][x1] * Z_x[y2][x2] + Z_y[y1][x1] * Z_y[y2][x2])
                corr_num[int(dist)] += 1

for i in range(size_R):
    if corr_num[i] > 0:
        correlations[i] = correlations[i] / corr_num[i]
        correlations[i] = correlations[i] / corr_dem
    if variable==1:
        print(corr_dist[i], correlations[i])



if variable==2:
    plt.figure(figsize=(5.452423529,4.089317647))
    plt.plot(corr_dist, correlations, '--o')
    plt.ylabel(r'$C_Q$', fontsize=18)
    plt.xlabel('R', fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
    #plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/mean_velocity_1.png", transparent=True)
    plt.show()
