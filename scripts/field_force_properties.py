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

tracer_particle = [75, 33]
time = []
fig = plt.figure(figsize=(6,6))

for i in range(start_frame, end_frame, 1):

    if i<10:
        file = sys.argv[2] + "force_field_00" + str(i) + ".txt"
    elif i<100:
        file = sys.argv[2] + "force_field_0" + str(i) + ".txt"
    else:
        file = sys.argv[2] + "force_field_" + str(i) + ".txt"

    cfile=open(file,"r")
    header=cfile.readline().split()
    t=int(header[2])
    time.append(i)

    header=cfile.readline().split()
    lx=int(float(header[2]))
    ly=int(float(header[3]))

    Z_x=[[0 for q in range(lx)] for k in range(ly)]
    Z_y=[[0 for q in range(lx)] for k in range(ly)]
    values_iso_stress = []
    start_value=0
    read_line = 0
    for line in cfile:

        if read_line==1:
            read_line+=1
            continue

        words=line.split()
        Z_x=[[0 for q in range(lx)] for k in range(ly)]
        Z_y=[[0 for q in range(lx)] for k in range(ly)]
        for i in range(start_value,len(words),8):
            xx=float(words[i])
            yy=float(words[i+1])
            value_x=float(words[i+2])
            value_y=float(words[i+3])

            Z_x[int(yy)][int(xx)]=value_x
            Z_y[int(yy)][int(xx)]=value_y

        tracer_particle[0] += Z_x[int(tracer_particle[1])][int(tracer_particle[0])] * 10 / 0.5 
        tracer_particle[1] += Z_y[int(tracer_particle[1])][int(tracer_particle[0])] * 10 / 0.5
        cset1 = plt.plot(tracer_particle[0], tracer_particle[1], 'go', markersize=10)

        read_line += 1



if variable==1:
    for i in range(end_frame):
        print(time[i])

if variable==2:
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    #ax.set_xlim([0, lx])
    #ax.set_ylim([0, ly])
    plt.show()
