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

avg_stress = []
variance_stress = []
avg_stress_var = []
void_area = []
time = []

for i in range(start_frame, end_frame, 1):

    #if i%3 != 0:
    #    continue

    if i<10:
        file = sys.argv[2] + "stress_field_00" + str(i) + ".txt"
    elif i<100:
        file = sys.argv[2] + "stress_field_0" + str(i) + ".txt"
    else:
        file = sys.argv[2] + "stress_field_" + str(i) + ".txt"
    #print(file)

    max_value = 0
    min_value = 0

    cfile=open(file,"r")
    header=cfile.readline().split()
    t=int(header[2])
    time.append(i)
    avg_stress.append(0.)
    variance_stress.append(0.)
    size_s = len(avg_stress)-1

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
        v_area=0
        for i in range(start_value,len(words),13):
            xx=float(words[i])
            yy=float(words[i+1])
            value_x=float(words[i+2])
            value_y=float(words[i+3])

            Z_x[int(yy)][int(xx)]=value_x
            Z_y[int(yy)][int(xx)]=value_y
            if value_x == 0 and value_y == 0:
                v_area+= 1 / (lx * ly)

            stress_iso = (value_x + value_y) / 2 
            values_iso_stress.append(stress_iso)
            avg_stress[size_s] += stress_iso / (lx * ly)
            if stress_iso > max_value:
                max_value = stress_iso
            if stress_iso < min_value:
                min_value = stress_iso

        for i in range(lx * ly):
            variance_stress[size_s] += (values_iso_stress[i] - avg_stress[size_s]) * (values_iso_stress[i] - avg_stress[size_s]) / (lx * ly - 1)
        avg_stress_var.append( (avg_stress[size_s] - min_value) / (max_value - min_value) )
        void_area.append(v_area)
        read_line += 1



if variable==1:
    for i in range(end_frame):
        print(time[i], avg_stress[i], variance_stress[i], avg_stress_var[i])

if variable==2:
    #plt.figure(figsize=(5.452423529,4.089317647))
    #fig1 = plt.figure(figsize=(5.452423529,4.089317647))
    fig1, ax1 = plt.subplots(figsize=(5.452423529,4.089317647))
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    #ax1 = fig1.add_subplot()
    ax1.plot(time, avg_stress, '--o', color='firebrick')
    ax2 = ax1.twinx()
    ax2.plot(time, void_area, '--o', color='forestgreen')
    ax1.set_xlabel('Time', fontsize=18)
    ax1.set_ylabel(r'$\langle \sigma_{iso} \rangle$', fontsize=18)
    ax2.set_ylabel('Voids area', fontsize=18)
    ax1.tick_params(axis='y', colors='firebrick')
    ax2.tick_params(axis='y', colors='forestgreen')
    ax1.yaxis.label.set_color('firebrick')
    ax2.yaxis.label.set_color('forestgreen')
    ax2.spines['right'].set_color('forestgreen')
    ax2.spines['left'].set_color('firebrick')
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    fig1.tight_layout()
    #plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)

    #plt.plot(time, avg_stress_var, '--o')
    #plt.plot(time, variance_stress, '--o')
    #plt.ylabel(r'$\sigma_{iso}$', fontsize=18)
    #plt.xlabel('Time', fontsize=18)
    #plt.xticks(fontsize=18)
    #plt.yticks(fontsize=18)
    #plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
    #plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/mean_velocity_1.png", transparent=True)
    plt.show()
