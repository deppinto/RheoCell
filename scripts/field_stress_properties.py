import sys
import numpy as np
from math import *
import os.path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgba, to_hex, Normalize
import numpy as np
import scipy.ndimage
from scipy.interpolate import interp1d

from matplotlib import cm
import matplotlib
#matplotlib.use('Agg')


if len(sys.argv)!=7:
    print(sys.argv[0]," [topology file] [frames location] [trajectory file] [start frame] [end frame] [1:save conf; 2:make plot]")
    sys.exit(1)

variable=int(float(sys.argv[6])) 
start_frame=int(float(sys.argv[4])) 
end_frame=int(float(sys.argv[5])) 


def rescale_array(x, M, method='linear'):
    """
    Rescales an array x of size N into an array y of size M, maintaining the same scale.
    
    Parameters:
        x (numpy array): Input array of size N.
        M (int): Desired output size (M > N).
        method (str): Interpolation method ('linear', 'cubic', 'nearest', etc.).
    
    Returns:
        numpy array: Interpolated array of size M.
    """
    N = len(x)
    original_indices = np.linspace(0, 1, N)  # Normalized indices of x
    new_indices = np.linspace(0, 1, M)  # New indices for interpolation
    interpolator = interp1d(original_indices, x, kind=method, fill_value="extrapolate")
    y = interpolator(new_indices)
    return y


tfile=open(sys.argv[1],"r")
N=int(tfile.readline().split()[0])
species=[]
line=tfile.readline()
lline=line.split()
for l in lline:
        species.append(int(l))
tfile.close()

numspecies=len(set(species))


file_t = sys.argv[3]
traj_file=open(file_t,"r")
header=traj_file.readline().split()
t=int(header[2])
header=traj_file.readline().split()
lx=int(float(header[2]))
ly=int(float(header[3]))
start_line=(N+2)*start_frame
for i in range(start_line):
    traj_file.readline()


R0 = 8
area0 = pi * R0 * R0
thresh = 0.5 * pi * R0 * R0

avg_stress_xy_time = [0. for i in range(end_frame-start_frame)]
avg_stress_xx_time = [0. for i in range(end_frame-start_frame)]
avg_stress_yy_time = [0. for i in range(end_frame-start_frame)]
avg_stress_xy = 0.
avg_stress_xx = 0.
avg_stress_yy = 0.
time = []


for i in range(start_frame, end_frame, 1):

    '''
    start_value = 11
    for j in range(N):
        words=traj_file.readline().split()
        for k in range(start_value,len(words),2):
            site=int(float(words[k]))
            value=float(words[k+1])
            yy=int(site/lx)
            xx=site-int(yy*lx)

            voids[yy][xx] += value
    '''


    #header=traj_file.readline().split()
    #t=int(header[2])
    #header=traj_file.readline().split()

    if i<10:
        file = sys.argv[2] + "stress_field_00" + str(i) + ".txt"
    elif i<100:
        file = sys.argv[2] + "stress_field_0" + str(i) + ".txt"
    else:
        file = sys.argv[2] + "stress_field_" + str(i) + ".txt"
    #print(file)

    cfile=open(file,"r")
    header=cfile.readline().split()
    t=int(header[2])
    time.append(i)

    header=cfile.readline().split()
    lx=int(float(header[2]))
    ly=int(float(header[3]))

    Z_xx=[[0 for q in range(lx)] for k in range(ly)]
    Z_yy=[[0 for q in range(lx)] for k in range(ly)]
    Z_xy=[[0 for q in range(lx)] for k in range(ly)]
    start_value = 0
    read_line = 0
    for line in cfile:

        if read_line==1:
            read_line+=1
            continue

        words=line.split()
        Z_xx=[[0 for q in range(lx)] for k in range(ly)]
        Z_yy=[[0 for q in range(lx)] for k in range(ly)]
        Z_xy=[[0 for q in range(lx)] for k in range(ly)]
        for k in range(start_value,len(words),5):
            xx=float(words[k])
            yy=float(words[k+1])
            value_xx=float(words[k+2])
            value_yy=float(words[k+3])
            value_xy=float(words[k+4])

            Z_xx[int(yy)][int(xx)]=value_xx
            Z_yy[int(yy)][int(xx)]=value_yy
            Z_xy[int(yy)][int(xx)]=value_xy

            if yy>50 and yy<ly-50:
                avg_stree_xy_time[i] += value_xy 
                avg_stree_xx_time[i] += value_xx 
                avg_stree_yy_time[i] += value_yy 

                avg_stress_xy += value_xy
                avg_stress_xx += value_xx
                avg_stress_yy += value_yy


        read_line += 1


if variable==1:
    #print("Stress:", avg_stress_xx, avg_stress_yy, avg_stress_xy)

    with open('stress_time_avg.txt', 'w') as f:
        print(avg_stress_xx, avg_stress_yy, avg_stress_xy)

    with open('stress_time.txt', 'w') as f:
        for i in range(len(avg_stress_xy_time)):
            print(time[i], avg_stress_xx[i], avg_stress_yy[i], avg_stress_xy[i], file=f)


if variable==2:
    #plt.figure(figsize=(5.452423529,4.089317647))
    #fig1 = plt.figure(figsize=(5.452423529,4.089317647))
    fig1, ax1 = plt.subplots(figsize=(5.452423529,4.089317647))
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    x = np.array(avg_stress)

    #ax1 = fig1.add_subplot()
    #ax1.plot(time, avg_stress, '--o', color='firebrick')
    #ax1.plot(time, max_stress, '--o', color='firebrick')
    #ax1.plot(time, variance_stress, '--o', color='firebrick')
    #ax1.plot(time, positive_iso, '--o', color='firebrick')
    ax1.plot(time, x, '--o', color='firebrick')
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
    

    #x = np.array(positive_iso)
    '''
    y = np.array(void_area)
    corr = np.correlate(x - np.mean(x), y - np.mean(y), mode="full") / (np.std(x) * np.std(y) * len(x))
    lags = np.arange(-len(x) + 1, len(x))
    fig2 = plt.subplots(figsize=(5.452423529,4.089317647))
    plt.plot(lags, corr)
    plt.xlabel("Lag")
    plt.ylabel("Cross-Correlation")
    plt.title("Cross-Correlation Function")
    '''

    plt.show()
