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

avg_stress = []
max_stress = []
variance_stress = []
avg_stress_var = []
void_area = []
positive_iso = []
time = []


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

start_hole_time = -1
R0 = 8
area0 = pi * R0 * R0
thresh = 0.5 * pi * R0 * R0
radius_stress = 2 * R0 
cell_value_thresh = 0.01
voids_time = []
number_holes = 0
time_window = 10
hole_stats = 0.

save_stress_lattice = np.zeros((time_window, lx*ly))
save_max_stress_lattice = np.zeros(time_window)
avg_iso_stress_tau = np.zeros(time_window)
max_iso_stress_tau = np.zeros(time_window)
max_stress_all_time = 0.

for i in range(start_frame+1, end_frame, 1):

    voids=[[0. for q in range(lx)] for k in range(ly)]
    voids_area = lx*ly
    start_value = 11
    for j in range(N):
        words=traj_file.readline().split()
        for k in range(start_value,len(words),2):
            site=int(float(words[k]))
            value=float(words[k+1])
            yy=int(site/lx)
            xx=site-int(yy*lx)

            if voids[yy][xx] < cell_value_thresh and voids[yy][xx] + value > cell_value_thresh:
                voids_area -= 1

            voids[yy][xx] += value


    voids_time.append(voids_area)
    if voids_area > thresh and start_hole_time < 0:
        start_hole_time = t
        #if i - 1 > time_window:
        if i > 1:
            hole_stats += 1.
            max_save_tau = np.zeros(time_window)
            avg_save_tau = np.zeros(time_window)
            temp = np.zeros(4)
            for k in range(ly):
                for l in range(lx):
                    if voids[k][l] < cell_value_thresh:
                        temp[0] += cos(2*pi*l/lx); 
                        temp[1] += sin(2*pi*l/lx);
                        temp[2] += cos(2*pi*k/ly);
                        temp[3] += sin(2*pi*k/ly);
            avg_x = lx * ( atan2(-temp[1]/voids_area, -temp[0]/voids_area) + pi ) / (2*pi); 
            avg_y = ly * ( atan2(-temp[3]/voids_area, -temp[2]/voids_area) + pi ) / (2*pi);
            print(avg_x, avg_y)
            
            available_window = time_window
            if i - 1 < time_window:
                available_window = i - 1 

            for k in range(0, 2 * radius_stress):
                yy = int(avg_y - radius_stress + k)
                if yy < 0:
                    yy += ly
                if yy >= ly:
                    yy -= ly

                for l in range(0, 2 * radius_stress):
                    xx = int(avg_x - radius_stress + l)
                    if xx < 0:
                        xx += lx
                    if xx >= lx:
                        xx -= lx

                    distx = xx - avg_x
                    if distx <= -lx/2:
                        distx += lx
                    disty = yy - avg_y
                    if disty <= -ly/2:
                        disty += ly

                    if distx * distx + disty * disty < radius_stress * radius_stress:
                        for q in range(available_window):
                            #tt_window = int( (int((i - 1)%time_window) - q + time_window) ) % time_window
                            tt_window = int( i - 1 - q + time_window ) % time_window
                            #print(q, int( (int((i - 1)%time_window) - q + time_window)) % time_window, l, k, save_stress_lattice[tt_window][l + k * lx])

                            avg_save_tau[q] += save_stress_lattice[tt_window][xx + yy * lx] / (pi * radius_stress * radius_stress) #/ voids_area
                            if save_stress_lattice[tt_window][xx + yy * lx] > max_save_tau[q]:
                                max_save_tau[q] = save_stress_lattice[tt_window][xx + yy * lx]

            x1 = max_save_tau[0:available_window]
            x2 = avg_save_tau[0:available_window]
            if i - 1 < time_window:
                y1 = rescale_array(x1, time_window)
                max_iso_stress_tau[:] += y1[:] / max(save_max_stress_lattice)
                y2 = rescale_array(x2, time_window)
                avg_iso_stress_tau[:] += y2[:] / y1[:]
            else:
                max_iso_stress_tau[:] += x1[:] / max(save_max_stress_lattice)
                avg_iso_stress_tau[:] += x2[:] / x1[:]

            #max_iso_stress_tau[:] += max_save_tau[:] / max(save_max_stress_lattice)
            #avg_iso_stress_tau[:] += avg_save_tau[:] / max_save_tau[:]

    elif voids_area < thresh and start_hole_time >= 0:
        number_holes += 1
        start_hole_time = -1

    header=traj_file.readline().split()
    t=int(header[2])
    header=traj_file.readline().split()
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
    positive_iso_stress = 0
    count_positive_iso = 0
    start_value = 0
    read_line = 0
    for line in cfile:

        if read_line==1:
            read_line+=1
            continue

        words=line.split()
        Z_x=[[0 for q in range(lx)] for k in range(ly)]
        Z_y=[[0 for q in range(lx)] for k in range(ly)]
        for k in range(start_value,len(words),13):
            xx=float(words[k])
            yy=float(words[k+1])
            value_x=float(words[k+2])
            value_y=float(words[k+3])

            Z_x[int(yy)][int(xx)]=value_x
            Z_y[int(yy)][int(xx)]=value_y

            stress_iso = (value_x + value_y) / 2 
            if voids[int(yy)][int(xx)] < cell_value_thresh:
                stress_iso = 0.

            save_stress_lattice[int(i%time_window)][int(xx) + int(yy) * lx] = stress_iso
            if stress_iso>0:
                positive_iso_stress += stress_iso
                count_positive_iso += 1
            values_iso_stress.append(stress_iso)
            avg_stress[size_s] += stress_iso / (lx * ly)
            if stress_iso > max_value:
                max_value = stress_iso
            if stress_iso < min_value:
                min_value = stress_iso
            if stress_iso > max_stress_all_time:
                max_stress_all_time = stress_iso

        save_max_stress_lattice[int(i%time_window)] = max_value
        for i in range(lx * ly):
            variance_stress[size_s] += (values_iso_stress[i] - avg_stress[size_s]) * (values_iso_stress[i] - avg_stress[size_s]) / (lx * ly - 1)

        avg_stress_var.append( (avg_stress[size_s] - min_value) / (max_value - min_value) )
        void_area.append(voids_area)
        max_stress.append(max_value)
        positive_iso.append(positive_iso_stress/count_positive_iso)
        read_line += 1


if start_hole_time >= 0:
    number_holes += 1

#avg_iso_stress_tau[:] = avg_iso_stress_tau[:] / (max_stress_all_time * float(number_holes))
if hole_stats > 0:
    avg_iso_stress_tau[:] = avg_iso_stress_tau[:] / (hole_stats)
    max_iso_stress_tau[:] = max_iso_stress_tau[:] / (hole_stats)

if variable==1:
    #for i in range(len(time)):
        #print(time[i], avg_stress[i], variance_stress[i], avg_stress_var[i])
    #print("Final:", hole_stats, time_window, avg_iso_stress_tau, max_iso_stress_tau)

    with open('voids_stress_time.txt', 'w') as f:
        for i in range(len(avg_stress)):
            print(i, void_area[i], avg_stress[i], max_stress[i], positive_iso[i], file=f)

    with open('voids_stress_histogram_tau10.txt', 'w') as f:
        for i in range(len(avg_iso_stress_tau)):
            print(i, avg_iso_stress_tau[i], max_iso_stress_tau[i], file=f)


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
