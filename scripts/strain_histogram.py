import sys
import numpy as np
from math import *
import os.path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgba, to_hex, Normalize
import numpy as np
import scipy.ndimage
from scipy.interpolate import interp1d
from itertools import combinations

from matplotlib import cm
import matplotlib
matplotlib.use('Agg')


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


start_hole_time = -1
R0 = 8
area0 = pi * R0 * R0
thresh = 0.5 * pi * R0 * R0
radius_stress = 2 * R0 
radius_stress_2 = radius_stress * radius_stress
cell_value_thresh = 0.01
number_holes = 0

total_frames = end_frame - start_frame

strain_bins = 20
strain_min = 0.
strain_max = 0.005
strain_delta = strain_max - strain_min
delta_bin = strain_delta / strain_bins
x_st = [(i * delta_bin + strain_min + delta_bin/2) for i in range(strain_bins)]
strain_hist = np.zeros(strain_bins)
strain_count = np.zeros(strain_bins)

strain_bins_30 = 30
delta_bin_30 = strain_delta / strain_bins_30
x_st_30 = [(i * delta_bin_30 + strain_min + delta_bin_30/2) for i in range(strain_bins_30)]
strain_hist_30 = np.zeros(strain_bins_30)
strain_count_30 = np.zeros(strain_bins_30)

strain_bins_40 = 40
delta_bin_40 = strain_delta / strain_bins_40
x_st_40 = [(i * delta_bin_40 + strain_min + delta_bin_40/2) for i in range(strain_bins_40)]
strain_hist_40 = np.zeros(strain_bins_40)
strain_count_40 = np.zeros(strain_bins_40)


strain=[[0. for q in range(lx)] for k in range(ly)]
voids=[[0. for q in range(lx)] for k in range(ly)]
for i in range(start_frame+1, end_frame+1, 1):

    voids_prev = voids
    voids=[[0. for q in range(lx)] for k in range(ly)]
    voids_area = lx*ly
    start_value = 11
    pt_num = 0
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

        pt_num += 1


    if i < end_frame:
        header=traj_file.readline().split()
        if len(header)<3:
            break
        t=int(header[2])
        header=traj_file.readline().split()

    if voids_area > thresh and start_hole_time < 0:
        start_hole_time = t
        for y1 in range(0, ly):
            for x1 in range(0, lx):
                if voids[y1][x1] < cell_value_thresh:
                    if voids_prev[y1][x1] > cell_value_thresh and i > start_frame + 1:
                        if strain[y1][x1] >= strain_min and strain[y1][x1] < strain_max:
                            index = int( (strain[y1][x1] - strain_min)/delta_bin )
                            strain_hist[index] += 1.
                            index = int( (strain[y1][x1] - strain_min)/delta_bin_30 )
                            strain_hist_30[index] += 1.
                            index = int( (strain[y1][x1] - strain_min)/delta_bin_40 )
                            strain_hist_40[index] += 1.

    elif voids_area < 1 and start_hole_time >= 0:
        number_holes += 1
        start_hole_time = -1

    if start_hole_time >= 0:
        continue

    if i<10:
        file = sys.argv[2] + "velocity_field_00" + str(i) + ".txt"
    elif i<100:
        file = sys.argv[2] + "velocity_field_0" + str(i) + ".txt"
    else:
        file = sys.argv[2] + "velocity_field_" + str(i) + ".txt"
    #print(file)

    cfile=open(file,"r")
    header=cfile.readline().split()
    t=int(header[2])
    if len(header)<3:
        break

    header=cfile.readline().split()
    lx=int(float(header[2]))
    ly=int(float(header[3]))

    Z_x=[[0 for q in range(lx)] for k in range(ly)]
    Z_y=[[0 for q in range(lx)] for k in range(ly)]
    start_value = 0
    read_line = 0
    for line in cfile:

        if read_line==0:
            read_line+=1
            continue

        words=line.split()
        Z_x=[[0 for q in range(lx)] for k in range(ly)]
        Z_y=[[0 for q in range(lx)] for k in range(ly)]
        for k in range(start_value,len(words), 4):
            xx=float(words[k])
            yy=float(words[k+1])

            value_x=float(words[k+2])
            value_y=float(words[k+3])
            Z_x[int(yy)][int(xx)]=value_x
            Z_y[int(yy)][int(xx)]=value_y

        read_line += 1


    #strain=[[0. for q in range(lx)] for k in range(ly)]
    for y1 in range(0, ly):
        for x1 in range(0, lx):

            if voids[y1][x1] < cell_value_thresh:
                if voids_prev[y1][x1] > cell_value_thresh and i > start_frame + 1:
                    #print(strain[y1][x1], y1, x1, voids_area, i)
                    if strain[y1][x1] >= strain_min and strain[y1][x1] < strain_max:
                        index = int( (strain[y1][x1] - strain_min)/delta_bin )
                        strain_hist[index] += 1.
                        index = int( (strain[y1][x1] - strain_min)/delta_bin_30 )
                        strain_hist_30[index] += 1.
                        index = int( (strain[y1][x1] - strain_min)/delta_bin_40 )
                        strain_hist_40[index] += 1.
                continue

            xnext = (x1 + 1) % lx
            xprev = (x1 - 1 + lx) % lx 
            ynext = (y1 + 1) % ly 
            yprev = (y1 - 1 + ly) % ly

            dvxdx = (Z_x[y1][xnext] - Z_x[y1][xprev])/2
            dvydy = (Z_y[ynext][x1] - Z_y[yprev][x1])/2
            strain_value = 0.5 * (dvxdx + dvydy)
            strain[y1][x1] = strain_value

            if strain_value >= strain_min and strain_value < strain_max:
                index = int( (strain_value - strain_min)/delta_bin )
                strain_count[index] += 1.
                index = int( (strain_value - strain_min)/delta_bin_30 )
                strain_count_30[index] += 1.
                index = int( (strain_value - strain_min)/delta_bin_40 )
                strain_count_40[index] += 1.

            if strain_value > strain_max:
                print("strain: ", value)

if start_hole_time >= 0:
    number_holes += 1

if variable==1:
    #for i in range(len(time)):
        #print(time[i], avg_stress[i], variance_stress[i], avg_stress_var[i])
    #print("final:", strain_hist, strain_count, number_holes)

    with open('voids_strain_histogram_20bins.txt', 'w') as f:
        for i in range(len(strain_hist)):
            print(x_st[i], strain_hist[i], strain_count[i], file=f)

    with open('voids_strain_histogram_30bins.txt', 'w') as f:
        for i in range(len(strain_hist_30)):
            print(x_st_30[i], strain_hist_30[i], strain_count_30[i], file=f)

    with open('voids_strain_histogram_40bins.txt', 'w') as f:
        for i in range(len(strain_hist_40)):
            print(x_st_40[i], strain_hist_40[i], strain_count_40[i], file=f)


if variable==2:
    for i in range(strain_bins):
        if strain_count[i] > 0:
            strain_hist[i] = (strain_hist[i] / strain_count[i])

    for i in range(strain_bins_30):
        if strain_count_30[i] > 0:
            strain_hist_30[i] = (strain_hist_30[i] / strain_count_30[i])

    for i in range(strain_bins_40):
        if strain_count_40[i] > 0:
            strain_hist_40[i] = (strain_hist_40[i] / strain_count_40[i])

    x_st = np.array(x_st)
    x = x_st[0:]
    y = strain_hist[0:]

    plt.figure(figsize=(5.452423529,4.089317647))
    plt.plot(x, y, '--o')
    #plt.plot(x_st_30, strain_hist_30 , '--o')
    #plt.plot(x_st_40, strain_hist_40 , '--o')
    #plt.yscale('log')
    #plt.xscale('log')
    plt.xlabel("Strain")
    plt.ylabel(r"P($\varepsilon$)")

    plt.show()
