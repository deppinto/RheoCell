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

def rotate(n,p):
    '''
    takes as arguments a vector n and an integer p
    rotates v by 2pi/p and returns the result
    '''

    t = 2*np.pi/p
    nx = cos(t)*n[0] - sin(t)*n[1]
    ny = sin(t)*n[0] + sin(t)*n[1]
    return [nx,ny]


def wang(a, b):
    """Infamous chinese function"""
    p = 1.
    ang = atan2(abs(a[0]*b[1]-a[1]*b[0]), a[0]*b[0]+a[1]*b[1])

    #if(ang > pi/2.):
        #b = [-i for i in b]

    while(abs(ang) > np.pi/p + 1e-3):
        b = rotate(b,p)
        ang = atan2(abs(a[0]*b[1]-a[1]*b[0]), a[0]*b[0]+a[1]*b[1])

    m = a[0]*b[1]-a[1]*b[0]
    return -np.sign(m)*atan2(abs(m), a[0]*b[0]+a[1]*b[1])


def collapse(i, j, s, LX, LY, w, x=0, y=0, n=0,rng = [0.4,0.6]):

    if (s*w[i][j] > rng[0]) and (s*w[i][j] < rng[1]):
        w[i][j] = 0
        x1,y1,n1 = collapse((i+1) % LY, j, s, LX, LY, w, x, y, n,rng)
        x2,y2,n2 = collapse((i-1+LY) % LY, j, s, LX, LY, w, x, y, n, rng)
        x3,y3,n3 = collapse(i, (j+1) % LX, s, LX, LY, w, x, y, n, rng)
        x4,y4,n4 = collapse(i, (j-1+LX) % LX, s, LX, LY, w, x, y, n, rng)
        x = j + x1 + x2 +x3 +x4
        y = i + y1 + y2 +y3 +y4
        n = 1 + n1 + n2 +n3 +n4
        return x,y,n
    else:
        return 0,0,0

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
radius_stress_2 = radius_stress * radius_stress
cell_value_thresh = 0.01
voids_time = []
number_holes = 0
time_window = 10
hole_stats = 0.

maximum_defect_number = 40
total_frames = end_frame - start_frame
number_defects_plusone = np.zeros(total_frames)
defect_plusone_x = np.zeros(total_frames * maximum_defect_number)
defect_plusone_y = np.zeros(total_frames * maximum_defect_number)
avg_distance = np.zeros(time_window)

xP = np.zeros((N , total_frames))
yP = np.zeros((N , total_frames))
hole_timestamp = []
hole_position_x = []
hole_position_y = []

for i in range(start_frame+1, end_frame+1, 1):

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

            CoMX=float(words[2])
            CoMY=float(words[3])

            xP[pt_num, i - 1] = CoMX
            yP[pt_num, i - 1] = CoMY

            if voids[yy][xx] < cell_value_thresh and voids[yy][xx] + value > cell_value_thresh:
                voids_area -= 1
            voids[yy][xx] += value

        pt_num += 1


    voids_time.append(voids_area)
    if voids_area > thresh and start_hole_time < 0:
        start_hole_time = t
        #if i - 1 > time_window:
        if i > 1:
            hole_stats += 1.
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

            hole_timestamp.append(i - 1)
            hole_position_x.append(avg_x)
            hole_position_y.append(avg_y)
            #print(avg_x, avg_y)
            
            available_window = time_window
            if i - 2 < time_window:
                available_window = i - 2 

            avg_distance_tau = np.zeros(time_window)
            for q in range(available_window):
                tt_window = int(i - 2 - q)
                min_value = lx * ly
                for k in range(int(number_defects_plusone[tt_window])):
                    distx = defect_plusone_x[k + tt_window * maximum_defect_number] - avg_x
                    if distx <= -lx/2:
                        distx += lx
                    disty = defect_plusone_y[k + tt_window * maximum_defect_number] - avg_y
                    if disty <= -ly/2:
                        disty += ly

                    dist_val = sqrt(distx * distx + disty * disty)
                    if min_value > dist_val:
                        min_value = dist_val

                avg_distance_tau[q] = min_value

            small_array = avg_distance_tau[0:available_window]
            if available_window == 1:
                small_array = avg_distance_tau[0:available_window+1]
                small_array[1] = small_array[0]
                #print(small_array)

            if i - 2 < time_window:
                new_array = rescale_array(small_array, time_window)
                avg_distance[:] += new_array[:]
            else:
                avg_distance[:] += small_array[:]


    elif voids_area < thresh and start_hole_time >= 0:
        number_holes += 1
        start_hole_time = -1

    if i < end_frame:
        header=traj_file.readline().split()
        t=int(header[2])
        header=traj_file.readline().split()

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
    time.append(i)

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


    LLX = lx
    LLY = ly
    vecfield_nx = [[0. for q in range(0, LLX)] for p in range(0, LLY)]
    vecfield_ny = [[0. for q in range(0, LLX)] for p in range(0, LLY)]
    vecfield_nx = Z_x
    vecfield_ny = Z_y
    winding_number = [[0. for q in range(0, LLX)] for p in range(0, LLY)]
    strain=[[0. for q in range(lx)] for k in range(ly)]
    for p in range(0, LLY):
        for q in range(0, LLX):
            y1 = p
            x1 = q

            xnext = (x1 + 1)
            xprev = (x1 - 1) 
            ynext = (y1 + 1) 
            yprev = (y1 - 1)

            if xnext>=lx:
                xnext-=lx
            if xprev<0:
                xprev+=lx
            if ynext>=lx:
                ynext-=ly
            if yprev<0:
                yprev+=ly

            dvxdx = (vecfield_nx[y1][xnext] - vecfield_nx[y1][xprev])/2
            dvydy = (vecfield_ny[ynext][x1] - vecfield_ny[yprev][x1])/2
            strain[y1][x1] = 0.5 * (dvxdx + dvydy)

            ax1 = [vecfield_nx[p][(q+1) % LLX], vecfield_ny[p][(q+1) % LLX]]
            ax2 = [vecfield_nx[p][(q-1+LLX) % LLX], vecfield_ny[p][(q-1+LLX) % LLX]]
            ax3 = [vecfield_nx[(p+1) % LLY][q], vecfield_ny[(p+1) % LLY][q]]
            ax4 = [vecfield_nx[(p-1+LLY) % LLY][q], vecfield_ny[(p-1+LLY) % LLY][q]]
            ax5 = [vecfield_nx[(p-1+LLY) % LLY][(q+1) % LLX], vecfield_ny[(p-1+LLY) % LLY][(q+1) % LLX]]
            ax6 = [vecfield_nx[(p-1+LLY) % LLY][(q-1+LLX) % LLX], vecfield_ny[(p-1+LLY)%LLY][(q-1+LLX)%LLX]]
            ax7 = [vecfield_nx[(p+1) % LLY][(q+1) % LLX], vecfield_ny[(p+1) % LLY][(q+1) % LLX]]
            ax8 = [vecfield_nx[(p+1) % LLY][(q-1+LLX) % LLX], vecfield_ny[(p+1) % LLY][(q-1+LLX) % LLX]]

            winding_number[p][q] = wang(ax1, ax5)
            winding_number[p][q] += wang(ax5, ax4)
            winding_number[p][q] += wang(ax4, ax6)
            winding_number[p][q] += wang(ax6, ax2)
            winding_number[p][q] += wang(ax2, ax8)
            winding_number[p][q] += wang(ax8, ax3)
            winding_number[p][q] += wang(ax3, ax7)
            winding_number[p][q] += wang(ax7, ax1)
            winding_number[p][q] /= 2. * pi

    defect_thresh = 0.05
    for p in range(0,LLY):
        for q in range(0,LLX):
            # keep this just in case our other symmetries give us integer defects
            if (abs(winding_number[p][q]) > 1-defect_thresh) and (abs(winding_number[p][q]) < 1+defect_thresh):
                # charge sign
                s = np.sign(winding_number[p][q])
                # bfs
                sum_x, sum_y, n = collapse(p, q, s, LLX, LLY, winding_number, rng = [1-defect_thresh, 1+defect_thresh])
                x,y = sum_x/n,sum_y/n
                # add defect to list
                if s==1:
                    if number_defects_plusone[i - 1] >= maximum_defect_number:
                        print("Too many defects", i)
                    defect_plusone_x[int(number_defects_plusone[i - 1]) + (i-1) * maximum_defect_number] = x
                    defect_plusone_y[int(number_defects_plusone[i - 1]) + (i-1) * maximum_defect_number] = y
                    number_defects_plusone[i - 1] += 1

                    if start_hole_time<0:
                        avg_x = x
                        avg_y = y
                        avg_strain = 0.
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
                                    avg_strain += strain[yy][xx] / (pi * radius_stress * radius_stress)
                        print(avg_strain, voids_area, i)


if hole_stats > 0:
    avg_distance[:] = avg_distance[:] / hole_stats


max_FTLE = 0.
avg_FTLE_area = 0.
FTLE_area_counts = 0.
LCS_window = 5
radius_stress_2 = radius_stress * (3/2) * radius_stress * (3/2)

# Number of particles
NP = xP.shape[0]
lam=1e-10

FTLE_save = []
time_save = []
idx_save = []
for k in range(int(hole_stats)):

    tf = 0
    t0 = 0
    if hole_timestamp[k] + LCS_window < end_frame and hole_timestamp[k] - LCS_window >= 0:
        tf = hole_timestamp[k] + LCS_window
        t0 = hole_timestamp[k] - LCS_window
    elif hole_timestamp[k] + LCS_window >= end_frame and hole_timestamp[k] - LCS_window < 0:
        print("LCS problem: time window too small")
        exit (1)
    elif hole_timestamp[k] + LCS_window >= end_frame and hole_timestamp[k] - LCS_window >= 0:
        tf = end_frame 
        t0 = end_frame - 2 * LCS_window 
    elif hole_timestamp[k] + LCS_window < end_frame and hole_timestamp[k] - LCS_window < 0:
        t0 = start_frame
        tf = 2 * LCS_window
    else:
        print("Forgot something?")

    T = (tf-t0) * 0.1 * 1000
    xP0, yP0 = xP[:,t0], yP[:,t0]
    xPf, yPf = xP[:,tf-1], yP[:,tf-1]

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

        distx = xP0[i] - hole_position_x[k] 
        disty = yP0[i] - hole_position_y[k] 
        if distx <= -lx/2:
            distx += lx
        if disty <= -ly/2:
            disty += ly

        if distx * distx + disty * disty < radius_stress_2:
            avg_FTLE_area += FTLE[i]
            FTLE_area_counts += 1.

            FTLE_save.append(FTLE[i])
            time_save.append(hole_timestamp[k])
            idx_save.append(i)
    max_FTLE += max(FTLE)




if variable==1:
    #for i in range(len(time)):
        #print(time[i], avg_stress[i], variance_stress[i], avg_stress_var[i])
    #print("final:", hole_stats, time_window, avg_distance, avg_FTLE_area/FTLE_area_counts, max_FTLE/hole_stats)

    with open('voids_velocity_stats.txt', 'w') as f:
        if hole_stats>0:
            print(hole_stats, time_window, avg_FTLE_area/FTLE_area_counts, max_FTLE/hole_stats, file=f)
        else:
            print(hole_stats, time_window, 0, 0, file=f)

    with open('voids_velocity_histogram_tau10.txt', 'w') as f:
        for i in range(len(avg_distance)):
            print(i, avg_distance[i], file=f)

    with open('voids_velocity_LCS_distribution_window5.txt', 'w') as f:
        for i in range(len(FTLE_save)):
            print(time_save[i], idx_save[i], FTLE_save[i], file=f)



if variable==2:

    print("final:", hole_stats, time_window, avg_distance, avg_FTLE_area/FTLE_area_counts, max_FTLE/hole_stats)
    '''
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


    x = np.array(positive_iso)
    y = np.array(void_area)
    corr = np.correlate(x - np.mean(x), y - np.mean(y), mode="full") / (np.std(x) * np.std(y) * len(x))
    lags = np.arange(-len(x) + 1, len(x))
    fig2 = plt.subplots(figsize=(5.452423529,4.089317647))
    plt.plot(lags, corr)
    plt.xlabel("Lag")
    plt.ylabel("Cross-Correlation")
    plt.title("Cross-Correlation Function")

    plt.show()
    '''
