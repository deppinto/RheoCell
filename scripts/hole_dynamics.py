import sys
import numpy as np
from math import *
import os.path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgba, to_hex, Normalize
import numpy as np
import scipy.ndimage
from skimage import measure
from scipy.interpolate import interp1d


from matplotlib import cm
import matplotlib
matplotlib.use('Agg')


if len(sys.argv)!=6:
    print(sys.argv[0]," [topology file] [frames location] [start frame] [end frame] [1:save conf; 2:make plot]")
    sys.exit(1)

variable=int(float(sys.argv[5])) 
start_frame=int(float(sys.argv[3])) 
end_frame=int(float(sys.argv[4])) 

def apply_periodic_boundary_conditions(binary_array):
    """Extend binary array to handle periodic boundaries."""
    extended_array = np.tile(binary_array, (3, 3))  # Create a 3x3 tiled version
    return extended_array


def compute_fourier_descriptors(contour, n_harmonics=10):
    """Compute Fourier descriptors of a 2D shape boundary."""
    # Center the contour
    contour_centered = contour - np.mean(contour, axis=0)

    # Convert contour to complex representation
    z = contour_centered[:, 0] + 1j * contour_centered[:, 1]
    
    # Compute FFT
    coeffs = np.fft.fft(z)
    
    # Normalize by first coefficient (scale invariance)
    coeffs /= np.abs(coeffs[1])
    
    # Compute power spectrum
    power_spectrum = np.abs(coeffs[:n_harmonics])**2
    
    return power_spectrum

def compute_anisotropy_ratio(power_spectrum):
    """Compute anisotropy ratio based on Fourier power spectrum."""
    return np.sum(power_spectrum[3:]) / np.sum(power_spectrum) - power_spectrum[2] / np.sum(power_spectrum)

def compute_quadrupole_ratio(power_spectrum):
    """Compute quadrupole ratio based on Fourier power spectrum."""
    return power_spectrum[2] / np.sum(power_spectrum)

def compute_isotropy_ratio(power_spectrum):
    """Compute anisotropy ratio based on Fourier power spectrum."""
    return (power_spectrum[0] + power_spectrum[1]) / np.sum(power_spectrum)


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

def downsample_interpolate(x, M, method='linear'):
    """
    Downsamples an array x (size N) to size M while preserving trend via interpolation.

    Parameters:
        x (numpy array): Input array of size N.
        M (int): Desired output size (M < N).
        method (str): Interpolation type ('linear', 'cubic', 'nearest', etc.).

    Returns:
        numpy array: Downsampled array of size M.
    """
    N = len(x)
    original_indices = np.linspace(0, 1, N)
    new_indices = np.linspace(0, 1, M)
    interpolator = interp1d(original_indices, x, kind=method)
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



cfile=open(sys.argv[2],"r")
header=cfile.readline().split()
t=int(header[2])

header=cfile.readline().split()
lx=int(float(header[2]))
ly=int(float(header[3]))


start_line=(N+2)*int(float(sys.argv[3])) 
for i in range(start_line):
    cfile.readline()

pt_num = 0
time_line = 0

start_hole_time = -1
R0 = 8
area0 = pi * R0 * R0

thresh = 0.5 * pi * R0 * R0
cell_value_thresh = 0.01

voids=[[0. for q in range(lx)] for k in range(ly)]
voids_area = lx*ly
circ=0
hole_lifetime = 0
stats_hole_lifetime = []

avg_max_hole_area = 0
max_hole_area = 0
stats_max_hole_area = []

number_holes = 0

tv = np.zeros(end_frame - start_frame)
vv = np.zeros(end_frame - start_frame)
cct = np.zeros(end_frame - start_frame)
circularity = np.zeros(end_frame - start_frame)
ani = np.zeros(end_frame - start_frame)
radius_speed = np.zeros(end_frame - start_frame)
radius_time = np.zeros(end_frame - start_frame)

lattice_max=[[-1. for q in range(lx)] for k in range(ly)]
lattice_cell=[[-1 for q in range(lx)] for k in range(ly)]

histogram_size = int(1. * (end_frame - start_frame))
#histogram_hole_counts = np.zeros(histogram_size)
histogram_hole_area = np.zeros(histogram_size)
histogram_hole_circularity = np.zeros(histogram_size)
histogram_hole_ani = np.zeros(histogram_size)
histogram_hole_radius_speed_time = np.zeros(histogram_size)

histogram_hole_radius_speed = np.zeros(lx)
histogram_hole_radius_speed_counts = np.zeros(lx)
deltat = 100
dt = 0.1

start_value = 11
nemX = np.zeros(N)
nemY = np.zeros(N)
CoMX = np.zeros(N)
CoMY = np.zeros(N)

nn_count_for_time_series = 0
area_time_series = []
ani_time_series = []
circularity_time_series = []
radius_speed_time_series = []

for line in cfile:
    words=line.split()

    if words[0]=='t':
        if voids_area > thresh and voids_area>max_hole_area:
            max_hole_area = voids_area
            
        if voids_area > thresh and start_hole_time < 0:
            start_hole_time = t
        elif voids_area < thresh and start_hole_time >= 0:
            hole_lifetime += (t - start_hole_time)
            stats_hole_lifetime.append( (t - start_hole_time) * dt)
            avg_max_hole_area += max_hole_area
            stats_max_hole_area.append(max_hole_area / area0)
            max_hole_area = 0
            number_holes += 1

            #pearson_coeff_1 = np.corrcoef(np.array(area_time_series), np.array(ani_time_series))
            #pearson_coeff_2 = np.corrcoef(np.array(area_time_series), np.array(circularity_time_series))
            #print("Corr coeficients:", pearson_coeff_1, pearson_coeff_2)
            
            if (t - start_hole_time)/deltat >= 5:
                nn_count_for_time_series += 1
                x1 = np.array(area_time_series)
                x2 = np.array(ani_time_series)
                x3 = np.array(circularity_time_series)
                x4 = np.array(radius_speed_time_series)
                N = len(x1)
                M = histogram_size

                if(N>M):
                    y1 = downsample_interpolate(x1, M)
                    histogram_hole_area[:] = histogram_hole_area[:] + y1[:]
                    y2 = downsample_interpolate(x2, M)
                    histogram_hole_ani[:] = histogram_hole_ani[:] + y2[:]
                    y3 = downsample_interpolate(x3, M)
                    histogram_hole_circularity[:] = histogram_hole_circularity[:] + y3[:]
                    y4 = downsample_interpolate(x4, M)
                    histogram_hole_radius_speed_time[:] = histogram_hole_radius_speed_time[:] + y4[:]
                elif(M>N):
                    y1 = rescale_array(x1, M)
                    histogram_hole_area[:] = histogram_hole_area[:] + y1[:]
                    y2 = rescale_array(x2, M)
                    histogram_hole_ani[:] = histogram_hole_ani[:] + y2[:]
                    y3 = rescale_array(x3, M)
                    histogram_hole_circularity[:] = histogram_hole_circularity[:] + y3[:]
                    y4 = rescale_array(x4, M)
                    histogram_hole_radius_speed_time[:] = histogram_hole_radius_speed_time[:] + y4[:]
                else:
                    histogram_hole_area[:] = histogram_hole_area[:] + area_time_series[:]
                    histogram_hole_ani[:] = histogram_hole_ani[:] + ani_time_series[:]
                    histogram_hole_circularity[:] = histogram_hole_circularity[:] + circularity_time_series[:]
                    histogram_hole_radius_speed_time[:] = histogram_hole_radius_speed_time[:] + radius_speed_time_series[:]

            area_time_series = []
            ani_time_series = []
            circularity_time_series = []
            radius_speed_time_series = []
            start_hole_time = -1


        if start_hole_time >= 0:
            #temp = np.zeros(4)
            avg_x = 0.
            avg_y = 0.
            cell_list = []
            save_surface_sites = []
            area_sites = np.zeros(shape=(lx, ly))
            normal_x = np.zeros(lx*ly)
            normal_y = np.zeros(lx*ly)
            for i in range(lx*ly):
                y = int(i/lx)
                x = i - y * lx

                xx = (x + 1) % lx
                yy = (y + 1) % ly
                xxx = (x - 1 + lx) % lx
                yyy = (y - 1 + ly) % ly

                if voids[y][x] < cell_value_thresh:
                    area_sites[y,x] = 1
                    #temp[0] += cos(2*pi*x/lx); 
                    #temp[1] += sin(2*pi*x/lx);
                    #temp[2] += cos(2*pi*y/ly);
                    #temp[3] += sin(2*pi*y/ly);
                    exx = 1.
                    exxx = 1.
                    eyy = 1.
                    eyyy = 1.
                    if voids[yy][x] > cell_value_thresh:
                        cell_list.append(lattice_cell[yy][x])
                        save_surface_sites.append(i)
                        eyy = 0.
                    if voids[yyy][x] > cell_value_thresh:
                        cell_list.append(lattice_cell[yyy][x])
                        save_surface_sites.append(i)
                        eyyy = 0.
                    if voids[y][xx] > cell_value_thresh:
                        cell_list.append(lattice_cell[y][xx])
                        save_surface_sites.append(i)
                        exx = 0.
                    if voids[y][xxx] > cell_value_thresh:
                        cell_list.append(lattice_cell[y][xxx])
                        save_surface_sites.append(i)
                        exxx = 0.

                    normal_x[i] = 0.5 * (exx - exxx)
                    normal_y[i] = 0.5 * (eyy - eyyy)
                    #if save_surface_sites[len(save_surface_sites)-1] == i and normal_x[i] == 0 and normal_y[i]==0:
                        #print("Problem", i, voids_area, exx, exxx, eyy, eyyy, t)
            #avg_x = lx * ( atan2(-temp[1]/voids_area, -temp[0]/voids_area) + pi ) / (2*pi); 
            #avg_y = ly * ( atan2(-temp[3]/voids_area, -temp[2]/voids_area) + pi ) / (2*pi);

            step = int((t - start_hole_time)/deltat)
            radius = sqrt(voids_area/pi) / R0
            neigh_hole = list(set(cell_list))
            surface_sites = list(set(save_surface_sites))
            radius_time[time_line] = radius

            extended_array = apply_periodic_boundary_conditions(area_sites)
            contours = measure.find_contours(extended_array, level=0.5)
            if len(contours) == 0:
                raise ValueError("No contours detected!")
            # Use the largest contour
            contour = max(contours, key=len)
            contour[:,0] %= area_sites.shape[0]  # Wrap coordinates into original domain
            contour[:,1] %= area_sites.shape[1]  # Wrap coordinates into original domain
            # Compute Fourier descriptors
            power_spectrum = compute_fourier_descriptors(contour)
            # Compute anisotropy ratio
            anisotropy = compute_anisotropy_ratio(power_spectrum)
            quadrupole = compute_quadrupole_ratio(power_spectrum)
            isotropy = compute_isotropy_ratio(power_spectrum)
            #histogram_hole_ani[step] += anisotropy
            ani_time_series.append(anisotropy)
            ani[time_line] += anisotropy

            '''
            print(t, isotropy, quadrupole, anisotropy)
            if t == 60000:
                # Plot results
                plt.figure(figsize=(10, 4))
                plt.subplot(1, 2, 1)
                plt.imshow(area_sites, cmap='gray')
                plt.plot(contour[:, 1], contour[:, 0], 'r', linewidth=2)
                plt.title("Extracted Contour")
                plt.subplot(1, 2, 2)
                plt.plot(power_spectrum, marker='o')
                plt.title(f"Fourier Power Spectrum (Anisotropy: {anisotropy:.2f})")
                plt.xlabel("Harmonic Order")
                plt.ylabel("Power")
                plt.show()
            '''


            circ = 0.
            for i in neigh_hole:
                tangent = []
                vec_x = 0.
                vec_y = 0.
                for j in surface_sites:
                    y = int(j/lx)
                    x = j - y * lx

                    xx = (x + 1) % lx
                    yy = (y + 1) % ly
                    xxx = (x - 1 + lx) % lx
                    yyy = (y - 1 + ly) % ly

                    if lattice_cell[yy][x] == i:
                        tangent.append(j)
                        vec_x += normal_x[j]
                        vec_y += normal_y[j]
                    elif lattice_cell[yyy][x] == i:
                        tangent.append(j)
                        vec_x += normal_x[j]
                        vec_y += normal_y[j]
                    elif lattice_cell[y][xx] == i:
                        tangent.append(j)
                        vec_x += normal_x[j]
                        vec_y += normal_y[j]
                    elif lattice_cell[y][xxx] == i:
                        tangent.append(j)
                        vec_x += normal_x[j]
                        vec_y += normal_y[j]

                vec_x /= len(tangent)
                vec_y /= len(tangent)
                norm = sqrt(vec_x * vec_x + vec_y * vec_y)
                if norm == 0:
                    norm = 1
                dot_1 = (vec_x * nemX[i] + vec_y * nemY[i]) / norm
                dot_2 = (vec_x * (-nemX[i]) + vec_y * (-nemY[i])) / norm
                dot = -1
                if dot_1 >= 0:
                    #histogram_hole_circularity[step] += dot_1/len(neigh_hole)
                    circularity[time_line] += dot_1 / len(neigh_hole)
                    circ += dot_1 / len(neigh_hole)
                    dot = dot_1
                    #if t==13000:
                        #print(i, CoMX[i], CoMY[i], dot_1)
                elif dot_2 >= 0:
                    #histogram_hole_circularity[step] += dot_2/len(neigh_hole)
                    circularity[time_line] += dot_2 / len(neigh_hole)
                    circ += dot_2 / len(neigh_hole)
                    dot = dot_2
                    #if t==13000:
                        #print(i, CoMX[i], CoMY[i], dot_2)
                else:
                    print("Problem circularity", dot_1, dot_2, vec_x, vec_y, tangent)

            circularity_time_series.append(circ)
            area_time_series.append(voids_area/area0)
            #histogram_hole_area[step] += voids_area/area0
            #histogram_hole_counts[step] += 1

            if t - start_hole_time > 0:
                histogram_hole_radius_speed[int(radius)] += (radius_time[time_line] - radius_time[time_line - 1]) / (deltat * dt)
                #histogram_hole_radius_speed_time[step] += (radius_time[time_line] - radius_time[time_line - 1]) / (deltat * dt)
                radius_speed_time_series.append((radius_time[time_line] - radius_time[time_line - 1]) / (deltat * dt))
                radius_speed[time_line] = (radius_time[time_line] - radius_time[time_line - 1]) / (deltat * dt)
            else:
                histogram_hole_radius_speed[int(radius)] += (radius) / (deltat * dt)
                #histogram_hole_radius_speed_time[step] += (radius) / (deltat * dt)
                radius_speed_time_series.append((radius) / (deltat * dt))
                radius_speed[time_line] = radius / (deltat * dt)

            histogram_hole_radius_speed_counts[int(radius)] += 1


        tv[time_line] = t
        vv[time_line] = voids_area / area0
        cct[time_line] = circ
        voids=[[0. for q in range(lx)] for k in range(ly)]
        voids_area = lx*ly
        lattice_max=[[-1. for q in range(lx)] for k in range(ly)]
        lattice_cell=[[-1 for q in range(lx)] for k in range(ly)]

        t=int(float(words[2]))
        time_line += 1
        if end_frame == start_frame + time_line:
            break

        pt_num = 0

    elif words[0]=='b':
        lx=int(float(words[2]))
        ly=int(float(words[3]))

    else:

        CoMX[pt_num]=float(words[2])
        CoMY[pt_num]=float(words[3])
        nemX[pt_num]=float(words[9])
        nemY[pt_num]=float(words[10])
        for i in range(start_value,len(words),2):
            site=int(float(words[i]))
            value=float(words[i+1])
            yy=int(site/lx)
            xx=site-int(yy*lx)

            if voids[yy][xx] < cell_value_thresh and voids[yy][xx] + value > cell_value_thresh:
                voids_area -= 1
            voids[yy][xx] += value

            if value > lattice_max[yy][xx]:
                lattice_max[yy][xx] = value
                lattice_cell[yy][xx] = pt_num

        pt_num += 1


if voids_area > thresh and voids_area>max_hole_area:
    max_hole_area = voids_area            

if start_hole_time >= 0:
    hole_lifetime += (t - start_hole_time)
    stats_hole_lifetime.append( (t - start_hole_time) * dt )
    avg_max_hole_area += max_hole_area
    stats_max_hole_area.append(max_hole_area / area0)
    number_holes += 1

    if (t - start_hole_time)/deltat >= 5:
        nn_count_for_time_series += 1
        x1 = np.array(area_time_series)
        x2 = np.array(ani_time_series)
        x3 = np.array(circularity_time_series)
        x4 = np.array(radius_speed_time_series)
        N = len(x1)
        M = histogram_size
        if(N>M):
            y = downsample_interpolate(x1, M)
            histogram_hole_area[:] = histogram_hole_area[:] + y[:]
            y = downsample_interpolate(x2, M)
            histogram_hole_ani[:] = histogram_hole_ani[:] + y[:]
            y = downsample_interpolate(x3, M)
            histogram_hole_circularity[:] = histogram_hole_circularity[:] + y[:]
            y = downsample_interpolate(x4, M)
            histogram_hole_radius_speed_time[:] = histogram_hole_radius_speed_time[:] + y[:]
        elif(M>N):
            y = rescale_array(x1, M)
            histogram_hole_area[:] = histogram_hole_area[:] + y[:]
            y = rescale_array(x2, M)
            histogram_hole_ani[:] = histogram_hole_ani[:] + y[:]
            y = rescale_array(x3, M)
            histogram_hole_circularity[:] = histogram_hole_circularity[:] + y[:]
            y = rescale_array(x4, M)
            histogram_hole_radius_speed_time[:] = histogram_hole_radius_speed_time[:] + y[:]
        else:
            histogram_hole_area[:] = histogram_hole_area[:] + area_time_series[:]
            histogram_hole_ani[:] = histogram_hole_ani[:] + ani_time_series[:]
            histogram_hole_circularity[:] = histogram_hole_circularity[:] + circularity_time_series[:]
            histogram_hole_radius_speed_time[:] = histogram_hole_radius_speed_time[:] + radius_speed_time_series[:]

    '''
    r1 = np.corrcoef(np.array(area_time_series), np.array(ani_time_series))
    r2 = np.corrcoef(np.array(area_time_series), np.array(circularity_time_series))
    print("Corr coefficient:", r1, r2)
    x = np.array(area_time_series)
    y = np.array(ani_time_series)
    corr1 = np.correlate(x - np.mean(x), y - np.mean(y), mode="full") / (np.std(x) * np.std(y) * len(x))
    y = np.array(circularity_time_series)
    corr2 = np.correlate(x - np.mean(x), y - np.mean(y), mode="full") / (np.std(x) * np.std(y) * len(x))
    lags = np.arange(-len(x) + 1, len(x))
    x = np.array(ani_time_series)
    y = np.array(circularity_time_series)
    corr3 = np.correlate(x - np.mean(x), y - np.mean(y), mode="full") / (np.std(x) * np.std(y) * len(x))
    '''


tv[time_line] = t
vv[time_line] = voids_area / area0
cct[time_line] = circ


for i in range(len(histogram_hole_area)):
    if nn_count_for_time_series>0:
        histogram_hole_area[i] /= nn_count_for_time_series
        histogram_hole_radius_speed_time[i] /= nn_count_for_time_series
        histogram_hole_circularity[i] /= nn_count_for_time_series
        histogram_hole_ani[i] /= nn_count_for_time_series

for i in range(len(histogram_hole_radius_speed_counts)):
    if histogram_hole_radius_speed_counts[i]>0:
        histogram_hole_radius_speed[i] /= histogram_hole_radius_speed_counts[i]


if variable==1:

    with open('voids_stats.txt', 'w') as f:
        if number_holes > 0:
            print('Avg stats:', number_holes, (hole_lifetime * dt)/number_holes, (avg_max_hole_area / area0)/number_holes, 'Lifetimes:', *stats_hole_lifetime, 'MaxAreas:', *stats_max_hole_area, file=f)
        else:
            print('Avg stats:', number_holes, 0, 0, 'Lifetimes:', 0, 'MaxAreas:', 0, file=f)

    with open('voids_histogram.txt', 'w') as f:
        for i in range(len(histogram_hole_area)):
            print(i * 1./histogram_size, histogram_hole_area[i], histogram_hole_circularity[i], histogram_hole_radius_speed_time[i], histogram_hole_ani[i], file=f)


    with open('voids_area_time.txt', 'w') as f:
        for i in range(len(vv)):
            print(i, vv[i] * area0, cct[i], file=f)

    '''
    with open('voids_histogram_area.txt', 'w') as f:
        for i in range(len(histogram_hole_area)):
            print(i * 1./histogram_size, histogram_hole_area[i], file=f)

    with open('voids_histogram_circularity.txt', 'w') as f:
        for i in range(len(histogram_hole_circularity)):
            print(i * 1./histogram_size, histogram_hole_circularity[i], file=f)

    with open('voids_histogram_radius_speed_time.txt', 'w') as f:
        for i in range(len(histogram_hole_radius_speed_time)):
            print(i * 1./histogram_size, histogram_hole_radius_speed_time[i], file=f)

    with open('voids_histogram_ani.txt', 'w') as f:
        for i in range(len(histogram_hole_ani)):
            print(i * 1./histogram_size, histogram_hole_ani[i], file=f)
    '''

    with open('voids_histogram_radius_speed.txt', 'w') as f:
        for i in range(len(histogram_hole_radius_speed)):
            print(i, histogram_hole_radius_speed[i], file=f)


if variable==2:

    #plt.figure(figsize=(6,6))
    fig = plt.figure(figsize=(5.452423529, 4.089317647))
    plt.plot(tv, vv, '-o')
    fig = plt.figure(figsize=(5.452423529, 4.089317647))
    plt.plot(tv, circularity, '-o')
    fig = plt.figure(figsize=(5.452423529, 4.089317647))
    plt.plot(tv, ani, '-o')
    #plt.plot(tv, radius_time, '-s')

    #fig2 = plt.figure(figsize=(5.452423529, 4.089317647))
    #plt.plot(tv, radius_speed, '-o')
    #plt.show()

    fig = plt.figure(figsize=(5.452423529, 4.089317647))
    #plt.plot(histogram_hole_radius, '-o')
    #plt.plot(histogram_hole_area, '-v')
    #plt.plot(histogram_hole_circularity, '-s')
    #plt.plot(histogram_hole_ellipse, '-p')
    plt.plot(histogram_hole_ani, '-s')

    fig = plt.figure(figsize=(5.452423529, 4.089317647))
    plt.plot(histogram_hole_area, '-o')

    fig = plt.figure(figsize=(5.452423529, 4.089317647))
    plt.plot(histogram_hole_circularity, '-o')

    #fig = plt.figure(figsize=(5.452423529, 4.089317647))
    #plt.plot(histogram_hole_radius_speed_time, '-o')
    plt.show()
    #fig = plt.figure(figsize=(5.452423529, 4.089317647))
    #plt.plot(histogram_hole_radius_speed, '-o')
    #plt.show()

    #plt.ylabel(r'$C_v$', fontsize=18)
    #plt.xlabel('R', fontsize=18)
    #plt.xticks(fontsize=18)
    #plt.yticks(fontsize=18)
    #plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
    #plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/mean_velocity_1.png", transparent=True)
    #ax = plt.gca()
    #ax.set_aspect('equal', adjustable='box')
    #ax.set_xlim([0, lx])
    #ax.set_ylim([0, ly])
