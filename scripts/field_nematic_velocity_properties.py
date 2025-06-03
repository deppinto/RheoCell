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
    print(sys.argv[0]," [topology file] [frames location] [start frame] [end frame] [1 or 2:save conf; 3 or 4:make plot]")
    sys.exit(1)

variable=int(float(sys.argv[5])) 
start_frame=int(float(sys.argv[3])) 
end_frame=int(float(sys.argv[4])) 


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
    p = 2.
    ang = atan2(abs(a[0]*b[1]-a[1]*b[0]), a[0]*b[0]+a[1]*b[1])

    if(ang > pi/2.):
        b = [-i for i in b]

    #while(abs(ang) > np.pi/p + 1e-3):
        #b = rotate(b,p)
        #ang = atan2(abs(a[0]*b[1]-a[1]*b[0]), a[0]*b[0]+a[1]*b[1])

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


tfile=open(sys.argv[1],"r")
N=int(tfile.readline().split()[0])
species=[]
line=tfile.readline()
lline=line.split()
for l in lline:
        species.append(int(l))
tfile.close()

numspecies=len(set(species))


plus_defects_x = []
minus_defects_x = []
plus_defects_y = []
minus_defects_y = []

theta_degree = []

for i in range(start_frame, end_frame, 1):

    #if i%10==0:
        #print(i)

    if variable == 1 or variable == 3:
        if i<10:
            fileS = sys.argv[2] + "shape_field_00" + str(i) + ".txt"
            fileQ = sys.argv[2] + "nematic_field_00" + str(i) + ".txt"
            fileV = sys.argv[2] + "velocity_field_00" + str(i) + ".txt"
        elif i<100:
            fileS = sys.argv[2] + "shape_field_0" + str(i) + ".txt"
            fileQ = sys.argv[2] + "nematic_field_0" + str(i) + ".txt"
            fileV = sys.argv[2] + "velocity_field_0" + str(i) + ".txt"
        else:
            fileS = sys.argv[2] + "shape_field_" + str(i) + ".txt"
            fileQ = sys.argv[2] + "nematic_field_" + str(i) + ".txt"
            fileV = sys.argv[2] + "velocity_field_" + str(i) + ".txt"


    Qfile=open(fileQ,"r")
    header=Qfile.readline().split()
    t=int(header[2])

    header=Qfile.readline().split()
    lx=int(float(header[2]))
    ly=int(float(header[3]))

    Z_Qx=[[0 for q in range(lx)] for k in range(ly)]
    Z_Qy=[[0 for q in range(lx)] for k in range(ly)]
    Z_Q00=[[0 for q in range(lx)] for k in range(ly)]
    Z_Q01=[[0 for q in range(lx)] for k in range(ly)]
    start_value=0
    read_line = 0
    for line in Qfile:

        if read_line==0:
            read_line+=1
            continue

        words=line.split()
        Z_Qx=[[0 for q in range(lx)] for k in range(ly)]
        Z_Qy=[[0 for q in range(lx)] for k in range(ly)]
        Z_Q00=[[0 for q in range(lx)] for k in range(ly)]
        Z_Q01=[[0 for q in range(lx)] for k in range(ly)]
        for i in range(start_value,len(words),4):
            xx=float(words[i])
            yy=float(words[i+1])
            value_x=float(words[i+2])
            value_y=float(words[i+3])

            Q00=float(words[i+2])
            Q01=float(words[i+3])
            if Q00 != 0 and Q01 != 0:
                S = sqrt(Q00**2 + Q01**2)
                nx = sqrt(2*S) * sqrt((1 + Q00/S)/2)
                ny = sqrt(2*S) * np.sign(Q01) * sqrt((1 - Q00/S)/2)
            else:
                nx = 0.
                ny = 0.

            Z_Qx[int(yy)][int(xx)]=nx
            Z_Qy[int(yy)][int(xx)]=ny
            Z_Q00[int(yy)][int(xx)]=Q00
            Z_Q01[int(yy)][int(xx)]=Q01

        read_line += 1
    Qfile.close()


    Sfile=open(fileS,"r")
    header=Sfile.readline().split()
    t=int(header[2])

    header=Sfile.readline().split()
    lx=int(float(header[2]))
    ly=int(float(header[3]))

    Z_Sx=[[0 for q in range(lx)] for k in range(ly)]
    Z_Sy=[[0 for q in range(lx)] for k in range(ly)]
    Z_S00=[[0 for q in range(lx)] for k in range(ly)]
    Z_S01=[[0 for q in range(lx)] for k in range(ly)]
    start_value=0
    read_line = 0
    for line in Sfile:

        if read_line==0:
            read_line+=1
            continue

        words=line.split()
        Z_Sx=[[0 for q in range(lx)] for k in range(ly)]
        Z_Sy=[[0 for q in range(lx)] for k in range(ly)]
        Z_S00=[[0 for q in range(lx)] for k in range(ly)]
        Z_S01=[[0 for q in range(lx)] for k in range(ly)]
        for i in range(start_value,len(words),4):
            xx=float(words[i])
            yy=float(words[i+1])
            value_x=float(words[i+2])
            value_y=float(words[i+3])

            S00=float(words[i+2])
            S01=float(words[i+3])
            if S00 != 0 and S01 != 0:
                S = sqrt(S00**2 + S01**2)
                nx = sqrt(2*S) * sqrt((1 + S00/S)/2)
                ny = sqrt(2*S) * np.sign(S01) * sqrt((1 - S00/S)/2)
            else:
                nx = 0.
                ny = 0.

            Z_Sx[int(yy)][int(xx)]=nx
            Z_Sy[int(yy)][int(xx)]=ny
            Z_S00[int(yy)][int(xx)]=S00

            Sx = Z_Qx[int(yy)][int(xx)]
            Sy = Z_Qy[int(yy)][int(xx)]
            value_x = nx
            value_y = ny

            norm1 = sqrt(value_x * value_x + value_y * value_y)
            norm2 = sqrt(Sx * Sx + Sy * Sy)

            value_x = value_x / norm1
            value_y = value_y / norm1
            Sx = Sx / norm2
            Sy = Sy / norm2
            angle1 = 0
            angle2 = 0

            if norm1<1e-8 or norm2<1e-8:
                continue

            angle1 = np.acos(value_x * Sx + value_y * Sy)
            if isnan(angle1):
                angle1 = 0

            Sx = -Sx
            Sy = -Sy

            angle2 = np.acos(value_x * Sx + value_y * Sy)
            if isnan(angle2):
                angle2 = 0

            if angle1 < angle2:
                theta_degree.append(angle1)
                angle = angle1
            else:
                theta_degree.append(angle2)
                angle = angle2
            Z_S01[int(yy)][int(xx)]=S01

        read_line += 1
    Sfile.close()


    Vfile=open(fileV,"r")
    header=Vfile.readline().split()
    t=int(header[2])

    header=Vfile.readline().split()
    lx=int(float(header[2]))
    ly=int(float(header[3]))

    Z_Vx=[[0 for q in range(lx)] for k in range(ly)]
    Z_Vy=[[0 for q in range(lx)] for k in range(ly)]
    start_value=0
    read_line = 0
    for line in Vfile:

        if read_line==0:
            read_line+=1
            continue

        words=line.split()
        Z_Vx=[[0 for q in range(lx)] for k in range(ly)]
        Z_Vy=[[0 for q in range(lx)] for k in range(ly)]
        for i in range(start_value,len(words),4):
            xx=float(words[i])
            yy=float(words[i+1])
            value_x=float(words[i+2])
            value_y=float(words[i+3])

            Z_Vx[int(yy)][int(xx)]=value_x
            Z_Vy[int(yy)][int(xx)]=value_y

            '''
            Sx = Z_Sx[int(yy)][int(xx)]
            Sy = Z_Sy[int(yy)][int(xx)]

            norm1 = sqrt(value_x * value_x + value_y * value_y)
            norm2 = sqrt(Sx * Sx + Sy * Sy)

            value_x = value_x / norm1
            value_y = value_y / norm1
            Sx = Sx / norm2
            Sy = Sy / norm2
            angle1 = 0
            angle2 = 0

            if norm1<1e-8 or norm2<1e-8:
                continue

            angle1 = np.acos(value_x * Sx + value_y * Sy)
            if isnan(angle1):
                angle1 = 0

            Sx = -Sx
            Sy = -Sy

            angle2 = np.acos(value_x * Sx + value_y * Sy)
            if isnan(angle2):
                angle2 = 0

            if angle1 < angle2:
                theta_degree.append(angle1)
                angle = angle1
            else:
                theta_degree.append(angle2)
                angle = angle2
            '''

        read_line += 1
    Vfile.close()



if variable == 1 or variable == 2:
    for j in range(len(theta_degree)):
        print(theta_degree[j])

if variable == 3 or variable == 4:
    '''
    plt.figure(figsize=(5.452423529,4.089317647))
    counts, bins = np.histogram(velocity_defects_plus, bins=20)
    bin_width = abs(bins[1] - bins[0]) / 2
    bin_length = len(bins)
    total_counts = sum(counts)
    print(total_counts)
    probability = []
    for i in range(len(counts)):
        probability.append(float(counts[i]) / float(total_counts))
        bins[i] = bins[i] + bin_width
    #plt.stairs(counts, bins)
    plt.plot(bins[0:bin_length-1], probability, '--o')
    plt.ylabel(r'$P(v)$', fontsize=18)
    plt.xlabel(r'Velocity $+1/2$', fontsize=18)
    #plt.xticks(fontsize=18)
    #plt.yticks(fontsize=18)
    plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
    #plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/mean_velocity_1.png", transparent=True)
    plt.show()
    '''

    plt.figure(figsize=(5.452423529,4.089317647))
    counts, bins = np.histogram(theta_degree, bins=20)
    bin_width = abs(bins[1] - bins[0]) / 2
    bin_length = len(bins)
    total_counts = sum(counts)
    print(total_counts)
    probability = []
    for i in range(len(counts)):
        probability.append(float(counts[i]) / float(total_counts))
        bins[i] = (bins[i] + bin_width) * 180 / pi
    #plt.stairs(counts, bins)
    plt.plot(bins[0:bin_length-1], probability, '--o')
    plt.ylabel(r'$P_{QS}(\theta)$', fontsize=18)
    plt.xlabel(r'$\theta$', fontsize=18)
    #plt.xticks(fontsize=18)
    #plt.yticks(fontsize=18)
    plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
    #plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/mean_velocity_1.png", transparent=True)
    plt.show()




