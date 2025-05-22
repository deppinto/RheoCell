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

r_thresh = 9
velocity_defects_plus = []
beta_tail_degree = []


for i in range(start_frame, end_frame, 1):

    #if i%10==0:
        #print(i)

    if variable == 1 or variable == 3:
        if i<10:
            fileS = sys.argv[2] + "shape_field_00" + str(i) + ".txt"
            fileQ = sys.argv[2] + "nematic_field_00" + str(i) + ".txt"
        elif i<100:
            fileS = sys.argv[2] + "shape_field_0" + str(i) + ".txt"
            fileQ = sys.argv[2] + "nematic_field_0" + str(i) + ".txt"
        else:
            fileS = sys.argv[2] + "shape_field_" + str(i) + ".txt"
            fileQ = sys.argv[2] + "nematic_field_" + str(i) + ".txt"
    '''
    elif variable == 2 or variable == 4:
        if i<10:
            file = sys.argv[2] + "nematic_field_00" + str(i) + ".txt"
        elif i<100:
            file = sys.argv[2] + "nematic_field_0" + str(i) + ".txt"
        else:
            file = sys.argv[2] + "nematic_field_" + str(i) + ".txt"
    '''

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
            Z_S01[int(yy)][int(xx)]=S01

        read_line += 1
    Sfile.close()


    LLX = lx
    LLY = ly

    vecfield_nx = [[0. for q in range(0, LLX)] for p in range(0, LLY)]
    vecfield_ny = [[0. for q in range(0, LLX)] for p in range(0, LLY)]
    vecfield_nx = Z_Sx
    vecfield_ny = Z_Sy

    vecfield_Q00 = [[0. for j in range(0, LLX)] for i in range(0, LLY)]
    vecfield_Q01 = [[0. for j in range(0, LLX)] for i in range(0, LLY)]
    vecfield_Q00 = Z_S00
    vecfield_Q01 = Z_S01

    winding_number = [[0. for q in range(0, LLX)] for p in range(0, LLY)]
    for p in range(0, LLY):
        for q in range(0, LLX):
            y1 = p
            x1 = q

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


    plus_defects_old_x = plus_defects_x
    minus_defects_old_x = minus_defects_x
    plus_defects_old_y = plus_defects_y
    minus_defects_old_y = minus_defects_y
    matched_defects = [0 for zz in range(len(plus_defects_x))]

    plus_defects_x = []
    minus_defects_x = []
    plus_defects_y = []
    minus_defects_y = []

    defect_thresh = 0.05
    for p in range(0,LLY):
        for q in range(0,LLX):
            # keep this just in case our other symmetries give us integer defects
            if (abs(winding_number[p][q]) > 0.5-defect_thresh) and (abs(winding_number[p][q]) < 0.5+defect_thresh):
                # charge sign
                s = np.sign(winding_number[p][q])
                # bfs
                sum_x, sum_y, n = collapse(p, q, s, LLX, LLY, winding_number, rng = [0.5-defect_thresh, 0.5+defect_thresh])
                x,y = sum_x/n,sum_y/n
                # add defect to list
                num = 0
                den = 0
                for (dx, dy) in [(0, 0), (0, 1), (1, 1), (1, 0)]:
                    # coordinates of nodes around the defect
                    kk = (int(x) + LLX + dx) % LLX
                    ll = (int(y) + LLY + dy) % LLY
                    # derivative at these points
                    dxQxx = .5*(vecfield_Q00[ll][(kk+1) % LLX] - vecfield_Q00[ll][(kk-1+LLX) % LLX])
                    dxQxy = .5*(vecfield_Q01[ll][(kk+1) % LLX] - vecfield_Q01[ll][(kk-1+LLX) % LLX])
                    dyQxx = .5*(vecfield_Q00[(ll+1) % LLY][kk] - vecfield_Q00[(ll-1+LLY) % LLY][kk])
                    dyQxy = .5*(vecfield_Q01[(ll+1) % LLY][kk] - vecfield_Q01[(ll-1+LLY) % LLY][kk])
                    # accumulate numerator and denominator
                    num += s*dxQxy - dyQxx
                    den += dxQxx + s*dyQxy
                psi = s/(2.-s)*atan2(num, den)
                if s==1:
                    plus_defects_x.append(x)
                    plus_defects_y.append(y)

                    x_direction = cos(psi)
                    y_direction = sin(psi)
                    v = np.array([x_direction, y_direction])

                    eps = 1.5
                    for nn in range(1, 9):
                        pp = np.array([x,y]) + nn * v
                        x_min, x_max = int(np.floor(pp[0] - eps)), int(np.ceil(pp[0] + eps))
                        y_min, y_max = int(np.floor(pp[1] - eps)), int(np.ceil(pp[1] + eps))
                        for xi in range(x_min, x_max + 1):
                            for yi in range(y_min, y_max + 1):
                                # Compute perpendicular distance from lattice point to line point
                                xxx = xi
                                yyy = yi
                                if xxx<0:
                                    xxx += LLX
                                if xxx>=LLX:
                                    xxx -= LLX
                                if yyy<0:
                                    yyy += LLY
                                if yyy>=LLY:
                                    yyy -= LLY

                                dd = np.linalg.norm(np.array([xxx, yyy]) - pp)
                                if dd <= eps:
                                    angle = np.acos((Z_Qx[yyy][xxx] * x_direction + Z_Qy[yyy][xxx] * y_direction))
                                    #angle = np.acos((Z_Sx[yyy][xxx] * x_direction + Z_Sy[yyy][xxx] * y_direction))
                                    #if angle>pi/2:
                                        #angle -= pi/2
                                    beta_tail_degree.append(angle)


                    if i > start_frame:
                        for zz in range(len(plus_defects_old_x)):
                            if matched_defects[zz] == 1:
                                continue

                            dist_x = x - plus_defects_old_x[zz]
                            if dist_x<-lx/2:
                                dist_x += lx

                            dist_y = y - plus_defects_old_y[zz]
                            if dist_y<-ly/2:
                                dist_y += ly

                            dist = dist_x * dist_x + dist_y * dist_y
                            if dist < r_thresh and dist>0:
                                matched_defects[zz] = 1
                                norm = sqrt(dist)
                                angle = np.acos((dist_x * cos(psi) + dist_y * sin(psi)) / norm)
                                velocity_defects_plus.append(norm * cos(angle))
                                break
                            

                    #cset1 = plt.plot(x, y, 'go', markersize=10)
                    #cset1 = plt.arrow(x, y, 4*cos(psi), 4*sin(psi), color='g', head_width=1.5, head_length=1.5, width=0.5)
                elif s==-1:
                    minus_defects_x.append(x)
                    minus_defects_y.append(y)

                    #cset1 = plt.plot(x, y, 'b^', markersize=10)

if variable == 1 or variable == 2:
    for j in range(len(velocity_defects_plus)):
        print(velocity_defects_plus[j])

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
    counts, bins = np.histogram(beta_tail_degree, bins=20)
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
    plt.ylabel(r'$P(\beta)$', fontsize=18)
    plt.xlabel(r'$\beta$', fontsize=18)
    #plt.xticks(fontsize=18)
    #plt.yticks(fontsize=18)
    plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
    #plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/mean_velocity_1.png", transparent=True)
    plt.show()




