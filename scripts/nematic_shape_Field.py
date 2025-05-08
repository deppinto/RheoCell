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
matplotlib.use('Agg')

if len(sys.argv)!=5:
    print(sys.argv[0]," [topology file] [nematic file] [shape file] [1:save conf; 2:make plot]")
    sys.exit(1)

variable=int(float(sys.argv[4])) 

def set_walls(lx,ly, walls):
    for y in range(ly):
        for x in range(lx):
            k = x + y * lx
            walls[k] = exp(-float(y)/5) + exp(-float(ly-y-1)/5)


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

walls = [0. for i in range(lx*ly)]
set_walls(lx,ly,walls)

x=np.arange(0,lx,1)
y=np.arange(0,ly,1)
Z_x=[[0 for q in range(lx)] for k in range(ly)]
Z_y=[[0 for q in range(lx)] for k in range(ly)]
Z_Q00=[[0 for q in range(lx)] for k in range(ly)]
Z_Q01=[[0 for q in range(lx)] for k in range(ly)]
fig = plt.figure(figsize=(6,6))
start_value=0
read_line=0 
for line in cfile:

    if read_line==0:
        read_line+=1
        continue

    words=line.split()
    Z_x=[[0 for q in range(lx)] for k in range(ly)]
    Z_y=[[0 for q in range(lx)] for k in range(ly)]
    Z_Q00=[[0 for q in range(lx)] for k in range(ly)]
    Z_Q01=[[0 for q in range(lx)] for k in range(ly)]
    for i in range(start_value,len(words),4):
        xx=float(words[i])
        yy=float(words[i+1])

        Q00=float(words[i+2])
        Q01=float(words[i+3])
        if Q00 != 0 and Q01 != 0:
            S = sqrt(Q00**2 + Q01**2)
            nx = sqrt(2*S) * sqrt((1 + Q00/S)/2)
            ny = sqrt(2*S) * np.sign(Q01) * sqrt((1 - Q00/S)/2)
        else:
            nx = 0.
            ny = 0.

        Z_x[int(yy)][int(xx)]=nx
        Z_y[int(yy)][int(xx)]=ny
        Z_Q00[int(yy)][int(xx)]=Q00
        Z_Q01[int(yy)][int(xx)]=Q01

        if int(xx)%4==0 and int(yy)%4==0:
            cset1 = plt.arrow(xx, yy, 1*nx, 1*ny, width=0.25, color="k", head_width=0)
            cset1 = plt.arrow(xx, yy, -1*nx, -1*ny, width=0.25, color="k", head_width=0)
    read_line+=1


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


LLX = lx
LLY = ly
vecfield_nx = [[0. for j in range(0, LLX)] for i in range(0, LLY)]
vecfield_ny = [[0. for j in range(0, LLX)] for i in range(0, LLY)]
vecfield_Q00 = [[0. for j in range(0, LLX)] for i in range(0, LLY)]
vecfield_Q01 = [[0. for j in range(0, LLX)] for i in range(0, LLY)]

vecfield_nx = Z_x
vecfield_ny = Z_y
vecfield_Q00 = Z_Q00
vecfield_Q01 = Z_Q01

winding_number = [[0. for j in range(0, LLX)] for i in range(0, LLY)]
norm_gradient_Q = [[0. for j in range(0, LLX)] for i in range(0, LLY)]
splay_field = [[0. for j in range(0, LLX)] for i in range(0, LLY)]
bend_field = [[0. for j in range(0, LLX)] for i in range(0, LLY)]
for p in range(0, LLY):
    for q in range(0, LLX):
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

        dxQxx = .5*(vecfield_Q00[p][(q+1) % LLX] - vecfield_Q00[p][(q-1+LLX) % LLX])
        dxQxy = .5*(vecfield_Q01[p][(q+1) % LLX] - vecfield_Q01[p][(q-1+LLX) % LLX])
        dyQxx = .5*(vecfield_Q00[(p+1) % LLY][q] - vecfield_Q00[(p-1+LLY) % LLY][q])
        dyQxy = .5*(vecfield_Q01[(p+1) % LLY][q] - vecfield_Q01[(p-1+LLY) % LLY][q])
        norm_gradient_Q[p][q] = dxQxx * dxQxx + dyQxx * dyQxx + 2 * (dxQxy * dxQxy + dyQxy * dyQxy)

        dxNx = .5*(vecfield_nx[p][(q+1) % LLX] - vecfield_nx[p][(q-1+LLX) % LLX])
        dxNy = .5*(vecfield_ny[p][(q+1) % LLX] - vecfield_ny[p][(q-1+LLX) % LLX])
        dyNx = .5*(vecfield_nx[(p+1) % LLY][q] - vecfield_nx[(p-1+LLY) % LLY][q])
        dyNy = .5*(vecfield_ny[(p+1) % LLY][q] - vecfield_ny[(p-1+LLY) % LLY][q])
        splay_field[p][q] = (dxNx + dyNy) * (dxNx + dyNy)
        bend_field[p][q] = (vecfield_ny[p][q] * dxNy - vecfield_nx[p][q] * dxNx)**2 + (vecfield_ny[p][q] * dyNy - vecfield_nx[p][q] * dyNx)**2

#z_min, z_max = 0., np.abs(norm_gradient_Q).max()
#cset1 = plt.imshow(norm_gradient_Q, cmap='RdBu', interpolation='nearest', vmin=0, vmax=z_max)
#z_min, z_max = 0., np.abs(splay_field).max()
#cset1 = plt.imshow(splay_field, cmap='RdBu', interpolation='nearest', vmin=0, vmax=z_max)
#z_min, z_max = 0., np.abs(bend_field).max()
#cset1 = plt.imshow(bend_field, cmap='RdBu', interpolation='nearest', vmin=0, vmax=z_max)



charge = 1.0/2.0
thresh = 0.05
for p in range(0,LLY):
    for q in range(0,LLX):
        # detect simplest charge 1/p defects
        if  (abs(winding_number[p][q]) > charge - thresh) and (abs(winding_number[p][q]) < charge + thresh):
            # charge sign
            s = np.sign(winding_number[p][q])
            # bfs
            sum_x, sum_y, n = collapse(p, q, s, LLX, LLY, winding_number, rng = [charge - thresh, charge + thresh])
            x,y = sum_x/n,sum_y/n
            # compute angle, see doi:10.1039/c6sm01146b
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
                cset1 = plt.plot(x, y, 'ko', markersize=10)
                cset1 = plt.arrow(x, y, 4*cos(psi), 4*sin(psi), color='k', head_width=1.5, head_length=1.5, width=0.5)
            elif s==-1:
                cset1 = plt.plot(x, y, 'ko', markersize=5)
                cset1 = plt.arrow(x, y, 4*cos(psi), 4*sin(psi), color='k', head_width=1.5, head_length=1.5, width=0.5)
                cset1 = plt.arrow(x, y, 4*cos(psi + 2.0944), 4*sin(psi + 2.0944), color='k', head_width=1.5, head_length=1.5, width=0.5)
                cset1 = plt.arrow(x, y, 4*cos(psi + 2 * 2.0944), 4*sin(psi + 2 * 2.0944), color='k', head_width=1.5, head_length=1.5, width=0.5)


        # keep this just in case our other symmetries give us integer defects
        elif (abs(winding_number[p][q]) > 1 - thresh) and (abs(winding_number[p][q])< 1 + thresh):
            # charge sign
            s = np.sign(winding_number[p][q])
            # bfs
            sum_x, sum_y, n = collapse(p, q, s, LLX, LLY, winding_number, rng = [1-thresh, 1+thresh])
            x,y = sum_x/n,sum_y/n
            # add defect to list
            if s==1:
                cset1 = plt.plot(x, y, 'r*', markersize=10)
            elif s==-1:
                cset1 = plt.plot(x, y, 'kX', markersize=10)




cfile=open(sys.argv[3],"r")
header=cfile.readline().split()
t=int(header[2])

header=cfile.readline().split()
lx=int(float(header[2]))
ly=int(float(header[3]))

x=np.arange(0,lx,1)
y=np.arange(0,ly,1)
Z_xx=[[0 for q in range(lx)] for k in range(ly)]
Z_yy=[[0 for q in range(lx)] for k in range(ly)]
Z_S00=[[0 for q in range(lx)] for k in range(ly)]
Z_S01=[[0 for q in range(lx)] for k in range(ly)]
alignment=[[0 for q in range(lx)] for k in range(ly)]
start_value=0
read_line=0 
for line in cfile:

    if read_line==0:
        read_line+=1
        continue

    words=line.split()
    Z_xx=[[0 for q in range(lx)] for k in range(ly)]
    Z_yy=[[0 for q in range(lx)] for k in range(ly)]
    Z_S00=[[0 for q in range(lx)] for k in range(ly)]
    Z_S01=[[0 for q in range(lx)] for k in range(ly)]
    for i in range(start_value,len(words),4):
        xx=float(words[i])
        yy=float(words[i+1])

        S00=float(words[i+2])
        S01=float(words[i+3])
        if S00 != 0 and S01 != 0:
            S = sqrt(S00**2 + S01**2)
            nx = sqrt(2*S) * sqrt((1 + S00/S)/2)
            ny = sqrt(2*S) * np.sign(S01) * sqrt((1 - S00/S)/2)
        else:
            nx = 0.
            ny = 0.

        Z_xx[int(yy)][int(xx)]=nx
        Z_yy[int(yy)][int(xx)]=ny
        Z_S00[int(yy)][int(xx)]=S00
        Z_S01[int(yy)][int(xx)]=S01

        dot_prod = Z_xx[int(yy)][int(xx)] * Z_x[int(yy)][int(xx)] + Z_yy[int(yy)][int(xx)] * Z_y[int(yy)][int(xx)]
        norm1 = sqrt(Z_xx[int(yy)][int(xx)] * Z_xx[int(yy)][int(xx)] + Z_yy[int(yy)][int(xx)] * Z_yy[int(yy)][int(xx)])
        norm2 = sqrt(Z_x[int(yy)][int(xx)] * Z_x[int(yy)][int(xx)] + Z_y[int(yy)][int(xx)] * Z_y[int(yy)][int(xx)])
        if norm1 == 0:
            norm1 = 1.
        if norm2 == 0:
            norm2 = 1.
        ang1 = np.acos(dot_prod / (norm1 * norm2))

        dot_prod = Z_xx[int(yy)][int(xx)] * (-1.) * Z_x[int(yy)][int(xx)] + Z_yy[int(yy)][int(xx)] * (-1.) * Z_y[int(yy)][int(xx)]
        ang2 = np.acos(dot_prod / (norm1 * norm2))


        if abs(ang1) < abs(ang2):
            alignment[int(yy)][int(xx)] = ang1
            #if ang1 > pi/2 - 0.1 or ang1 < 0 + 0.15:
                #alignment[int(yy)][int(xx)] = ang1
            #else:
                #alignment[int(yy)][int(xx)] = pi/4
        else:
            alignment[int(yy)][int(xx)] = ang2
            #if ang2 > pi/2 - 0.1 or ang1< 0 + 0.15:
                #alignment[int(yy)][int(xx)] = ang2
            #else:
                #alignment[int(yy)][int(xx)] = pi/4


    z_min, z_max = 0., np.abs(alignment).max()
    #print(z_max)
    cset1 = plt.imshow(alignment, cmap='RdBu', interpolation='nearest', vmin=0, vmax=z_max)
    read_line+=1

vecfield_nx = Z_xx
vecfield_ny = Z_yy
vecfield_Q00 = Z_S00
vecfield_Q01 = Z_S01

winding_number = [[0. for j in range(0, LLX)] for i in range(0, LLY)]
for p in range(0, LLY):
    for q in range(0, LLX):
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

charge = 1.0/2.0
thresh = 0.05
for p in range(0,LLY):
    for q in range(0,LLX):
        # detect simplest charge 1/p defects
        if  (abs(winding_number[p][q]) > charge - thresh) and (abs(winding_number[p][q]) < charge + thresh):
            # charge sign
            s = np.sign(winding_number[p][q])
            # bfs
            sum_x, sum_y, n = collapse(p, q, s, LLX, LLY, winding_number, rng = [charge - thresh, charge + thresh])
            x,y = sum_x/n,sum_y/n
            # compute angle, see doi:10.1039/c6sm01146b
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
                cset1 = plt.plot(x, y, 'ro', markersize=10)
                cset1 = plt.arrow(x, y, 4*cos(psi), 4*sin(psi), color='r', head_width=1.5, head_length=1.5, width=0.5)
            elif s==-1:
                cset1 = plt.plot(x, y, 'ro', markersize=5)
                cset1 = plt.arrow(x, y, 4*cos(psi), 4*sin(psi), color='r', head_width=1.5, head_length=1.5, width=0.5)
                cset1 = plt.arrow(x, y, 4*cos(psi + 2.0944), 4*sin(psi + 2.0944), color='r', head_width=1.5, head_length=1.5, width=0.5)
                cset1 = plt.arrow(x, y, 4*cos(psi + 2 * 2.0944), 4*sin(psi + 2 * 2.0944), color='r', head_width=1.5, head_length=1.5, width=0.5)


        # keep this just in case our other symmetries give us integer defects
        elif (abs(winding_number[p][q]) > 1 - thresh) and (abs(winding_number[p][q])< 1 + thresh):
            # charge sign
            s = np.sign(winding_number[p][q])
            # bfs
            sum_x, sum_y, n = collapse(p, q, s, LLX, LLY, winding_number, rng = [1-thresh, 1+thresh])
            x,y = sum_x/n,sum_y/n
            # add defect to list
            if s==1:
                cset1 = plt.plot(x, y, 'r*', markersize=10)
            elif s==-1:
                cset1 = plt.plot(x, y, 'kX', markersize=10)



ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
ax.set_xlim([0, lx])
ax.set_ylim([0, ly])
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.set_xticks([])
ax.set_yticks([])
#fig.tight_layout()
if variable==1:
    plt.savefig('frame.png')
if variable==2:
    plt.show()
