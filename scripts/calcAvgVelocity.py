import sys
import numpy as np
from math import *
import os.path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgba, to_hex, Normalize
import numpy as np
import scipy.ndimage

from matplotlib import cm

if len(sys.argv)!=3:
    print(sys.argv[0]," [input] [variable]")
    sys.exit(1)


variable=int(float(sys.argv[2])) 

def set_walls(lx,ly, walls):
    for y in range(ly):
        for x in range(lx):
            k = x + y * lx
            walls[k] = exp(-float(y)/5) + exp(-float(ly-y-1)/5)

seed=4982
eq_steps=0
steps = 1000000
print_conf_interval = 2000
dt = 1
R = 8
J0 = 0.005
gamma = 1.4
mu = 120
llambda = 2.0
kappa = 1.5
friction = 3
omega = 0.4
friction_cell = 3
J_Q = 0.1
zetaQ_self = 0.5
zetaQ_inter = 0
anchoring = 0
lambda_wall = 3

external_forces_file = 'external.conf '
topology = 'test.top'
conf_file = 'start.conf'
trajectory_file = 'trajectory.dat'

ifile=open(sys.argv[1],"r")
for line in ifile:
    words=line.split()
    if len(words)==0:
        continue
    if words[0]=='topology':
        topology = words[2]
    elif words[0]=='conf_file':
        conf_file = words[2]
    elif words[0]=='trajectory_file':
        trajectory_file = words[2]
    elif words[0]=='seed':
        seed=int(float(words[2]))
    elif words[0]=='equilibration_steps':
        eq_steps=int(float(words[2]))
    elif words[0]=='print_conf_interval':
        print_conf_interval=int(float(words[2]))
    elif words[0]=='dt':
        dt=float(words[2])
    elif words[0]=='R':
        R=int(float(words[2]))
    elif words[0]=='J0':
        J0=float(words[2])
    elif words[0]=='gamma':
        gamma=float(words[2])
    elif words[0]=='mu':
        mu=float(words[2])
    elif words[0]=='lambda':
        llambda=float(words[2])
    elif words[0]=='kappa':
        kappa=float(words[2])
    elif words[0]=='friction':
        friction=float(words[2])
    elif words[0]=='omega':
        omega=float(words[2])
    elif words[0]=='friction_cell':
        friction_cell=float(words[2])
    elif words[0]=='J_Q':
        J_Q=float(words[2])
    elif words[0]=='zetaQ_self':
        zetaQ_self=float(words[2])
    elif words[0]=='zetaQ_inter':
        zetaQ_inter=float(words[2])
    elif words[0]=='anchoring':
        anchoring=int(float(words[2]))
    elif words[0]=='lambda_wall':
        lambda_wall=float(words[2])


tfile=open(topology,"r")
N=int(tfile.readline().split()[0])
species=[]
line=tfile.readline()
lline=line.split()
for l in lline:
        species.append(int(l))
tfile.close()

numspecies=len(set(species))


cfile=open(trajectory_file,"r")
v_com_x_t = []
v_com_y_t = []
com_x_t = []
com_y_t = []
time_conf = []
start_value = 9

CoMX=[0. for i in range(0,N)]
CoMY=[0. for i in range(0,N)]
theta_nem=[0. for i in range(0,N)]
area=[0. for i in range(0,N)]
Q00=[0. for i in range(0,N)]
Q01=[0. for i in range(0,N)]
LsubX=[0 for i in range(0,N)]
LsubY=[0 for i in range(0,N)]
offsetX=[0 for i in range(0,N)]
offsetY=[0 for i in range(0,N)]
cornerSite=[0 for i in range(0,N)]
com_velocity_x = [0. for i in range(0,N)]
com_velocity_y = [0. for i in range(0,N)]
lx=0
ly=0
x=np.arange(0,lx,1)
y=np.arange(0,ly,1)

cont_line=0
vmax=0.005

fig = plt.figure(figsize=(6,6))
for line in cfile:
    cont_line+=1
    words=line.split()

    if words[0]=='t':
        t=int(float(words[2]))
        time_conf.append(t)

        pt_num=0
        CoMX=[0. for i in range(0,N)]
        CoMY=[0. for i in range(0,N)]
        theta_nem=[0. for i in range(0,N)]
        area=[0. for i in range(0,N)]
        Q00=[0. for i in range(0,N)]
        Q01=[0. for i in range(0,N)]
        LsubX=[0 for i in range(0,N)]
        LsubY=[0 for i in range(0,N)]
        offsetX=[0 for i in range(0,N)]
        offsetY=[0 for i in range(0,N)]
        cornerSite=[0 for i in range(0,N)]
        com_velocity_x = [0. for i in range(0,N)]
        com_velocity_x = [0. for i in range(0,N)]

    elif words[0]=='b':
        lx=int(float(words[2]))
        ly=int(float(words[3]))
        x=np.arange(0,lx,1)
        y=np.arange(0,ly,1)
        #walls = [0. for i in range(lx*ly)]
        #set_walls(lx,ly,walls)

    else:
        Z=[[0 for q in range(lx)] for k in range(ly)]
        LsubX[pt_num]=int(float(words[0]))
        LsubY[pt_num]=int(float(words[1]))
        CoMX[pt_num]=float(words[2])
        CoMY[pt_num]=float(words[3])
        offsetX[pt_num]=int(float(words[4]))
        offsetY[pt_num]=int(float(words[5]))
        cornerSite[pt_num]=int(float(words[6]))

        if t==print_conf_interval:
            com_velocity_x[pt_num]=0.
            com_velocity_y[pt_num]=0.
        else:
            t1=len(time_conf)-1
            t2=len(time_conf)-2
            com_velocity_x[pt_num]=(CoMX[pt_num]-com_x_t[t2][pt_num])/((time_conf[t1]-time_conf[t2])*dt)
            com_velocity_y[pt_num]=(CoMY[pt_num]-com_y_t[t2][pt_num])/((time_conf[t1]-time_conf[t2])*dt)
            #print(com_velocity_x[pt_num], com_velocity_y[pt_num])

        nemX=float(words[7])
        nemY=float(words[8])
        theta_nem[pt_num]=asin((nemX*nemY)/0.5)/2
        Q00[pt_num]= 0.5 * (nemX * nemX - nemY * nemY)
        Q01[pt_num]= nemX * nemY

        for i in range(start_value,len(words),2):
            site=int(float(words[i]))
            value=float(words[i+1])
            yy=int(site/lx)
            xx=site-int(yy*lx)
            Z[yy][xx]=value
            area[pt_num]+=value*value

        X, Y = np.meshgrid(x, y)
        if pt_num<0:
            cset1 = plt.contour(X, Y, Z, levels, cmap=cm.winter, alpha=0.5)
        else:
            cset1 = plt.contour(X, Y, Z, levels=[0.5], cmap=cm.winter, alpha=0.5)

        norm=sqrt(com_velocity_x[pt_num]*com_velocity_x[pt_num]+com_velocity_y[pt_num]*com_velocity_y[pt_num])
        color_val=norm/vmax
        if norm==0:
            norm=1
        cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], com_velocity_x[pt_num]/norm, com_velocity_y[pt_num]/norm, width=0.5, color=cm.hot(color_val))

        #print(pt_num, area)
        pt_num+=1


    if cont_line%(N+2)==0:
            com_x_t.append(CoMX)
            com_y_t.append(CoMY)
            v_com_x_t.append(com_velocity_x)
            v_com_y_t.append(com_velocity_y)
            ax = plt.gca()
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlim([0, lx])
            ax.set_ylim([0, ly])
            frame_num=int(t/print_conf_interval)-1
            print(frame_num)
            if variable==1:
                plt.savefig('./Video_velocity/frame_'+str(frame_num)+'.png')
            if variable==2:
                plt.show()

            if t==3*print_conf_interval:
                exit(1)


print('done')
exit (1)


def getElementX(site, distX, sidex, sidey):
    x = int( (site-(int(site/sidex)*sidex))+distX )
    if x<0:
        x+=sidex
    if x>=sidex:
        x-=sidex
    return x

def getElementY(site, distY, sidex, sidey):
    y = int( int(site/sidex)+distY )
    if y<0:
        y+=sidey
    if y>=sidey:
        y-=sidey
    return y

def getElement(site, LsubX, LsubY, offsetX, offsetY, sub_corner_bottom_left, boxXsize, boxYsize):
    return getElementX(sub_corner_bottom_left, ((site - int(site/LsubX) * LsubX)+offsetX)%LsubX, boxXsize, boxYsize ) + getElementY(sub_corner_bottom_left, ((site/LsubX)+offsetY)%LsubY , boxXsize, boxYsize) * boxXsize

def getX(site, LsubX, LsubY, offsetX, offsetY, sub_corner_bottom_left, boxXsize, boxYsize):
    return getElementX(sub_corner_bottom_left, ((site - int(site/LsubX) * LsubX)+offsetX)%LsubX , boxXsize, boxYsize)

def getY(site, LsubX, LsubY, offsetX, offsetY, sub_corner_bottom_left, boxXsize, boxYsize):
    return getElementY(sub_corner_bottom_left, ((site/LsubX)+offsetY)%LsubY , boxXsize, boxYsize)


all_phi=row_phi
for i in range(len(row_phi)):
    for j in range(LsubY[i]):
        for k in range(LsubX[i]):
            row_site=k+j*LsubX[i]
            column_site=j+LsubY[i]*k
            all_phi[i][row_site]=row_phi[i][row_site]
            if j==1 and i==0:
                print("this: ", k,j,row_phi[i][row_site])
            if i==0 and getY(row_site, LsubX[i], LsubY[i], offsetX[i], offsetY[i], cornerSite[i], lx, ly)==40:
                print( getX(row_site, LsubX[i], LsubY[i], offsetX[i], offsetY[i], cornerSite[i], lx, ly), getY(row_site, LsubX[i], LsubY[i], offsetX[i], offsetY[i], cornerSite[i], lx, ly), getElement(row_site, LsubX[i], LsubY[i], offsetX[i], offsetY[i], cornerSite[i], lx, ly), row_phi[i][row_site])


print("end")
