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

if len(sys.argv)!=4:
    print(sys.argv[0]," [input] [variable] [start line]")
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
lastconf_file = 'last_conf.dat'

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
    elif words[0]=='lastconf_file':
        lastconf_file = words[2]
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

com_x_t = []
com_y_t = []
time_conf = []
start_value = 11


unrap_comx=[0. for i in range(0,N)]
unrap_comx_save_1=[0. for i in range(0,N)]
unrap_comx_save_2=[0. for i in range(0,N)]
unrap_comy=[0. for i in range(0,N)]

CoMX=[0. for i in range(0,N)]
CoMY=[0. for i in range(0,N)]
CoMX_old=[0. for i in range(0,N)]
CoMY_old=[0. for i in range(0,N)]
theta_nem=[0. for i in range(0,N)]
area=[0. for i in range(0,N)]
Q00=[0. for i in range(0,N)]
Q01=[0. for i in range(0,N)]
LsubX=[0 for i in range(0,N)]
LsubY=[0 for i in range(0,N)]
offsetX=[0 for i in range(0,N)]
offsetY=[0 for i in range(0,N)]
cornerSite=[0 for i in range(0,N)]
cornerSite_x=[0. for i in range(0,N)]
cornerSite_y=[0. for i in range(0,N)]

lfile=open(lastconf_file,"r")
header=lfile.readline().split()
t=int(header[2])
header=lfile.readline().split()
lx=int(float(header[2]))
ly=int(float(header[3]))
x=np.arange(0,lx,1)
y=np.arange(0,ly,1)

sizex_coarse=int(lx/(R))+1
sizey_coarse=int(ly/(R))+1

deltax_coarse=float(lx)/float(sizex_coarse)
deltay_coarse=float(ly)/float(sizey_coarse)

nematic_grid_coarse_00=[[0. for j in range(0, sizex_coarse)] for i in range(0,sizey_coarse)]
nematic_grid_coarse_01=[[0. for j in range(0, sizex_coarse)] for i in range(0,sizey_coarse)]
n_coarse=[[0. for j in range(0, sizex_coarse)] for i in range(0,sizey_coarse)]


timeavg_nematic_grid_coarse_00=[[0. for j in range(0, sizex_coarse)] for i in range(0,sizey_coarse)]
timeavg_nematic_grid_coarse_01=[[0. for j in range(0, sizex_coarse)] for i in range(0,sizey_coarse)]
timeavg_n_coarse=[[0. for j in range(0, sizex_coarse)] for i in range(0,sizey_coarse)]


lambda_wall+=2
lambda_wall=0


cont_line=0
vmax=0.1
start_line=(N+2)*int(float(sys.argv[3])) 
cfile=open(trajectory_file,"r")
for i in range(start_line):
    cfile.readline()


fig = plt.figure(figsize=(6,6))
for line in cfile:
    cont_line+=1
    words=line.split()

    if words[0]=='t':
        t=int(float(words[2]))
        time_conf.append(t)

        pt_num=0
        CoMX_old=CoMX
        CoMY_old=CoMY
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
        cornerSite_x=[0. for i in range(0,N)]
        cornerSite_y=[0. for i in range(0,N)]

    elif words[0]=='b':
        lx=int(float(words[2]))
        ly=int(float(words[3]))
        x=np.arange(0,lx,1)
        y=np.arange(0,ly,1)
        nematic_grid_coarse_00=[[0. for j in range(0, sizex_coarse)] for i in range(0,sizey_coarse)]
        nematic_grid_coarse_01=[[0. for j in range(0, sizex_coarse)] for i in range(0,sizey_coarse)]
        n_coarse=[[0. for q in range(0, sizex_coarse)] for k in range(0,sizey_coarse)]
        #walls = [0. for i in range(lx*ly)]
        #set_walls(lx,ly,walls)

    else:
        Z=[[0. for q in range(lx)] for k in range(ly)]
        LsubX[pt_num]=int(float(words[0]))
        LsubY[pt_num]=int(float(words[1]))
        CoMX[pt_num]=float(words[2])
        CoMY[pt_num]=float(words[3])
        offsetX[pt_num]=int(float(words[4]))
        offsetY[pt_num]=int(float(words[5]))
        cornerSite[pt_num]=int(float(words[6]))
        cornerSite_x[pt_num]=int(float(words[7]))
        cornerSite_y[pt_num]=int(float(words[8]))

        nemX=float(words[9])
        nemY=float(words[10])
        normNem = sqrt(nemX * nemX + nemY * nemY)
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

            if cont_line>N+2 and yy>=int(ceil(2*lambda_wall)/2) and yy<ly-int(ceil(2*lambda_wall)/2):
                n_coarse[int(yy/deltay_coarse)][int(xx/deltax_coarse)]+=value
                nematic_grid_coarse_00[int(yy/deltay_coarse)][int(xx/deltax_coarse)]+=value*Q00[pt_num]
                nematic_grid_coarse_01[int(yy/deltay_coarse)][int(xx/deltax_coarse)]+=value*Q01[pt_num]

        X, Y = np.meshgrid(x, y)
        step = 0.01
        m = np.amax(Z)
        #if m<0.000001:
            #continue

        levels = np.arange(0.0, m, step) + step

        if variable==1 or variable==2:
            if pt_num==0:
                cset1 = plt.contour(X, Y, Z, levels, cmap=cm.winter, alpha=0.5)
            else:
                cset1 = plt.contour(X, Y, Z, levels=[0.5], cmap=cm.winter, alpha=0.5)

            cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], 3*nemX, 3*nemY, width=0.5, head_width=0, color='k')

        #increment phase field index
        pt_num+=1
        

    if cont_line%(N+2)==0:


            for p in range(0,sizey_coarse):
                for q in range(0,sizex_coarse):
                    S = 0.
                    nx = 0.
                    ny = 0.
                    if n_coarse[p][q]>0:
                        Q00=nematic_grid_coarse_00[p][q]/n_coarse[p][q]
                        Q01=nematic_grid_coarse_01[p][q]/n_coarse[p][q]
                        S = sqrt(Q00**2 + Q01**2)
                        nx = sqrt((1 + Q00/S)/2)
                        ny = np.sign(Q01)*sqrt((1 - Q00/S)/2)

                        timeavg_nematic_grid_coarse_00[p][q]+=Q00
                        timeavg_nematic_grid_coarse_01[p][q]+=Q01
                        timeavg_n_coarse[p][q]+=1
                        
                    else:
                        Q00=0.
                        Q01=0.

                    if n_coarse[p][q]>0:
                        if variable==3 or variable==4:
                            cset1 = plt.arrow(q*deltax_coarse+deltax_coarse/2, p*deltay_coarse+deltay_coarse/2, 0.4*deltax_coarse*nx, 0.4*deltay_coarse*ny, width=deltax_coarse/15, color="k", head_width=0)
                            cset1 = plt.arrow(q*deltax_coarse+deltax_coarse/2, p*deltay_coarse+deltay_coarse/2, -0.4*deltax_coarse*nx, -0.4*deltay_coarse*ny, width=deltax_coarse/15, color="k", head_width=0)
                    else:
                        if variable==3 or variable==4:
                            cset1 = plt.arrow(q*deltax_coarse+deltax_coarse/2, p*deltay_coarse+deltay_coarse/2, 0.4*deltax_coarse*nx, 0.4*deltay_coarse*ny, width=deltax_coarse/15, color="k", head_width=0)
                            cset1 = plt.arrow(q*deltax_coarse+deltax_coarse/2, p*deltay_coarse+deltay_coarse/2, -0.4*deltax_coarse*nx, -0.4*deltay_coarse*ny, width=deltax_coarse/15, color="k", head_width=0)
                    

            frame_num=int(t/print_conf_interval)-1
            #if frame_num%1==0:
                #print(frame_num, cont_line, t)
            #if cont_line>N+2:
                #cset1 = plt.imshow(velocity_grid, vmin=velmin, vmax=velmax, cmap=cm.Reds)
            com_x_t.append(CoMX)
            com_y_t.append(CoMY)
            ax = plt.gca()
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlim([0, lx])
            ax.set_ylim([0, ly])
            if variable==1 or variable==3:
                if frame_num<10:
                    plt.savefig('./Video/frame_00'+str(frame_num)+'.png', transparent=True)
                elif frame_num<100:
                    plt.savefig('./Video/frame_0'+str(frame_num)+'.png')
                elif frame_num<1000:
                    plt.savefig('./Video/frame_'+str(frame_num)+'.png')
            if variable==2 or variable==4:
                plt.show()
                #plt.savefig('./newfig_'+str(frame_num)+'.png', transparent=True)
            if variable<=4:
                plt.clf()

plt.close()

y=np.arange(0., ly, deltay_coarse)
Q00_width=[0. for i in range(0, sizey_coarse)]
Q01_width=[0. for i in range(0, sizey_coarse)]
theta_width=[0. for i in range(0, sizey_coarse)]
avg_val=0
for p in range(0, sizey_coarse):
    avg_val=0
    for q in range(0, sizex_coarse):
        Q00_width[p]+=timeavg_nematic_grid_coarse_00[p][q]/float(timeavg_n_coarse[p][q])
        Q01_width[p]+=timeavg_nematic_grid_coarse_01[p][q]/float(timeavg_n_coarse[p][q])
        avg_val+=1

    Q00_width[p] = Q00_width[p]/float(avg_val)
    Q01_width[p] = Q01_width[p]/float(avg_val)
    S = sqrt(Q00_width[p]**2 + Q01_width[p]**2)
    nx = sqrt((1 + Q00_width[p]/S)/2)
    ny = np.sign(Q01_width[p])*sqrt((1 - Q00_width[p]/S)/2)
    theta_width[p]=asin((nx*ny)/0.5)/2

if variable==5:

    fig = plt.figure(figsize=(6,6))
    for p in range(0, sizey_coarse):
        for q in range(0, sizex_coarse):
            Q00=timeavg_nematic_grid_coarse_00[p][q]/float(timeavg_n_coarse[p][q])
            Q01=timeavg_nematic_grid_coarse_01[p][q]/float(timeavg_n_coarse[p][q])
            S = sqrt(Q00**2 + Q01**2)
            nx = sqrt((1 + Q00/S)/2)
            ny = np.sign(Q01)*sqrt((1 - Q00/S)/2)

            cset1 = plt.arrow(q*deltax_coarse+deltax_coarse/2, p*deltay_coarse+deltay_coarse/2, 0.4*deltax_coarse*nx, 0.4*deltay_coarse*ny, width=deltax_coarse/15, color="k", head_width=0)
            cset1 = plt.arrow(q*deltax_coarse+deltax_coarse/2, p*deltay_coarse+deltay_coarse/2, -0.4*deltax_coarse*nx, -0.4*deltay_coarse*ny, width=deltax_coarse/15, color="k", head_width=0)

    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim([0, lx])
    ax.set_ylim([0, ly])
    #plt.show()
    plt.savefig('./theta_avg_coarse.png')
    plt.close()

    #fig = plt.figure(figsize=(8,6))
    fig = plt.figure(figsize=(5.452423529, 4.089317647))
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    plt.plot(theta_width, y, '-o' , color='darkred')
    #plt.xlabel('Channel width', fontsize=18, fontname='Times New Roman')
    #plt.ylabel(r'Velocity $v_y$', fontsize=18, fontname='Times New Roman')
    #plt.xticks(fontsize=18, fontname='Times New Roman')
    #plt.yticks(fontsize=18, fontname='Times New Roman')
    plt.ylabel('Channel width', fontsize=18)
    plt.xlabel(r'$\theta$', fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    #plt.xlim(velmin_x,velmax_x)
    #plt.xlim(-3*1e-5,2.5*1e-5)
    plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
    plt.locator_params(axis='x', nbins=6)
    #plt.show()
    plt.savefig('./theta_width_coarse.png')
    plt.close()

if variable==6:

    '''
    with open('MSD.txt', 'w') as f:
        for i in range(len(MSD)):
            print(time_conf[i]*dt,MSD[i], file=f)  

    with open('mean_velocity.txt', 'w') as f:
        print(abs(avg_mean_velocity), file=f)  

    y=np.arange(int(ceil(2*lambda_wall)/2), ly-int(ceil(2*lambda_wall)/2), 1)
    with open('v_width.txt', 'w') as f:
        for i in range(len(y)):
            print(avg_velocity_y[i]/counter_for_avg[i], avg_velocity_x[i]/counter_for_avg[i], y[i], file=f)  
    '''

print('done')
