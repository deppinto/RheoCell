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

if len(sys.argv)!=4:
    print(sys.argv[0]," [input] [variable] [start line]")
    sys.exit(1)


variable=int(float(sys.argv[2])) 

def set_walls(lx,ly, walls):
    for y in range(ly):
        for x in range(lx):
            k = x + y * lx
            walls[k] = exp(-float(y)/5) + exp(-float(ly-y-1)/5)

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

v_com_x_t = []
v_com_y_t = []
com_x_t = []
com_y_t = []
MSD = []
time_conf = []
start_value = 11
avg_mean_velocity=0.

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
com_velocity_x = [0. for i in range(0,N)]
com_velocity_y = [0. for i in range(0,N)]

lfile=open(lastconf_file,"r")
header=lfile.readline().split()
t=int(header[2])
header=lfile.readline().split()
lx=int(float(header[2]))
ly=int(float(header[3]))
x=np.arange(0,lx,1)
y=np.arange(0,ly,1)
velocity_grid=[[0. for q in range(lx)] for k in range(ly)]
velocity_grid_x=[[0. for q in range(lx)] for k in range(ly)]
velocity_grid_y=[[0. for q in range(lx)] for k in range(ly)]
sum_phi=[[0. for q in range(lx)] for k in range(ly)]



cont_line=0
start_line=(N+2)*int(float(sys.argv[3])) 
cfile=open(trajectory_file,"r")
for i in range(start_line):
    cfile.readline()


fig = plt.figure(figsize=(6,6))

com_avg_velocity_calc = 0.
com_all_x = 0.
com_all_y = 0.
Gamma_rot = []
for line in cfile:
    cont_line+=1
    words=line.split()

    if words[0]=='t':
        t=int(float(words[2]))
        time_conf.append(t)
        com_all_x = 0.
        com_all_y = 0.
        MSD.append(0.)

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
        com_velocity_x = [0. for i in range(0,N)]
        com_velocity_y = [0. for i in range(0,N)]

        #print(com_avg_velocity_calc/N)
        com_avg_velocity_calc = 0.

    elif words[0]=='b':
        lx=int(float(words[2]))
        ly=int(float(words[3]))
        x=np.arange(0,lx,1)
        y=np.arange(0,ly,1)

    else:
        Z=[[0. for q in range(lx)] for k in range(ly)]
        #velocity_grid=[[0. for q in range(lx)] for k in range(ly)]
        LsubX[pt_num]=int(float(words[0]))
        LsubY[pt_num]=int(float(words[1]))
        CoMX[pt_num]=float(words[2])
        CoMY[pt_num]=float(words[3])
        offsetX[pt_num]=int(float(words[4]))
        offsetY[pt_num]=int(float(words[5]))
        cornerSite[pt_num]=int(float(words[6]))
        cornerSite_x[pt_num]=int(float(words[7]))
        cornerSite_y[pt_num]=int(float(words[8]))

        if cont_line<=N+2:
            com_velocity_x[pt_num]=0.
            com_velocity_y[pt_num]=0.
        else:
            t1=len(time_conf)-1
            t2=len(time_conf)-2
            dist_com_x=(CoMX[pt_num]-CoMX_old[pt_num])
            if dist_com_x>lx/2:
                dist_com_x-=lx
            if dist_com_x<-lx/2:
                dist_com_x+=lx
            dist_com_y=(CoMY[pt_num]-CoMY_old[pt_num])
            if dist_com_y>ly/2:
                dist_com_y-=ly
            if dist_com_y<-ly/2:
                dist_com_y+=ly
            unrap_comx[pt_num]+=dist_com_x

            #if t==10000:
                #unrap_comx_save_1[pt_num]=unrap_comx[pt_num]
            #if t==30000:
                #unrap_comx_save_2[pt_num]=unrap_comx[pt_num]
                #avg_mean_velocity+=(unrap_comx_save_2[pt_num]-unrap_comx_save_1[pt_num])/(20000*dt*N)

            unrap_comy[pt_num]+=dist_com_y
            MSD[len(MSD)-1]+=(unrap_comx[pt_num]*unrap_comx[pt_num]+unrap_comy[pt_num]*unrap_comy[pt_num])/N
            com_velocity_x[pt_num]=dist_com_x/((time_conf[t1]-time_conf[t2])*dt)
            com_velocity_y[pt_num]=dist_com_y/((time_conf[t1]-time_conf[t2])*dt)
            com_avg_velocity_calc += sqrt(com_velocity_x[pt_num]*com_velocity_x[pt_num] + com_velocity_y[pt_num]*com_velocity_y[pt_num])

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

            cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], com_velocity_x[pt_num]/norm, com_velocity_y[pt_num]/norm, width=0.5, color="k")
            #cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], com_velocity_x[pt_num]/norm, com_velocity_y[pt_num]/norm, width=0.5, color=cm.hot(color_val))

        com_all_x += CoMX[pt_num]
        com_all_y += CoMY[pt_num]
        #increment phase field index
        pt_num+=1
        

    if cont_line%(N+2)==0:

            com_all_x = com_all_x / N
            com_all_y = com_all_y / N
            vector_sum = 0.
            for p in range(N):
                distx = CoMX[p] - com_all_x
                if distx>lx/2:
                    distx-=lx
                if distx<-lx/2:
                    distx+=lx
                disty = CoMY[p] - com_all_y
                if disty>ly/2:
                    disty-=ly
                if disty<-ly/2:
                    disty+=ly

                dist_com_x=(CoMX[p]-CoMX_old[p])
                if dist_com_x>lx/2:
                    dist_com_x-=lx
                if dist_com_x<-lx/2:
                    dist_com_x+=lx
                dist_com_y=(CoMY[p]-CoMY_old[p])
                if dist_com_y>ly/2:
                    dist_com_y-=ly
                if dist_com_y<-ly/2:
                    dist_com_y+=ly

                vector_x = distx * dist_com_y
                vector_y = - disty * dist_com_x
                norm = sqrt(vector_x**2 + vector_y**2)
                if norm>0:
                    vector_x = vector_x / norm
                    vector_y = vector_y / norm
                else:
                    vector_x = 0.
                    vector_y = 0.
                vector_sum += vector_x

            Gamma_rot.append(vector_sum/N) 

            frame_num=int(t/print_conf_interval)-1
            #print(frame_num)
            #if frame_num%1==0:
                #print(frame_num, cont_line, t)
            #if cont_line>N+2:
                #cset1 = plt.imshow(velocity_grid, vmin=velmin, vmax=velmax, cmap=cm.Reds)
            com_x_t.append(CoMX)
            com_y_t.append(CoMY)
            v_com_x_t.append(com_velocity_x)
            v_com_y_t.append(com_velocity_y)
            ax = plt.gca()
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlim([0, lx])
            ax.set_ylim([0, ly])
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            ax.set_xticks([])
            ax.set_yticks([])
            if variable==1 or variable==3:
                if frame_num<10:
                    plt.savefig('./Video/frame_00'+str(frame_num)+'.png', transparent=True)
                elif frame_num<100:
                    plt.savefig('./Video/frame_0'+str(frame_num)+'.png')
                elif frame_num<1000:
                    plt.savefig('./Video/frame_'+str(frame_num)+'.png')
                else:
                    plt.savefig('./Video/frame_'+str(frame_num)+'.png')
            if variable==2 or variable==4:
                plt.show()
                #plt.savefig('./newfig_'+str(frame_num)+'.png', transparent=True)

            if variable<=4:
                plt.clf()

plt.close()


if variable==5:

    '''
    plot_vel=[]
    for i in range(len(time_conf)):
        norm=sqrt(v_com_x_t[i][0]*v_com_x_t[i][0])/0.01
        plot_vel.append(norm)
    fig = plt.figure(figsize=(5.452423529, 4.089317647))
    plt.plot(time_conf, plot_vel, '-o' , color='darkgreen')
    plt.ylabel(r'$v_a/v_0$', fontsize=18)
    plt.xlabel('time', fontsize=18)
    plt.xlim(500,7000)
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
    plt.show()
    plt.close()
    exit(1)
    '''

    vx_max = 0.3
    vy_max = 0.3
    fig = plt.figure(figsize=(6,6))
    for p in range(0, sizey_coarse):
        for q in range(0, sizex_coarse):
            vx=timeavg_velocity_grid_coarse_x[p][q]/float(timeavg_n_coarse[p][q])
            vy=timeavg_velocity_grid_coarse_y[p][q]/float(timeavg_n_coarse[p][q])

            cset1 = plt.arrow(q*deltax_coarse+deltax_coarse/2, p*deltay_coarse+deltay_coarse/2, 0.8*deltax_coarse*vx/(vx_max), 0.8*deltay_coarse*vy/(vy_max), width=deltax_coarse/15, color="k")

    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim([0, lx])
    ax.set_ylim([0, ly])
    #plt.show()
    plt.savefig('./v_avg_coarse.png')
    plt.close()

    #fig = plt.figure(figsize=(8,6))
    fig = plt.figure(figsize=(5.452423529, 4.089317647))
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    plt.plot(vx_width, y, '-o' , color='darkred')
    #plt.xlabel('Channel width', fontsize=18, fontname='Times New Roman')
    #plt.ylabel(r'Velocity $v_y$', fontsize=18, fontname='Times New Roman')
    #plt.xticks(fontsize=18, fontname='Times New Roman')
    #plt.yticks(fontsize=18, fontname='Times New Roman')
    plt.ylabel('Channel width', fontsize=18)
    plt.xlabel(r'Velocity $v_x$', fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    #plt.xlim(velmin_x,velmax_x)
    #plt.xlim(-3*1e-5,2.5*1e-5)
    plt.locator_params(axis='x', nbins=6)
    plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
    plt.savefig('./vx_width_coarse.png')
    plt.show()
    #plt.savefig('./vx_width.png')
    plt.clf()

    #fig = plt.figure(figsize=(8,6))
    #fig = plt.figure(figsize=(5.452423529, 4.089317647))
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    plt.plot(vy_width, y, '-o' , color='darkgreen')
    #plt.xlabel('Channel width', fontsize=18, fontname='Times New Roman')
    #plt.ylabel(r'Velocity $v_x$', fontsize=18, fontname='Times New Roman')
    #plt.xticks(fontsize=18, fontname='Times New Roman')
    #plt.yticks(fontsize=18, fontname='Times New Roman')
    plt.ylabel('Channel width', fontsize=18)
    plt.xlabel(r'Velocity $v_y$', fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    #plt.xlim(velmin_x,velmax_x)
    #plt.xlim(-3*1e-5,2.5*1e-5)
    plt.locator_params(axis='x', nbins=6)
    plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
    #plt.show()
    plt.savefig('./vy_width_coarse.png')
    plt.close()


if variable==6:

    #fig = plt.figure(figsize=(8,6))
    fig = plt.figure(figsize=(5.452423529, 4.089317647))
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    plt.plot(time_conf, MSD, '-o' , color='darkgreen')
    #plt.plot(Gamma_rot, '-o' , color='darkgreen')
    #plt.plot(averages, '-o' , color='royalblue')
    #plt.ylabel(r'$\Gamma$', fontsize=18)
    plt.ylabel(r'MSD', fontsize=18)
    plt.xlabel('time', fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xscale('log')
    plt.yscale('log')
    #plt.ylim(-1,1)
    #plt.axhline(y=-0.5, color='green', linestyle='--', linewidth=1)
    #plt.axhline(y=0.5, color='green', linestyle='--', linewidth=1)
    plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
    #plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/Slides/Results11/gamma_time_circle_Velocity_NN.svg", transparent=True)
    #plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/Slides/Results11/gamma_time_circle_Velocity_NN.png", transparent=True)
    plt.show()
    #plt.savefig('./MSD_time.png')
    plt.close()

    MSD = np.array(MSD)
    time = np.array(time_conf)
    print(len(MSD))

    # 1. Slice the last 10 points
    time_last_10 = time[-15:]
    MSD_last_10 = MSD[-15:]

    # 2. Transform the sliced data into log space (Base e or Base 10 both work)
    log_time = np.log(time_last_10)
    log_MSD = np.log(MSD_last_10)

    # 3. Fit a 1st-degree polynomial (linear regression: y = mx + c)
    # The first returned value [0] is the slope, which is your power-law exponent.
    slope, intercept = np.polyfit(log_time, log_MSD, 1)
    print(f"Power-law exponent (slope): {slope:.4f}")



    '''
    with open('MSD.txt', 'w') as f:
        for i in range(len(MSD)):
            print(time_conf[i]*dt,MSD[i], file=f)  
    '''


    '''
    with open('mean_velocity.txt', 'w') as f:
        print(abs(avg_mean_velocity), file=f)  
    '''


    '''
    with open('v_width.txt', 'w') as f:
        for i in range(len(y)):
            print(vy_width[i], vx_width[i], y[i], file=f)  
    '''


print('done')
