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


sizex_coarse=int(lx/(R))+1
sizey_coarse=int(ly/(R))+1

deltax_coarse=float(lx)/float(sizex_coarse)
deltay_coarse=float(ly)/float(sizey_coarse)

velocity_grid_coarse_x=[[0. for j in range(0, sizex_coarse)] for i in range(0,sizey_coarse)]
velocity_grid_coarse_y=[[0. for j in range(0, sizex_coarse)] for i in range(0,sizey_coarse)]
n_coarse=[[0. for j in range(0, sizex_coarse)] for i in range(0,sizey_coarse)]

timeavg_velocity_grid_coarse_x=[[0. for j in range(0, sizex_coarse)] for i in range(0,sizey_coarse)]
timeavg_velocity_grid_coarse_y=[[0. for j in range(0, sizex_coarse)] for i in range(0,sizey_coarse)]
timeavg_n_coarse=[[0. for j in range(0, sizex_coarse)] for i in range(0,sizey_coarse)]

lambda_wall+=2
lambda_wall=0
avg_velocity_x=np.zeros(ly-ceil(2*lambda_wall))
avg_velocity_y=np.zeros(ly-ceil(2*lambda_wall))
counter_for_avg=np.zeros(ly-ceil(2*lambda_wall))

cont_line=0
vmax=0.1
start_line=(N+2)*int(float(sys.argv[3])) 
cfile=open(trajectory_file,"r")
for i in range(start_line):
    cfile.readline()


fig = plt.figure(figsize=(6,6))
velmax=0.
velmin=10000.

velmax_x=0.
velmin_x=10000.
for line in cfile:
    cont_line+=1
    words=line.split()

    if words[0]=='t':
        t=int(float(words[2]))
        time_conf.append(t)
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
        com_velocity_x = [0. for i in range(0,N)]

    elif words[0]=='b':
        lx=int(float(words[2]))
        ly=int(float(words[3]))
        x=np.arange(0,lx,1)
        y=np.arange(0,ly,1)
        velocity_grid=[[0. for q in range(lx)] for k in range(ly)]
        velocity_grid_x=[[0. for q in range(lx)] for k in range(ly)]
        velocity_grid_y=[[0. for q in range(lx)] for k in range(ly)]
        velocity_grid_coarse_x=[[0. for q in range(0, sizex_coarse)] for k in range(0,sizey_coarse)]
        velocity_grid_coarse_y=[[0. for q in range(0, sizex_coarse)] for k in range(0,sizey_coarse)]
        n_coarse=[[0. for q in range(0, sizex_coarse)] for k in range(0,sizey_coarse)]
        sum_phi=[[0. for q in range(lx)] for k in range(ly)]
        velmax=0.
        velmin=10000.
        #walls = [0. for i in range(lx*ly)]
        #set_walls(lx,ly,walls)

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
            if t==10000:
                unrap_comx_save_1[pt_num]=unrap_comx[pt_num]
            if t==30000:
                unrap_comx_save_2[pt_num]=unrap_comx[pt_num]
                avg_mean_velocity+=(unrap_comx_save_2[pt_num]-unrap_comx_save_1[pt_num])/(20000*dt*N)

            unrap_comy[pt_num]+=dist_com_y
            MSD[len(MSD)-1]+=sqrt(unrap_comx[pt_num]*unrap_comx[pt_num]+unrap_comy[pt_num]*unrap_comy[pt_num])/N
            com_velocity_x[pt_num]=dist_com_x/((time_conf[t1]-time_conf[t2])*dt)
            com_velocity_y[pt_num]=dist_com_y/((time_conf[t1]-time_conf[t2])*dt)

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
                velocity_grid[yy][xx]=value*sqrt(com_velocity_x[pt_num]*com_velocity_x[pt_num]+com_velocity_y[pt_num]*com_velocity_y[pt_num])
                velocity_grid_x[yy][xx]+=value*com_velocity_x[pt_num]
                velocity_grid_y[yy][xx]+=value*com_velocity_y[pt_num]
                sum_phi[yy][xx]+=value
                velocity_grid_coarse_x[int(yy/deltay_coarse)][int(xx/deltax_coarse)]+=value*com_velocity_x[pt_num]
                velocity_grid_coarse_y[int(yy/deltay_coarse)][int(xx/deltax_coarse)]+=value*com_velocity_y[pt_num]
                n_coarse[int(yy/deltay_coarse)][int(xx/deltax_coarse)]+=value
                avg_velocity_x[yy-int(ceil(2*lambda_wall)/2)]+=value*com_velocity_x[pt_num]
                avg_velocity_y[yy-int(ceil(2*lambda_wall)/2)]+=value*com_velocity_y[pt_num]

                if value*com_velocity_x[pt_num]<velmin_x:
                    velmin_x=value*com_velocity_x[pt_num]
                if value*com_velocity_y[pt_num]<velmin_x:
                    velmin_x=value*com_velocity_y[pt_num]

                if value*com_velocity_x[pt_num]>velmax_x:
                    velmax_x=value*com_velocity_x[pt_num]
                if value*com_velocity_y[pt_num]>velmax_x:
                    velmax_x=value*com_velocity_y[pt_num]

                counter_for_avg[yy-int(ceil(2*lambda_wall)/2)]+=1


        X, Y = np.meshgrid(x, y)
        norm=sqrt(com_velocity_x[pt_num]*com_velocity_x[pt_num]+com_velocity_y[pt_num]*com_velocity_y[pt_num])
        color_val=norm/vmax
        if norm>velmax:
            velmax=norm
        if norm<velmin:
            velmin=norm
        if norm==0:
            norm=1

        step = 0.01
        m = np.amax(Z)
        #if m<0.000001:
            #continue

        levels = np.arange(0.0, m, step) + step

        if variable==1 or variable==2:
            if pt_num==1:
                cset1 = plt.contour(X, Y, Z, levels, cmap=cm.winter, alpha=0.5)
            else:
                cset1 = plt.contour(X, Y, Z, levels=[0.5], cmap=cm.winter, alpha=0.5)

            cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], com_velocity_x[pt_num]/norm, com_velocity_y[pt_num]/norm, width=0.5, color="k")
            #cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], com_velocity_x[pt_num]/norm, com_velocity_y[pt_num]/norm, width=0.5, color=cm.hot(color_val))

        #increment phase field index
        pt_num+=1
        

    if cont_line%(N+2)==0:

            #vx_max = max(map(max, velocity_grid_coarse_x))
            #vy_max = max(map(max, velocity_grid_coarse_y))
            vx_max = 0.2
            vy_max = 0.2

            vorticity=[[0. for q in range(0, sizex_coarse)] for k in range(0,sizey_coarse)]
            Q_criterion=[[0. for q in range(0, sizex_coarse)] for k in range(0,sizey_coarse)]

            for p in range(0,sizey_coarse):
                ynext = p + 1
                yprev = p - 1
                if ynext>=sizey_coarse:
                    ynext-=sizey_coarse
                if yprev<0:
                    yprev+=sizey_coarse

                for q in range(0,sizex_coarse):
                    if n_coarse[p][q]>0:
                        timeavg_velocity_grid_coarse_x[p][q]+=velocity_grid_coarse_x[p][q]/n_coarse[p][q]
                        timeavg_velocity_grid_coarse_y[p][q]+=velocity_grid_coarse_y[p][q]/n_coarse[p][q]
                        timeavg_n_coarse[p][q]+=1

                    xnext = q + 1
                    xprev = q - 1
                    if xnext>=sizex_coarse:
                        xnext-=sizex_coarse
                    if xprev<0:
                        xprev+=sizex_coarse

                    if n_coarse[ynext][q]>0:
                        vx = velocity_grid_coarse_x[ynext][q]/n_coarse[ynext][q]
                        vy = velocity_grid_coarse_y[ynext][q]/n_coarse[ynext][q]
                        norm = sqrt(vx*vx+vy*vy)
                        vx_ynext = vx
                        vy_ynext = vy
                    else:
                        vx_ynext = 0
                        vy_ynext = 0

                    if n_coarse[yprev][q]>0:
                        vx = velocity_grid_coarse_x[yprev][q]/n_coarse[yprev][q]
                        vy = velocity_grid_coarse_y[yprev][q]/n_coarse[yprev][q]
                        norm = sqrt(vx*vx+vy*vy)
                        vx_yprev = vx
                        vy_yprev = vy
                    else:
                        vx_yprev = 0
                        vy_yprev = 0

                    if n_coarse[p][xnext]>0:
                        vx = velocity_grid_coarse_x[p][xnext]/n_coarse[p][xnext]
                        vy = velocity_grid_coarse_y[p][xnext]/n_coarse[p][xnext]
                        norm = sqrt(vx*vx+vy*vy)
                        vx_xnext = vx
                        vy_xnext = vy
                    else:
                        vx_xnext = 0
                        vy_xnext = 0

                    if n_coarse[p][xprev]>0:
                        vx = velocity_grid_coarse_x[p][xprev]/n_coarse[p][xprev]
                        vy = velocity_grid_coarse_y[p][xprev]/n_coarse[p][xprev]
                        norm = sqrt(vx*vx+vy*vy)
                        vx_xprev = vx
                        vy_xprev = vy
                    else:
                        vx_xprev = 0
                        vy_xprev = 0

                    dvxdx = (vx_xnext-vx_xprev)/(2*deltax_coarse)
                    dvydx = (vy_xnext-vy_xprev)/(2*deltax_coarse)
                    dvxdy = (vx_ynext-vx_yprev)/(2*deltay_coarse)
                    dvydy = (vy_ynext-vy_yprev)/(2*deltay_coarse)

                    vorticity[p][q] = dvydx - dvxdy
                    Q_criterion[p][q] = dvxdx * dvydy - dvxdy * dvydx

                    '''
                    if Q_criterion[p][q]>0:
                        Q_criterion[p][q]=1
                    elif Q_criterion[p][q]<0:
                        Q_criterion[p][q]=-1
                    else:
                        Q_criterion[p][q]=0
                    '''

                    '''
                    if abs(Q_criterion[p][q])>2e-5:
                       print("max: ",Q_criterion[p][q])
                    '''

                    if n_coarse[p][q]>0:
                        if variable==3 or variable==4:
                            cset1 = plt.arrow(q*deltax_coarse+deltax_coarse/2, p*deltay_coarse+deltay_coarse/2, 0.8*deltax_coarse*velocity_grid_coarse_x[p][q]/(n_coarse[p][q]*vx_max), 0.8*deltay_coarse*velocity_grid_coarse_y[p][q]/(n_coarse[p][q]*vy_max), width=deltax_coarse/15, color="k")
                    else:
                        if variable==3 or variable==4:
                            cset1 = plt.arrow(q*deltax_coarse+deltax_coarse/2, p*deltay_coarse+deltay_coarse/2, 0.8*deltax_coarse*velocity_grid_coarse_x[p][q], 0.8*deltay_coarse*velocity_grid_coarse_y[p][q], width=deltax_coarse/15, color="k")
                    

            #X, Y = np.meshgrid( np.arange(0+deltax_coarse/2, lx+deltax_coarse/2, deltax_coarse) , np.arange(0+deltay_coarse/2, ly+deltay_coarse/2, deltay_coarse) )
            #z_min, z_max = -np.abs(vorticity).max(), np.abs(vorticity).max()
            z_min, z_max = -np.abs(Q_criterion).max(), np.abs(Q_criterion).max()
            #print(z_min, z_max)

            #z_min = -0.015
            #z_max = 0.015 
            #z_min = -0.009
            #z_max = 0.009 
            #cset1 = plt.pcolormesh(X, Y, vorticity, cmap='RdBu', vmin=z_min, vmax=z_max, alpha=0.25)
            #cset1 = plt.imshow(vorticity, cmap='RdBu', vmin=z_min, vmax=z_max, alpha=0.6, interpolation='lanczos', extent=[0,lx,0,ly])

            #z_min = -3.5e-5
            #z_max = 3.5e-5
            #z_min = -2e-5
            #z_max = 2e-5
            if variable==3 or variable==4:
                cset1 = plt.imshow(Q_criterion, cmap='RdBu', vmin=z_min, vmax=z_max, alpha=0.6, interpolation='lanczos', extent=[0,lx,0,ly])

            frame_num=int(t/print_conf_interval)-1
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
vx_width=[0. for i in range(0, sizey_coarse)]
vy_width=[0. for i in range(0, sizey_coarse)]
avg_val=0
for p in range(0, sizey_coarse):
    avg_val=0
    for q in range(0, sizex_coarse):
        vx_width[p]+=timeavg_velocity_grid_coarse_x[p][q]/float(timeavg_n_coarse[p][q])
        vy_width[p]+=timeavg_velocity_grid_coarse_y[p][q]/float(timeavg_n_coarse[p][q])
        avg_val+=1

    vx_width[p] = vx_width[p]/float(avg_val)
    vy_width[p] = vy_width[p]/float(avg_val)


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

    '''
    #fig = plt.figure(figsize=(8,6))
    fig = plt.figure(figsize=(5.452423529, 4.089317647))
    plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    plt.plot(time_conf, MSD, '-o' , color='darkgreen')
    plt.ylabel('MSD', fontsize=18)
    plt.xlabel('time', fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
    plt.show()
    #plt.savefig('./MSD_time.png')
    plt.close()
    '''

    with open('MSD.txt', 'w') as f:
        for i in range(len(MSD)):
            print(time_conf[i]*dt,MSD[i], file=f)  

    with open('mean_velocity.txt', 'w') as f:
        print(abs(avg_mean_velocity), file=f)  

    with open('v_width.txt', 'w') as f:
        for i in range(len(y)):
            print(vy_width[i], vx_width[i], y[i], file=f)  


print('done')
