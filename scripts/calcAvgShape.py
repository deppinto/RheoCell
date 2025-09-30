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
    elif words[0]=='steps':
        steps=int(float(words[2]))


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


total_time_frames = int(steps / print_conf_interval)
n_rows = 20
n_columns = 5
theta_time = [[0. for j in range(total_time_frames)] for i in range(n_rows)]
elongation_time = [[0. for j in range(total_time_frames)] for i in range(n_rows)]
minor_axis_time = [[0. for j in range(total_time_frames)] for i in range(n_rows)]
aspect_ratio = [[0. for j in range(total_time_frames)] for i in range(n_rows)]
S_time = []


cont_line=0
vmax=0.1
start_line=(N+2)*int(float(sys.argv[3])) 
cfile=open(trajectory_file,"r")
for i in range(start_line):
    cfile.readline()


fig = plt.figure(figsize=(6,6))
frame_num=int(t/print_conf_interval)-1
theta_all = []
for line in cfile:
    cont_line+=1
    words=line.split()

    if words[0]=='t':
        t=int(float(words[2]))
        frame_num=int(t/print_conf_interval)-1
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
        theta_all = []

    elif words[0]=='b':
        lx=int(float(words[2]))
        ly=int(float(words[3]))
        x=np.arange(0,lx,1)
        y=np.arange(0,ly,1)

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

        nemQ_mod = sqrt(float(words[9])**2 + float(words[10])**2)
        nemX = sqrt((1 + float(words[9])/nemQ_mod)/2)
        nemY = np.sign(float(words[10]))*sqrt((1 - float(words[9])/nemQ_mod)/2)
        Q00[pt_num]=float(words[9])
        Q01[pt_num]=float(words[10])

        for i in range(start_value,len(words),2):
            site=int(float(words[i]))
            value=float(words[i+1])
            yy=int(site/lx)
            xx=site-int(yy*lx)

            Z[yy][xx]=value
            area[pt_num]+=value*value


        S00 = 0
        S01 = 0
        test_s00 = 0.
        test_s11 = 0.
        test_s01 = 0.
        for k in range(LsubX[pt_num]*LsubY[pt_num]):
            yy = int(k/LsubX[pt_num])
            xx = int(k - yy * LsubX[pt_num])
            xleft = int((xx - 1 + LsubX[pt_num]) % LsubX[pt_num])
            ybottom = int((yy - 1 + LsubY[pt_num]) % LsubY[pt_num])
            xright = int((xx + 1) % LsubX[pt_num])
            ytop = int((yy + 1) % LsubY[pt_num])

        for k in range(lx*ly):
            yy = int(k/lx)
            xx = int(k - yy * lx)
            xleft = int((xx - 1 + lx) % lx)
            ybottom = int((yy - 1 + ly) % ly)
            xright = int((xx + 1) % lx)
            ytop = int((yy + 1) % ly)

            field_dx = 0.5*(Z[yy][xright] - Z[yy][xleft])
            field_dy = 0.5*(Z[ytop][xx] - Z[ybottom][xx])

            S00 += -0.5*(field_dx * field_dx - field_dy * field_dy)
            S01 += -field_dx * field_dy

            test_s00 += field_dx * field_dx
            test_s11 += field_dy * field_dy
            test_s01 += field_dx * field_dy

        D_major_axis = 0.5 * np.atan2(S01, S00)
        D_major_axis_vec_x = np.cos(D_major_axis)
        D_major_axis_vec_y = np.sin(D_major_axis)
        D_i = np.sqrt(S00 * S00 + S01 * S01)

        D_i = np.sqrt(S00 * S00 + S01 * S01)
        if D_i > 0.000000001:
            D_major_axis_vec_x = D_i * sqrt((1 + S00/D_i)/2)
            D_major_axis_vec_y = D_i * np.sign(S01) * sqrt((1 - S00/D_i)/2)

            D_minor_axis_vec_x = - D_i * np.sign(S01) * np.sqrt((1 - S00/D_i)/2)
            D_minor_axis_vec_y =   D_i * np.sqrt((1 + S00/D_i)/2)
        else:
            D_major_axis_vec_x = D_major_axis_vec_y = 0.
            D_minor_axis_vec_x = D_minor_axis_vec_y = 0.

        #F = np.array([[S00, S01],[S01, -S00]], dtype=float)
        F = np.array([[test_s00, test_s01],[test_s01, test_s11]], dtype=float)

        U, S, Vt = np.linalg.svd(F)   # S[0] >= S[1]
        V = Vt.T                      # columns are principal directions in the reference frame

        D_major_axis_vec = S[0] * V[:, 0]   # scaled major axis (length = σ1)
        D_minor_axis_vec = S[1] * V[:, 1]   # scaled minor axis (length = σ2)

        
        theta_i = (0.5 * np.atan2(S01, S00) * 180 / pi)
        if pt_num - int(pt_num / n_columns) * n_columns == 0:
            theta_time[int(pt_num / n_columns)][frame_num] = (0.5 * np.atan2(S01, S00) * 180 / pi)
            elongation_time[int(pt_num / n_columns)][frame_num] = sqrt(D_major_axis_vec_x**2 + D_major_axis_vec_y**2)
            minor_axis_time[int(pt_num / n_columns)][frame_num] = sqrt(D_minor_axis_vec_x**2 + D_minor_axis_vec_y**2)
            aspect_ratio[int(pt_num / n_columns)][frame_num] = sqrt(S[0]/S[1]) - 1
            #elongation_time[int(pt_num / n_columns)][frame_num] = sqrt(D_major_axis_vec[0]**2 + D_major_axis_vec[1]**2)
            #minor_axis_time[int(pt_num / n_columns)][frame_num] = sqrt(D_minor_axis_vec[0]**2 + D_minor_axis_vec[1]**2)
            #print(elongation_time[int(pt_num / n_columns)][frame_num], minor_axis_time[int(pt_num / n_columns)][frame_num])
            if int(pt_num/n_columns)==9:
                print(S[0], S[1], (0.5 * np.atan2(S01, S00) * 180 / pi), sqrt(S[0]/S[1]) - 1)

        if CoMY[pt_num] > 50 and CoMY[pt_num] < ly - 50: 
            theta_all.append(theta_i)


        X, Y = np.meshgrid(x, y)
        step = 0.01
        m = np.amax(Z)
        #if m<0.000001:
            #continue

        levels = np.arange(0.0, m, step) + step

        if variable==1 or variable==2 or variable==3:
            if pt_num==-1:
                cset1 = plt.contour(X, Y, Z, levels, cmap=cm.winter, alpha=0.5)
            else:
                cset1 = plt.contour(X, Y, Z, levels=[0.5], cmap=cm.winter, alpha=0.5)

            #cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], 2*D_major_axis_vec[0], 2*D_major_axis_vec[1], width=0.5, head_width=0, color='r')
            #cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], -2*D_major_axis_vec[0], -2*D_major_axis_vec[1], width=0.5, head_width=0, color='r')
            cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], 2*D_major_axis_vec_x, 2*D_major_axis_vec_y, width=0.5, head_width=0, color='r')
            cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], -2*D_major_axis_vec_x, -2*D_major_axis_vec_y, width=0.5, head_width=0, color='r')

        #increment phase field index
        pt_num+=1
        

    if cont_line%(N+2)==0:

            sin_sum = 0.
            cos_sum = 0.
            for q in theta_all:
                sin_sum += sin(2 * q)
                cos_sum += cos(2 * q)
            mean_phi = 0.5 * atan2(sin_sum, cos_sum)

            avg_value_S = 0.
            for q in theta_all:
                avg_value_S += cos(2*( q - mean_phi))
            S_time.append(avg_value_S / len(theta_all))
            #print(avg_value_S / len(theta_all))


            frame_num=int(t/print_conf_interval)-1
            #print(frame_num)
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
            if variable==2 or variable==4:
                plt.show()
                #plt.savefig('./newfig_'+str(frame_num)+'.png', transparent=True)
            if variable<=4:
                plt.clf()

plt.close()



if variable==5:

    #fig = plt.figure(figsize=(8,6))
    fig = plt.figure(figsize=(5.452423529, 4.089317647))
    plt.plot(theta_time, '-o' , color='firebrick', label='Row 9')
    #plt.plot(theta_time_2, '-s' , color='green', label='Row 7')
    #plt.plot(theta_time_3, '-^' , color='royalblue', label='Row 5')
    #plt.plot(theta_time_4, '-p' , color='goldenrod', label='Row 3')
    plt.ylabel(r'$\theta_i$', fontsize=18)
    plt.xlabel('Time', fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=12, frameon=False)
    plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
    plt.show()
    #plt.savefig('./theta_width_coarse.png')
    plt.close()

if variable==6:

    with open('theta_shape.txt', 'w') as f:
        for i in range(total_time_frames):
            print_str = ''
            for j in range(n_rows):
                print_str += str(theta_time[j][i])
                print_str += ' '

            print(print_str , file=f)  

    with open('elongation_shape.txt', 'w') as f:
        for i in range(total_time_frames):
            print_str = ''
            for j in range(n_rows):
                print_str += str(elongation_time[j][i])
                print_str += ' '

            print(print_str , file=f)  

    with open('elongation_minor_shape.txt', 'w') as f:
        for i in range(total_time_frames):
            print_str = ''
            for j in range(n_rows):
                print_str += str(minor_axis_time[j][i])
                print_str += ' '

            print(print_str , file=f)  


    with open('aspect_ratio_shape.txt', 'w') as f:
        for i in range(total_time_frames):
            print_str = ''
            for j in range(n_rows):
                print_str += str(aspect_ratio[j][i])
                print_str += ' '

            print(print_str , file=f)  

print('done')
