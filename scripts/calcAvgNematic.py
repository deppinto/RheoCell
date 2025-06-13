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


Q00_grid =[[0. for j in range(0, lx)] for i in range(0, ly)]
Q01_grid =[[0. for j in range(0, lx)] for i in range(0, ly)]
n_grid=[[0. for j in range(0, lx)] for i in range(0, ly)]
Q00_avg_grid =[[0. for j in range(0, lx)] for i in range(0, ly)]
Q01_avg_grid =[[0. for j in range(0, lx)] for i in range(0, ly)]
size_grid = 3*R+1


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
        Q00_grid =[[0. for j in range(0, lx)] for i in range(0, ly)]
        Q01_grid =[[0. for j in range(0, lx)] for i in range(0, ly)]
        Q00_avg_grid =[[0. for j in range(0, lx)] for i in range(0, ly)]
        Q01_avg_grid =[[0. for j in range(0, lx)] for i in range(0, ly)]
        n_grid=[[0. for j in range(0, lx)] for i in range(0, ly)]
        #walls = [0. for i in range(lx*ly)]
        #set_walls(lx,ly,walls)

    else:
        Z=[[0. for q in range(lx)] for k in range(ly)]
        LsubX[pt_num]=int(float(words[0]))
        LsubY[pt_num]=int(float(words[1]))
        CoMX[pt_num]=float(words[2])
        CoMY[pt_num]=float(words[3])

        #CoMX[pt_num] += 35
        #CoMY[pt_num] -= 25
        #if CoMY[pt_num] < 0:
            #CoMY[pt_num] += ly
        #if CoMX[pt_num] >= lx:
            #CoMX[pt_num] -= lx

        offsetX[pt_num]=int(float(words[4]))
        offsetY[pt_num]=int(float(words[5]))
        cornerSite[pt_num]=int(float(words[6]))
        cornerSite_x[pt_num]=int(float(words[7]))
        cornerSite_y[pt_num]=int(float(words[8]))

        nemX=float(words[9])
        nemY=float(words[10])
        normNem = sqrt(nemX * nemX + nemY * nemY)
        #theta_nem[pt_num]=asin((nemX*nemY)/0.5)/2
        Q00[pt_num]= 0.5 * (nemX * nemX - nemY * nemY)
        Q01[pt_num]= nemX * nemY
        #print(normNem)

        for i in range(start_value,len(words),2):
            site=int(float(words[i]))
            value=float(words[i+1])
            yy=int(site/lx)
            xx=site-int(yy*lx)

            #xx = xx + 35
            #yy = yy - 25
            #if yy < 0:
                #yy += ly
            #if xx >= lx:
                #xx -= lx

            Z[yy][xx]=value
            area[pt_num]+=value*value

            Q00_grid[yy][xx] += value * Q00[pt_num]
            Q01_grid[yy][xx] += value * Q01[pt_num]
            n_grid[yy][xx] += value

            n_coarse[int(yy/deltay_coarse)][int(xx/deltax_coarse)]+=value
            nematic_grid_coarse_00[int(yy/deltay_coarse)][int(xx/deltax_coarse)]+=value*Q00[pt_num]
            nematic_grid_coarse_01[int(yy/deltay_coarse)][int(xx/deltax_coarse)]+=value*Q01[pt_num]

        S00 = 0
        S01 = 0
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

        #D_major_axis = 0.5 * np.atan2(S01, S00)
        #D_major_axis_vec_x = np.cos(D_major_axis)
        #D_major_axis_vec_y = np.sin(D_major_axis)
        #D_i = np.sqrt(S00 * S00 + S01 * S01)

        D_i = np.sqrt(S00 * S00 + S01 * S01)
        if D_i > 0.000000001:
            D_major_axis_vec_x = D_i * sqrt((1 + S00/D_i)/2)
            D_major_axis_vec_y = D_i * np.sign(S01) * sqrt((1 - S00/D_i)/2)
        else:
            D_major_axis_vec_x = 0
            D_major_axis_vec_y = 0

        #print(2 * D_i)

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

            cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], 3*nemX, 3*nemY, width=0.5, head_width=0, color='k')
            cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], -3*nemX, -3*nemY, width=0.5, head_width=0, color='k')

            cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], 1*D_major_axis_vec_x, 1*D_major_axis_vec_y, width=0.5, head_width=0, color='r')
            cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], -1*D_major_axis_vec_x, -1*D_major_axis_vec_y, width=0.5, head_width=0, color='r')

        #increment phase field index
        pt_num+=1
        

    if cont_line%(N+2)==0:

            '''
            LLX = sizex_coarse
            LLY = sizey_coarse
            vecfield_nx = [[0. for j in range(0, LLX)] for i in range(0, LLY)]
            vecfield_ny = [[0. for j in range(0, LLX)] for i in range(0, LLY)]
            vecfield_Q00 = [[0. for j in range(0, LLX)] for i in range(0, LLY)]
            vecfield_Q01 = [[0. for j in range(0, LLX)] for i in range(0, LLY)]
            for p in range(0, LLY):
                for q in range(0, LLX):
                    S = 0.
                    nx = 0.
                    ny = 0.
                    if n_coarse[p][q]>0:
                        Q00=nematic_grid_coarse_00[p][q]/n_coarse[p][q]
                        Q01=nematic_grid_coarse_01[p][q]/n_coarse[p][q]
                        S = sqrt(Q00**2 + Q01**2)
                        nx = sqrt(2*S)*sqrt((1 + Q00/S)/2)
                        ny = sqrt(2*S)*np.sign(Q01)*sqrt((1 - Q00/S)/2)

                        vecfield_nx[p][q] = nx
                        vecfield_ny[p][q] = ny
                        vecfield_Q00[p][q] = Q00
                        vecfield_Q01[p][q] = Q01

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
                            cset1 = plt.arrow(q*deltax_coarse+deltax_coarse/2, p*deltay_coarse+deltay_coarse/2, 0.4*deltax_coarse*nx, 0.4*deltay_coarse*ny, width=deltax_coarse/15, color="w", head_width=0)
                            cset1 = plt.arrow(q*deltax_coarse+deltax_coarse/2, p*deltay_coarse+deltay_coarse/2, -0.4*deltax_coarse*nx, -0.4*deltay_coarse*ny, width=deltax_coarse/15, color="w", head_width=0)
            '''

            LLX = lx
            LLY = ly
            vecfield_nx = [[0. for j in range(0, LLX)] for i in range(0, LLY)]
            vecfield_ny = [[0. for j in range(0, LLX)] for i in range(0, LLY)]
            vecfield_Q00 = [[0. for j in range(0, LLX)] for i in range(0, LLY)]
            vecfield_Q01 = [[0. for j in range(0, LLX)] for i in range(0, LLY)]

            if variable==3 or variable==4:
                for p in range(0, LLY):
                    for q in range(0, LLX):
                        if n_grid[p][q]>1e-1:
                            for k in range(0, size_grid):
                                for l in range(0, size_grid):
                                    pp = (((p - int(size_grid/2) + LLY) % LLY) + k) % LLY
                                    qq = (((q - int(size_grid/2) + LLX) % LLX) + l) % LLX
                                    if n_grid[pp][qq]>1e-1:
                                        Q00_avg_grid[p][q] += (Q00_grid[pp][qq] / n_grid[pp][qq]) / (size_grid * size_grid)
                                        Q01_avg_grid[p][q] += (Q01_grid[pp][qq] / n_grid[pp][qq]) / (size_grid * size_grid)

                        S = 0.
                        nx = 0.
                        ny = 0.
                        if n_grid[p][q]>1e-1:
                            Q00=Q00_avg_grid[p][q]
                            Q01=Q01_avg_grid[p][q]
                            S = sqrt(Q00**2 + Q01**2)
                            nx = sqrt(2*S) * sqrt((1 + Q00/S)/2)
                            ny = sqrt(2*S) * np.sign(Q01) * sqrt((1 - Q00/S)/2)

                            vecfield_nx[p][q] = nx
                            vecfield_ny[p][q] = ny
                            vecfield_Q00[p][q] = Q00
                            vecfield_Q01[p][q] = Q01

                            if variable==3 or variable==4:
                                if q%2==0 and p%2==0:
                                    cset1 = plt.arrow(q+1/2, p+1/2, 0.5*nx, 0.5*ny, width=1/15, color="k", head_width=0)
                                    cset1 = plt.arrow(q+1/2, p+1/2, -0.5*nx, -0.5*ny, width=1/15, color="k", head_width=0)
                        else:
                            if variable==3 or variable==4:
                                if q%2==0 and p%2==0:
                                    cset1 = plt.arrow(q+1/2, p+1/2, 0.5*nx, 0.5*ny, width=1/15, color="w", head_width=0)
                                    cset1 = plt.arrow(q+1/2, p+1/2, -0.5*nx, -0.5*ny, width=1/15, color="w", head_width=0)

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
                        plotted_flag=0
                        # detect simplest charge 1/p defects
                        if  (abs(winding_number[p][q]) > charge - thresh) and (abs(winding_number[p][q]) < charge + thresh):
                            # charge sign
                            s = np.sign(winding_number[p][q])
                            # bfs
                            sum_x, sum_y, n = collapse(p, q, s, LLX, LLY, winding_number, rng = [charge - thresh, charge + thresh])
                            x,y = sum_x/n,sum_y/n

                            # compute angle, see doi:10.1039/c6sm01146b
                            '''
                            num = 0
                            den = 0
                            for (dx, dy) in [(0, 0), (0, 1), (1, 1), (1, 0)]:
                                # coordinates of nodes around the defect
                                kk = (int(x) + sizex_coarse + dx) % sizex_coarse
                                ll = (int(y) + sizey_coarse + dy) % sizey_coarse
                                # derivative at these points
                                dxQxx = .5*(vecfield_Q00[(kk+1) % sizex_coarse, ll] - vecfield_Q00[(kk-1+sizex_coarse) % sizex_coarse, ll])
                                dxQxy = .5*(vecfield_Q01[(kk+1) % sizex_coarse, ll] - vecfield_Q01[(kk-1+sizex_coarse) % sizex_coarse, ll])
                                dyQxx = .5*(vecfield_Q00[kk, (ll+1) % sizey_coarse] - vecfield_Q00[kk, (ll-1+sizey_coarse) % sizey_coarse])
                                dyQxy = .5*(vecfield_Q01[kk, (ll+1) % sizey_coarse] - vecfield_Q01[kk, (ll-1+sizey_coarse) % sizey_coarse])
                                # accumulate numerator and denominator
                                num += s*dxQxy - dyQxx
                                den += dxQxx + s*dyQxy
                            psi = s/(2.-s)*atan2(num, den)
                            '''
                            if variable==3 or variable==4:
                                if s==1:
                                    if n_grid[int(x)][int(y)]>1e-1:
                                        cset1 = plt.plot(x, y, 'go', markersize=10)
                                        #cset1 = plt.plot(x*sizex_coarse+sizex_coarse/2, y*sizey_coarse+sizey_coarse/2, 'go', markersize=10)
                                        #cset1 = plt.arrow(x, y, -0.4*deltax_coarse*cos(psi), -0.4*deltay_coarse*sin(psi), width=deltax_coarse/15, color="r")
                                        plotted_flag=1
                                elif s==-1:
                                    if n_grid[int(x)][int(y)]>1e-1:
                                        cset1 = plt.plot(x, y, 'b^', markersize=10)
                                        #cset1 = plt.plot(x*sizex_coarse+sizex_coarse/2, y*sizey_coarse+sizey_coarse/2, 'b^', markersize=10)
                                        plotted_flag=1

                        # keep this just in case our other symmetries give us integer defects
                        elif (abs(winding_number[p][q]) > 1 - thresh) and (abs(winding_number[p][q])< 1 + thresh):
                            # charge sign
                            s = np.sign(winding_number[p][q])
                            # bfs
                            sum_x, sum_y, n = collapse(p, q, s, sizex_coarse, sizey_coarse, winding_number, rng = [0.9,1.1])
                            x,y = sum_x/n,sum_y/n
                            # add defect to list
                            if variable==3 or variable==4:
                                if s==1:
                                    if n_grid[int(x)][int(y)]>1e-1:
                                        cset1 = plt.plot(x, y, 'r*', markersize=10)
                                        #cset1 = plt.plot(x*sizex_coarse, y*sizey_coarse, 'r*', markersize=10)
                                        plotted_flag=1
                                elif s==-1:
                                    if n_grid[int(x)][int(y)]>1e-1:
                                        cset1 = plt.plot(x, y, 'kX', markersize=10)
                                        #cset1 = plt.plot(x*sizex_coarse, y*sizey_coarse, 'kX', markersize=10)
                                        plotted_flag=1

                        #if plotted_flag==0 and abs(winding_number[p][q]) > 0.1:
                            #print("Defect topology: ", q, p, winding_number[p][q])

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
