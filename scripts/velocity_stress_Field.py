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

from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import font_manager, rcParams

'''
# Load font from file
#font_path = "/home/p/pinto/Fonts/Times_New_Roman_Normal.ttf"
font_path = "/home/p/pinto/Fonts/times.ttf"
italic_font_path = "/home/p/pinto/Fonts/timesi.ttf"
bold_font_path = "/home/p/pinto/Fonts/timesbd.ttf"

custom_font = FontProperties(fname=font_path)
legend_font = FontProperties(fname=font_path, size=12)
font = font_manager.FontProperties(fname=font_path)
fonti = font_manager.FontProperties(fname=italic_font_path)
fontbd = font_manager.FontProperties(fname=bold_font_path)

# Register font with a name
font_manager.fontManager.addfont(font_path)
font_manager.fontManager.addfont(italic_font_path)
font_manager.fontManager.addfont(bold_font_path)

# Set custom mathtext font to match your font
rcParams['mathtext.fontset'] = 'custom'

# Set roman (upright), italic, and bold versions (all Times New Roman if needed)
#rcParams['mathtext.rm'] = font.get_name()  # e.g. "Times New Roman"
#rcParams['mathtext.it'] = fonti.get_name()
#rcParams['mathtext.bf'] = fontbd.get_name()
rcParams['mathtext.it'] = 'Times New Roman:italic'
rcParams['mathtext.rm'] = 'Times New Roman'
rcParams['mathtext.bf'] = 'Times New Roman:bold'
'''


if len(sys.argv)!=5:
    print(sys.argv[0]," [topology file] [velocity file] [stress file] [1:save conf; 2:make plot]")
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
    for i in range(start_value,len(words),4):
        xx=float(words[i])
        yy=float(words[i+1])

        value_x=float(words[i+2])
        value_y=float(words[i+3])
        Z_x[int(yy)][int(xx)]=value_x
        Z_y[int(yy)][int(xx)]=value_y

        if int(xx)%4==0 and int(yy)%4==0:
            cset1 = plt.arrow(xx, yy, 60*value_x, 60*value_y, width=0.5, color='k')
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
    p = 1.
    ang = atan2(abs(a[0]*b[1]-a[1]*b[0]), a[0]*b[0]+a[1]*b[1])

    #if(ang > pi/2.):
        #b = [-i for i in b]

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


LLX = lx
LLY = ly
vecfield_nx = [[0. for j in range(0, LLX)] for i in range(0, LLY)]
vecfield_ny = [[0. for j in range(0, LLX)] for i in range(0, LLY)]

vecfield_nx = Z_x
vecfield_ny = Z_y

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
            '''
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
            '''
            if s==1:
                cset1 = plt.plot(x, y, 'go', markersize=10)
                #cset1 = plt.arrow(x, y, 4*cos(psi), 4*sin(psi), color='g', head_width=1.5, head_length=1.5, width=0.5)
            elif s==-1:
                cset1 = plt.plot(x, y, 'b^', markersize=10)


        # keep this just in case our other symmetries give us integer defects
        elif (abs(winding_number[p][q]) > 1 - thresh) and (abs(winding_number[p][q])< 1 + thresh):
            # charge sign
            s = np.sign(winding_number[p][q])
            # bfs
            sum_x, sum_y, n = collapse(p, q, s, LLX, LLY, winding_number, rng = [1-thresh, 1+thresh])
            x,y = sum_x/n,sum_y/n
            # add defect to list
            if s==1:
                #if y<130:
                    #print(x,y)
                cset1 = plt.plot(x, y, '*', color='#00FFFF', markersize=10)
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
Z_xy=[[0 for q in range(lx)] for k in range(ly)]
Z_iso=[[0 for q in range(lx)] for k in range(ly)]
start_value=0

for line in cfile:
    words=line.split()
    Z_xx=[[0 for q in range(lx)] for k in range(ly)]
    Z_yy=[[0 for q in range(lx)] for k in range(ly)]
    Z_xy=[[0 for q in range(lx)] for k in range(ly)]
    Z_iso=[[0 for q in range(lx)] for k in range(ly)]
    for i in range(start_value,len(words),13):
        xx=float(words[i])
        yy=float(words[i+1])

        value_field_xx=float(words[i+2])
        value_field_yy=float(words[i+3])
        value_field_xy=float(words[i+4])
        value_pass_xx=float(words[i+5])
        value_pass_yy=float(words[i+6])
        value_pass_xy=float(words[i+7])
        value_act_xx=float(words[i+8])
        value_act_yy=float(words[i+9])
        value_act_xy=float(words[i+10])
        value_pre_xx=float(words[i+11])
        value_pre_yy=float(words[i+12])

        Z_xx[int(yy)][int(xx)] = value_field_xx
        Z_yy[int(yy)][int(xx)] = value_field_yy
        Z_xy[int(yy)][int(xx)] = value_field_xy
        Z_iso[int(yy)][int(xx)] = (Z_xx[int(yy)][int(xx)] + Z_yy[int(yy)][int(xx)]) / 2
        #Z_iso[int(yy)][int(xx)] = Z_xy[int(yy)][int(xx)]

    z_min, z_max = 0., np.abs(Z_iso).max()
    for i in range(ly):
        for j in range(lx):
            Z_iso[i][j] = Z_iso[i][j] / z_max
    cset1 = plt.imshow(Z_iso, cmap='RdBu_r', interpolation='nearest', vmin=-1, vmax=1)


'''
div_stress = 0.
total_stress = 0.
for p in range(ly):
    ynext = (p + 1) % ly
    yprev = (p - 1 + ly) % ly
    for q in range(lx):
        total_stress += Z_iso[p][q]

        xnext = (q + 1) % lx
        xprev = (q - 1 + lx) % lx

        dvxdx = (Z_iso[p][xnext] - Z_iso[p][xprev])/2
        dvydy = (Z_iso[ynext][q] - Z_iso[yprev][q])/2

        div_stress += dvxdx + dvydy


print("Stresses: ", total_stress, div_stress)
'''


ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
ax.set_xlim([0, lx])
ax.set_ylim([0, ly])
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.set_xticks([])
ax.set_yticks([])

#fig = plt.gcf()
# Add colorbar
#cbar = fig.colorbar(cset1, ax=ax)
# Change tick font name and size
#cbar.ax.tick_params(labelsize=18)  # Font size
#for label in cbar.ax.get_yticklabels():
#    label.set_fontname('Times New Roman')  # Replace with desired font
#cbar.set_ticks([-1, -0.5, 0, 0.5, 1])

# Set label above the colorbar
#cbar.set_label(r'$\dot{\varepsilon}$,$\nabla f^a$', rotation=0, pad=10, loc='top', fontsize=18, fontname='Times New Roman')
#cbar.ax.xaxis.set_label_position('top')  # Move label to top
#cbar.ax.set_title(r'$\dot{\varepsilon}$ , $\nabla \cdot \mathbf{f}^\mathit{a}$', pad=10, fontsize=18, fontname='Times New Roman')

#fig.tight_layout()
if variable==1:
    plt.savefig('frame.png')
if variable==2:
    plt.show()
