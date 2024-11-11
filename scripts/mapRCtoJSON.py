import sys
import numpy as np
from math import *
import os.path
import matplotlib.pyplot as plt
import json

if len(sys.argv)!=3:
    print(sys.argv[0]," [input] [conf file]")
    sys.exit(1)

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

external_forces_file = 'external.conf '
lambda_wall = 3

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



tfile=open(topology,"r")
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

pt_num=0
x=np.arange(0,lx,1)
y=np.arange(0,ly,1)
Z=[[0 for q in range(lx)] for k in range(ly)]
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
cornerX=[0 for i in range(0,N)]
cornerY=[0 for i in range(0,N)]
patchMaxX=[0 for i in range(0,N)]
patchMaxY=[0 for i in range(0,N)]
row_phi = []

for line in cfile:
    phi_val=[]
    out_area=0
    words=line.split()
    Z=[[0 for q in range(lx)] for k in range(ly)]
    track_problem = 0
    LsubX[pt_num]=int(float(words[0]))
    LsubY[pt_num]=int(float(words[1]))
    CoMX[pt_num]=float(words[2])
    CoMY[pt_num]=float(words[3])
    offsetX[pt_num]=int(float(words[4]))
    offsetY[pt_num]=int(float(words[5]))
    corner=int(float(words[6]))
    cornerSite[pt_num]=corner
    cornerY[pt_num]=int(corner/lx)
    cornerX[pt_num]=corner-cornerY[pt_num]*lx

    patchMaxX[pt_num]=cornerX[pt_num]+LsubX[pt_num]
    if patchMaxX[pt_num]>=lx:
        patchMaxX[pt_num]=patchMaxX[pt_num]-lx
    patchMaxY[pt_num]=cornerY[pt_num]+LsubY[pt_num]
    if patchMaxY[pt_num]>=ly:
        patchMaxY[pt_num]=patchMaxY[pt_num]-ly
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
        if yy==40 and pt_num==0:
            print(xx, yy, value)
        phi_val.append(value)

        if value>0:
            Z[yy][xx]=value
        else:
            out_area=value*value

        area[pt_num]+=value*value
        if value>1.5 or value<-0.5:
            print("phase field is not in [0,1]!: ", pt_num, xx, yy, value)
            track_problem+=1

    if pt_num==0:
        print(Z[40][:])


    if abs(1-area[pt_num]/(pi*8*8)>0.5):
        print("area is not conserved: ", pt_num, area)
        cmap=cm.winter

    if out_area/area[pt_num] > 0.9:
        print("cell is leaking: ", pt_num, area, out_area)
        cmap=cm.autumn

    X, Y = np.meshgrid(x, y)
    row_phi.append(phi_val)

    #print(pt_num, area)
    pt_num+=1

print('done')


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


# a Python object (dict):
frame = {
    "id": "frame",
    "version": 1,
    "data": {
    "nphases": {"type" : "unsigned", "value" : N},
    "phi": {"type" : "array(array(double))" , "value" : all_phi},
    "offset": {"type" : "array(array(unsigned))" , "value" : [[offsetX[i], offsetY[i]] for i in range(N)]},
    "com": {"type" : "array(array(double))" , "value" : [[CoMX[i], CoMY[i]] for i in range(N)]},
    "Q00": {"type" : "array(double)", "value" : [Q00[i] for i in range(N)]},
    "Q01": {"type" : "array(double)", "value" : [Q01[i] for i in range(N)]},
    "theta_nem": {"type" : "array(double)", "value" : [theta_nem[i] for i in range(N)]},
    "patch_min": {"type" : "array(array(unsigned))" , "value" : [[cornerX[i], cornerY[i]] for i in range(N)]},
    "patch_max": {"type" : "array(array(unsigned))" , "value" : [[patchMaxX[i], patchMaxY[i]] for i in range(N)]}
    }
}

init_config = {
    "id": "init_config",
    "version": 1,
    "data": {
    "nphases": {"type" : "unsigned", "value" : N},
    "phi": {"type" : "array(array(double))" , "value" : all_phi},
    "offset": {"type" : "array(array(unsigned))" , "value" : [[offsetX[i], offsetY[i]] for i in range(N)]},
    "com": {"type" : "array(array(double))" , "value" : [[CoMX[i], CoMY[i]] for i in range(N)]},
    "Q00": {"type" : "array(double)", "value" : [Q00[i] for i in range(N)]},
    "Q01": {"type" : "array(double)", "value" : [Q01[i] for i in range(N)]},
    "theta_nem": {"type" : "array(double)", "value" : [theta_nem[i] for i in range(N)]},
    "patch_min": {"type" : "array(array(unsigned))" , "value" : [[cornerX[i], cornerY[i]] for i in range(N)]},
    "patch_max": {"type" : "array(array(unsigned))" , "value" : [[patchMaxX[i], patchMaxY[i]] for i in range(N)]}
    }
}

params = {
    "id": "parameters",
    "version": 1,
    "data": {
        "Size": {"type" : "array(unsigned)", "value" : [lx, ly]},
        "BC": {"type" : "unsigned", "value" : 0},
        "nsteps": {"type" : "unsigned", "value" : int(steps/dt)},
        "nsubsteps": {"type" : "unsigned", "value" : int(1/dt)},
        "ninfo": {"type" : "unsigned", "value" : int(print_conf_interval/dt)},
        "nstart": {"type" : "unsigned", "value" : 0},
        "gam": {"type" : "array(double)", "value" : [gamma for i in range(N)]},
        "mu": {"type" : "array(double)", "value" : [mu for i in range(N)]},
        "J0": {"type" : "double", "value" : J0},
        "lambda": {"type" : "double", "value" : llambda},
        "nphases": {"type" : "unsigned", "value" : N},
        "init_config": {"type" : "string", "value" : conf_file},
        "kappa": {"type" : "double", "value" : kappa},
        "types": {"type" : "array(string)", "value" : ["cell" for i in range(N)]},
        "xi": {"type" : "array(double)", "value" : [friction for i in range(N)]},
        "R": {"type" : "array(double)", "value" : [R for i in range(N)]},
        "omega": {"type" : "double", "value" : omega},
        "xi_cell": {"type" : "double", "value" : friction_cell},
        "wall_type": {"type" : "string", "value" : "nonadhesive"},
        "wall_thickness": {"type" : "double", "value" : lambda_wall},
        "wall_kappa": {"type" : "double", "value" : 4},
        "xi_wall": {"type" : "double", "value" : 0},
        "walls": {"type" : "array(double)", "value" : [walls[i] for i in range(lx*ly)]},
        "relax_time": {"type" : "unsigned", "value" : int(eq_steps/dt)},
        "relax_nsubsteps": {"type" : "unsigned", "value" : int(1/dt)},
        "npc": {"type" : "unsigned", "value" : 1},
        "seed": {"type" : "unsigned long", "value" : seed},
        "Jnem": {"type" : "array(double)", "value" : [J_Q for i in range(N)]},
        "patch_size": {"type" : "array(array(unsigned))" , "value" : [[LsubX[i], LsubY[i]] for i in range(N)]},
        "anchoring": {"type" : "bool", "value" : anchoring}
        }
}

# convert into JSON:
ff = json.dumps(frame)
pp = json.dumps(params)
ini = json.dumps(init_config)

# the result is a JSON string:
#print(ff)

pfile = open("parameters.json",'w')
ffile = open("frame0.json",'w')
ifile = open("init_config0.json",'w')
print(ff, file=ffile)
print(pp, file=pfile)
print(init_config, file=ifile)
pfile.close()
ffile.close()
ifile.close()

print("end")
