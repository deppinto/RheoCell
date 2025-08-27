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


from shapely.geometry import Polygon, LineString, Point, box
from scipy.spatial import Voronoi
from scipy.spatial import Delaunay

bounds = (0, 14, 70, 294)
def bounded_voronoi(points, bounds):
    vor = Voronoi(points)
    min_x, min_y, max_x, max_y = bounds
    bbox = box(min_x, min_y, max_x, max_y)

    regions = []
    for point_idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if -1 in region or len(region) == 0:
            # Skip infinite or degenerate regions
            continue
        poly_coords = [vor.vertices[i] for i in region]
        poly = Polygon(poly_coords)
        clipped_poly = poly.intersection(bbox)
        if not clipped_poly.is_empty:
            regions.append((point_idx, clipped_poly))
    return vor, regions


def plot_periodic_delaunay(points, box_lengths):
    """
    Plot points and their Delaunay triangulation under periodic BCs
    in a rectangular 2D box.
    """
    N = len(points)
    Lx, Ly = box_lengths

    # Tile the system
    shifts = np.array([[dx * Lx, dy * Ly] 
                       for dx in (-1,0,1) 
                       for dy in (-1,0,1)])
    tiled_points = np.vstack([points + shift for shift in shifts])
    tiled_indices = np.repeat(np.arange(N), len(shifts))

    # Triangulate
    tri = Delaunay(tiled_points)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.set_aspect("equal")

    # Draw box
    ax.plot([0, Lx, Lx, 0, 0], [0, 0, Ly, Ly, 0], 'k-', lw=1)

    # Draw triangulation edges
    for simplex in tri.simplices:
        for i in range(3):
            a, b = simplex[i], simplex[(i+1)%3]
            ia, ib = tiled_indices[a], tiled_indices[b]

            # Only draw edges if at least one point is in the central box
            if (0 <= tiled_points[a,0] < Lx and 0 <= tiled_points[a,1] < Ly) or \
               (0 <= tiled_points[b,0] < Lx and 0 <= tiled_points[b,1] < Ly):
                x = [tiled_points[a,0] % Lx, tiled_points[b,0] % Lx]
                y = [tiled_points[a,1] % Ly, tiled_points[b,1] % Ly]
                if ia==50 or ib==50:
                    ax.plot(x, y, 'b-', lw=0.8, alpha=0.7)

    # Plot points
    ax.scatter(points[:,0], points[:,1], c="red", zorder=5)

    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)
    plt.show()


def minimum_image(vec, box_lengths):
    """Apply minimum image convention for a displacement vector."""
    return vec - box_lengths * np.round(vec / box_lengths)

def periodic_delaunay_neighbors(points, box_lengths):
    """
    Compute nearest neighbors in 2D with periodic boundary conditions
    using Delaunay triangulation in a rectangular box.

    Parameters
    ----------
    points : (N, 2) ndarray
        Array of point coordinates inside [0, Lx) x [0, Ly).
    box_lengths : (2,) array_like
        (Lx, Ly) lengths of the periodic box.

    Returns
    -------
    neighbors : list of sets
        neighbors[i] is a set of indices of nearest neighbors of point i.
    """
    N = len(points)
    Lx, Ly = box_lengths

    # Generate shifts for 3x3 tiling (center + neighbors in x and y)
    shifts = np.array([
        [dx * Lx, dy * Ly]
        for dx in (-1, 0, 1)
        for dy in (-1, 0, 1)
    ])

    # Tile the system
    tiled_points = np.vstack([points + shift for shift in shifts])
    tiled_indices = np.tile(np.arange(N), len(shifts))

    # Delaunay triangulation
    tri = Delaunay(tiled_points)


    # Build neighbor sets
    neighbors = [set() for _ in range(N)]
    for simplex in tri.simplices:
        for i in range(3):
            for j in range(i+1, 3):
                a, b = simplex[i], simplex[j]
                ia, ib = tiled_indices[a], tiled_indices[b]
                pa, pb = tiled_points[a], tiled_points[b]

                # Only keep edges if at least one endpoint is in the central box
                if (0 <= pa[0] < Lx and 0 <= pa[1] < Ly) or (0 <= pb[0] < Lx and 0 <= pb[1] < Ly):
                    if ia != ib:  # skip self-links from periodic images
                        neighbors[ia].add(ib)
                        neighbors[ib].add(ia)

    return neighbors

def bond_orientational_order(points, neighbors, box_lengths, n=6):
    """
    Compute the n-fold bond orientational order parameter for each point.

    Parameters
    ----------
    points : (N,2) ndarray
        Particle positions
    neighbors : list of lists
        neighbors[i] is a list of indices of neighbors of particle i
    box_lengths : (2,) tuple
        Box dimensions (Lx, Ly)
    n : int
        Order parameter symmetry (6 for hexatic)

    Returns
    -------
    psi_n : (N,) complex ndarray
        n-fold bond orientational order parameter for each particle
    """
    N = len(points)
    psi_n = np.zeros(N, dtype=complex)

    for i in range(N):
        psi = 0.0 + 0.0j
        Ni = len(neighbors[i])
        if Ni == 0:
            continue
        for j in neighbors[i]:
            vec = minimum_image(points[j] - points[i], box_lengths)
            theta = np.arctan2(vec[1], vec[0])
            psi += np.exp(1j * n * theta)
        psi_n[i] = psi / Ni
    return psi_n


def nematic_order(deformations, neighbors, box_lengths, n=2):

    N = len(points)
    psi_n = np.zeros(N, dtype=complex)

    for i in range(N):
        psi = 0.0 + 0.0j
        Ni = len(neighbors[i])
        if Ni == 0:
            continue
        for j in neighbors[i]:
            theta1 = np.arctan2(deformations[i][0]*deformations[j][1] - deformations[i][1]*deformations[j][0], np.dot(deformations[i], deformations[j]))
            theta2 = np.arctan2(-deformations[i][0]*deformations[j][1] + deformations[i][1]*deformations[j][0], np.dot(-deformations[i], deformations[j]))
            theta=theta2
            if theta1 < theta2:
                theta = theta1
            psi += np.exp(1j * n * theta)
        psi_n[i] = psi / Ni
    return psi_n


def lane_order_local(points, neighbors, box_lengths, Ly):
    """
    Compute lane order with direction switching at Ly/2.

    points: (N,2)
    neighbors: list of lists
    box_lengths: (Lx, Ly)
    Ly: box height (for midline)

    Returns:
        psi_lane: (N,) alignment magnitude along local lane direction
    """
    N = len(points)
    psi_lane = np.zeros(N)

    for i in range(N):
        # choose local lane direction
        lane_dir = np.array([1.0, 0.0]) if points[i,1] > Ly/2 else np.array([-1.0, 0.0])
        lane_dir /= np.linalg.norm(lane_dir)

        Ni = len(neighbors[i])
        if Ni == 0:
            continue
        proj_sum = 0.0
        for j in neighbors[i]:
            vec = points[j] - points[i]
            # minimum image
            vec -= box_lengths * np.round(vec / box_lengths)
            proj_sum += np.abs(np.dot(vec, lane_dir)) / np.linalg.norm(vec)
        psi_lane[i] = proj_sum / Ni
    return psi_lane



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
D_X=[0. for i in range(0,N)]
D_Y=[0. for i in range(0,N)]

lfile=open(lastconf_file,"r")
header=lfile.readline().split()
t=int(header[2])
header=lfile.readline().split()
lx=int(float(header[2]))
ly=int(float(header[3]))
x=np.arange(0,lx,1)
y=np.arange(0,ly,1)


order_psi6 = []
order_psiN = []
order_psiL = []


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
        D_X=[0. for i in range(0,N)]
        D_Y=[0. for i in range(0,N)]

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
        #normNem = sqrt(nemX * nemX + nemY * nemY)
        #theta_nem[pt_num]=asin((nemX*nemY)/0.5)/2
        #Q00[pt_num]= 0.5 * (nemX * nemX - nemY * nemY)
        #Q01[pt_num]= nemX * nemY
        #print(normNem)

        for i in range(start_value,len(words),2):
            site=int(float(words[i]))
            value=float(words[i+1])
            yy=int(site/lx)
            xx=site-int(yy*lx)

            Z[yy][xx]=value
            area[pt_num]+=value*value

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

        D_X[pt_num] = D_major_axis_vec_x
        D_Y[pt_num] = D_major_axis_vec_y

        X, Y = np.meshgrid(x, y)
        step = 0.01
        m = np.amax(Z)
        #if m<0.000001:
            #continue

        levels = np.arange(0.0, m, step) + step

        if variable==1 or variable==2 or variable==3:
            #if pt_num==52:
                #cset1 = plt.contour(X, Y, Z, levels, cmap=cm.winter, alpha=0.5)
            #else:
                #cset1 = plt.contour(X, Y, Z, levels=[0.5], cmap=cm.winter, alpha=0.5)

            #cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], 3*nemX, 3*nemY, width=0.5, head_width=0, color='k')
            #cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], -3*nemX, -3*nemY, width=0.5, head_width=0, color='k')

            cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], 2*D_major_axis_vec_x, 2*D_major_axis_vec_y, width=0.5, head_width=0, color='r')
            cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], -2*D_major_axis_vec_x, -2*D_major_axis_vec_y, width=0.5, head_width=0, color='r')

        #increment phase field index
        pt_num+=1
        

    if cont_line%(N+2)==0:

            deformations = np.array([[D_X[qq], D_Y[qq]] for qq in range(0, N)])
            points = np.array([[CoMX[qq], CoMY[qq]] for qq in range(0, N)])
            vor = Voronoi(points)

            neigh = periodic_delaunay_neighbors(points, (lx, ly))
            #plot_periodic_delaunay(points, (lx, ly))
            #for i, nset in enumerate(neigh):
                #print(f"Point {i} neighbors: {sorted(list(nset))}")

            psi6 = bond_orientational_order(points, neigh, (lx, ly), n=6)
            psiN = nematic_order(deformations, neigh, (lx, ly), n=2)
            psiL = lane_order_local(points, neigh, (lx, ly), ly)
            cont_bulk = 0
            Psi6_global = 0.
            PsiN_global = 0.
            PsiL_global = 0.
            for pp in range(N):
                if CoMY[pp]>50 and CoMY[pp]<ly-50:
                    cont_bulk+=1
                    Psi6_global += psi6[pp]
                    PsiN_global += psiN[pp]
                    PsiL_global += psiL[pp]
            Psi6_global = np.abs(Psi6_global/cont_bulk)
            PsiN_global = np.abs(PsiN_global/cont_bulk)
            PsiL_global = PsiL_global/cont_bulk
            order_psi6.append(Psi6_global)
            order_psiN.append(PsiN_global)
            order_psiL.append(PsiL_global)
            #print("Global hexatic order:", Psi6_global)
            #print("Global nematic order:", PsiN_global)
            #print("Global lane order:", PsiL_global)


            frame_num=int(t/print_conf_interval)-1
            #print(frame_num)
            #if frame_num%1==0:
                #print(frame_num, cont_line, t)
            #if cont_line>N+2:
                #cset1 = plt.imshow(velocity_grid, vmin=velmin, vmax=velmax, cmap=cm.Reds)
            ax = plt.gca()
            for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
                simplex = np.asarray(simplex)
                if np.all(simplex >= 0):
                    # Finite edge
                    ax.plot(vor.vertices[simplex, 0], vor.vertices[simplex, 1], 'b-')
            else:
                # Infinite edge — extend in direction
                i = simplex[simplex >= 0][0]  # Finite vertex
                t = vor.points[pointidx[1]] - vor.points[pointidx[0]]  # tangent
                t = t / np.linalg.norm(t)
                n = np.array([-t[1], t[0]])  # normal

                midpoint = vor.points[pointidx].mean(axis=0)
                far_point = vor.vertices[i] + n * 10  # Extend far for visibility

                ax.plot([vor.vertices[i, 0], far_point[0]],[vor.vertices[i, 1], far_point[1]], 'b--')



            ax.set_aspect('equal', adjustable='box')
            ax.set_xlim([0, lx])
            ax.set_ylim([165, 200])
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

    fig = plt.figure(figsize=(6,6))
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

    with open('order_parameters.txt', 'w') as f:
        for i in range(len(order_psi6)):
            print(time_conf[i]*dt, order_psi6[i], order_psiN[i], order_psiL[i], file=f)  


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
