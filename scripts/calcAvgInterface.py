import sys
import numpy as np
from math import *
import os.path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgba, to_hex, Normalize
import numpy as np
from scipy import ndimage
from skimage import measure
import pandas as pd


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



def wrap_angle(a):
    return (a + pi) % (2*pi) - pi

def unwrap_series(a):
    return np.unwrap(a)

def moving_average(x, w=5):
    if w <= 1:
        return x
    # pad edges with reflection to keep length
    pad = w//2
    xp = np.pad(x, pad, mode='reflect')
    y = np.convolve(xp, np.ones(w)/w, mode='valid')
    return y

# ----------------------------- Contour & orientation -----------------------------


def extract_zero_contours(f, smooth_sigma=0.8):
    """Compute contours of f at level 0. Returns list of contours in (row,col) coords.
    We optionally smooth f first to reduce grid noise.
    """
    if smooth_sigma and smooth_sigma > 0:
        f_s = ndimage.gaussian_filter(f, sigma=smooth_sigma)
    else:
        f_s = f
    contours = measure.find_contours(f_s, level=0.0)
    return contours

def pick_main_contour(contours, image_center=None):
    """Pick the main contour when multiple are present: prioritize longest, then closest to center."""
    if not contours:
        return None
    # sort by length descending
    contours_sorted = sorted(contours, key=lambda c: c.shape[0], reverse=True)
    if image_center is None:
        return contours_sorted[0]
    # check if top two are similar size; choose one closer to center
    c0 = contours_sorted[0]
    if len(contours_sorted) == 1:
        return c0
    c1 = contours_sorted[1]
    if abs(c0.shape[0] - c1.shape[0]) / max(1, c0.shape[0]) > 0.2:
        return c0
    # otherwise choose the one whose centroid is closer to image_center
    def centroid_dist(c):
        # convert to (x,y)
        xs = c[:,1]
        ys = c[:,0]
        cx = xs.mean()
        cy = ys.mean()
        return np.hypot(cx - image_center[0], cy - image_center[1])
    chosen = min(contours_sorted, key=centroid_dist)
    return chosen

def contour_to_xy(contour):
    """Convert skimage contour (row,col) to array Nx2 of (x,y) floats."""
    # contour is (row, col) == (y,x)
    return np.vstack([contour[:,1], contour[:,0]]).T


def principal_axis_angle(contour_xy):
    pts = contour_xy - contour_xy.mean(axis=0)
    # SVD
    U, S, Vt = np.linalg.svd(pts, full_matrices=False)
    v = Vt[0] # principal direction vector (x,y)
    theta = atan2(v[1], v[0])
    return wrap_angle(theta)

def tangent_average_angle(contour_xy):
    # compute segment tangents
    dx = np.diff(contour_xy[:,0])
    dy = np.diff(contour_xy[:,1])
    thetas = np.arctan2(dy, dx)
    # average on unit circle
    M = np.mean(np.exp(1j * thetas))
    theta = np.angle(M)
    return wrap_angle(theta)

def boundary_intersections(contour_xy, circle_center, circle_radius, tol_px=1.0):
    """Return intersection points (x,y) between contour and approximate circle boundary.
    Because contour points are continuous in general, we'll select contour points within tol_px of circle radius.
    """

    circle_center = contour_xy.mean(axis=0)
    dists = np.linalg.norm(contour_xy - circle_center, axis=1)
    circle_radius = dists.max()
    #if circle_center is None or circle_radius is None:
        #return np.empty((0,2))

    d = np.hypot(contour_xy[:,0] - circle_center[0], contour_xy[:,1] - circle_center[1])
    mask = np.abs(d - circle_radius) <= tol_px
    pts = contour_xy[mask]
    # if sparse sampling, cluster and choose unique intersection pairs
    if pts.shape[0] == 0:
        return pts
    # cluster by angular coordinate to merge nearby points
    angs = np.arctan2(pts[:,1] - circle_center[1], pts[:,0] - circle_center[0])
    # sort by angle and pick representative points for groups separated by >0.2 rad
    order = np.argsort(angs)
    angs_s = angs[order]
    pts_s = pts[order]
    clusters = [pts_s[0:1]]
    for i in range(1, len(pts_s)):
        if abs(angs_s[i] - angs_s[i-1]) > 0.2:
            clusters.append(pts_s[i:i+1])
        else:
            clusters[-1] = np.vstack([clusters[-1], pts_s[i:i+1]])
    reps = np.array([c.mean(axis=0) for c in clusters])
    return reps

def endpoints_line_angle(p1, p2):
    theta = atan2(p2[1] - p1[1], p2[0] - p1[0])
    return wrap_angle(theta)


# ----------------------------- Analysis pipeline -----------------------------

def extract_theta_timeseries(phi1_ts, phi2_ts, method='pca', dt=1.0,
                            smooth_phi_sigma=1.0, contour_smooth_sigma=0.8,
                            circle_center=None, circle_radius=None,
                            select_largest=True):
    nt = phi1_ts.shape[0]
    H = phi1_ts.shape[1]
    W = phi1_ts.shape[2]
    theta_wrapped = np.zeros(nt)
    contours_all = [None] * nt
    image_center = (W/2.0, H/2.0)
    eps=0.0001
    for t in range(nt):
        phi1 = phi1_ts[t]
        phi2 = phi2_ts[t]
        f = phi1 - phi2
        mask = (phi1 > eps) | (phi2 > eps)
        f_masked = np.where(mask, f, np.nan)
        #f_s = ndimage.gaussian_filter(f_masked, sigma=sigma)
        contours = extract_zero_contours(f_masked, smooth_sigma=contour_smooth_sigma)
        if not contours:
            theta_wrapped[t] = np.nan
            contours_all[t] = []
            continue
        # pick main
        if select_largest:
            main = pick_main_contour(contours, image_center=image_center)
        else:
            main = contours[0]
        contour_xy = contour_to_xy(main)
        theta = None
        if method == 'endpoints' and circle_center is not None and circle_radius is not None:
            pts = boundary_intersections(contour_xy, circle_center, circle_radius)
            if pts.shape[0] >= 2:
                # choose the two farthest apart cluster reps
                # if more than 2, pick pair with largest Euclidean distance
                if pts.shape[0] > 2:
                    pdist = np.sqrt(((pts[:,None,:] - pts[None,:,:])**2).sum(axis=2))
                    i,j = np.unravel_index(np.argmax(pdist), pdist.shape)
                    p1 = pts[i]; p2 = pts[j]
                else:
                    p1, p2 = pts[0], pts[1]
                theta = endpoints_line_angle(p1, p2)
                #print(theta*180/pi)
            else:
                # fallback to PCA
                theta = principal_axis_angle(contour_xy)
        elif method == 'tangent':
            try:
                theta = tangent_average_angle(contour_xy)
            except Exception:
                theta = principal_axis_angle(contour_xy)
        else:
            # default PCA
            theta = principal_axis_angle(contour_xy)
        theta_wrapped[t] = wrap_angle(theta)
        contours_all[t] = contour_xy
    # postprocess: interpolate nans
    t_idx = np.arange(nt)
    nans = np.isnan(theta_wrapped)
    if nans.any():
        good = ~nans
        if good.sum() >= 2:
            theta_wrapped[nans] = np.interp(t_idx[nans], t_idx[good], theta_wrapped[good])
        else:
            # too few good values
            pass
    theta_unwrapped = unwrap_series(theta_wrapped)
    return theta_wrapped, theta_unwrapped, contours_all

def angular_velocity(theta_unwrapped, dt=1.0):
    v = np.zeros_like(theta_unwrapped)
    # central differences
    v[1:-1] = (theta_unwrapped[2:] - theta_unwrapped[:-2]) / (2*dt)
    v[0] = (theta_unwrapped[1] - theta_unwrapped[0]) / dt
    v[-1] = (theta_unwrapped[-1] - theta_unwrapped[-2]) / dt
    return v


def autocorr(x, maxlag=None):
    x = np.asarray(x)
    n = len(x)
    if maxlag is None:
        maxlag = n//2
    x = x - np.nanmean(x)
    res = np.empty(maxlag+1)
    for lag in range(maxlag+1):
        if lag == 0:
            res[lag] = 1.0
            continue
        a = x[:-lag]
        b = x[lag:]
        valid = ~np.isnan(a) & ~np.isnan(b)
        if valid.sum() < 2:
            res[lag] = np.nan
        else:
            res[lag] = np.corrcoef(a[valid], b[valid])[0,1]
    return res


def fit_exponential_to_acf(acf, dt=1.0, fit_range=(1,100)):
    lags = np.arange(len(acf))
    lo, hi = fit_range
    mask = (lags >= lo) & (lags <= hi) & (~np.isnan(acf))
    if mask.sum() < 3:
        return np.nan, np.nan
    x = lags[mask] * dt
    y = acf[mask]
    # keep positive y only
    pos = y > 0
    if pos.sum() < 3:
        return np.nan, np.nan
    x = x[pos]
    y = y[pos]
    coef = np.polyfit(x, np.log(y), 1)
    slope, intercept = coef[0], coef[1]
    tau = -1.0/slope if slope != 0 else np.nan
    A = np.exp(intercept)
    return A, tau

def extract_runs_from_sign(sign_series):
    # sign_series values in {-1,0,1}
    runs = []
    s = np.asarray(sign_series)
    n = len(s)
    i = 0
    while i < n:
        if s[i] == 0:
            i += 1
            continue
        val = int(s[i])
        start = i
        i += 1
        while i < n and int(s[i]) == val:
            i += 1
        end = i-1
        runs.append({'start': start, 'end': end, 'length': end-start+1, 'direction': val})
    return pd.DataFrame(runs)


# -----------------------------
# Main function: visualize interface + orientation
# -----------------------------

def contour_principal_angle(contour):
    pts = contour - contour.mean(axis=0)
    U, S, Vt = np.linalg.svd(pts, full_matrices=False)
    v = Vt[0]
    return atan2(v[1], v[0]), v

def contour_tangent_average_angle(contour):
    dx = np.diff(contour[:,0])
    dy = np.diff(contour[:,1])
    thetas = np.arctan2(dy, dx)
    M = np.mean(np.exp(1j*thetas))
    return np.angle(M)


def endpoints_angle(contour, tol=1.0):
    # contour: Nx2 array (x,y)
    # center: (cx, cy)
    # radius: circle radius

    center = contour.mean(axis=0)
    dists = np.linalg.norm(contour - center, axis=1)
    radius = dists.max()
    p3=center

    cx, cy = center
    d = np.sqrt((contour[:,0]-cx)**2 + (contour[:,1]-cy)**2)
    pts = contour[np.abs(d - radius) < tol]
    if pts.shape[0] < 2:
        return None

    # shift coords relative to center
    rel = pts - np.array(center)

    # angle of each point relative to center
    angs = np.arctan2(rel[:,1], rel[:,0])

    # find two points roughly opposite (max angular separation)
    #i, j = np.unravel_index(np.argmax(np.abs(np.subtract.outer(angs, angs))), (len(pts), len(pts)))
    i, j = np.unravel_index(np.argmax(np.sum((pts[:,None,:]-pts[None,:,:])**2,axis=2)), (len(pts),len(pts)))
    if i == j:
        return None
    p1, p2 = pts[i], pts[j]

    theta = atan2(p2[1]-p1[1], p2[0]-p1[0])
    return wrap_angle(theta), radius, (p1, p2, p3)


def visualize_interface(phi1, phi2, method='pca', circle_center=None, circle_radius=None, sigma=1.0, eps=0.0001):
    f = phi1 - phi2
    mask = (phi1 > eps) | (phi2 > eps)
    f_masked = np.where(mask, f, np.nan)
    f_s = ndimage.gaussian_filter(f_masked, sigma=sigma)
    contours = measure.find_contours(f_s, 0.0)
    if not contours:
        raise RuntimeError("No interface found")
    contour = max(contours, key=len)  # largest
    contour = np.vstack([contour[:,1], contour[:,0]]).T  # (x,y)
    
    fig, ax = plt.subplots(figsize=(6,6))
    ax.imshow(f, cmap='RdBu', origin='lower')
    ax.plot(contour[:,0], contour[:,1], 'k-', lw=2)
    
    if method == 'pca':
        theta, v = contour_principal_angle(contour)
        center = contour.mean(axis=0)
        ax.quiver(center[0], center[1], np.cos(theta), np.sin(theta),
                  angles='xy', scale_units='xy', scale=5, color='lime', lw=2)
        ax.set_title(f"PCA orientation: {theta*180/pi:.1f}°")
    
    elif method == 'tangent':
        theta = contour_tangent_average_angle(contour)
        center = contour.mean(axis=0)
        ax.quiver(center[0], center[1], np.cos(theta), np.sin(theta),
                  angles='xy', scale_units='xy', scale=5, color='orange', lw=2)
        ax.set_title(f"Tangent-average orientation: {theta*180/pi:.1f}°")
    
    elif method == 'endpoints':
        if circle_center is None or circle_radius is None:
            raise ValueError("circle_center and circle_radius must be provided for endpoints method")
        res = endpoints_angle(contour)
        if res is None:
            ax.set_title("Endpoints not found on boundary")
        else:
            theta, circle_radius, (p1,p2, circle_center) = res
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], 'm--', lw=2)
            ax.set_title(f"Endpoints orientation: {theta*180/pi:.1f}°")
            ax.plot(*circle_center, 'mo')
            circ = plt.Circle(circle_center, circle_radius, fill=False, color='gray', ls=':')
            ax.add_patch(circ)
    
    ax.set_aspect('equal')
    plt.show()


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


cont_line=0
vmax=0.1
start_line=(N+2)*int(float(sys.argv[3])) 
cfile=open(trajectory_file,"r")
for i in range(start_line):
    cfile.readline()


fig = plt.figure(figsize=(6,6))
frame_num=int(t/print_conf_interval)-1

time_total = int(steps / print_conf_interval)
phi1 = np.zeros((time_total, ly, lx))
phi2 = np.zeros((time_total, ly, lx))


for line in cfile:
    cont_line+=1
    words=line.split()

    if words[0]=='t':
        t=int(float(words[2]))
        frame_num=int(t/print_conf_interval)-1

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
            if pt_num==0:
                phi1[frame_num][yy][xx] = value
            if pt_num==1:
                phi2[frame_num][yy][xx] = value

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

            cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], 2*D_major_axis_vec[0], 2*D_major_axis_vec[1], width=0.5, head_width=0, color='r')
            cset1 = plt.arrow(CoMX[pt_num], CoMY[pt_num], -2*D_major_axis_vec[0], -2*D_major_axis_vec[1], width=0.5, head_width=0, color='r')

        #increment phase field index
        pt_num+=1
        

    if cont_line%(N+2)==0:


            frame_num=int(t/print_conf_interval)-1
            #print(frame_num)
            #if frame_num%1==0:
                #print(frame_num, cont_line, t)
            #if cont_line>N+2:
                #cset1 = plt.imshow(velocity_grid, vmin=velmin, vmax=velmax, cmap=cm.Reds)
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

T = time_total 
H = ly
W = lx


method = 'endpoints'
omega_threshold = 0.1
theta_smooth_window = 100
omega_smooth_window = 100
contour_smooth_sigma = smooth_phi_sigma = 1


# default: circle centered at image center, radius min(H,W)/2 - 1
circle_center = (W/2.0, H/2.0)
#circle_radius = min(H, W)/2.0 - 1.0
circle_radius = (W - 2 * 6) / 2 -2 
print('Auto circle_center=', circle_center, 'radius=', circle_radius)

#for i in range(time_total):
pphi1 = phi1[4]
pphi2 = phi2[4]
visualize_interface(pphi1, pphi2, method=method, circle_radius = circle_radius, circle_center = circle_center)

print('Extracting theta(t) from phase fields using method=', method)
theta_wrapped, theta_unwrapped, contours = extract_theta_timeseries(phi1, phi2, method=method, dt=dt, 
                                                                    smooth_phi_sigma=smooth_phi_sigma,
                                                                    contour_smooth_sigma=contour_smooth_sigma,
                                                                    circle_center=circle_center,
                                                                    circle_radius=circle_radius)


if theta_smooth_window and theta_smooth_window > 1:
    theta_unwrapped_sm = moving_average(theta_unwrapped, w=theta_smooth_window)
else:
    theta_unwrapped_sm = theta_unwrapped

omega = angular_velocity(theta_unwrapped_sm, dt=dt)
if omega_smooth_window and omega_smooth_window > 1:
    omega_sm = moving_average(omega, w=omega_smooth_window)
else:
    omega_sm = omega

# sign series
sign = np.sign(omega_sm)
sign[np.abs(omega_sm) < omega_threshold] = 0

runs_df = extract_runs_from_sign(sign)
#print(runs_df['length'])
total_runs = 0
time_rot = []
for i in runs_df['length']:
    total_runs += i
    time_rot.append(i/time_total)
#print(total_runs/time_total)
print(time_rot)

if not runs_df.empty:
    runs_df['time_start'] = runs_df['start'] * dt
    runs_df['time_end'] = runs_df['end'] * dt
    runs_df['duration'] = runs_df['length'] * dt


cumrot = theta_unwrapped_sm - theta_unwrapped_sm[0]
net_rot_rad = cumrot[-1]
net_rot_turns = net_rot_rad / (2*pi)


frac_ccw = np.mean(sign == 1)
frac_cw = np.mean(sign == -1)
n_switches = int(np.sum(np.abs(np.diff(sign)) > 0))
mean_abs_speed = np.nanmean(np.abs(omega_sm))
median_abs_speed = np.nanmedian(np.abs(omega_sm))


# autocorrelation and fit
maxlag = min(500, len(omega_sm)//2)
acf = autocorr(omega_sm, maxlag=maxlag)
A_fit, tau_fit = fit_exponential_to_acf(acf, dt=dt, fit_range=(1, min(200, maxlag)))

# Plots
times = np.arange(T) * dt
times1 = np.arange(T+1) * dt
times2 = np.arange(T+2) * dt
plt.figure(figsize=(10, 6))
plt.subplot(3,1,1)
plt.plot(times, wrap_angle(theta_wrapped), label='theta wrapped')
plt.plot(times1, np.mod(theta_unwrapped_sm, 2*pi)-2*pi, label='theta sm (mod 2pi)')
plt.ylabel('angle (rad)')
plt.legend()


plt.subplot(3,1,2)
plt.plot(times1, omega, label='omega (raw)')
plt.plot(times2, omega_sm, label='omega sm')
plt.axhline(omega_threshold, linestyle='--', alpha=0.6)
plt.axhline(-omega_threshold, linestyle='--', alpha=0.6)
plt.ylabel('angular velocity')
plt.legend()

plt.subplot(3,1,3)
plt.plot(times1, cumrot/(2*pi), label='cumulative turns')
plt.ylabel('turns')
plt.xlabel('time')
plt.legend()
plt.tight_layout()
#fig1 = os.path.join(outdir, 'theta_and_omega.png')
plt.show()
plt.close()
#print('Wrote', fig1)


# ACF
lags = np.arange(len(acf)) * dt
plt.figure(figsize=(6,4))
plt.plot(lags, acf, label='ACF(omega)')
if not np.isnan(tau_fit):
    plt.plot(lags, A_fit * np.exp(-lags / tau_fit), '--', label=f'exp fit tau={tau_fit:.2f}')
plt.xlabel('lag')
plt.ylabel('autocorrelation')
plt.legend()
#fig2 = os.path.join(outdir, 'acf_omega.png')
plt.tight_layout()
plt.show()
plt.close()
#print('Wrote', fig2)


# histogram of run durations
if not runs_df.empty:
    plt.figure(figsize=(6,3))
    plt.hist(runs_df['duration'], bins=30)
    plt.xlabel('run duration')
    plt.ylabel('count')
    #fig3 = os.path.join(outdir, 'run_duration_hist.png')
    plt.tight_layout()
    plt.show()
    plt.close()
    #print('Wrote', fig3)



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

    with open('time_interface_rotation.txt', 'w') as f:
        for i in range(len(time_rot)):
            print(time_rot[i], file=f)  

    '''
    with open('theta_shape.txt', 'w') as f:
        for i in range(total_time_frames):
            print_str = ''
            for j in range(n_rows):
                print_str += str(theta_time[j][i])
                print_str += ' '

            print(print_str , file=f)  
    '''

print('done')
