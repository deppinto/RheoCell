# This file is part of CELADRO, Copyright (C) 2016-17, Romain Mueller
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib import colors
from math import sqrt, pi, atan2, cos, sin
from matplotlib.colors import LinearSegmentedColormap, to_rgba, to_hex, Normalize
from scipy import ndimage
from itertools import product
from scipy.spatial import Voronoi,voronoi_plot_2d
from scipy.spatial import ConvexHull
import matplotlib.ticker as ticker
# import seaborn

def change_object_to_float(q, Lx, Ly):

    p = np.zeros(q.shape)
    #print(q.shape)
    #print(q)
    #exit (1)
    
    for i in range(Lx):
        for j in range(Ly):
            p[i,j] = float(q[i,j])

    return p


def grad(phi):
    # returns the gradient of phase field grad_phi

    fx = np.zeros(phi.shape)
    fy = np.zeros(phi.shape)

    fx1 = np.gradient(phi,axis=0)
    fy1 = np.gradient(phi,axis=1)

    fx = np.where(abs(fx)>abs(fx1),abs(fx),abs(fx1))
    fy = np.where(abs(fy)>abs(fy1),abs(fy),abs(fy1))

    return fx,fy

def grad_phi(frame, engine=plt,cbar = True,avg=1,width = 1,scale = None):
    '''
    plot grad_phi for individual cells
    '''
    fx = np.zeros(frame.phi[0].shape)
    fy = np.zeros(frame.phi[0].shape)
    for i in range(len(frame.phi)):
        fx1 = np.gradient(frame.phi[i],axis=0)
        fy1 = np.gradient(frame.phi[i],axis=1)
        fx = np.where(abs(fx)>abs(fx1),abs(fx),abs(fx1))
        fy = np.where(abs(fy)>abs(fy1),abs(fy),abs(fy1))

    m = np.sqrt(fx**2 + fy**2)
    im = engine.imshow(m.T, interpolation='lanczos', cmap='GnBu',
                            origin='lower')
    if cbar:
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(engine)
        ax_cb = divider.new_horizontal(size="3%", pad=0.03)
        fig = engine.get_figure()
        fig.add_axes(ax_cb)
        plt.colorbar(im, cax=ax_cb)

def _get_field(phases, vals, size=1, mode='wrap'):
    """
    Compute the coarse grained field from a collection of phase-fields and
    associated values: ret[i] = sum_j phases[j]*values[i, j].

    Args:
        phases: List of phase-fields.
        vals: List of lists of size (None, len(phases)) of values to be
            associated with each phase-field.
        size: Coarse-graining size.
        mode: How to treat the boundaries, see
            scipy.ndimage.filters.uniform_filter.

    Returns:
        A list of fields, each corresponding to the individual values.
    """
    ret = []

    for vlist in vals:
        # print(vlist)
        assert len(vlist) == len(phases)
        field = np.zeros(phases[0].shape)
        for n in range(len(phases)):
            field += vlist[n]*phases[n]
        #field = ndimage.filters.uniform_filter(field, size=size, mode=mode)
        # field = ndimage.filters.gaussian_filter(field, sigma=size, mode=mode)
        ret.append(field)
    return ret

def _get_density_sort(phases, zetas, activity, size=1, mode='wrap'):
    """
    Compute the coarse grained field from a collection of phase-fields and
    associated values: ret[i] = sum_j phases[j], separated according to activity

    Args:
        phases: List of phase-fields.
        vals: List of lists of size (None, len(phases)) of values to be
            associated with each phase-field.
        size: Coarse-graining size.
        mode: How to treat the boundaries, see
            scipy.ndimage.filters.uniform_filter.

    Returns:
        A list of fields, each corresponding to the individual values.
    """
    ret = []

    for z in zetas:
        field = np.zeros(phases[0].shape)
        for n in range(len(phases)):
            # print(activity[n])
            if activity[n]==z:
                field += phases[n]
        #field = ndimage.filters.uniform_filter(field, size=size, mode=mode)
        # field = ndimage.filters.gaussian_filter(field, sigma=size, mode=mode)
        ret.append(field)
    return ret

def get_grad_phi(phase, norm, size=1, mode='wrap'):

    (LX,LY) = phase.shape
    ret = np.zeros(phase.shape)

    for i in range(LX):
        start = i*LY
        for j in range(LY):
            ret[i,j] = norm[start+j]

    return ret

def get_chemical_field(phases, chem, size=1, mode='wrap'):

    (LX,LY) = phases[0].shape
    # print(LX)
    # print(LY)
    ret = np.zeros(phases[0].shape)

    for i in range(LX):
        start = i*LY
        for j in range(LY):
            ret[i,j] = chem[start+j]

    return ret

def get_repulsive_field(c, phases, area, rep_x, rep_y, Fdp_tot, size=1, mode='wrap'):
    """
    Compute coarse-grained nematic field.

    Args:
        phases: List of phase fields.
        vel: List of 2d velocities associated with each phase field.
        size: Coarse-graining size.
        mode: How to treat the boundaries, see
            scipy.ndimage.filters.uniform_filter.
    """
    (LX,LY) = phases[0].shape
    frep_x = np.zeros(phases[0].shape)
    frep_y = np.zeros(phases[0].shape)

    charge = 0

    for i in range(LX):
        start = i*LY
        for j in range(LY):
            disp_x = i - c[0]
            disp_y = j - c[1]
            com_angle = atan2(disp_y,disp_x)
            r2 = disp_x**2 + disp_y**2
            #component of the active force in direction of com
            magn = rep_x[start+j]*cos(com_angle) + rep_y[start+j]*sin(com_angle)
            #magn = sqrt((rep_x[start+j])**2+(rep_y[start+j])**2)
            charge += magn/area[0]

    #print(charge)

    for i in range(LX):
        start = i*LY
        for j in range(LY):
            disp_x = i - c[0]
            disp_y = j - c[1]
            com_angle = atan2(disp_y,disp_x)
            r2 = disp_x**2 + disp_y**2

            frep_x[i,j] = rep_x[start+j]# - Fdp_tot[0][0]*(1-phases[0][i][j])
            frep_y[i,j] = rep_y[start+j]# - Fdp_tot[0][1]*(1-phases[0][i][j])
            #if(r2>25 and r2<81):
            #    frep_x[i,j] += -charge*cos(com_angle)#/(r2)
            #    frep_y[i,j] += -charge*sin(com_angle)#/(r2)

    return frep_x,frep_y

def get_avg_force_field(phases, vel, xi, area, size=1, mode='wrap'):
    """
    Compute coarse-grained velocity field = sum_j v_com_j*xi_j*phi_j

    Args:
        phases: List of phase fields.
        vel: List of 2d velocities associated with each phase field.
        size: Coarse-graining size.
        mode: How to treat the boundaries, see
            scipy.ndimage.filters.uniform_filter.
    """

    v0 = []
    v1 = []

    # 'distribute' in a sense the COM velocity across the entire cell

    for i in range(len(xi)):
        v0.append(vel[i][0]*xi[i]/area[i])
        v1.append(vel[i][1]*xi[i]/area[i])

    return _get_field(phases, [v0, v1], size, mode)

def get_velocity_field(phases, vel, size=1, mode='wrap'):
    """
    Compute coarse-grained velocity field.

    Args:
        phases: List of phase fields.
        vel: List of 2d velocities associated with each phase field.
        size: Coarse-graining size.
        mode: How to treat the boundaries, see
            scipy.ndimage.filters.uniform_filter.
    """
    v0 = [v[0] for v in vel]
    v1 = [v[1] for v in vel]
    return _get_field(phases, [v0, v1], size, mode)

def get_polarity_field(phases, pol, size=1, mode='wrap'):
    """
    Compute coarse-grained nematic field.

    Args:
        phases: List of phase fields.
        vel: List of 2d velocities associated with each phase field.
        size: Coarse-graining size.
        mode: How to treat the boundaries, see
            scipy.ndimage.filters.uniform_filter.
    """
    p0 = [p[0] for p in pol]
    p1 = [p[1] for p in pol]
    return _get_field(phases, [p0, p1], size, mode)

def get_p_atic_field(phases, Gx, Gy, size=1, mode='wrap'):
    '''
    no different to get_nematic_field since you'll be choosing which p to use according to the order parameter you pass as an argument
    '''

    return _get_field(phases, [Gx, Gy], size, mode)

def get_nematic_field(phases, qxx, qxy, size=1, mode='wrap'):
    """
    Compute coarse-grained nematic field.

    Args:
        phases: List of phase fields.
        qxx, qxy: Components of the nematic field of the individual phase
            fields.
        size: Coarse-graining size.
        mode: How to treat the boundaries, see
            scipy.ndimage.filters.uniform_filter.
    """
    return _get_field(phases, [qxx, qxy], size, mode)

def get_vorticity_field(ux, uy, pbc=True):
    """
    Compute vorticity field from velocity field

    Args:
        ux, uy: the individual components of the velocity field.
        pbc: How to treat boundaries, set to true if using pbc.

    Returns:
        Vorticity field.
    """
    if pbc:
        dxuy = np.gradient(np.pad(uy, 2, mode='wrap'), axis=0)[2:-2, 2:-2]
        dyux = np.gradient(np.pad(ux, 2, mode='wrap'), axis=1)[2:-2, 2:-2]
    else:
        dxuy = np.gradient(uy, axis=0)
        dyux = np.gradient(ux, axis=1)
    return dxuy - dyux

def get_gradient_field(ux, uy, pbc=True):
    """
    Compute gradient field from velocity field

    Args:
        ux, uy: the individual components of the velocity field.
        pbc: How to treat boundaries, set to true if using pbc.

    Returns:
        Gradient field.
    """
    if pbc:
        dxux = np.gradient(np.pad(ux, 2, mode='wrap'), axis=0)[2:-2, 2:-2]
        dyuy = np.gradient(np.pad(uy, 2, mode='wrap'), axis=1)[2:-2, 2:-2]
    else:
        dxux = np.gradient(ux, axis=0)
        dyuy = np.gradient(uy, axis=1)
    return dxux + dyuy

def autocorr(u):
    '''
    computes the autocorrelation of a scalar time series.

    Arguments:
        u: the scalar time series

    Returns:
        the autocorrelation of u as a time s
    '''

    T = len(u)
    corr = np.zeros(T)

    corr[0] = np.sum(np.multiply(u,u))/T

    for i in range(1,T):
        # can't have i=0 because we get a shape error
        # print(u[i:].size)
        # print(u[:-i].size)
        total = np.sum(np.multiply(u[i:],u[:-1*i]))
        corr[i] = total/(T-i)

    corr = corr/corr[0]

    return corr

def get_corr(u):
    """
    Compute the cross-correlation (as a function of distance) of a real two-
    dimensional scalar field.

    Arguments:
        u: The scalar field.

    Returns:
        The cross-correlation of u as an array.
    """
    # get 2d correlations
    c = np.fft.rfft2(u)
    c = np.fft.irfft2(np.multiply(c, np.conj(c)))
    # go to polar coords
    s = int(sqrt(c.size)/2)
    r = np.zeros(s)
    n = np.zeros(s)
    k = 0
    for (i, j), v in np.ndenumerate(c):
        k = int(sqrt(i**2 + j**2))
        if k >= s:
            continue
        r[k] += v
        n[k] += 1
    r = np.divide(r, n)
    r /= r[0]
    return r

def get_corr2(ux, uy):
    """
    Compute the correlation (as a function of distance) of two real two-
    dimensional scalar fields.

    Arguments:
        ux, uy: The scalar fields.

    Returns:
        The correlation of ux and uy as an array.
    """
    # get 2d correlations
    cx = np.fft.rfft2(ux)
    cx = np.fft.irfft2(np.multiply(cx, np.conj(cx)))
    cy = np.fft.rfft2(uy)
    cy = np.fft.irfft2(np.multiply(cy, np.conj(cy)))
    c = cx + cy
    # go to polar coords
    s = int(sqrt(c.size)/2)
    r = np.zeros(s)
    n = np.zeros(s)
    k = 0
    for (i, j), v in np.ndenumerate(c):
        k = int(sqrt(i**2 + j**2))
        if k >= s:
            continue
        r[k] += v
        n[k] += 1
    r = np.divide(r, n)
    r /= r[0]
    return r

def first_min(u):
    """
    compute the location and value of the first minimum of a set of data
    """

    index = 0
    value = u[0]
    
    for i in range(len(u)):
        if u[i] < value:
            index = i
            value = u[i]
        elif u[i] >= value:
            continue

    return index,value

def eqt_corr(u):
    """
    Compute the equal-time correlation function of u
    Typically density

    <u(R)u(R+r)>-<u(R)><u(R+r)>

    Arguments:
        u: The scalar field.

    Returns:
        The equal-time correlation of u with itself as an array.
    """

    Lx,Ly = u.shape
    N = Lx*Ly
    tot = u.sum()
    avg = tot/N

    # get 2d correlations
    c = np.fft.rfft2(u)
    d = np.zeros(c.shape,dtype=np.complex_)
    e = np.fft.rfft2(np.ones(u.shape))
    d[0,0] = e[0,0]

    d = np.fft.irfft2(np.multiply(e,np.conj(c)))
    c = np.fft.irfft2(np.multiply(c, np.conj(c)))
    #

    # go to polar coords
    s = int(sqrt(c.size)/2)
    rc = np.zeros(s)
    nc = np.zeros(s)
    rd = np.zeros(s)
    nd = np.zeros(s)

    k = 0

    for (i, j), v in np.ndenumerate(c):
        k = int(sqrt(i**2 + j**2))
        if k >= s:
            continue
        rc[k] += v
        nc[k] += 1

    for (i, j), v in np.ndenumerate(d):
        k = int(sqrt(i**2 + j**2))
        if k >= s:
            continue
        rd[k] += avg*v
        nd[k] += 1

    rc = np.divide(rc, nc)
    rd = np.divide(rd, nd)

    r = np.subtract(rc,rd)

    r /= r[0]
    return r

def get_stress(Tx,Ty,nu=0):

    fTx = np.fft.rfft2(Tx)
    fTy = np.fft.rfft2(Ty)

    nx,ny = fTx.shape

    # print(fTx.shape)
    # print(fTy.shape)

    kx,ky = np.mgrid[:nx,:ny]
    # ky,kx = np.mgrid[:nx,:ny]

    ky = np.where(ky==0,1e-8,ky)
    kx = np.where(kx==0,1e-8,kx)

    # print(kx)
    # print(ky)

    k2 = kx*kx + ky*ky
    k4 = k2*k2

    fsigma_xy = (-1.j)*(ky*fTx+kx*fTy)/k2 + (1.j)*(1+nu)*kx*ky*(kx*fTx+ky*fTy)/k4
    # fsigma_xy = np.nan_to_num(fsigma_xy)

    fsigma_xx = (-1.j)*fTx/kx - ky*fsigma_xy/kx
    # fsigma_xx = np.nan_to_num(fsigma_xy)
    fsigma_yy = (-1.j)*fTy/ky - kx*fsigma_xy/ky
    # fsigma_yy = np.nan_to_num(fsigma_xy)

    sigma_xx = np.fft.irfft2(fsigma_xx)
    sigma_xy = np.fft.irfft2(fsigma_xy)
    sigma_yy = np.fft.irfft2(fsigma_yy)

    return sigma_xx, sigma_yy, sigma_xy

def nematic_corr(frame,engine = plt,size = 2,show = True):
    # get field
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    (Qxx, Qxy) = get_nematic_field(frame.phi, frame.Q00, frame.Q01,
                                   size=size, mode=mode)
    Qxx *= (1.-frame.parameters['walls'])
    Qxy *= (1.-frame.parameters['walls'])
    (LX,LY) = Qxx.shape
    corr = get_corr2(Qxx, Qxy)
    radius = int(sqrt(LX*LY)/2)
    if show:
        engine.plot(range(radius), corr,'r')
    return corr

def shape_corr(frame,engine = plt,size = 2,show = True):
    # get field
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    (Qxx, Qxy) = get_nematic_field(frame.phi, frame.S00, frame.S01,
                                   size=size, mode=mode)
    Qxx *= (1.-frame.parameters['walls'])
    Qxy *= (1.-frame.parameters['walls'])
    (LX,LY) = Qxx.shape
    corr = get_corr2(Qxx, Qxy)
    radius = int(sqrt(LX*LY)/2)
    if show:
        engine.plot(range(radius), corr,'g')
    return corr

def velocity_corr(frame,engine = plt,size =1,show = True):
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    vx, vy = get_velocity_field(frame.phi, frame.com_velocity, size, mode=mode)
    vx *= (1.-frame.parameters['walls'])
    vy *= (1.-frame.parameters['walls'])
    corr = get_corr2(vx, vy)
    (LX,LY) = vx.shape
    radius = int(sqrt(LX*LY)/2)
    if show:
        engine.plot(range(radius), corr,'k')
    return corr

def vorticity_corr(frame,engine = plt,size = 1,show = True):
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    vx, vy = get_velocity_field(frame.phi, frame.com_velocity, size, mode=mode)
    vx *= (1.-frame.parameters['walls'])
    vy *= (1.-frame.parameters['walls'])
    w = get_vorticity_field(vx, vy)
    
    w = ndimage.filters.uniform_filter(w, size=size, mode=mode)
    corr = get_corr(w)
    (LX,LY) = vx.shape
    radius = int(sqrt(LX*LY)/2)
    if show:
        engine.plot(range(radius), corr,'r')
    return corr

def _charge_array(vx,vy,p=2):
    """
    Compute the charge array of vector field (vx,vy)

    Args:
        vx: x component of vector field
        vy: y component of vector field
        Type: nematic -- (vx,vy) is head-tail symmetrical
              polar-- (vx,vy) is not head-tail symmetical
    Returns:
        Field of the charge distribution with the same shape as vx and vy
    """
    # compute angle
    # def wang(a, b):
    #     """Infamous chinese function"""
    #     '''Type = 'nematic' or 'polar' '''
    #     ang = atan2(abs(a[0]*b[1]-a[1]*b[0]), a[0]*b[0]+a[1]*b[1])
    #     if (Type == 'nematic') and (ang > pi/2.):
    #         b = [-i for i in b]
    #     m = a[0]*b[1]-a[1]*b[0]
    #     return -np.sign(m)*atan2(abs(m), a[0]*b[0]+a[1]*b[1])

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
        ang = atan2(abs(a[0]*b[1]-a[1]*b[0]), a[0]*b[0]+a[1]*b[1])

        while(abs(ang) > np.pi/p + 1e-3):
            b = rotate(b,p)
            ang = atan2(abs(a[0]*b[1]-a[1]*b[0]), a[0]*b[0]+a[1]*b[1])

        m = a[0]*b[1]-a[1]*b[0]
        return -np.sign(m)*atan2(abs(m), a[0]*b[0]+a[1]*b[1])

    (LX, LY) = vx.shape
    w = np.zeros((LX,LY))
    # This mysterious part was stolen from Amin's code.
    # (calculate the winding number)
    for i in range(LX):
        for j in range(LY):
            ax1 = [vx[(i+1) % LX, j],
                   vy[(i+1) % LX, j]]
            ax2 = [vx[(i-1+LX) % LX, j],
                   vy[(i-1+LX) % LX, j]]
            ax3 = [vx[i, (j-1+LY) % LY],
                   vy[i, (j-1+LY) % LY]]
            ax4 = [vx[i, (j+1) % LY],
                   vy[i, (j+1) % LY]]
            ax5 = [vx[(i+1) % LX, (j-1+LY) % LY],
                   vy[(i+1) % LX, (j-1+LY) % LY]]
            ax6 = [vx[(i-1+LX) % LX, (j-1+LY) % LY],
                   vy[(i-1+LX) % LX, (j-1+LY) % LY]]
            ax7 = [vx[(i+1) % LX, (j+1) % LY],
                   vy[(i+1) % LX, (j+1) % LY]]
            ax8 = [vx[(i-1+LX) % LX, (j+1) % LY],
                   vy[(i-1+LX) % LX, (j+1) % LY]]

            w[i, j] = wang(ax1, ax5)
            w[i, j] += wang(ax5, ax3)
            w[i, j] += wang(ax3, ax6)
            w[i, j] += wang(ax6, ax2)
            w[i, j] += wang(ax2, ax8)
            w[i, j] += wang(ax8, ax4)
            w[i, j] += wang(ax4, ax7)
            w[i, j] += wang(ax7, ax1)
            w[i, j] /= 2.*pi

    return w

def G_charge_array(G00,G01,p=2):
    '''
    Compute the charge array associated with a Q-tensor field. The defects
    then show up as small regions of non-zero charge (typically 3x3).
    Args:
        G00: xx component of  G-tensor field
        G01: xy component of  G-tensor field
        p: symmetry associated wit the xx, xy components of the tensor field
    Returns:
        Field of the charge distribution with the same shape as G00 and G01
    '''
    (LX, LY) = G00.shape
    w = np.zeros((LX, LY))
    # get shape and init charge array
    # we use the director field

    theta = np.arctan2(G01,G00)/p
    theta[np.where(theta<1e-3)] += 2*np.pi/p

    nx = np.vectorize(cos)(theta)
    ny = np.vectorize(sin)(theta)

    w = _charge_array(nx,ny,p)
    return w

def Q_charge_array(Q00,Q01):
    '''
    Compute the charge array associated with a Q-tensor field. The defects
    then show up as small regions of non-zero charge (typically 3x3).
    Args:
        Q00: xx component of  Q-tensor field
        Q01: xy component of  Q-tensor field
    Returns:
        Field of the charge distribution with the same shape as Qxx and Qxy

    deprecated by G_charge_array
    '''
    (LX, LY) = Q00.shape
    w = np.zeros((LX, LY))
    # get shape and init charge array
    # we use the director field instead of Q
    S = np.vectorize(sqrt)(Q00**2 + Q01**2)
    # here vx and vy represents nematic director
    nx = np.sqrt(2.0*S)*np.vectorize(sqrt)((1 + Q00/S)/2)
    ny = np.sqrt(2.0*S)*np.sign(Q01)*np.vectorize(sqrt)((1 - Q00/S)/2)
    w = _charge_array(nx,ny,p=2)
    return w

def polar_charge_array(px,py):
    '''
    Compute the charge array associated with a polarity field. The defects
    then show up as small regions of non-zero charge (typically 3x3).
    Args:
        px: x component of polarisation field
        py: y component of polarisation field
    Returns:
        Field of the charge distribution with the same shape as px and py

    deprecated by G_charge_array
    '''
    return _charge_array(px,py,p=1)

def get_p_atic_defects(w, G00, G01, p=2, cal_angle = False):
    """
    Returns list of defects from charge array.

    Args:
        w: Charge array.
        G00: xx component of p-atic tensor
        G01: xy component of p-atic tensor

    Returns:
        List of the form [ [ (x, y), charge] ].
    """
    # defects show up as 2x2 regions in the charge array w and must be
    # collapsed to a single point by taking the average position of
    # neighbouring points with the same charge (by breath first search).

    # bfs recursive function
    def collapse(i, j, s, x=0, y=0, n=0,rng = [0.4,0.6]):
        (LX,LY) = w.shape
        if (s*w[i, j] > rng[0]) and (s*w[i, j] < rng[1]):
            w[i, j] = 0
            x1,y1,n1 = collapse((i+1) % LX, j, s, x, y, n,rng)
            x2,y2,n2 = collapse((i-1+LX) % LX, j, s, x, y, n, rng)
            x3,y3,n3 = collapse(i, (j+1) % LY, s, x, y, n, rng)
            x4,y4,n4 = collapse(i, (j-1+LY) % LY, s, x, y, n, rng)
            x = i + x1 + x2 +x3 +x4
            y = j + y1 + y2 +y3 +y4
            n = 1 + n1 + n2 +n3 +n4
            return x,y,n
        else:
            return 0,0,0

    (LX, LY) = w.shape
    d = []

    charge = 1.0/p
    thresh = 0.1

    for i in range(LX):
        for j in range(LY):
            # detect simplest charge 1/p defects
            if  (abs(w[i, j]) > charge - thresh) and (abs(w[i,j])< charge + thresh):
                # charge sign
                s = np.sign(w[i, j])
                # bfs
                sum_x, sum_y, n = collapse(i, j, s,rng = [charge - thresh, charge + thresh])
                x,y = sum_x/n,sum_y/n
                # compute angle, see doi:10.1039/c6sm01146b
                if cal_angle:
                    num = 0
                    den = 0
                    for (dx, dy) in [(0, 0), (0, 1), (1, 1), (1, 0)]:
                        # coordinates of nodes around the defect
                        kk = (int(x) + LX + dx) % LX
                        ll = (int(y) + LY + dy) % LY
                        # derivative at these points
                        dxQxx = .5*(G00[(kk+1) % LX, ll] - G00[(kk-1+LX) % LX, ll])
                        dxQxy = .5*(G01[(kk+1) % LX, ll] - G01[(kk-1+LX) % LX, ll])
                        dyQxx = .5*(G00[kk, (ll+1) % LY] - G00[kk, (ll-1+LY) % LY])
                        dyQxy = .5*(G01[kk, (ll+1) % LY] - G01[kk, (ll-1+LY) % LY])
                        # accumulate numerator and denominator
                        num += s*dxQxy - dyQxx
                        den += dxQxx + s*dyQxy
                    psi = s/(2.-s)*atan2(num, den)
                else:
                    psi = 0

                # add defect to list
                d.append({"pos": np.array([x, y]),
                          "charge": charge*s,
                          "angle": psi})

            # keep this just in case our other symmetries give us integer defects
            elif (abs(w[i, j]) > 0.9) and (abs(w[i,j])<1.1):
    		# charge sign
                s = np.sign(w[i, j])
                # bfs
                sum_x, sum_y, n = collapse(i, j, s,rng = [0.9,1.1])
                x,y = sum_x/n,sum_y/n
                # add defect to list
                d.append({"pos": np.array([x, y]),
                          "charge": 1*s,
                          "angle": 0})

    # print(d)
    return d

def get_defects(w, Qxx = 0, Qxy = 0, cal_angle = False):
    """
    Returns list of defects from charge array.

    Args:
        w: Charge array.
        if Type == 'nematic' then Qxx,Qxy is useful
        if Type == 'polar' then Qxx,Qxy can be omitted

    Returns:
        List of the form [ [ (x, y), charge] ].

    deprecated by get_p_atic_defects
    """
    # defects show up as 2x2 regions in the charge array w and must be
    # collapsed to a single point by taking the average position of
    # neighbouring points with the same charge (by breath first search).

    # bfs recursive function
    def collapse(i, j, s, x=0, y=0, n=0,rng = [0.4,0.6]):
        (LX,LY) = w.shape
        if (s*w[i, j] > rng[0]) and (s*w[i, j] < rng[1]):
            w[i, j] = 0
            x1,y1,n1 = collapse((i+1) % LX, j, s, x, y, n,rng)
            x2,y2,n2 = collapse((i-1+LX) % LX, j, s, x, y, n, rng)
            x3,y3,n3 = collapse(i, (j+1) % LY, s, x, y, n, rng)
            x4,y4,n4 = collapse(i, (j-1+LY) % LY, s, x, y, n, rng)
            x = i + x1 + x2 +x3 +x4
            y = j + y1 + y2 +y3 +y4
            n = 1 + n1 + n2 +n3 +n4
            return x,y,n
        else:
            return 0,0,0

    (LX, LY) = w.shape
    d = []

    for i in range(LX):
        for j in range(LY):
            # detect
            if  (abs(w[i, j]) > 0.4) and (abs(w[i,j])<0.6):
                # charge sign
                s = np.sign(w[i, j])
                # bfs
                sum_x, sum_y, n = collapse(i, j, s,rng = [0.4,0.6])
                x,y = sum_x/n,sum_y/n
                # compute angle, see doi:10.1039/c6sm01146b
                if cal_angle:
                    num = 0
                    den = 0
                    for (dx, dy) in [(0, 0), (0, 1), (1, 1), (1, 0)]:
                        # coordinates of nodes around the defect
                        kk = (int(x) + LX + dx) % LX
                        ll = (int(y) + LY + dy) % LY
                        # derivative at these points
                        dxQxx = .5*(Qxx[(kk+1) % LX, ll] - Qxx[(kk-1+LX) % LX, ll])
                        dxQxy = .5*(Qxy[(kk+1) % LX, ll] - Qxy[(kk-1+LX) % LX, ll])
                        dyQxx = .5*(Qxx[kk, (ll+1) % LY] - Qxx[kk, (ll-1+LY) % LY])
                        dyQxy = .5*(Qxy[kk, (ll+1) % LY] - Qxy[kk, (ll-1+LY) % LY])
                        # accumulate numerator and denominator
                        num += s*dxQxy - dyQxx
                        den += dxQxx + s*dyQxy
                    psi = s/(2.-s)*atan2(num, den)
                else:
                    psi = 0

                # add defect to list
                d.append({"pos": np.array([x, y]),
                          "charge": .5*s,
                          "angle": psi})
            elif (abs(w[i, j]) > 0.9) and (abs(w[i,j])<1.1):
		# charge sign
                s = np.sign(w[i, j])
                # bfs
                sum_x, sum_y, n = collapse(i, j, s,rng = [0.9,1.1])
                x,y = sum_x/n,sum_y/n
                # add defect to list
                d.append({"pos": np.array([x, y]),
                          "charge": 1*s,
                          "angle": 0})

    return d

def G_defects(G00, G01, p=2, engine=plt, arrow_len=0):
    """
    Plot defects of the p-atic field G.

    Args:
        G00, G01: Components of the p-atic field.
        p: symmetry (typically nematic, triatic, tetratic, hexatic)
        engine: Plotting engine or axis.
        arrow_len: If non-zero plot speed of defects as well.
    """
    w = G_charge_array(G00, G01, p)
    defects = get_p_atic_defects(w, G00, G01, p)
    for d in defects:
        if d['charge'] == 1/p:
            engine.plot(d["pos"][0], d["pos"][1], 'go', markersize=4)
            # plot direction of pos defects
            if not arrow_len == 0:
                engine.arrow(d['pos'][0], d['pos'][1],
                             -arrow_len*cos(d['angle']),
                             -arrow_len*sin(d['angle']),
                             color='r', head_width=3, head_length=3)
        elif d['charge'] == -1/p:
            engine.plot(d["pos"][0], d["pos"][1], 'b^', markersize=4)
        elif d['charge'] == 1:
            engine.plot(d["pos"][0], d["pos"][1], 'r*', markersize=4)
        elif d['charge'] == -1:
            engine.plot(d["pos"][0], d["pos"][1], 'kX', markersize=4)

def calculate_defects(G00, G01, p=2, arrow_len=0, cal_angle=False):
    """
    Calculate defects of the p-atic field G.

    Args:
        G00, G01: Components of the p-atic field.
        p: symmetry (typically nematic, triatic, tetratic, hexatic)
        engine: Plotting engine or axis.
        arrow_len: If non-zero plot speed of defects as well.

    Returns:
        list of defects with charge, position, and angle (if specified)
    """

    w = G_charge_array(G00, G01, p)
    return get_p_atic_defects(w, G00, G01, p, cal_angle)

def plot_defects(defects, engine=plt, p=2, arrow_len=0):
    """
    Plot defects of p-fold orientational symmetry, calculated elsewhere.

    Args:
        G00, G01: Components of the p-atic field.
        p: symmetry (typically nematic, triatic, tetratic, hexatic)
        engine: Plotting engine or axis.
        arrow_len: If non-zero plot speed of defects as well.
    """
    for d in defects:
        if d['charge'] == 1/p:
            engine.plot(d["pos"][0], d["pos"][1], 'go', markersize=4)
            # plot direction of pos defects
            if not arrow_len == 0:
                engine.arrow(d['pos'][0], d['pos'][1],
                             -arrow_len*cos(d['angle']),
                             -arrow_len*sin(d['angle']),
                             color='r', head_width=3, head_length=3)
        elif d['charge'] == -1/p:
            engine.plot(d["pos"][0], d["pos"][1], 'b^', markersize=4)
        elif d['charge'] == 1:
            engine.plot(d["pos"][0], d["pos"][1], 'r*', markersize=4)
        elif d['charge'] == -1:
            engine.plot(d["pos"][0], d["pos"][1], 'kX', markersize=4)

def Q_defects(Q00, Q01, engine=plt, arrow_len=0):
    """
    Plot defects of the nematic field Q.

    Args:
        Q00, Q01: Components of the nematic field.
        engine: Plotting engine or axis.
        arrow_len: If non-zero plot speed of defects as well.

    deprecated by G_defects
    """
    w = Q_charge_array(Q00, Q01)
    defects = get_defects(w, Q00, Q01)
    for d in defects:
        if d['charge'] == 0.5:
            engine.plot(d["pos"][0], d["pos"][1], 'go')
            # plot direction of pos defects
            if not arrow_len == 0:
                engine.arrow(d['pos'][0], d['pos'][1],
                             -arrow_len*cos(d['angle']),
                             -arrow_len*sin(d['angle']),
                             color='r', head_width=3, head_length=3)
        elif d['charge'] == -0.5:
            engine.plot(d["pos"][0], d["pos"][1], 'b^')
        elif d['charge'] == 1:
            engine.plot(d["pos"][0], d["pos"][1], 'r*')
        elif d['charge'] == -1:
            engine.plot(d["pos"][0], d["pos"][1], 'kX')

def polar_defects(px, py, engine=plt, arrow_len=0):
    """
    Plot defects of the polar field (px,py)

    Args:
        px, py: Components of the nematic field.
        engine: Plotting engine or axis.
        arrow_len: If non-zero plot speed of defects as well.
    """
    w = polar_charge_array(px, py)
    defects = get_defects(w)
    for d in defects:
        if d['charge'] == 0.5:
            engine.plot(d["pos"][0], d["pos"][1], 'go')
            # plot direction of pos defects
            if not arrow_len == 0:
                engine.arrow(d['pos'][0], d['pos'][1],
                             -arrow_len*cos(d['angle']),
                             -arrow_len*sin(d['angle']),
                             color='r', head_width=3, head_length=3)
        elif d['charge'] == -0.5:
            engine.plot(d["pos"][0], d["pos"][1], 'b^')
        elif d['charge'] == 1:
            engine.plot(d["pos"][0], d["pos"][1], 'r*')
        elif d['charge'] == -1:
            engine.plot(d["pos"][0], d["pos"][1], 'kX')

def local_field(frame,i,name,engine=plt,size=2,magn=False,car=True,avg=2,width=0.2):
    """
    Plot the local field of individual cell i
    Args:
        frame:Frame to plot,from archive module
        i:index of the cell to plot
        name:name of the local field
             options are 'passive','polar','nematic' ,'shape'
                     or 'active','total'
        engine:Plotting engine or axis
        size:average filter size
        magn:plot the magnitude plot or not
        car:plot the color bar or not
        avg:average over avg*avg lattices
        width:width of the arrow

    """
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    if name =='passive':
        fx,fy = frame.fp_x[i],frame.fp_y[i]
    elif name =='polar':
        fx,fy = frame.fpol_x[i],frame.fpol_y[i]
    elif name =='nematic':
        fx,fy = frame.fnem_x[i],frame.fnem_y[i]
    elif name =='shape':
        fx,fy = frame.fshape_x[i],frame.shape_y[i]
    elif name =='active':
        fx = frame.fpol_x[i] + frame.fnem_x[i] + frame.fshape_x[i]
        fy = frame.fpol_y[i] + frame.fnem_y[i] + frame.fshape_y[i]
    elif name =='total':
        fx = frame.fp_x[i]+frame.fpol_x[i] + frame.fnem_x[i] + frame.fshape_x[i]
        fy = frame.fp_y[i]+frame.fpol_y[i] + frame.fnem_y[i] + frame.fshape_y[i]

    if magn:
        m = np.sqrt(fx**2 + fy**2)
        im = engine.imshow(m.T, interpolation='lanczos', cmap='plasma',
                            origin='lower')
        if cbar:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(engine)
            ax_cb = divider.new_horizontal(size="3%", pad=0.03)
            fig = engine.get_figure()
            fig.add_axes(ax_cb)
            plt.colorbar(im, cax=ax_cb)

    fx = fx.reshape((fx.shape[0]//avg, avg, fx.shape[1]//avg, avg))
    fx = np.mean(fx, axis=(1, 3))
    fy = fy.reshape((fy.shape[0]//avg, avg, fy.shape[1]//avg, avg))
    fy = np.mean(fy, axis=(1, 3))

    cax = engine.quiver(np.arange(0, frame.parameters['Size'][0], step=avg),
                        np.arange(0, frame.parameters['Size'][1], step=avg),
                        fx.T, fy.T,width = width,
                        pivot='tail',units='dots',scale=0.05,scale_units='inches')

def cell(frame, i, engine=plt, color='k',contour_levels=0.5):
    """
    Plot a single phase field as a contour.

    Args:
        frame: Frame to plot, from archive module.
        i: Index of the cell to plot.
        engine: Plotting engine or axis.
        color: Color to use for the contour.
    """
    # if i == 1:
    #     color='r'
    # elif i == 16:
    #     color = 'r'
    # elif i == 18:
    #     color = 'b'
    # elif i == 17:
    #     color = 'b'

    #     # 16 nbrs: [ 1, 2, 12, 17, 18, 19 ]
    #     # 1 nbrs: [ 6, 16, 17, 18, 20, 21 ]
    
    q = frame.phi[i]
    #p = change_object_to_float(q, frame.parameters['Size'][0], frame.parameters['Size'][1])
    p = q

    if isinstance(contour_levels,list):
        lvs = contour_levels
    else:
        lvs = [contour_levels]

    lx=frame.parameters['Size'][0]
    ly=frame.parameters['Size'][1]

    engine.contour(np.arange(0, frame.parameters['Size'][0], 1),
                   np.arange(0, frame.parameters['Size'][1], 1),
                   p,
                   cmap=cm.winter,
                   levels=[0.5],
                   #colors=color,
                   alpha=0.5)

def cells(frame, engine=plt, colors='k',contour_levels=0.5):
    """
    Plot all cells as contours.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        colors: Colors to use for the contour. Can also be a list of colors,
            one for each cell.
    """
    if not isinstance(colors, list):
        colors = len(frame.phi)*[colors]

    for i in range(len(frame.phi)):
        cell(frame, i, engine, color=colors[i],contour_levels=contour_levels)

    # cb = plt.colorbar(matplotlib.cm.ScalarMappable(cmap='viridis'),ax=engine)
    # cb.remove()

def nbr_cells(frame, j, engine=plt, colors='k',contour_levels=0.5):
    """
    Plot all cells as contours.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        colors: Colors to use for the contour. Can also be a list of colors,
            one for each cell.
    """
    if not isinstance(colors, list):
        colors = len(frame.phi)*[colors]

    cell(frame, j, engine, color='r',contour_levels=contour_levels)
    for i in frame.nbr_cells[j]:
        cell(frame, i, engine, color=colors[i],contour_levels=contour_levels)

def interfaces(frame, engine=plt):
    """
    Plot the overlap between cells as heatmap using beautiful shades of gray
    for an absolutely photorealistic effect that will impress all your friends.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    totphi = np.zeros(frame.parameters['Size'])
    for i in range(len(frame.phi)):
        totphi += frame.phi[i]*frame.parameters['walls']
        for j in range(i+1, len(frame.phi)):
            totphi += frame.phi[i]*frame.phi[j]

    cmap = LinearSegmentedColormap.from_list('mycmap', ['grey', 'white'])
    engine.imshow(totphi.T, interpolation='lanczos', cmap=cmap, origin='lower')

def interfaces2(frame, engine=plt):
    """
    Plot the overlap between cells as heatmap in a different but also beatiful
    way.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    totphi = [np.zeros(frame.parameters['Size']),
              np.zeros(frame.parameters['Size'])]
    for i in range(len(frame.phi)):
        k = 0 if i < 64 else 1
        totphi[k] += frame.phi[i]*frame.parameters['walls']
        for j in range(i+1, len(frame.phi)):
            totphi[k] += frame.phi[i]*frame.phi[j]

    cmap0 = LinearSegmentedColormap.from_list('mycmap', ['grey', 'white'])
    cmap1 = LinearSegmentedColormap.from_list('mycmap', ['blue', 'white'])
    engine.imshow(totphi[0].T, interpolation='lanczos',
                  cmap=cmap0, origin='lower')
    engine.imshow(totphi[1].T, interpolation='lanczos',
                  cmap=cmap1, origin='lower')

def QSfill(frame, engine=plt, colors='gray', colormap = 'RdBu'):
    """
    Plot all phase fields with solid colours corresponding to individual areas.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        colors: Colors to use for the contour. Can also be a list of colors,
            one for each cell.
    """

    cmap = matplotlib.cm.get_cmap(colormap)

    if not isinstance(colors, list):
        colors = len(frame.phi)*[colors]
    for i in range(len(frame.phi)):
        # Q:S
        tmp = 2*frame.Q00[i]*frame.S00[i] + 2*frame.Q01[i]*frame.S01[i];

        # cos(2\theta), \theta the difference in angle between Q, S directors
        val = 2*tmp

        alpha = 1
        colors[i] = to_hex(cmap(val))

        engine.contourf(np.arange(0, frame.parameters['Size'][0]),
                        np.arange(0, frame.parameters['Size'][1]),
                        frame.phi[i].T,
                        levels=[.5, 10.],
                        # cmap = colormap,
                        colors=colors[i],
                        alpha=alpha)
        # engine.add_colorbar()

    # engine.set_title(label)

    # plt.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap),ax=engine)
    # cmap.set_array([])
    # plt.colorbar(cmap,ax=engine)

def zetafill(frame, engine=plt, colors='gray', colormap = 'coolwarm'):
    """
    Plot all phase fields with solid colours corresponding to individual areas.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        colors: Colors to use for the contour. Can also be a list of colors,
            one for each cell.
    """

    def map(val,min,max):
        if min == max:
            return 0.5
        else:
            return (val-min)/(max-min)
        
    cmap = matplotlib.cm.get_cmap(colormap)

    zetamax = frame.parameters['zetaQ'][0]
    zetamin = frame.parameters['zetaQ'][0]

    for i in range(1,len(frame.phi)):
        val = frame.parameters['zetaQ'][i]
        if val > zetamax:
            zetamax = val
        if val < zetamin:
            zetamin = val

    # norm = Normalize(vmin=-1,vmax=1)
    # norm = Normalize(vmin=zetamin,vmax=zetamax)

    if not isinstance(colors, list):
        colors = len(frame.phi)*[colors]

    for i in range(len(frame.phi)):
        val = frame.parameters['zetaQ'][i]

        # alpha = 1-min(1, frame.area[i]/(pi*frame.parameters['R'][i]**2))
        alpha = 1
        # print(len(frame.nbr_cells[i])-5)
        colors[i] = to_hex(cmap(map(val,zetamin,zetamax)))
        # colors[i] = to_hex(norm(val))
        # rgba = cmap(val)
        # print(colors)
        engine.contourf(np.arange(0, frame.parameters['Size'][0]),
                        np.arange(0, frame.parameters['Size'][1]),
                        frame.phi[i].T,
                        levels=[.5, 10.],
                        # cmap = colormap,
                        colors=colors[i],
                        alpha=alpha)

def psi6(frame, i):
    """
    Calculates bond-orientational order parameter for cell n
    """

    Lx = frame.parameters['Size'][0]
    Ly = frame.parameters['Size'][1]

    nnbrs = len(frame.nbr_cells[i])

    psi = [0,0]
    c = frame.com[i]

    for j in frame.nbr_cells[i]:

        cj = frame.com[j]
        d=[(cj[0]-c[0])%Lx,(cj[1]-c[1])%Ly]
        if d[0]>Lx/2:
            d[0]-=Lx
        if d[1]>Ly/2:
            d[1]-=Ly

        theta = atan2(d[1],d[0])
        psi[0]+=cos(6*theta)
        psi[1]+=sin(6*theta)

    if nnbrs>0:
        psi[0]/=nnbrs
        psi[1]/=nnbrs

    return psi

def psifill(frame, engine=plt, colors='gray', colormap='viridis'):
    """"
    Plot phase field cells with solid colours corresponding to individual valuse of the 6-fold bond orientational order parameter.

    Args:
    frame: Frame to plot, from archive module.
    engine: Plotting engine or axis.
    colors: Colors to use for the contour. Can also be a list of colors,
        one for each cell.
    """

    cmap = matplotlib.cm.get_cmap(colormap)

    if not isinstance(colors, list):
        colors = len(frame.phi)*[colors]

    for i in range(len(frame.phi)):

        # print(i)
        psi = psi6(frame,i)

        val = sqrt(psi[0]**2+psi[1]**2)

        alpha = 1
        colors[i] = to_hex(cmap(val))
        engine.contourf(np.arange(0, frame.parameters['Size'][0]),
                        np.arange(0, frame.parameters['Size'][1]),
                        frame.phi[i].T,
                        levels=[.5, 10.],
                        # cmap = colormap,
                        colors=colors[i],
                        alpha=alpha)
        # engine.set_title(label)
    
def valfill(frame, engine=plt, colors='gray',p = 2, colormap = 'viridis'):
    """
    Plot all phase fields with solid colours corresponding to individual areas.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        colors: Colors to use for the contour. Can also be a list of colors,
            one for each cell.
    """

    cmap = matplotlib.cm.get_cmap(colormap)

    if not isinstance(colors, list):
        colors = len(frame.phi)*[colors]
    for i in range(len(frame.phi)):
        val = 0
        if(p == 2):
            val = np.sqrt(frame.G00[i]**2+frame.G01[i]**2)
            label = 'nematic order'
        elif(p == 3):
            val = np.sqrt(frame.G000[i]**2+frame.G001[i]**2)
            label = 'triatic order'
        elif(p == 4):
            val = np.sqrt(frame.G0000[i]**2+frame.G0001[i]**2)
            label = 'tetratic order'
        elif(p == 6):
            val = np.sqrt(frame.G000000[i]**2+frame.G000001[i]**2)
            label = 'hexatic order'

        # alpha = 1-min(1, frame.area[i]/(pi*frame.parameters['R'][i]**2))
        alpha = 1
        # print(len(frame.nbr_cells[i])-5)
        # colors[i] = cmap(val)
        colors[i] = to_hex(cmap(val))
        # rgba = cmap(val)
        # print(colors)
        engine.contourf(np.arange(0, frame.parameters['Size'][0]),
                        np.arange(0, frame.parameters['Size'][1]),
                        frame.phi[i].T,
                        levels=[.5, 10.],
                        # cmap = colormap,
                        colors=colors[i],
                        alpha=alpha)
        engine.set_title(label)
        # engine.add_colorbar()

    plt.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap),ax=engine)
    # cmap.set_array([])
    # plt.colorbar(cmap,ax=engine)

def solidarea(frame, engine=plt, colors='gray',nbrs = False):
    """
    Plot all phase fields with solid colours corresponding to individual areas.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        colors: Colors to use for the contour. Can also be a list of colors,
            one for each cell.
        nbr: the colors of cells are determined by the number of neighbours of each cells
             #neihbours= 5(blue) 6(grey) 7(red)
    """
    if not isinstance(colors, list):
        colors = len(frame.phi)*[colors]
    for i in range(len(frame.phi)):
        alpha = 1-min(1, frame.area[i]/(pi*frame.parameters['R'][i]**2))
        if nbrs:
            c = ['blue','white','red','green','purple','olive','cyan','pink','yellow']
            if len(frame.nbr_cells[i]) > 4:
                alpha = 1
                # print(len(frame.nbr_cells[i])-5)
                colors[i] = c[len(frame.nbr_cells[i])-5]
        engine.contourf(np.arange(0, frame.parameters['Size'][0]),
                        np.arange(0, frame.parameters['Size'][1]),
                        frame.phi[i].T,
                        levels=[.5, 10.],
                        colors=colors[i],
                        alpha=alpha)

def segregationindex(frame, engine = plt):
    # computes segregation index for each cell in a frame, returns dictionary with indices the values of activity and 'entire' for the entire system

    # raw = {}
    # result = {}
    raw = []
    # result = []

    def remove_duplicates(x):
        out = []
        for i in x:
            if i not in out:
                out.append(i)

        out.sort()
        return out

    # activity = remove_duplicates(frame.parameters['zetaQ'])
    # N_activity = len(activity)

    #raw['entire'] = []

    # for i in activity:
        # raw[i] = []

    for i in range(len(frame.phi)):
        tot = 0.0
        zeta = frame.parameters['zetaQ'][i]
        nbrs = len(frame.nbr_cells[i])
        for j in frame.nbr_cells[i]:
            if zeta == frame.parameters['zetaQ'][j]:
                tot += 1

        # raw[zeta].append(tot/nbrs)
        if nbrs==0:
            raw.append(0)
        else:
            raw.append(tot/nbrs)
        # raw['entire'].append(tot/nbrs)

    # for i in raw.keys():
        # result[i] = (np.average(raw[i]),np.std(raw[i]))
    result = (np.average(raw),np.std(raw))

    return result

def nbr_hist(frame, engine = plt):
    max = 11
    nbr_vals = range(max)
    nbr_nums = np.zeros(max)

    for i in range(len(frame.phi)):
        nbr_nums[len(frame.nbr_cells[i])] += 1

    engine.bar(nbr_vals,nbr_nums)


def show_neighbour(frame,engine = plt,cell_index = 0):
    for n in frame.nbr_cells[cell_index]:
        engine.text(5+2*n,25,str(n),fontsize = 15)

def com(frame, engine=plt,plotIndex = False):
    """
    Plot the center-of-mass of each cell as a red dot. Not really
    photorealistic.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    for c in frame.com:
        engine.plot(c[0], c[1], 'ro')
    if plotIndex:
        for i in range(len(frame.com)):
            engine.text(frame.com[i][0],frame.com[i][1],str(i),fontsize = 15)

def shape(frame, engine=plt, color='k'):
    """
    Print shape tensor of each cell as the director of a nematic tensor.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    for i in range(frame.nphases):
        #Q00 = frame.S00[i]
        #Q01 = frame.S01[i]
        Q00 = frame.Q00[i]
        Q01 = frame.Q01[i]
        S = sqrt(Q00**2 + Q01**2)
        w = atan2(Q01, Q00)/2
        nx = cos(w)
        ny = sin(w)
        c = frame.com[i]
        engine.arrow(c[0], c[1],  S*nx,  S*ny, color=color)
        engine.arrow(c[0], c[1], -S*nx, -S*ny, color=color)

def shapeQ(frame, engine=plt, color='k'):
    """
    Print shape tensor of each cell as the director of a nematic tensor.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    for i in range(frame.nphases):
        Q00 = frame.Q00[i]
        Q01 = frame.Q01[i]
        S = 4*sqrt(Q00**2 + Q01**2)
        w = atan2(Q01, Q00)/2
        nx = cos(w)
        ny = sin(w)
        c = frame.com[i]
        engine.arrow(c[0], c[1],  S*nx,  S*ny, color=color)
        engine.arrow(c[0], c[1], -S*nx, -S*ny, color=color)

def director(Qxx, Qxy, avg=1, scale=False, engine=plt):
    """
    Plot director field associated with a given nematic field.

    Args:
        avg: Coarse-graining size.
        scale: Scale factor that controls the size of the director.
        engine: Plotting engine or axis.
    """

    # obtain S, nx, and ny
    S = np.vectorize(sqrt)(Qxy**2 + Qxx**2)
    nx = np.vectorize(sqrt)((1 + Qxx/S)/2)
    ny = np.sign(Qxy)*np.vectorize(sqrt)((1 - Qxx/S)/2)

    (LX, LY) = S.shape

    # construct nematic lines
    x = []
    y = []
    Smax = np.amax(S)
    for i, j in product(np.arange(LX, step=avg),
                        np.arange(LY, step=avg)):
        f = avg*(S[i, j]/Smax if scale else 1.)
        x.append(i + .5 - f*nx[i, j]/2.)
        x.append(i + .5 + f*nx[i, j]/2.)
        x.append(None)
        y.append(j + .5 - f*ny[i, j]/2.)
        y.append(j + .5 + f*ny[i, j]/2.)
        y.append(None)

    engine.plot(x, y, color='k', linestyle='-', linewidth=1)

def p_angle(frame, p, size=1, show_def=False, engine=plt, cmap='hsv'):
    """
    Plot angle of director field associated with a given p-atic field.

    Args:
        size: size of uniform filter
        engine: Plotting engine or axis.
        show_def: show defect or not
    """

    # get field
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    fx,fy = coarse_grain_p(frame, p, size)

    fx *= (1.-frame.parameters['walls'])
    fy *= (1.-frame.parameters['walls'])
    # defects
    if show_def:
        G_defects(fx, fy, p=p, engine=engine)

    theta = np.arctan2(fy,fx)/p
    theta[np.where(theta<1e-3)] += 2*np.pi/p

    # nx = np.vectorize(cos)(theta)
    # ny = np.vectorize(sin)(theta)

    (LX, LY) = fx.shape
    x = np.arange(0,LX)
    y = np.arange(0,LY)
    xv,yv = np.meshgrid(x,y)
    cb = engine.pcolor(theta.T,cmap=cmap,norm=colors.Normalize(vmin=0,vmax=2*np.pi/p))
    fig = engine.get_figure()
    cbar = fig.colorbar(cb, ax=engine, ticks=[0,np.pi/p,2*np.pi/p])

    if(p==2):
        engine.set_title('coarse-grained nematic phase')
    elif(p==3):
        engine.set_title('coarse-grained triatic phase')
    elif(p==4):
        engine.set_title('coarse-grained tetratic phase')
    elif(p==6):
        engine.set_title('coarse-grained hexatic phase')

def shape_angle(frame, size=1,show_def=False, engine=plt):
    """
    Plot angle of director field associated with a given nematic field.

    Args:
        Qxx, Qxy: Components of the nematic field.
        size: size of uniform filter
        engine: Plotting engine or axis.
        show_def: show defect or not
    """
    # get field
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    (Qxx, Qxy) = get_nematic_field(frame.phi, frame.S00, frame.S01,
                                   size=size, mode=mode)
    Qxx *= (1.-frame.parameters['walls'])
    Qxy *= (1.-frame.parameters['walls'])
    # defects
    if show_def:
        Q_defects(Qxx, Qxy, engine=engine)
    # obtain S, nx, and ny
    # does a sneaky half-angle formula for us
    S = np.vectorize(sqrt)(Qxy**2 + Qxx**2)
    nx = np.vectorize(sqrt)((1 + Qxx/S)/2)
    ny = np.sign(Qxy)*np.vectorize(sqrt)((1 - Qxx/S)/2)
    # get angle
    angle = np.arctan2(ny,nx) # -pi < angle < pi
    angle[np.where(angle<1e-3)] += np.pi   # make sure angle is between 0,pi
    (LX, LY) = S.shape
    x = np.arange(0,LX)
    y = np.arange(0,LY)
    xv,yv = np.meshgrid(x,y)
    cb = engine.pcolor(angle.T)
    fig = engine.get_figure()
    fig.colorbar(cb, ax=engine)

def nematic_field(frame, size=1, avg=1, show_def=False, arrow_len=0,
                  engine=plt):
    """
    Plot nematic field associated with the internal degree of freedom

    Args:
        frame: Frame to plot, from archive module.
        size: Coarse-graining size.
        avg: Average size (reduces the number of points plotted)
        show_def: If true, show defects.
        arrow_len: If non-zero, prints defect speed.
        engine: Plotting engine or axis.
    """
    # get field
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    (Qxx, Qxy) = get_nematic_field(frame.phi, frame.Q00, frame.Q01,
                                   size=size, mode=mode)
    Qxx *= (1.-frame.parameters['walls'])
    Qxy *= (1.-frame.parameters['walls'])
    # defects
    if show_def:
        Q_defects(Qxx, Qxy, engine=engine, arrow_len=arrow_len)
    # plot
    director(Qxx, Qxy, avg=avg, engine=engine,scale = True)

def shape_field(frame, size=1, avg=1, show_def=False, arrow_len=0,
                engine=plt):
    """
    Plot nematic field associated with the shape tensor of each cell.

    Args:
        frame: Frame to plot, from archive module.
        size: Coarse-graining size.
        avg: Average size (reduces the number of points plotted)
        show_def: If true, show defects.
        arrow_len: If non-zero, prints defect speed.
        engine: Plotting engine or axis.
    """
    # get field
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    (Qxx, Qxy) = get_nematic_field(frame.phi, frame.S00, frame.S01,
                                   size=size, mode=mode)
    Qxx *= (1.-frame.parameters['walls'])
    Qxy *= (1.-frame.parameters['walls'])
    # defects
    if show_def:
        Q_defects(Qxx, Qxy, engine=engine, arrow_len=arrow_len)
    # plot
    director(Qxx, Qxy, avg=avg, engine=engine,scale = True)

def shape_alignment(frame,engine=plt,size=1,avg=2,cbar = True):
    """
    Plot the alignment of different field <cos 2theta >

    Args:
        frame: Frame to plot, from archive module.
        size: Coarse-graining size.
        avg: Average size (reduces the number of points plotted)
        show_def: If true, show defects.
        engine: Plotting engine or axis.
    """
    o = []
    for i in range(frame.nphases):
        cos2t = cos(atan2(frame.S01[i],frame.S00[i]))
        for j in frame.nbr_cells[i]:
            cos2t += cos(atan2(frame.S01[j],frame.S00[j]))
        cos2t /= (len(frame.nbr_cells[i]) + 1)
        o.append(cos2t)
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    S = _get_field(frame.phi,[o],size=size, mode=mode)
    '''
    # get the angle
    if field_type == 'polar':
        theta = np.arctan2(frame.polarization[:,1],frame.polarization[:,0])
    elif field_type == 'nematic':
        theta = 0.5*np.arctan2(frame.Q01,frame.Q00)
    elif field_type == 'shape':
        theta = 0.5*np.arctan2(frame.S01,frame.S00)
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    S = _get_field(frame.phi,[np.cos(2.0*theta)],size=size, mode=mode)
    '''
    im = engine.imshow(S[0].T, interpolation='lanczos', cmap='plasma',
                            origin='lower')
    if cbar:
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(engine)
        ax_cb = divider.new_horizontal(size="3%", pad=0.03)
        fig = engine.get_figure()
        fig.add_axes(ax_cb)
        ticks = [-1+i*0.2 for i in range(11)]
        plt.colorbar(im, cax=ax_cb,ticks=ticks)

'''
def grad_phi(frame, i, size=15, engine=plt, cbar=True):
    """
    Plot the norm of the phase-field gradient of a single cell as a heat map.

    Args:
        frame: Frame to plot, from archive module.
        size: Coarse-graining size.
        engine: Plotting engine or axis.
        cbar: Show color bar?
    """
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'nearest'
    # vx, vy = get_velocity_field(frame.phi, frame.com_velocity, size, mode=mode)
    # vx *= (1.-frame.parameters['walls'])
    # vy *= (1.-frame.parameters['walls'])

    p = frame.phi[i]
    print(p.shape)
    g = frame.norm_grad_phi[i]
    print(g.shape)
    w = get_grad_phi(p, g, size, mode)

    mode = 'wrap' if frame.parameters['BC'] == 0 else 'nearest'
    w = ndimage.filters.uniform_filter(w, size=size, mode=mode)
    # w is only one-dimensional
    # print(w.ndim)
    im = engine.imshow(w.T, interpolation='lanczos', cmap='coolwarm',
                        origin='lower')
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(engine)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)

    if cbar:
        c =plt.colorbar(im,cax=cax,format=ticker.FuncFormatter(fmt))
        c.ax.tick_params(labelsize=12)
    return w
'''

def p_directors(frame, g00, g01, engine=plt, p=2, dx=10, dy=10):
    '''
    plots unit directors from fields g00, g01 spaced according to dx, dy

    takes as arguments the coarse-grained order parameters, p, spacing dx, dy
    '''

    lx,ly = frame.phi[0].shape

    for i in np.arange(0,lx,dx):
        for j in np.arange(0,ly,dy):
            theta = atan2(g01[i,j],g00[i,j])/p

            for k in range(p):
                engine.arrow(i, j, 2*cos(theta + k*2*np.pi/p), 2*sin(theta + k*2*np.pi/p), color='k', width = 0.0001, head_width = None, head_length = None)

def coarse_grain_p(frame, p=2, size=1):
    '''
    coarse-grains the p-atic order parameter. more or less a wrapper for get_p_atic_field and _get_field
    so we don't have to think to hard when writing our scripts

    takes as arguments the frame, vector order parameter assigned to each cell, and p

    it's untidy to make two heat maps in different subfigures from within plot.py, so we need two functions:
    coarse_grain_theta and coarse_grain_magnitude
    '''

    lx,ly = frame.phi[0].shape

    fx = []
    fy = []

    if(p==2):
        fx, fy = get_p_atic_field(frame.phi, frame.G00, frame.G01, size=16)
    elif(p==3):
        fx, fy = get_p_atic_field(frame.phi, frame.G000, frame.G001, size=16)
    elif(p==4):
        fx, fy = get_p_atic_field(frame.phi, frame.G0000, frame.G0001, size=16)
    elif(p==6):
        fx, fy = get_p_atic_field(frame.phi, frame.G000000, frame.G000001, size=16)

    return fx,fy

def coarse_grain_magnitude(frame, g00, g01, engine=plt, p=2, colormap = 'viridis'):
    '''
    plots as heat map the magnitude of a coarse-grained vector field
    plots also the directors

    takes as arguments coarse-grained vector order parameter, p
    '''

    magn = np.sqrt(np.square(g00)+np.square(g01))
    magn.reshape(frame.phi[0].shape)

    # print(frame.phi[0].shape)
    # print(magn.shape)

    engine.matshow(magn.transpose(), cmap=colormap)

    p_directors(frame, g00, g01, engine, p)

    # if(p==2):
    #     engine.set_title('coarse-grained nematic order')
    # elif(p==3):
    #     engine.set_title('coarse-grained triatic order')
    # elif(p==4):
    #     engine.set_title('coarse-grained tetratic order')
    # elif(p==6):
    #     engine.set_title('coarse-grained hexatic order')

    cb = engine.pcolor(magn.T,cmap=colormap,norm=colors.Normalize(vmin=0,vmax=0.5))
    # fig = engine.get_figure()
    # cbar = fig.colorbar(cb, ax=engine)
    # plt.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap),ax=engine)
    # plt.colorbar()

    return cb

def coarse_grain_theta(frame, g00, g01, show_def=False, engine=plt, p=2, colormap = 'hsv'):
    """
    Plot angle of director field associated with a given p-atic field.

    Args:
        g00, g01 the coarse-grained x and y components of the p-atic shape tensor
        size: size of uniform filter
        engine: Plotting engine or axis.
        show_def: show defect or not

    Previously wrote this function from scratch but its output was a bit of a mess. Co-opting the code from p_angle.
    Only difference is that you pass in the coarse-grained fields to this function while p_angle calls coarse_grain_p itself
    """

    # get field
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    # fx,fy = coarse_grain_p(frame, p, size)

    g00 *= (1.-frame.parameters['walls'])
    g01 *= (1.-frame.parameters['walls'])
    # defects
    if show_def:
        G_defects(g00, g01, p=p, engine=engine)

    theta = np.arctan2(g00,g01)/p
    theta[np.where(theta<1e-3)] += 2*np.pi/p

    # nx = np.vectorize(cos)(theta)
    # ny = np.vectorize(sin)(theta)

    (LX, LY) = g00.shape
    x = np.arange(0,LX)
    y = np.arange(0,LY)
    xv,yv = np.meshgrid(x,y)
    cb = engine.pcolor(theta.T,cmap=colormap,norm=colors.Normalize(vmin=0,vmax=2*np.pi/p))
    # fig = engine.get_figure()
    # cbar = fig.colorbar(cb, ax=engine, ticks=[0,np.pi/p,2*np.pi/p])

    # if(p==2):
    #     engine.set_title('coarse-grained nematic phase')
    # elif(p==3):
    #     engine.set_title('coarse-grained triatic phase')
    # elif(p==4):
    #     engine.set_title('coarse-grained tetratic phase')
    # elif(p==6):
    #     engine.set_title('coarse-grained hexatic phase')

    return cb

def traction(frame, engine=plt):
    '''
    Return traction force map.

    Traction is the mismatch between velocity & active forces, i.e. the passive forces

    Args:
        frame: Frame to extract force data
        cg: Whether to coarse-grain by cell or
        engine: Plotting engine or axis
    '''

    # mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'

    # get passive force density

    fp_x,fp_y = frame.fp_field_x,frame.fp_field_y
    # ffric_x,ffric_y = frame.ffric_field_x,frame.ffric_field_y

    # fx,fy = get_avg_force_field(frame.phi,frame.com_velocity,frame.parameters['xi'],frame.area)

    fpol_x, fpol_y = frame.fpol_field_x,frame.fpol_field_y
    fdipole_x,fdipole_y = frame.fdipole_field_x,frame.fdipole_field_y

    traction_x = fp_x.T + fpol_x.T + fdipole_x.T
    traction_y = fp_y.T + fpol_y.T + fdipole_y.T

    # active force density
    # fpol_x,fpol_y = frame.fpol_field_x,frame.fpol_field_y
    # fnem_x,fnem_y = frame.fdipole_field_x,frame.fdipole_field_y

    # traction_x = fp_x.T #np.add(fp_x,ffric_x)
    # traction_y = fp_y.T #np.add(fp_y,ffric_y)

    # traction_x = fx.T - fpol_x.T - fnem_x.T
    # traction_y = fy.T - fpol_y.T - fnem_y.T

    # traction_x = -fpol_x.T - fnem_x.T
    # traction_y = -fpol_y.T - fnem_y.T

    return traction_x,traction_y

def cg_stress(frame, size=3, engine=plt):
    '''
    Return coarse-grained monolayer stress.

    Stress calculated using formula
        originated by: Christoffersen, Mehrabadi, Nemat-Nasser (1981)
        dug up by: Monfared, Ravichandran, Andrade, Doostmohammadi (2022)

    Args:
        frame: Frame to calculate stress map
        size: Coarse-graining size (must be odd number)
        engine: Plotting engine or axis.
    '''

    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'

    V_cg = np.power(size,2)
    d = int(np.ceil((size-1)/2))

    dims = frame.phi[0].shape

    stress_xx = np.zeros(dims)
    stress_yy = np.zeros(dims)
    stress_xy = np.zeros(dims)

    traction_x,traction_y = traction(frame)

    for i in np.arange(-d,d+1):
        for j in np.arange(-d,d+1):
            r = np.sqrt(i**2+j**2)

            # the contribution from x_m=0 is zero anyway
            if r!=0:
                ex = -i/r
                ey = -j/r

                tempx = np.roll(np.roll(traction_x,-i,axis=1),-j,axis=0)
                tempy = np.roll(np.roll(traction_y,-i,axis=1),-j,axis=0)

                stress_xx += tempx*ex
                stress_yy += tempy*ey
                stress_xy += tempx*ey + tempy*ex

    stress_xx /= V_cg
    stress_yy /= V_cg
    stress_xy /= (2*V_cg)

    return (stress_xx, stress_yy, stress_xy), (traction_x,traction_y)

def chemical_field(frame, size=15, engine=plt, cbar=True):
    """
    Plot the value of the chemical concentration as a heat map.

    Args:
        frame: Frame to plot, from archive module.
        size: Coarse-graining size.
        engine: Plotting engine or axis.
        cbar: Show color bar?
    """
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'nearest'
    # vx, vy = get_velocity_field(frame.phi, frame.com_velocity, size, mode=mode)
    # vx *= (1.-frame.parameters['walls'])
    # vy *= (1.-frame.parameters['walls'])
    w = get_chemical_field(frame.phi, frame.chem, size, mode)
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'nearest'
    w = ndimage.filters.uniform_filter(w, size=size, mode=mode)
    # w is only one-dimensional
    # print(w.ndim)
    im = engine.imshow(w.T, interpolation='lanczos', cmap='coolwarm',
                        origin='lower')
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(engine)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)

    if cbar:
        c =plt.colorbar(im,cax=cax,format=ticker.FuncFormatter(fmt))
        c.ax.tick_params(labelsize=12)
    return w

def polarity_field(frame, size=2, show_def = False,engine=plt,magn=True, cbar=True, avg=2,width = 0.2):
    """
    Plot polarity field associated with the internal degree of freedom

    Args:
        frame: Frame to plot, from archive module.
        size: Coarse-graining size.
        avg: Average size (reduces the number of points plotted)
        show_def: If true, show defects.
        arrow_len: If non-zero, prints defect speed.
        engine: Plotting engine or axis.
    """
    # get field
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    (px, py) = get_polarity_field(frame.phi, frame.polarization, size=size, mode=mode)
    px *= (1.-frame.parameters['walls'])
    py *= (1.-frame.parameters['walls'])
    # defects
    if show_def:
        polar_defects(px, py, engine=engine)
    # plot
    if magn:
        m = np.sqrt(px**2 + py**2)
        im = engine.imshow(m.T, interpolation='lanczos', cmap='plasma',
                            origin='lower')
        if cbar:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(engine)
            ax_cb = divider.new_horizontal(size="3%", pad=0.03)
            fig = engine.get_figure()
            fig.add_axes(ax_cb)
            plt.colorbar(im, cax=ax_cb)

    px = px.reshape((px.shape[0]//avg, avg, px.shape[1]//avg, avg))
    px = size*np.mean(px, axis=(1, 3))
    py = py.reshape((py.shape[0]//avg, avg, py.shape[1]//avg, avg))
    py = size*np.mean(py, axis=(1, 3))

    cax = engine.quiver(np.arange(0, frame.parameters['Size'][0], step=avg),
                        np.arange(0, frame.parameters['Size'][1], step=avg),
                        px.T, py.T,
                        pivot='tail', units='xy',scale_units='xy',width = width)

def repulsive_field(frame, size=2, engine=plt,magn=False, cbar=True, avg=2,width = 0.2,cg = 'cell'):
    """
    Plot repulsive field associated with the overlap of each cell junction as deviation from symmetric circular repulsion.

    Args:
        frame: Frame to plot, from archive module.
        size: Coarse-graining size.
        engine: Plotting engine or axis.
        magn: Plot velocity magnitude as a heatmap?
        cbar: Show color bar?
        avg: Size of the averaging (drops points)
        cg: Coarse-graining by cell or lattice
    """

    mode = 'wrap' if frame.parameters['BC'] == 0 else 'nearest'
    frepx,frepy = get_repulsive_field(frame.com[0], frame.phi, frame.area, frame.frep_field_x, frame.frep_field_y, frame.Fdp_tot, size=1, mode='wrap')
    #frepx,frepy = get_repulsive_field(frame.com[0], frame.phi, frame.area, frame.fov_field_x, frame.fov_field_y, frame.Frep, size=1, mode='wrap')

    if magn:
        m = np.sqrt(frepx**2 + frepy**2)
        im = engine.imshow(m.T, interpolation='lanczos', cmap='plasma',
                            origin='lower')
        if cbar:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(engine)
            ax_cb = divider.new_horizontal(size="3%", pad=0.03)
            fig = engine.get_figure()
            fig.add_axes(ax_cb)
            plt.colorbar(im, cax=ax_cb)

    cax = engine.quiver(np.arange(0, frame.parameters['Size'][0], step=1),
                        np.arange(0, frame.parameters['Size'][1], step=1),
                        frepx.T, frepy.T,width = 2*width,
                        pivot='tail', units='dots', scale_units='dots')

def average_velocity_polar(frame, axis, engine=plt):
    '''
    Returns:
        axis = 0: radial velocity averaged over angle
        axis = 1: angular velocity averaged over angle

    Arguments:
        frame: frame to plot
        axis: 0 or 1, axis over which to average the velocity field
        engine: Plotting engine or axis
    '''

    bc = frame.parameters['BC']
    inner_radius = 0
    outer_radius = frame.parameters['Size'][0]
    v = []
    counter = []

    if(bc==2 and frame.parameters['Size'][0]==frame.parameters['Size'][1]):
        inner_radius = 0
        outer_radius = np.floor(frame.parameters['Size'][0]*frame.parameters['confinement_ratio'])
        v = np.zeros(outer_radius)
        counter = np.zeros(outer_radius)
    elif(bc==10):
        inner_radius = frame.parameters['annulus_inner']
        outer_radius = frame.parameters['annulus_outer']
        v = np.zeros(outer_radius - inner_radius)
        counter = np.zeros(outer_radius - inner_radius)
    else:
        return None

    Lx = frame.parameters['Size'][0]
    Ly = frame.parameters['Size'][1]

    for i in np.arange(Lx):
        for j in np.arange(Ly):
            rx = i - Lx/2
            ry = j - Ly/2
            r = np.sqrt(rx**2 + ry**2)

            if (r > inner_radius and r < outer_radius):
                index = int(np.floor(r) - inner_radius)
                if(axis==0):
                    proj = (rx*frame.velocity_field_x[i,j] + ry*frame.velocity_field_y[i,j])/r
                elif(axis==1):
                    proj = (-ry*frame.velocity_field_x[i,j] + rx*frame.velocity_field_y[i,j])/r
                v[index] += proj
                counter[index] +=1

    result = np.divide(v,counter)
    domain = np.arange(inner_radius,outer_radius)

    engine.plot(domain,result,'-g')

def average_velocity_field(frame, engine=plt, component=0, axis=0, size=2, cg='cell'):
    '''
    Returns:
        component = 0: an averaged x velocity
        component = 1: an averaged y velocity

        axis = 0: averaged over the x axis
        axis = 1: averaged over the y axis

    Arguments:
        frame: frame to plot
        axis: 0 or 1, axis over which to average the velocity field
        engine: Plotting engine or axis
    '''

    mode = 'wrap' if frame.parameters['BC'] == 0 else 'nearest'
    if cg =='cell':
        vx, vy = get_velocity_field(frame.phi, frame.com_velocity, size, mode=mode)
    else:
        vx,vy = frame.velocity_field_x,frame.velocity_field_y
        # print(vx.shape)
    # vx *= (1.-frame.parameters['walls'])
    # vy *= (1.-frame.parameters['walls'])

    # domain = []

    domain = np.arange(frame.parameters['Size'][1-axis],step = 1)

    if component == 0:
        result = np.average(vx,axis)
    elif component == 1:
        result = np.average(vy,axis)

    # if axis == 0:
    #     domain = np.arange(frame.parameters['Size'][1],step=1)
    #     if component == 0:
    #         result = np.average(vx,axis)
    #     if component == 1:
    #         result = np.average(vy, axis)
    #     # result = avg(v,axis)
    # elif axis == 1:
    #     domain = np.arange(frame.parameters['Size'][0],step=1)
    #     result = np.average(vy,axis)
    #     # result = avg(v,axis)

    engine.plot(domain,result,'-g')

def velocity_slice(frame, axis, size=2, coord=0, engine=plt, cg='cell'):
    '''
    Returns:
        axis = 0: x velocity at x=coord as function of y
        axis = 1: y velocity at y=coord as function of x

    Arguments:
        frame: frame to plot
        axis: 0 or 1, axis over which to view the velocity field
        engine: Plotting engine or axis
    '''

    mode = 'wrap' if frame.parameters['BC'] == 0 else 'nearest'
    if cg =='cell':
        vx, vy = get_velocity_field(frame.phi, frame.com_velocity, size, mode=mode)
    else:
        vx,vy = frame.velocity_field_x,frame.velocity_field_y
        # print(vx.shape)
    vx *= (1.-frame.parameters['walls'])
    vy *= (1.-frame.parameters['walls'])

    domain = np.arange(frame.parameters['Size'][1-axis],step=1)

    if axis == 0:
        result = vx.T[:,coord]
        # print(v.T[:,coord].shape)
    elif axis == 1:
        result = vy.T[coord,:]

    # print(result)

    engine.plot(domain,result,'-g')

def velocity_field(frame, size=2, engine=plt,magn=False, cbar=True, avg=2,width = 0.2,cg = 'cell',step=1):
    """

    Args:
        frame: Frame to plot, from archive module.
        size: Coarse-graining size.
        engine: Plotting engine or axis.
        magn: Plot velocity magnitude as a heatmap?
        cbar: Show color bar?
        avg: Size of the averaging (drops points)
        cg: Coarse-graining by cell or lattice
    """

    mode = 'wrap' if frame.parameters['BC'] == 0 else 'nearest'
    if cg =='cell':
        vx, vy = get_velocity_field(frame.phi, frame.com_velocity, size, mode=mode)
    else:
        vx,vy = frame.velocity_field_x,frame.velocity_field_y
        # print(vx.shape)
    vx *= (1.-frame.parameters['walls'])
    vy *= (1.-frame.parameters['walls'])
    # print("field")
    # print(vx.T)

    if magn:
        m = np.sqrt(vx**2 + vy**2)
        im = engine.imshow(m.T, interpolation='lanczos', cmap='plasma',
                            origin='lower')
        if cbar:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(engine)
            ax_cb = divider.new_horizontal(size="3%", pad=0.03)
            fig = engine.get_figure()
            fig.add_axes(ax_cb)
            plt.colorbar(im, cax=ax_cb)

    # print(vx.shape)
    # print('\n')
    # vx = vx.reshape((vx.shape[0]//avg, avg, vx.shape[1]//avg, avg))
    # vx = np.mean(vx, axis=(1, 3))
    # vy = vy.reshape((vy.shape[0]//avg, avg, vy.shape[1]//avg, avg))
    # vy = np.mean(vy, axis=(1, 3))

    cax = engine.quiver(np.arange(0, frame.parameters['Size'][0], step=step),
                        np.arange(0, frame.parameters['Size'][1], step=step),
                        100*vx.T[::step,::step], 100*vy.T[::step,::step],width = width,
                        pivot='tail', units='dots', scale=0.05, scale_units='dots')

    # cax = engine.quiver(np.arange(0, frame.parameters['Size'][0], step=1),
    #                     np.arange(0, frame.parameters['Size'][1], step=1),
    #                     1000*vx.T, 0*vy.T,width = width,
    #                     pivot='tail', units='dots', scale=0.05, scale_units='dots')

    # print(vx.T[:,120])

def strain_rate(frame,label = 'xx',size=2, engine=plt,cbar=True, avg=2,width = 0.2,cg = 'cell'):
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    if cg =='cell':
        vx, vy = get_velocity_field(frame.phi, frame.com_velocity, size, mode=mode)
    else:
        vx,vy = frame.velocity_field_x,frame.velocity_field_y
    vx *= (1.-frame.parameters['walls'])
    vy *= (1.-frame.parameters['walls'])

    vx = vx.reshape((vx.shape[0]//avg, avg, vx.shape[1]//avg, avg))
    vx = np.mean(vx, axis=(1, 3))
    vy = vy.reshape((vy.shape[0]//avg, avg, vy.shape[1]//avg, avg))
    vy = np.mean(vy, axis=(1, 3))
    if label == 'xx':
        e = np.gradient(vx,axis = 0)
    elif label == 'xy':
        e = np.gradient(vx,axis = 1)
    elif label =='yx':
        e = np.gradient(vy,axis = 0)
    elif label == 'yy':
        e = np.gradient(vy,axis = 1)
    im = engine.imshow(e.T, interpolation='lanczos', cmap='plasma',
                            origin='lower')
    if cbar:
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(engine)
        ax_cb = divider.new_horizontal(size="3%", pad=0.03)
        fig = engine.get_figure()
        fig.add_axes(ax_cb)
        plt.colorbar(im, cax=ax_cb)

def __modu(a,b):
    while a >= b:
        a -= b
    while a < 0:
        a += b
    return a

def vorticity_field(frame, size=15, engine=plt, cbar=True):
    """
    Plot nematic field associated with the shape tensor of each cell.

    Args:
        frame: Frame to plot, from archive module.
        size: Coarse-graining size.
        engine: Plotting engine or axis.
        cbar: Show color bar?
    """
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    vx, vy = get_velocity_field(frame.phi, frame.com_velocity, size, mode=mode)
    vx *= (1.-frame.parameters['walls'])
    vy *= (1.-frame.parameters['walls'])
    w = get_vorticity_field(vx, vy)
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    w = ndimage.filters.uniform_filter(w, size=size, mode=mode)
    im = engine.imshow(w.T, interpolation='lanczos', cmap='viridis',
                        origin='lower')
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(engine)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)

    if cbar:
        c =plt.colorbar(im,cax=cax,format=ticker.FuncFormatter(fmt))
        c.ax.tick_params(labelsize=12)
    return w

def force_density(frame, engine=plt,force_type = 'all',magn=True, cbar=True, avg=1,size = 8,
                  width = 1,scale = None,cg = 'cell',step=1):
    """
    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        type: all, polar,dipole,passive,friction
        magn: Plot forcedistribution  magnitude as a heatmap?
        cbar: Show color bar?
        avg: Size of the averaging (drops points)
        cg: coarse-grain by cell or lattice
    """
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    # get passive and active force density    '
    if cg == 'cell':
        fpol_x, fpol_y = get_velocity_field(frame.phi, frame.Fpol, size, mode=mode)
        fdipole_x,fdipole_y = get_velocity_field(frame.phi, frame.Fnem + frame.Fshape, size, mode=mode)
        fp_x,fp_y = get_velocity_field(frame.phi, frame.Fpassive, size, mode=mode)
        ffric_x,ffric_y = get_velocity_field(frame.phi, frame.Ffric, size, mode=mode)
    else:
        fpol_x, fpol_y = frame.fpol_field_x,frame.fpol_field_y
        fdipole_x,fdipole_y = frame.fdipole_field_x,frame.fdipole_field_y
        fp_x,fp_y = frame.fp_field_x,frame.fp_field_y

    if force_type == 'all':
        fx = fpol_x + fdipole_x + fp_x
        fy = fpol_y + fdipole_y + fp_y
    elif force_type == 'polar':
        fx,fy = fpol_x,fpol_y
    elif force_type == 'dipole':
        fx,fy = fdipole_x,fdipole_y
    elif force_type == 'passive':
        fx,fy = fp_x,fp_y

    if magn:
        m = np.sqrt(fx**2 + fy**2)
        im = engine.imshow(m.T, interpolation='lanczos', cmap='plasma',
                            origin='lower')
    fx = fx.reshape((fx.shape[0]//avg, avg, fx.shape[1]//avg, avg))
    fx = np.mean(fx, axis=(1, 3))
    fy = fy.reshape((fy.shape[0]//avg, avg, fy.shape[1]//avg, avg))
    fy = np.mean(fy, axis=(1, 3))

    if cbar:
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(engine)
        ax_cb = divider.new_horizontal(size="3%", pad=0.03)
        fig = engine.get_figure()
        fig.add_axes(ax_cb)
        plt.colorbar(im, cax=ax_cb)

    # cax = engine.quiver(np.arange(0, frame.parameters['Size'][0], step=step),
    #                     np.arange(0, frame.parameters['Size'][1], step=avg),
    #                     fx.T, fy.T,width = width,scale = scale,
    #                     pivot='tail', units='dots', scale_units='dots')

    cax = engine.quiver(np.arange(0, frame.parameters['Size'][0], step=step),
                        np.arange(0, frame.parameters['Size'][1], step=step),
                        100*fx.T[::step,::step], 100*fy.T[::step,::step],width = width,
                        pivot='tail', units='dots', scale=0.05, scale_units='dots')

    # if force_type=='passive':    print(fx.max())

def _force(frame, i, v, engine=plt, **kwargs):
    """
    Helper function to plot forces.
    """
    c = frame.com[i]
    l = sqrt(v[0]*v[0]+v[1]*v[1])
    engine.arrow(c[0], c[1], v[0], v[1],head_width=l/6, head_length=l/3,**kwargs)

def get_v_com(phase, area, vx, vy):
    """
    Calculate centre-of-mass velocity for a given phase field using the velocity field data
    """

    v_com_x = 0
    v_com_y = 0

    (LX,LY) = phase.shape
    print(LX)

    for j in range(LX-1):
        for k in range(LY):
            index = j*LY+k
            # print(j,k,index)
            v_com_x += phase[j][k]*vx[index]
            v_com_y += phase[j][k]*vy[index]

    v_com_x /= area
    v_com_y /= area

    # print("[",v_com_x,",",v_com_y,"]\n")

    return [v_com_x,v_com_y]

def v_com(frame, engine=plt, color='r'):
    """
    Plot the centre-of-mass velocity for phase fields, calculated in this script rather than in the simulation
    """

    scale = 100*frame.parameters['nsubsteps']
    for i in range(frame.nphases):
        _force(frame, i,
               1*scale*get_v_com(frame.phi[i], frame.area[i], frame.vx[i], frame.vy[i]),
               engine=engine,
               color=color)

def velocity(frame, engine=plt, color='r'):
    """
    Plot total velocity of each cell.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        color: Color of the arrow.
    """
    scale = 100*frame.parameters['nsubsteps']
    for i in range(frame.nphases):
        # print(frame.com_velocity[i])
        _force(frame, i,
               1*scale*frame.com_velocity[i],
               engine=engine,
               color=color)

def avg_velocity(frame,engine = plt, color = 'r'):
    """
    Plot averged CoM velocity of each cell.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        color: Color of the arrow.
    """
    scale = 100*frame.parameters['nsubsteps']
    for i in range(frame.nphases):
        _force(frame, i,
               1*scale*frame.avg_velocity[i],
               engine=engine,
               color=color)

def polar_force(frame, engine=plt, color='b'):
    """
    Plot polar force of each cell.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        color: Color of the arrow.
    """
    scale = 0.2*frame.parameters['ninfo']*frame.parameters['nsubsteps']
    for i in range(frame.nphases):
        _force(frame, i,
               scale*frame.Fpol[i],
               engine=engine,
               color=color)

def shape_force(frame, engine=plt, color='b'):
    """
    Plot shape force of each cell stemming from deformation tensor

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        color: Color of the arrow.
    """
    scale = frame.parameters['ninfo']*frame.parameters['nsubsteps']
    for i in range(frame.nphases):
        _force(frame, i,
               scale*frame.Fshape[i],
               engine=engine,
               color=color)

def nematic_force(frame, engine=plt, color='b'):
    """
    Plot shape nematic force of each cell stemming from deformation tensor

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        color: Color of the arrow.
    """
    scale = frame.parameters['ninfo']*frame.parameters['nsubsteps']
    for i in range(frame.nphases):
        _force(frame, i,
               scale*frame.Fnem[i],
               engine=engine,
               color=color)

def passive_force(frame, engine=plt, color='b'):
    """
    Plot passive force of each cell.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        color: Color of the arrow.
    """
    scale = 0.1*frame.parameters['ninfo']*frame.parameters['nsubsteps']
    for i in range(frame.nphases):
        _force(frame, i,
               scale*frame.Fpassive[i],
               engine=engine,
               color=color)

def interaction_force(frame, engine=plt, color='b'):
    """
    Plot interaction force of each cell.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        color: Color of the arrow.
    """
    scale = frame.parameters['ninfo']*frame.parameters['nsubsteps']*0.01
    for i in range(frame.nphases):
        _force(frame, i,
               scale*frame.Fint[i],
               engine=engine,
               color=color)

def nematic(frame, engine=plt):
    """
    Print director of each cell as a line at their center.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    for i in range(frame.nphases):
        # Q00 = frame.S00[i]
        # Q01 = frame.S01[i]
        Q00 = frame.Q00[i]
        Q01 = frame.Q01[i]
        S = sqrt(Q00**2 + Q01**2)
        nx = sqrt((1 + Q00/S)/2)
        ny = np.sign(Q01)*sqrt((1 - Q00/S)/2)
        c = frame.com[i]
        a = frame.parameters['R'][i]/2.0*S
        engine.arrow(c[0], c[1],  a*nx,  a*ny, color='k')
        engine.arrow(c[0], c[1], -a*nx, -a*ny, color='k')

def p_atic(p, frame, engine=plt):
    """
    Print p-atic director of each cell as an aster at their center.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        p: p-fold symmetry parameter \in (2,3,4,6)
    """

    px = []
    py = []

    if(p == 2):
        px = frame.G00
        py = frame.G01
    elif(p == 3):
        px = frame.G000
        py = frame.G001
    elif(p == 4):
        px = frame.G0000
        py = frame.G0001
    elif(p == 6):
        px = frame.G000000
        py = frame.G000001

    for i in range(frame.nphases):

        theta = atan2(py[i],px[i])/p
        # print(theta)

        c = frame.com[i]

        for j in range(p):
            engine.arrow(c[0], c[1], 2*cos(theta + j*2*np.pi/p), 2*sin(theta + j*2*np.pi/p), color='k', width = 0.0001, head_width = None, head_length = None)

def polarization(frame, engine=plt):
    """
    Print poalrization of each cell as a line at their center.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    for i in range(frame.nphases):
        #px = frame.polarization[i][0]*3
        #py = frame.polarization[i][1]*3
        px = 3.0*cos(frame.theta_pol[i])
        py = 3.0*sin(frame.theta_pol[i])
        S = sqrt(px**2 + py**2)
        c = frame.com[i]
        engine.arrow(c[0], c[1],  px,  py, head_width=2/3, head_length=0.8,color='g')

def phase(frame, n, engine=plt, cbar=False):
    """
    Plot single phase as a density plot.

    Args:
        frame: Frame to plot, from archive module.
        n: Index of the cell to plot.
        engine: Plotting engine or axis.
        cbar: Display cbar?
    """
    cax = engine.imshow(frame.phi[n].T, interpolation='lanczos', cmap='Greys',
                        origin='lower')
    if cbar:
        plt.colorbar(cax)

def walls(frame, engine=plt, cbar=False):
    """
    Plot walls.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        cbar: Display cbar?
    """
    q = frame.parameters['walls']
    p = change_object_to_float(q, frame.parameters['Size'][0], frame.parameters['Size'][1])
    cax = engine.imshow(p.T, cmap='Greys',
                        origin='lower', clim=(0., 1.))
    if cbar:
        plt.colorbar(cax)

def patch(frame, n, engine=plt):
    """Plot the restricted patch of a single cell

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    def plot(m, M): engine.fill([m[0], M[0], M[0], m[0], m[0], None],
                                [m[1], m[1], M[1], M[1], m[1], None],
                                color='b', alpha=0.04)
    LX, LY = frame.parameters['Size']
    m = frame.patch_min[n]
    M = frame.patch_max[n]

    if(m[0] == M[0]):
        m[0] += 1e-1
        M[0] -= 1e-1
    if(m[1] == M[1]):
        m[1] += 1e-1
        M[1] -= 1e-1

    if(m[0] > M[0] and m[1] > M[1]):
        plot(m, [LX, LY])
        plot([0, 0], M)
        plot([m[0], 0], [LX, M[1]])
        plot([0, m[1]], [M[0], LY])
    elif(m[0] > M[0]):
        plot(m, [LX, M[1]])
        plot([0, m[1]], M)
    elif(m[1] > M[1]):
        plot(m, [M[0], LY])
        plot([m[0], 0], M)
    else:
        plot(m, M)

def patches(frame, engine=plt):
    """
    Plot the subdomain patches of each cell.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    for n in range(frame.nphases):
        patch(frame, n, engine)

def masks(frame, engine=plt):
    """
    Plot division/death masks.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    m1 = np.array([1 if i else 0 for i in frame.division_mask])
    m2 = np.array([1 if i else 0 for i in frame.death_mask])

    engine.contour(np.arange(0, frame.parameters['LX']),
                   np.arange(0, frame.parameters['LY']),
                   m1.reshape(frame.parameters['Size']).T,
                   levels=[.5], colors=['b'])

    engine.contour(np.arange(0, frame.parameters['LX']),
                   np.arange(0, frame.parameters['LY']),
                   m2.reshape(frame.parameters['Size']).T,
                   levels=[.5], colors=['r'])

def trajectories(frame,engine = plt,color = 'm'):
    import animation
    for pos in animation.position: # pos is the postition for each time step until now
        x = pos[:,0]
        y = pos[:,1]
        engine.scatter(x,y,s = 0.3,c = color)

def voronoi_lattice(frame,engine = plt):
    """
    plot voronoi lattice polygon
    """
    LX,LY = frame.parameters['Size'][0],frame.parameters['Size'][1]
    com = frame.com
    if frame.parameters['BC'] == 0:
        com = np.vstack((com,frame.com + [LX,0.0]))
        com = np.vstack((com,frame.com + [-LX,0.0]))
        com = np.vstack((com,frame.com + [0.0,LY]))
        com = np.vstack((com,frame.com + [0.0,-LY]))
        com = np.vstack((com,frame.com + [LX,LY]))
        com = np.vstack((com,frame.com + [LX,-LY]))
        com = np.vstack((com,frame.com + [-LX,LY]))
        com = np.vstack((com,frame.com + [-LX,-LY]))
    else:
        reflected = []
        for c in frame.com:
            reflected.append([-c[0],c[1]])
            reflected.append([c[0],-c[1]])
            reflected.append([2.0*LX-c[0],c[1]])
            reflected.append([c[0],2.0*LY-c[1]])
        com = np.vstack((com,reflected))
    vor = Voronoi(com)
    voronoi_plot_2d(vor,ax = engine)

def convex_hull(frame,engine = plt,size=0.3):
    """
    plot the convex hull based on the center of mass
    """
    points = frame.com   # 30 random points in 2-D
    hull = ConvexHull(points)
    engine.plot(points[hull.vertices,0], points[hull.vertices,1], 'r--', lw=2)
    ps = points[hull.vertices[0]]
    pf = points[hull.vertices[-1]]
    engine.plot([ps[0],pf[0]],[ps[1],pf[1]], 'r--',lw=2)
    from numpy import linalg as LA
    S = np.zeros((2,2))
    p = points[hull.vertices[-1]]
    for idx in hull.vertices:
        v = (points[idx] - p).reshape(2,1)
        S += v.dot(v.T)
        p = points[idx]
    w,v = LA.eig(S)
    c = np.mean(frame.com,axis=0)
    for i in range(2):
        p1 = c- size*np.sqrt(w[i])*v[:,i]
        p2 = c+ size*np.sqrt(w[i])*v[:,i]
        engine.plot([p1[0],p2[0]],[p1[1],p2[1]],'r--')

def elliptical_contour(frame,engine = plt,axis = True):
    """
    fit the convex hull from center of mass into an ellipse
    fitting method: http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
    """
    def fitEllipse(x,y):
        # return the coeffcients of the fitted ellipse
        x = x[:,np.newaxis]
        y = y[:,np.newaxis]
        D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
        S = np.dot(D.T,D)
        C = np.zeros([6,6])
        C[0,2] = C[2,0] = 2; C[1,1] = -1
        E, V =  np.linalg.eig(np.dot(np.linalg.inv(S), C))
        n = np.argmax(np.abs(E))
        a = V[:,n]
        return a

    def angle(a):
        # return the orientational angle of this ellipse
        b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
        if b == 0:
            if a > c:
                return 0
            else:
                return np.pi/2
        else:
            if a > c:
                return np.arctan(2*b/(a-c))/2
            else:
                return np.pi/2 + np.arctan(2*b/(a-c))/2
    def semi_axis_length(a):
        # return the length of semi-axis [major,minor]
        b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
        up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
        sqt = np.sqrt((a-c)*(a-c)+4*b*b)
        down1=(b*b-a*c)*(sqt-(c+a))
        down2=(b*b-a*c)*(-sqt-(c+a))
        major=np.sqrt(up/down1)
        minor=np.sqrt(up/down2)
        return[major,minor]

    def ellipse_center(a):
        b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
        num = b*b-a*c
        x0=(c*d-b*f)/num
        y0=(a*f-b*d)/num
        return np.array([x0,y0])

    points = frame.com   # 30 random points in 2-D
    hull = ConvexHull(points)
    x = points[hull.vertices,0]
    y = points[hull.vertices,1]
    a = fitEllipse(x,y)
    theta = angle(a)
    center = ellipse_center(a)
    axis = semi_axis_length(a)

    r = np.arange(0,2.0*np.pi, 0.2)
    xx = center[0] + axis[0]*np.cos(r)*np.cos(theta) - axis[1]*np.sin(r)*np.sin(theta)
    yy = center[1] + axis[0]*np.cos(r)*np.sin(theta) + axis[1]*np.sin(r)*np.cos(theta)
    engine.plot(xx,yy,'r--')
    if axis:
        engine.plot([center[0]-axis[0]*cos(theta),center[0]+axis[0]*cos(theta)],
                  [center[1]-axis[0]*sin(theta),center[1]+axis[0]*sin(theta)],'r--')
        engine.plot([center[0]+axis[1]*sin(theta),center[0]-axis[1]*sin(theta)],
		  [center[1]-axis[1]*cos(theta),center[1]+axis[1]*cos(theta)],'b--')
        
def omegaij(frame):

    N = frame.nphases
    omega_ij = frame.omegaij

    result = np.zeros(int(N*(N-1)/2))

    k = 0

    for i in range(1,N):
        for j in range(i):
            result[k] = omega_ij[i][j]
            k+=1

    return result

def E_CH(frame, i):

    phi = frame.phi[i]
    i_dx,i_dy = grad(frame.phi[i])

    gamma = frame.parameters['gam'][i]

    l = frame.parameters['lambda']

    ones = np.ones(phi.shape)

    phi_squared = np.multiply(phi,phi)
    diff_squared = np.multiply(ones-phi,ones-phi)

    bulk = (gamma/l)*np.multiply(phi_squared,diff_squared)
    gradient = gamma*l*(np.multiply(i_dx,i_dx)+np.multiply(i_dy,i_dy))

    return np.sum(bulk+gradient)

def E_adh(frame, i):
    # returns the adhesive energy of cell i

    E_adh = 0

    walls = frame.parameters['walls']

    walls_dx,walls_dy = grad(walls)
    i_dx,i_dy = grad(frame.phi[i])

    omega = frame.parameters['omega']
    wall_omega = frame.parameters['wall_omega']

    l = frame.parameters['lambda']

    E_density_wall = wall_omega*l*(np.multiply(walls_dx,i_dx)+np.multiply(walls_dy,i_dy))

    E_adh += np.sum(E_density_wall)

    for j in frame.nbr_cells[i]:
        if i != j:
            phi_dx,phi_dy = grad(frame.phi[j])
            E_density_phi = omega*l*(np.multiply(phi_dx,i_dx)+np.multiply(phi_dy,i_dy))
            E_adh += np.sum(E_density_phi)

    return E_adh

def E_area(frame, i):
    # returns area energy of cell i

    phi = frame.phi[i]

    mu = frame.parameters['mu'][i]
    R = frame.parameters['R'][i]

    area = np.sum(np.multiply(phi,phi))
    a0 = np.pi*R*R

    return mu*np.power(1-area/a0,2)

def area(frame, i):
    # returns area energy of cell i

    phi = frame.phi[i]

    mu = frame.parameters['mu'][i]
    R = frame.parameters['R'][i]

    area = np.sum(np.multiply(phi,phi))

    return area

def E_rep(frame, i):
    # returns the repulsive energy of cell i

    E_rep = 0

    phi = frame.phi[i]
    phi_squared = np.multiply(phi,phi)

    walls = frame.parameters['walls']
    walls_squared = np.multiply(walls,walls)

    kappa = frame.parameters['kappa']
    wall_kappa = frame.parameters['wall_kappa']

    l = frame.parameters['lambda']

    E_density_wall = (wall_kappa/l)*np.multiply(phi_squared,walls_squared)

    E_rep += np.sum(E_density_wall)

    for j in frame.nbr_cells[i]:
        if i != j:
            phi_j = frame.phi[j]
            phi_j_squared = np.multiply(phi_j,phi_j)

            E_density_phi = (kappa/l)*np.multiply(phi_squared,phi_j_squared)
            E_rep += np.sum(E_density_phi)

    return E_rep
