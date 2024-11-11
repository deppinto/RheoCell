#
# This is a simple example file to show the plotting capabilities of the
# program. Uses python2.
#
# Usage:
#
#   python2 plot-cells.py input [output]
#
#  where
#
#    intput -- the input file or directory
#    output -- (optional) if present saves the animation as a video to show to
#              your mom.
import sys
import matplotlib
import numpy as np
from numpy import sign
from mpl_toolkits.axes_grid1 import AxesGrid

if len(sys.argv) == 3:
    matplotlib.use('Agg')


sys.path.insert(0, "../plot/")

import plot
import archive
import animation
from math import sqrt

##################################################
# Init

if len(sys.argv) == 1:
    print("Please provide an input file.")
    exit(1)

# load archive from file
ar = archive.loadarchive(sys.argv[1])

oname = ""
if len(sys.argv) == 3:
    oname = sys.argv[2]
    print("Output name is", sys.argv[2])


#################################################
# plot simple animation of phases

def myplot(frame, fig):
    matplotlib.rcParams.update({'font.size': 8})
    #fig.set_size_inches(18,9)
    # radius = 4

    # stresses, tractions = plot.cg_stress(frame,size=32)
    traction_x,traction_y = plot.traction(frame)
    stresses = plot.get_stress(traction_x,traction_y,nu=0)
    # print(stresses[0])

    Tx = np.gradient(stresses[0],axis=0) + np.gradient(stresses[2],axis=1)
    Ty = np.gradient(stresses[2],axis=0) + np.gradient(stresses[1],axis=1)

    Q00, Q01 = plot.get_nematic_field(frame.phi, frame.Q00, frame.Q01, size=8)
    defects = plot.get_defects(plot.Q_charge_array(Q00, Q01), Q00, Q01)

    # print(stresses[0])

    Lx = frame.parameters['Size'][0]
    Ly = frame.parameters['Size'][1]

    x = np.arange(0,Lx)
    y = np.arange(0,Ly)
    X, Y = np.meshgrid(x,y)

    grid = AxesGrid(fig, 111,
                nrows_ncols=(2, 3),
                axes_pad=0.05,
                share_all=True,
                label_mode="L",
                cbar_location="right",
                cbar_mode="single",
                )

    for stress, ax in zip(stresses,grid[0:3]):
        im = ax.imshow(stress)
        # plot.cells(frame,ax)
        plot.plot_defects(defects,ax)
        # ax.axes.set_aspect('equal', adjustable='box')
        # ax.set_xlim([25,75])
        # ax.set_ylim([40,60])
        ax.set_xlim([0, frame.parameters['Size'][0]-1])
        ax.set_ylim([0, frame.parameters['Size'][1]-1])

    # plot.cells(frame,engine=grid[3])
    # plot.nematic(frame,engine=grid[3])
    # plot.plot_defects(defects,engine=grid[4])
    # plot.nematic_field(frame,engine=grid[4],size=36,avg=5,show_def=True,arrow_len=0)

    grid[3].quiver(X,Y,Tx,Ty)
    grid[5].quiver(X,Y,0.1*traction_x,0.1*traction_y)
    # plot.force_density(frame,engine=grid[5],force_type='dipole')

    grid[0].set_title('xx stress')
    grid[1].set_title('yy stress')
    grid[2].set_title('xy stress')


    grid.cbar_axes[0].colorbar(im)

if len(oname) == 0:
    animation.animate(ar,myplot,inter = 1,show=True)
else:
    an = animation.animate(ar, myplot,show=False)
    animation.save(an, oname+'.mp4', 10, dpi = 360)
