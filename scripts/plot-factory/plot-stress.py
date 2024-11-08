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

    stresses, tractions = plot.cg_stress(frame,size=7)
    # print(stresses[0])

    Lx = frame.parameters['Size'][0]
    Ly = frame.parameters['Size'][1]

    # x = np.arange(0,Lx)
    # y = np.arange(0,Ly)
    # X, Y = np.meshgrid(x,y)

    grid = AxesGrid(fig, 111,
                nrows_ncols=(1, 3),
                axes_pad=0.05,
                share_all=True,
                label_mode="L",
                cbar_location="right",
                cbar_mode="single",
                )

    for stress, ax in zip(stresses,grid):
        im = ax.imshow(stress, vmin=-0.002, vmax=0.002)
        plot.cells(frame,ax)
        # ax.quiver(X,Y,tractions[0],tractions[1])
        # ax.axes.set_aspect('equal', adjustable='box')
        ax.set_xlim([25,75])
        ax.set_ylim([40,60])
        # ax.set_xlim([0, frame.parameters['Size'][0]-1])
        # ax.set_ylim([0, frame.parameters['Size'][1]-1])

    grid.cbar_axes[0].colorbar(im)

if len(oname) == 0:
    animation.animate(ar,myplot,inter = 1,show=True)
else:
    an = animation.animate(ar, myplot,show=False)
    animation.save(an, oname+'.mp4', 10, dpi = 360)
