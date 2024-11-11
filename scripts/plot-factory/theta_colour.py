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
from numpy import sign

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

    cmap = 'viridis'

    #fig.set_size_inches(18,9)
    radius = 4

    ax1 = fig.add_subplot(312)
    ax2 = fig.add_subplot(313)
    ax3 = fig.add_subplot(311)

    size=16
    p=4

    fx,fy = plot.coarse_grain_p(frame,p=p,size=size)
    defects = plot.calculate_defects(fx,fy,p=p)

    plot.p_atic(p, frame,ax3)
    plot.cells(frame,ax3)
    plot.coarse_grain_magnitude(frame,fx,fy,ax1,p=p,colormap=cmap)
    # plot.G_defects(fx,fy,p=4,engine=ax1)
    plot.p_angle(frame,p=p,size=size,engine=ax2,show_def=True)

    for ax in [ax1,ax2,ax3]:
        plot.plot_defects(defects,ax,p)
        ax.set_xlim([0, frame.parameters['Size'][0]-1])
        ax.set_ylim([0, frame.parameters['Size'][1]-1])
        ax.axis('off')

    for ax in [ax1,ax2,ax3]:
        ax.axes.set_aspect('equal', adjustable='box')


if len(oname) == 0:
    animation.animate(ar,myplot,inter = 10,show=True)
else:
    an = animation.animate(ar, myplot,show=False)
    animation.save(an, oname+'.mp4', fps = 10, dpi = 360)
