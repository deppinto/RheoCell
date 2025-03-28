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
    radius = 2

    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(223)
    ax3 = fig.add_subplot(222)
    ax4 = fig.add_subplot(224)

    p_old = 2
    fx,fy = plot.coarse_grain_p(frame,p=p_old,size=radius)
    plot.coarse_grain_magnitude(frame,fx,fy,ax1,p=p_old,colormap=cmap)
    plot.coarse_grain_theta(frame,fx,fy,show_def=True,engine=ax2,p=p_old)

    p_new = int(5-sign(frame.parameters['zetaQ'][0]))
    fx,fy = plot.coarse_grain_p(frame,p=p_new,size=radius)
    plot.coarse_grain_magnitude(frame,fx,fy,ax3,p=p_new,colormap=cmap)
    plot.coarse_grain_theta(frame,fx,fy,show_def=True,engine=ax4,p=p_new)


    for ax in [ax1,ax2,ax3,ax4]:
        ax.set_xlim([0, frame.parameters['Size'][0]-1])
        ax.set_ylim([0, frame.parameters['Size'][1]-1])
        ax.axis('off')
        ax.axes.set_aspect('equal', adjustable='box')

    # for ax in [ax1,ax2]:
    #     ax.axes.set_aspect('equal', adjustable='box')


if len(oname) == 0:
    animation.animate(ar,myplot,inter = 10,show=True)
else:
    an = animation.animate(ar, myplot,show=False)
    animation.save(an, oname+'.mp4', fps = 10, dpi = 360)
