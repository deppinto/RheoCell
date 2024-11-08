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
import matplotlib.gridspec as gridspec
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
    #fig.set_size_inches(18,9)
    radius = 4
    gs = gridspec.GridSpec(3,2)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,:])
    ax4 = fig.add_subplot(gs[2,:])

    p = int(5-sign(frame.parameters['zetaQ'][0]))
    plot.cells(frame,ax1)
    plot.p_atic(p, frame, ax1)
    # plot.velocity(frame,ax1,color = 'b')

    plot.cells(frame,ax2)
    plot.velocity_field(frame,radius,ax2)
    # plot.p_atic(2, frame,ax2)

    # plot.average_velocity_field(frame,0,ax3)
    # ax3.set_xlabel('y coordinate')
    # ax3.set_ylabel('average x velocity')

    # plot.average_velocity_field(frame,1,ax4)
    # ax4.set_xlabel('x coordinate')
    # ax4.set_ylabel('average y velocity')

    plot.average_velocity_polar(frame,0,ax3)
    ax3.set_xlabel('radial coordinate')
    ax3.set_ylabel('average radial velocity')

    plot.average_velocity_polar(frame,1,ax4)
    ax4.set_xlabel('radial coordiante')
    ax4.set_ylabel('average azimuthal velocity')

    # plot.nematic(frame,ax2)
    # plot.v_com(frame, ax1, color = 'r')
    #plot.avg_velocity(frame,ax1)
    #plot.polarization(frame,ax1)
    #plot.shape(frame,ax1)
    #plot.nematic_force(frame,engine = ax1)
    #plot.trajectories(frame,ax1,color = 'r')
    #plot.elliptical_contour(frame,ax1)
    #plot.convex_hull(frame,ax1)

    for ax in [ax1,ax2]:
        ax.axes.set_aspect('equal', adjustable='box')
        ax.set_xlim([0, frame.parameters['Size'][0]-1])
        ax.set_ylim([0, frame.parameters['Size'][1]-1])
        # ax.axis('off')

        x0,y0,width,height = ax.get_position().bounds
        # print(width)

        x = 0.5*(1-width)

        # ax.set_position([x,y0,width,height])

if len(oname) == 0:
    animation.animate(ar,myplot,inter = 1,show=True)
else:
    an = animation.animate(ar, myplot,show=False)
    animation.save(an, oname+'.mp4', 10, dpi = 360)
