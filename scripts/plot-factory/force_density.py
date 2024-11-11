#
# This is a simple example file to show the plotting capabilities of the
# program. Uses python3.
#
# Usage:
#
#   python3 dipolar.py input [output]
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
    #fig.set_size_inches(18,9)
    radius = 4
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(224)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(222)

    plot.force_density(frame,engine=ax1,force_type='all',magn=True,cbar=True,width=1,step=1,cg='lattice')
    ax1.set_title('all forces')
    plot.force_density(frame,engine=ax2,force_type='passive',magn=True,cbar=True,width=1,step=1,cg='lattice')
    ax2.set_title('passive forces')
    plot.force_density(frame,engine=ax3,force_type='dipole',magn=True,cbar=True,width=1,step=1,cg='lattice')
    ax3.set_title('dipolar forces')

    plot.velocity_field(frame,engine=ax4,cg=None,step=1,width=1)
    ax4.set_title('velocity field')

    # plot.velocity_field(frame,engine=ax3,avg=7,size=7)

    for ax in [ax1,ax2,ax3,ax4]:
        ax.axes.set_aspect('equal', adjustable='box')
        # ax.set_xlim([0, frame.parameters['Size'][0]-1])
        ax.set_xlim([0,50])
        ax.set_ylim([0, frame.parameters['Size'][1]-1])
        # ax.axis('off')

        # x0,y0,width,height = ax.get_position().bounds

        # print(width)
        # print(width)
        #
        # x = 0.5*(1-width)
        #
        # ax.set_position([x,y0,width,height])

    # for ax in [ax1,ax2,ax3]:
        # ax.axes.set_aspect('equal', adjustable='box')
        # ax.set_xlim([0, frame.parameters['Size'][0]-1])
        # ax.set_xlim([0,50])
        # ax.set_ylim([0, frame.parameters['Size'][1]-1])
        # ax.axis('off')

        # x0,y0,width,height = ax.get_position().bounds

        # print(width)
        # print(width)
        #
        # x = 0.5*(1-width)
        #
        # ax.set_position([x,y0,width*1.2,height*1.2])

if len(oname) == 0:
    animation.animate(ar,myplot,inter = 1,show=True)
else:
    an = animation.animate(ar, myplot,show=False)
    animation.save(an, oname+'.mp4', 10, dpi = 360)
