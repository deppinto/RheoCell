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
import matplotlib
import sys
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
    fig.set_size_inches(18,9)
    radius = 8
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    ax1.title.set_text('red-velocity\nblue-force')
    plot.cells(frame,ax1)
    plot.velocity(frame,ax1)
    plot.nematic(frame,ax1)

    ax2.title.set_text('director-deformation')
    plot.shape_field(frame,avg = 4,size = radius*4,show_def = True, engine = ax2)
    # plot.nematic_field(frame,avg = 4,size = radius*4,show_def = True, engine = ax2)


    for ax in [ax1,ax2]:
        ax.axes.set_aspect('equal', adjustable='box')
        ax.set_xlim([0, frame.parameters['Size'][0]-1])
        ax.set_ylim([0, frame.parameters['Size'][1]-1])
        ax.axis('off')


if len(oname) == 0:
    animation.animate(ar,myplot,inter = 1,show=True)
else:
    an = animation.animate(ar, myplot, show=False)
    animation.save(an, oname+'.mp4', 10, dpi = 360)
