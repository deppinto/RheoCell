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
    ax1 = fig.add_subplot(111)
    
    
    ax1.title.set_text('polar')
    plot.force_density(frame, engine=ax1,avg =1,size = 1,force_type = 'polar',cbar=True,cg='lattice')
    
    for ax in [ax1]:
        ax.axes.set_aspect('equal', adjustable='box')
        ax.set_xlim([0, frame.parameters['Size'][0]-1])
        ax.set_ylim([0, frame.parameters['Size'][1]-1])
        ax.axis('off')


if len(oname) == 0:
    animation.animate(ar,myplot,inter = 1,show=True)
else:
    an = animation.animate(ar, myplot,show=False)
    animation.save(an, oname+'.mp4', 10, dpi = 360)

