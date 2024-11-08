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

sys.path.insert(0, "../plot/")

import plot
import archive
import animation

##################################################
# Init

if len(sys.argv) == 1:
    print("Please provide an input file.")
    exit(1)

# load archive from file
ar = archive.loadarchive(sys.argv[1])

oname = ""
if len(sys.argv) == 3:
    oname = "movie_"+sys.argv[2]
    print("Output name is", sys.argv[2])


##################################################
# plot simple animation of phases


def myplot(frame, fig):

    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    plot.cells(frame, ax1)
    #plot.vorticity_field(frame, size=8, engine=ax1, cbar=True,show_center = True)
    #plot.velocity(frame, ax1)
    #plot.nematic_field(frame, engine = ax2,avg = 1,show_def = True)
    plot.voronoi_lattice(frame,ax2)
    plot.cells(frame,ax3)
    plot.voronoi_lattice(frame,ax3)
    plot.com(frame, ax1,plotIndex = True)
    plot.com(frame, ax2,plotIndex = True)
    #plot.shape_field(frame, engine=ax2, size = 16,avg=1, show_def=True)
    plot.show_neighbour(frame,engine = ax1,cell_index=6)
    #plot.force_density(frame, engine=ax1, avg =1,type = 'active',magn = False)
    #plot.force_density(frame, engine=ax2, type = 'passive',magn = False)
    #plot.velocity_field(frame, engine=ax1,size =0.1, magn = False, scale = False)
    for ax in [ax1, ax2, ax3]:
        ax.axes.set_aspect('equal', adjustable='box')
        ax.set_xlim([0, frame.parameters['Size'][0]-1])
        ax.set_ylim([0, frame.parameters['Size'][1]-1])
        ax.axis('off')


if len(oname) == 0:
    animation.animate(ar, myplot, show=True)
else:
    an = animation.animate(ar, myplot, show=False)
    animation.save(an, oname+'.mp4', 5)
