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
import os 
workspace = '/media/richard/Elements1/contractile/'
py_file = 'plot.py
'
ndir = 25
for i in range(ndir):
    input_dir = workspace + str(i+1) + '/data'
    output_name = str(i+1)
    os.system('addqueue -q long ./make-movie-bash.sh {} {} {}'.format(py_file,input_dir,output_name))
