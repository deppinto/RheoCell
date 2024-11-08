#!/bin/bash

# DIR=${1?Error:No Input Directory}

FILES=../../movies/channel/omega=0.4/LY=104
# cd $FILES

for f  in 22-07-16_19-27-58 22-07-17_03-49-55 22-07-19_14-34-38 22-07-19_22-54-26;do# LY=87/22-07-10_20-53-34 LY=87/22-07-11_04-24-51 LY=87/22-07-13_09-21-27 LY=87/22-07-13_16-53-11;do
echo $f

# python3 force_density.py $FILES/$f/data $FILES/$f/density_$f
python3 strip_velocity.py $FILES/$f/data $FILES/$f/strip_$f
# python3 plot-cell.py $FILES/$f/data $FILES/$f/cell_$f

cd ../scripts
# python3 center-of-mass-trajectory.py $FILES/$f/data $FILES/$f/com_$f
cd ../plot-factory

done

cd ..
