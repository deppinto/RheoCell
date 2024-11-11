#!/bin/bash

DIR=${1?Error: No Input Directory}

FILES=../../movies/$DIR

cd $FILES

# for f in *;
for f in 21-11-12_11-17-29 21-11-12_11-34-09 21-11-12_11-52-46 21-11-12_12-09-58 21-11-12_12-26-39 21-11-12_12-43-35 21-11-12_12-59-51 21-11-12_13-16-20 21-11-12_13-32-49 21-11-12_13-49-48 21-11-12_14-06-07 21-11-12_14-22-48;
do

  # echo $f\

  # pwd
  echo $f
  cd ../../celadro/plot-factory
  python3 plot-cells-directors.py $FILES/$f/data $FILES/$f/directors_$f
  # python3 autocorr.py ../../movies/phases//data ../../movies/phases//autocorr$f ../../models/phases/persistence_time.csv
  # python3 center-of-mass-trajectory.py $FILES/$f/data $FILES/$f/com_$f
  cd $FILES

done

cd ..
