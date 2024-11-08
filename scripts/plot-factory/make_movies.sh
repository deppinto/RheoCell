#!/bin/bash
#firts input is path
DIR=$1
echo "directory:$1"
for DIR0 in $(ls $DIR); do
    WD="$DIR/$DIR0"
    if [ -d "$WD" ]; then
        output_name="${WD}/movie.mp4"
        data_dir="${WD}/data"
        echo $data_dir
        addqueue -q long -m 1 /usr/bin/python3 movie-6.py $data_dir $output_name
    fi
done
