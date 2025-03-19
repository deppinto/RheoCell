#!/bin/bash   

START=0
END=100
ITER=1
CONT=0

mkdir Video

for ((ii=$START; ii<=$END; ii=ii+$ITER))
do
        echo "$ii $CONT"
	if (( CONT < 10 ));then
		/home/p/pinto/PythonPackages/bin/python3 /home/p/pinto/Phase_Field/RheoCell/scripts/velocityField.py "../test.top" "velocity_field_00$CONT.txt" "1"
		mv "frame.png" Video/"frame_00$CONT.png"
	elif (( CONT < 100 ));then
		/home/p/pinto/PythonPackages/bin/python3 /home/p/pinto/Phase_Field/RheoCell/scripts/velocityField.py "../test.top" "velocity_field_0$CONT.txt" "1"
		mv "frame.png" Video/"frame_0$CONT.png"
	else
		/home/p/pinto/PythonPackages/bin/python3 /home/p/pinto/Phase_Field/RheoCell/scripts/velocityField.py "../test.top" "velocity_field_$CONT.txt" "1"
		mv "frame.png" Video/"frame_$CONT.png"
	fi
        CONT=$((CONT+1))
done

cd Video
convert -dispose previous -delay 10 *.png animation.gif
#/home/p/pinto/LocalPackages/ffmpeg-7.0.2-amd64-static/ffmpeg -f image2 -framerate 1 -i frame_%003d.png out.mp4
cd ..
