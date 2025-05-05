#!/bin/bash   

START=1
END=$(cat trajectory.dat | wc -l)
ITER=$(cat test.top | head -1 | awk '{print $1;}')
ITER=$((ITER+2))
iii=$((START+ITER-1))
CONT=1

#mkdir Analysis
#cp test.top Analysis
#cp input Analysis
#cp external.conf Analysis
#cd Analysis
mkdir Video

cp '../equilibrated_conf.dat' 'start.conf'
FILENAME='start.conf'
#/home/p/pinto/Phase_Field/RheoCell/build/bin/analyser input
/home/p/pinto/PythonPackages/bin/python3 /home/p/pinto/Phase_Field/RheoCell/scripts/stressField.py "test.top" $FILENAME "1"
mv "frame.png" Video/"frame_000.png" 

for ((ii=$START; ii<=$END; ii=ii+$ITER))
do
        echo "$ii $iii $CONT"
        head -$iii '../trajectory.dat' | tail +$ii >> tst.dump
	mv 'tst.dump' 'start.conf'
	/home/p/pinto/Phase_Field/RheoCell/build/bin/analyser input
	/home/p/pinto/PythonPackages/bin/python3 /home/p/pinto/Phase_Field/RheoCell/scripts/stressField.py "test.top" $FILENAME "1"
	if (( CONT < 10 ));then
		mv "frame.png" Video/"frame_00$CONT.png"
	elif (( CONT < 100 ));then
		mv "frame.png" Video/"frame_0$CONT.png"
	else
		mv "frame.png" Video/"frame_$CONT.png"
	fi
        iii=$((iii+ITER))
        CONT=$((CONT+1))
        #rm tst.dump
done

cd Video
convert -dispose previous -delay 10 *.png animation.gif
#/home/p/pinto/LocalPackages/ffmpeg-7.0.2-amd64-static/ffmpeg -f image2 -framerate 1 -i frame_%003d.png out.mp4
cd ..
