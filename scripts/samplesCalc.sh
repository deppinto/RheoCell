#!/bin/bash

SCRIPT=70
START=1
END=12
INIJOB=1
SEED_VAR=1
#FILENAME="stress_time.txt"
FILENAME="nematic_correlations.txt"

mkdir /home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"$SCRIPT"
cp /scratch/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/params /home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"$SCRIPT"

for ((ii=$START; ii<=$END; ii++))
do

JOBS=`cat /scratch/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/params | sed -n ${ii}','${ii}'p' | awk '{print $1}'`

for ((CURRJOB=$INIJOB; CURRJOB<=JOBS; CURRJOB++))
#for ((CURRJOB=$JOBS; CURRJOB<JOBS; CURRJOB++))
do

if [ $(find /scratch/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/Job_"$SEED_VAR"/ -maxdepth 1 -name "last_conf.dat") ]; then
  echo $SEED_VAR

  mkdir /home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"$SCRIPT"/Job_"$SEED_VAR"/
  cd /scratch/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/Job_"$SEED_VAR"/
  #/home/p/pinto/PythonPackages/bin/python3 /home/p/pinto/Phase_Field/RheoCell/scripts/field_velocity_properties.py "./Analysis/test.top" "./Analysis/Results/" "100" "101" "1" >> $FILENAME
  /home/p/pinto/PythonPackages/bin/python3 /home/p/pinto/Phase_Field/RheoCell/scripts/field_nematic_properties.py "./Analysis/test.top" "./Analysis/Results/" "100" "101" "1" >> $FILENAME
  #/home/p/pinto/PythonPackages/bin/python3 /home/p/pinto/Phase_Field/RheoCell/scripts/field_stress_properties.py "./Analysis/test.top" "./Analysis/Results/" "0" "101" "1" >> $FILENAME
  mv /scratch/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/Job_"$SEED_VAR"/$FILENAME /home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"$SCRIPT"/Job_"$SEED_VAR"/

  #/home/p/pinto/PythonPackages/bin/python3 /home/p/pinto/Phase_Field/RheoCell/scripts/calcAvgVelocity.py "input" "5" "500"
  #/home/p/pinto/PythonPackages/bin/python3 /home/p/pinto/Phase_Field/RheoCell/scripts/calcAvgNematic.py "input" "5" "500"
  #mv /scratch/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/Job_"$SEED_VAR"/*.png /home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"$SCRIPT"/Job_"$SEED_VAR"/
  #mv /home/p/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/Job_"$SEED_VAR"/*.txt /home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"$SCRIPT"/Job_"$SEED_VAR"/
else
  echo "File is not found in Job_"$SEED_VAR""
fi

((SEED_VAR++))

done
done

