#!/bin/bash

SCRIPT=74
START=11
END=15
INIJOB=1
SEED_VAR=101
#FILENAME="stress_time.txt"
FILENAME="hole_avg.txt"

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
  #/home/p/pinto/PythonPackages/bin/python3 /home/p/pinto/Phase_Field/RheoCell/scripts/hole_dynamics.py "./test.top" "./trajectory.dat" "0" "100" "1"
  #/home/p/pinto/PythonPackages/bin/python3 /home/p/pinto/Phase_Field/RheoCell/scripts/field_stress_properties.py "./test.top" "./Analysis/Results/" "./trajectory.dat" "0" "100" "1"
  /home/p/pinto/PythonPackages/bin/python3 /home/p/pinto/Phase_Field/RheoCell/scripts/field_velocity_trajectory_properties.py "./test.top" "./Analysis/Results/" "./trajectory.dat" "0" "100" "4" >> "voids_strain_defectsplusone.txt"
  #/home/p/pinto/PythonPackages/bin/python3 /home/p/pinto/Phase_Field/RheoCell/scripts/field_velocity_properties.py "./Analysis/test.top" "./Analysis/Results/" "100" "101" "1" >> $FILENAME
  #/home/p/pinto/PythonPackages/bin/python3 /home/p/pinto/Phase_Field/RheoCell/scripts/field_nematic_properties.py "./Analysis/test.top" "./Analysis/Results/" "100" "101" "1" >> $FILENAME
  #/home/p/pinto/PythonPackages/bin/python3 /home/p/pinto/Phase_Field/RheoCell/scripts/field_stress_properties.py "./Analysis/test.top" "./Analysis/Results/" "0" "101" "1" >> $FILENAME
  #mv /scratch/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/Job_"$SEED_VAR"/$FILENAME /home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"$SCRIPT"/Job_"$SEED_VAR"/
  mv /scratch/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/Job_"$SEED_VAR"/voids_* /home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"$SCRIPT"/Job_"$SEED_VAR"/

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

