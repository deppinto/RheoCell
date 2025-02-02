#!/bin/bash

SCRIPT=43
START=1
END=5
INIJOB=1
SEED_VAR=1

mkdir /home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"$SCRIPT"
cp /scratch/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/params /home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"$SCRIPT"

for ((ii=$START; ii<=$END; ii++))
do

JOBS=`cat /scratch/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/params | sed -n ${ii}','${ii}'p' | awk '{print $1}'`

for ((CURRJOB=$INIJOB; CURRJOB<=JOBS; CURRJOB++))
#for ((CURRJOB=$JOBS; CURRJOB<JOBS; CURRJOB++))
do

if [ $(find /scratch/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/Job_"$SEED_VAR"/ -name "last_conf.dat") ]; then
  echo $SEED_VAR
  cd /scratch/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/Job_"$SEED_VAR"/
  /home/p/pinto/PythonPackages/bin/python3 /home/p/pinto/Phase_Field/RheoCell/scripts/calcAvgVelocity.py "input" "5" "25"
  #/home/p/pinto/PythonPackages/bin/python3 /home/p/pinto/Phase_Field/RheoCell/scripts/calcAvgNematic.py "input" "5" "25"
  mkdir /home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"$SCRIPT"/Job_"$SEED_VAR"/
  mv /scratch/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/Job_"$SEED_VAR"/*.png /home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"$SCRIPT"/Job_"$SEED_VAR"/
  #mv /home/p/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/Job_"$SEED_VAR"/*.txt /home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"$SCRIPT"/Job_"$SEED_VAR"/
else
  echo "File is not found in Job_"$SEED_VAR""
fi

((SEED_VAR++))

done
done

