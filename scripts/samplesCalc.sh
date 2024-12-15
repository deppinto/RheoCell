#!/bin/bash

SCRIPT=19
START=1
END=8
INIJOB=1
SEED_VAR=1

mkdir /home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"$SCRIPT"
cp /home/p/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/params /home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"$SCRIPT"

for ((ii=$START; ii<=$END; ii++))
do

JOBS=`cat /home/p/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/params | sed -n ${ii}','${ii}'p' | awk '{print $1}'`

for ((CURRJOB=$INIJOB; CURRJOB<=JOBS; CURRJOB++))
#for ((CURRJOB=$JOBS; CURRJOB<JOBS; CURRJOB++))
do

if [ $(find /home/p/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/Job_"$SEED_VAR"/ -name "last_conf.dat") ]; then
  echo $SEED_VAR
  cd /home/p/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/Job_"$SEED_VAR"/
  /home/p/pinto/PythonPackages/bin/python3 /home/p/pinto/Phase_Field/RheoCell/scripts/calcAvgVelocity.py "input" "6" "20"
  mkdir /home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"$SCRIPT"/Job_"$SEED_VAR"/
  #mv /home/p/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/Job_"$SEED_VAR"/*.png /home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"$SCRIPT"/Job_"$SEED_VAR"/
  mv /home/p/pinto/Phase_Field/RheoCell/Work/Results/scripts"$SCRIPT"/Job_"$SEED_VAR"/*.txt /home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"$SCRIPT"/Job_"$SEED_VAR"/
else
  echo "File is not found in Job_"$SEED_VAR""
fi

((SEED_VAR++))

done
done

