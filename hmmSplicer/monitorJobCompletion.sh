#!/bin/bash
#Written by Malachi Griffith

# This script counts the jobs on the cluster for a user every N seconds
E_BADARGS=65
if [ ! -n "$1" ]
then
  echo "Usage: `basename $0` /projects/alexa2/hmmSplicer/ "AML,DLBCL,HumanBodyMap,Oligo,REMC,SA_TN_Breast" (i.e. analysis dir; list of project names)"
  exit $E_BADARGS
fi

ANALYSISPATH=$1
LIST=$2

IFS=","
PROJECTS=( $LIST )

echo Examining each of the following projects: ${PROJECTS[@]}
echo 

for PROJECT in ${PROJECTS[@]}
do
  COMPLETE=$(tail $ANALYSISPATH/$PROJECT/logs/* | grep COMPLETE | wc -l)
  FAILED=$(tail $ANALYSISPATH/$PROJECT/logs/* | grep FAILED | wc -l)
  echo $PROJECT = $COMPLETE COMPLETE and $FAILED FAILED jobs
done

