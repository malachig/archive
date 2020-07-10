#!/bin/bash

#!!!!!!!!!!!!!!!!!!!!!!!!!! NOTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# cromwell.config location is hardcoded below
# you must also run from a docker image with python, and the other sys dependencies ...
# docker(registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:20) should work
#!!!!!!!!!!!!!!!!!!!!!!!!!! NOTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

# EXAMPLE COMMAND:
#bsub -q oncology -G compute-oncology -g /mgriffit/cwl -M 10000000 -R 'select[mem>10000] rusage[mem=10000]' -oo /storage1/fs1/mgriffit/Active/MyCancerDb/cwl_runs/JLF001_NewTumor_vs_NewNormal/logs/attempt-XX.log -a 'docker0(registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:20)' /bin/bash /storage1/fs1/mgriffit/Active/common/GMS_lite_noclean.sh JLF001_NewTumor_vs_NewNormal /storage1/fs1/mgriffit/Active/MyCancerDb/git/analysis-workflows/definitions/pipelines/immuno.cwl /storage1/fs1/mgriffit/Active/MyCancerDb/cwl_runs/JLF001_NewTumor_vs_NewNormal/JLF001_NewTumor_vs_NewNormal.yaml /storage1/fs1/mgriffit/Active/MyCancerDb/cwl_runs/JLF001_NewTumor_vs_NewNormal/result /storage1/fs1/mgriffit/Active/common/cromwell-jars/cromwell-51.jar NOCLEANUP


# grab the command line arguments
SAMPLE=$1
CWL=$2
CWL_YAML=$3
RESULTS_DIR=$4
CROMWELL_JAR_LOC=$5
CLEAN=$6


###########################################################################################
############################# pre-setup ###################################################
###########################################################################################

# get base name of the CWL pipeline being run
CWL_BASE=$(basename $CWL)

###########################################################################################
###################start up the cromwell server/start job #################################
###########################################################################################


# start cromwell server and give it time to setup
/usr/bin/java -Dconfig.file=/storage1/fs1/mgriffit/Active/common/cromwell.config -jar $CROMWELL_JAR_LOC server &
sleep 60

# create a yaml to label the cromwell job
echo -e "{\n\"model\":\"$SAMPLE\"\n}" >| $SAMPLE.label

# submit the cromwell job
/usr/bin/java -Dconfig.file=/storage1/fs1/mgriffit/Active/common/cromwell.config -jar $CROMWELL_JAR_LOC submit -h http://localhost:8000 -l $SAMPLE.label -t cwl -i $CWL_YAML $CWL

############################################################################################
################# query the cromwell server for status #####################################
############################################################################################


# infinity loop to check job status
x=0
while [ $x -le 1 ]
do
    curl -SL http://localhost:8000/api/workflows/v1/query?label=model:$SAMPLE >| $SAMPLE.status
    sleep 600
    if cat $SAMPLE.status | python -m json.tool | grep -q "Succeeded"; then
        break
    elif cat $SAMPLE.status | python -m json.tool | grep -q "Failed"; then
        exit 1
    else
        continue
    fi
done

############################################################################################
################ grab the final outputs and put them in the correct place ##################
############################################################################################

# with the cromwell job complete get the id so we can query the outputs
CROMWELL_ID="$(cat $SAMPLE.status | python -m json.tool | grep "\"id\":" | sed 's@.*\"id\": \"\(.*\)\".*@\1@')"

# when done continue and grab the outputs
curl -SL http://localhost:8000/api/workflows/v1/$CROMWELL_ID/outputs >| $SAMPLE.output

# loop through the output and put everything in the right place
cat $SAMPLE.output | python -m json.tool | grep location | sed 's@.*\"location\": \"\(.*\)\".*@\1@' >| $SAMPLE.final_results

for line in $(cat $SAMPLE.final_results)
do
    if [ $CLEAN == "NOCLEANUP" ]; then
        echo "Copying results file to results dir:" $line
        cp -r $line $RESULTS_DIR
    else
        echo "Moving results file to results dir:" $line
	mv $line $RESULTS_DIR
    fi
done

#############################################################################################
################ with everything now done clean up after yourself ###########################
#############################################################################################

rm -f $SAMPLE.final_results
rm -f $SAMPLE.output
rm -f $SAMPLE.status
rm -f $SAMPLE.label

if [ $CLEAN == "NOCLEANUP" ]; then
    echo "Leaving full cromwell-executions run dir in place"
else
    echo "Removing cromwell-executions dir"
    rm -rf cromwell-executions/$CWL_BASE/$CROMWELL_ID
fi
