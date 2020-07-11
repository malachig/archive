#!/bin/bash

# EXAMPLE COMMAND:
#bsub -q oncology -G compute-oncology -g /mgriffit/cwl -M 10000000 -R 'select[mem>10000] rusage[mem=10000]' -oo /storage1/fs1/mgriffit/Active/MyCancerDb/cwl_runs/JLF001_NewTumor_vs_NewNormal/logs/attempt-XX.log -a 'docker0(registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:20)' /bin/bash /storage1/fs1/mgriffit/Active/common/runCromwellCWL.sh --cromwell_config /storage1/fs1/mgriffit/Active/common/cromwell.config --sample JLF001_NewTumor_vs_NewNormal --cwl /storage1/fs1/mgriffit/Active/MyCancerDb/git/analysis-workflows/definitions/pipelines/immuno.cwl --yaml /storage1/fs1/mgriffit/Active/MyCancerDb/cwl_runs/JLF001_NewTumor_vs_NewNormal/JLF001_NewTumor_vs_NewNormal.yaml --results /storage1/fs1/mgriffit/Active/MyCancerDb/cwl_runs/JLF001_NewTumor_vs_NewNormal/result --cromwell_jar /storage1/fs1/mgriffit/Active/common/cromwell-jars/cromwell-51.jar

function usage
{
    echo ""
    echo "usage: runCromwellCWL.sh -q <queue> -m <memory> -a <docker image> -h"
    echo ""
    echo "  -g | --cromwell_config    Path to cromwell config file"
    echo "  -s | --sample             Sample name"
    echo "  -c | --cwl                Path to CWL pipeline file"
    echo "  -y | --yaml               Path to input YAML file"
    echo "  -r | --results            Path to final results dir where names outputs of the pipeline will be placed"
    echo "  -j | --cromwell_jar       Path to cromwell jar file (default: /storage1/fs1/mgriffit/Active/common/cromwell-jars/cromwell-51.jar)"
    echo "  -n | --clean              Whether to clean up or not (default: YES - things will be cleaned up unless you say --clean NO)"
    echo "  -h | --help               Show this message"
    echo ""
    exit 1
}

if [[ $1 == "" ]];then
    usage
    exit;
fi

while [[ "$1" != "" ]]; do
    case $1 in
         -g | --cromwell_config ) shift
                                  cromwell_config=$1
				  ;;
         -s | --sample )           shift
                                  sample=$1
                                  ;;
        -c | --cwl )              shift
                                  cwl=$1
                                  ;; 
        -y | --yaml )             shift
                                  yaml=$1
                                  ;;
        -r | --results )          shift
                                  results=$1
                                  ;;
        -j | --cromwell_jar )     shift
                                  cromwell_jar=$1
				  ;;
        -n | --clean )            shift
                                  clean=$1
				  ;;
       -h | --help )              usage
                                  exit
                                  ;;
        * )                       usage
                                  exit 1
    esac
    shift
done

if [[ $cromwell_config == "" ]];then
    echo "--cromwell_config must be specified"
    exit;
fi
if [[ $sample == "" ]];then
    echo "--sample must be specified"
    exit;
fi
if [[ $cwl == "" ]];then
    echo "--cwl must be specified"
    exit;
fi
if [[ $yaml == "" ]];then
    echo "--yaml must be specified"
    exit;
fi
if [[ $results == "" ]];then
    echo "--results must be specified"
    exit;
fi
if [[ $cromwell_jar = "" ]];then
    echo "using cromwell jar: /storage1/fs1/mgriffit/Active/common/cromwell-jars/cromwell-51.jar"
    cromwell_jar=/storage1/fs1/mgriffit/Active/common/cromwell-jars/cromwell-51.jar;
fi
if [[ $clean == "" ]];then
    echo "Temp files will be cleaned up"
    clean="YES";
fi

###########################################################################################
############################# pre-setup ###################################################
###########################################################################################

# get base name of the CWL pipeline being run
CWL_BASE=$(basename $cwl)

exit

###########################################################################################
###################start up the cromwell server/start job #################################
###########################################################################################

# start cromwell server and give it time to setup
/usr/bin/java -Dconfig.file=$cromwell_config -jar $cromwell_jar server &
sleep 60

# create a yaml to label the cromwell job
echo -e "{\n\"model\":\"$sample\"\n}" >| $sample.label

# submit the cromwell job
/usr/bin/java -Dconfig.file=$cromwell_config -jar $cromwell_jar submit -h http://localhost:8000 -l $sample.label -t cwl -i $yaml $cwl

############################################################################################
################# query the cromwell server for status #####################################
############################################################################################

# infinity loop to check job status
x=0
while [ $x -le 1 ]
do
    curl -SL http://localhost:8000/api/workflows/v1/query?label=model:$sample >| $sample.status
    sleep 600
    if cat $sample.status | python -m json.tool | grep -q "Succeeded"; then
        break
    elif cat $sample.status | python -m json.tool | grep -q "Failed"; then
        exit 1
    else
        continue
    fi
done

############################################################################################
################ grab the final outputs and put them in the correct place ##################
############################################################################################

# with the cromwell job complete get the id so we can query the outputs
CROMWELL_ID="$(cat $sample.status | python -m json.tool | grep "\"id\":" | sed 's@.*\"id\": \"\(.*\)\".*@\1@')"

# when done continue and grab the outputs
curl -SL http://localhost:8000/api/workflows/v1/$CROMWELL_ID/outputs >| $sample.output

# loop through the output and put everything in the right place
cat $sample.output | python -m json.tool | grep location | sed 's@.*\"location\": \"\(.*\)\".*@\1@' >| $sample.final_results

for line in $(cat $sample.final_results)
do
    if [ $clean == "NO" ]; then
        echo "Copying results file to results dir:" $line
        cp -r $line $results
    else
        echo "Moving results file to results dir:" $line
	mv $line $results
    fi
done

#############################################################################################
################ with everything now done clean up after yourself ###########################
#############################################################################################

rm -f $sample.final_results
rm -f $sample.output
rm -f $sample.status
rm -f $sample.label

if [ $clean == "NO" ]; then
    echo "Leaving full cromwell-executions run dir in place"
else
    echo "Removing cromwell-executions dir"
    rm -rf cromwell-executions/$CWL_BASE/$CROMWELL_ID
fi

