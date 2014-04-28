source $HOME/.bashrc

E_BADARGS=65
if [ ! -n "$1" ]
then
  echo "Usage: `basename $0` /home/malachig/svn/alexa_seq /projects/malachig/sequence_databases/temp jango.bcgsc.ca hs_53_36o 60 53"
  echo "Parameters needed: [alexa_seq code path] [working dir] [mysql hostname] [annotation build name] [junction seq length] [ensembl version]"
  exit $E_BADARGS
fi

echo
echo "Using args: $1 $2 $3 $4 $5 $6"

SCRIPT_DIR="$1"
WORKING_DIR="$2"
MYSQL_HOST="$3"
DB="$4"
SIZE="$5"
VERSION="$6"

DIR=$(echo $WORKING_DIR/$DB/s$SIZE/exonJunctions)
J1=$(echo $DIR/exonJunctions_$SIZE\mers_annotated.txt.gz)
J2=$(echo $DIR/exonJunctions_$SIZE\mers.fa.gz)
J3=$(echo $DIR/exonJunctions_$SIZE\mers_annotated.txt)

echo 
echo cd $DIR
cd $DIR

echo $SCRIPT_DIR/alternativeExpressionDatabase/identifyJunctionFrames.pl  --database=ALEXA_$DB  --server=$MYSQL_HOST  --user=viewer  --password=viewer  --junction_db=$J1  --junction_seq_file=$J2  --outfile=$J3   --ensembl_version=$VERSION
$SCRIPT_DIR/alternativeExpressionDatabase/identifyJunctionFrames.pl  --database=ALEXA_$DB  --server=$MYSQL_HOST  --user=viewer  --password=viewer  --junction_db=$J1  --junction_seq_file=$J2  --outfile=$J3   --ensembl_version=$VERSION

