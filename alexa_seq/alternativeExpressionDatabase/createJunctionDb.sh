source $HOME/.bashrc

E_BADARGS=65
if [ ! -n "$1" ]
then
  echo "Usage: `basename $0` /home/malachig/svn/alexa_seq /projects/malachig/sequence_databases jango.bcgsc.ca hs_53_36o 62"
  echo "Parameters are: [alexa_seq code path] [annotation database path] [mysql hostname with ALEXA dbs installed] [annotation build name] [boundary length]"
  exit $E_BADARGS
fi

echo
echo "Using args: $1 $2 $3 $4 $5"

SCRIPT_DIR="$1"
WORKING_DIR="$2/temp"
BASE_DIR="$2"
MYSQL_HOST="$3"
DB="$4"
SIZE="$5"

#Make a dir for this DB. e.g. hs_53_36o
DB_DIR="$WORKING_DIR/$DB"
echo
echo "mkdir -p $DB_DIR"
mkdir -p $DB_DIR

#Make a dir for this size. e.g. s62
SDIR="$DB_DIR/s$SIZE"
echo "mkdir -p $SDIR"
mkdir -p $SDIR

#Make junction dir and subdirs: 'logs' and 'temp'
JDIR="$SDIR/exonJunctions"
LOGDIR="$JDIR/logs"
TEMPDIR="$JDIR/temp"
ANNDIR="$LOGDIR/annotateExonJunctions"
echo "mkdir -p $JDIR"
mkdir -p $JDIR
echo "mkdir -p $LOGDIR"
mkdir -p $LOGDIR
echo "mkdir -p $TEMPDIR"
mkdir -p $TEMPDIR
echo "mkdir -p $ANNDIR"
mkdir -p $ANNDIR

#Generate the actual junctions
echo
echo "$SCRIPT_DIR/alternativeExpressionDatabase/createExonJunctionDatabase.pl  --database=ALEXA_$DB  --server=$MYSQL_HOST  --user=viewer  --password=viewer  --target_length=$SIZE  --outdir=$JDIR  --logfile=$LOGDIR/createExonJunctionDatabase_LOG.txt"
$SCRIPT_DIR/alternativeExpressionDatabase/createExonJunctionDatabase.pl  --database=ALEXA_$DB  --server=$MYSQL_HOST  --user=viewer  --password=viewer  --target_length=$SIZE  --outdir=$JDIR  --logfile=$LOGDIR/createExonJunctionDatabase_LOG.txt 1>/dev/null

#Then create jobs annotate them
echo
echo "Make annotation jobs and run"
ANNFILE="$JDIR/annotateExonJunctions.sh"

export DB
export JDIR
export SIZE
export TEMPDIR
export ANNDIR
export MYSQL_HOST
export BASE_DIR
export SCRIPT_DIR

perl -ne '$DB=$ENV{"DB"}; $JDIR=$ENV{"JDIR"}; $SIZE=$ENV{"SIZE"}; $TEMPDIR=$ENV{"TEMPDIR"}; $ANNDIR=$ENV{"ANNDIR"}; $MYSQL_HOST=$ENV{"MYSQL_HOST"}; $BASE_DIR=$ENV{"BASE_DIR"}; $SCRIPT_DIR=$ENV{"SCRIPT_DIR"}; unless($_ =~ /Chromosome/){chomp($_); @line=split("\t", $_); print "$SCRIPT_DIR/alternativeExpressionDatabase/annotateExonJunctions.pl  --database=ALEXA_$DB  --server=$MYSQL_HOST  --user=viewer  --password=viewer  --chr_filter=\"$line[0]:$line[1]:$line[2]-$line[3]\"  --junction_file=$JDIR/exonJunctions\_$SIZE\mers.txt  --ucsc_align_dir=$BASE_DIR/$DB/mrna_est/partitions/   --genbank_mapfile=$BASE_DIR/$DB/mrna_est/GenBankToOrganism.btree  --wiggle=3  --outfile=$TEMPDIR/exonJunctions_annotated_$line[0]_$line[1].txt  --logfile=$ANNDIR/annotateExonJunctions_$line[0]_$line[1]_LOG.txt\n"}' $BASE_DIR/$DB/Regions_250_Genes.txt > $ANNFILE

#Then run these jobs
echo "bash $ANNFILE"
bash $ANNFILE

#Now join the annotation and fasta files together
echo
echo "Now join the annotation and fasta files together"
echo "cd $JDIR" 
cd $JDIR

echo "grep -h -v Junction_ID $JDIR/temp/*.txt | sort -n > $JDIR/tmp.txt"
grep -h -v Junction_ID $JDIR/temp/*.txt | sort -n > $JDIR/tmp.txt

J_FILE_FA=$(echo $JDIR/exonJunctions_$SIZE\mers.fa)
J_FILE=$(echo $JDIR/exonJunctions_$SIZE\mers.txt)
FINAL_FILE=$(echo $JDIR/exonJunctions_$SIZE\mers_annotated.txt)
FINAL_FILE_GZ=$(echo $JDIR/exonJunctions_$SIZE\mers_annotated.txt.gz)

echo "head -q -n 1 $JDIR/temp/*.txt | sort | uniq | cat - $JDIR/tmp.txt > $FINAL_FILE"
head -q -n 1 $JDIR/temp/*.txt | sort | uniq | cat - $JDIR/tmp.txt > $FINAL_FILE

echo "rm -f $JDIR/tmp.txt"
rm -f $JDIR/tmp.txt

echo
echo "Compress results file"
echo "gzip $FINAL_FILE"
gzip $FINAL_FILE
echo "gzip $J_FILE_FA"
gzip $J_FILE_FA

echo
echo "Clean-up temp files..."
echo "cd $JDIR"
cd $JDIR

echo "rm -f $J_FILE"
rm -f $J_FILE

echo "rm -f $JDIR/temp/*.txt"
rm -f $JDIR/temp/*.txt

echo
echo "Add the FID column"
echo "$SCRIPT_DIR/alternativeExpressionDatabase/appendFeatureId.pl  --infile=$FINAL_FILE_GZ  --outfile=$JDIR/tmp.txt  --prefix='EJ'"
$SCRIPT_DIR/alternativeExpressionDatabase/appendFeatureId.pl  --infile=$FINAL_FILE_GZ  --outfile=$JDIR/tmp.txt  --prefix='EJ'




