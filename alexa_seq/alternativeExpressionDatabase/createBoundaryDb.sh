#!/bin/bash
#source $HOME/.bashrc

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

#Make boundary dir and subdirs: 'logs' and 'temp'
BDIR="$SDIR/exonBoundaries"
LOGDIR="$BDIR/logs"
TEMPDIR="$BDIR/temp"
ANNDIR="$LOGDIR/annotateExonBoundaries"
echo "mkdir -p $BDIR"
mkdir -p $BDIR
echo "mkdir -p $LOGDIR"
mkdir -p $LOGDIR
echo "mkdir -p $TEMPDIR"
mkdir -p $TEMPDIR
echo "mkdir -p $ANNDIR"
mkdir -p $ANNDIR

#Generate the actual boundaries
echo
echo "$SCRIPT_DIR/alternativeExpressionDatabase/createExonBoundaryDatabase.pl  --database=ALEXA_$DB  --server=$MYSQL_HOST  --user=viewer  --password=viewer  --target_length=$SIZE  --outdir=$BDIR  --logfile=$LOGDIR/createExonBoundaryDatabase_LOG.txt"
$SCRIPT_DIR/alternativeExpressionDatabase/createExonBoundaryDatabase.pl  --database=ALEXA_$DB  --server=$MYSQL_HOST  --user=viewer  --password=viewer  --target_length=$SIZE  --outdir=$BDIR  --logfile=$LOGDIR/createExonBoundaryDatabase_LOG.txt 1>/dev/null

#Then create jobs annotate them
echo
echo "Make annotation jobs and run"
ANNFILE="$BDIR/annotateExonBoundaries.sh"

export DB
export BDIR
export SIZE
export TEMPDIR
export ANNDIR
export MYSQL_HOST
export BASE_DIR
export SCRIPT_DIR

perl -ne '$DB=$ENV{"DB"}; $BDIR=$ENV{"BDIR"}; $SIZE=$ENV{"SIZE"}; $TEMPDIR=$ENV{"TEMPDIR"}; $ANNDIR=$ENV{"ANNDIR"}; $MYSQL_HOST=$ENV{"MYSQL_HOST"}; $BASE_DIR=$ENV{"BASE_DIR"}; $SCRIPT_DIR=$ENV{"SCRIPT_DIR"}; unless($_ =~ /Chromosome/){chomp($_); @line=split("\t", $_); print "$SCRIPT_DIR/alternativeExpressionDatabase/annotateExonBoundaries.pl  --database=ALEXA_$DB  --server=$MYSQL_HOST  --user=viewer  --password=viewer  --chr_filter=\"$line[0]:$line[1]:$line[2]-$line[3]\"  --target_size=$SIZE  --boundary_file=$BDIR/exonBoundaries\_$SIZE\mers.txt  --ucsc_align_dir=$BASE_DIR/$DB/mrna_est/partitions/   --genbank_mapfile=$BASE_DIR/$DB/mrna_est/GenBankToOrganism.btree  --wiggle=3  --outfile=$TEMPDIR/exonBoundaries_annotated_$line[0]_$line[1].txt  --logfile=$ANNDIR/annotateExonBoundaries_$line[0]_$line[1]_LOG.txt\n"}' $BASE_DIR/$DB/Regions_250_Genes.txt > $ANNFILE

#Then run these jobs
echo "bash $ANNFILE"
bash $ANNFILE

#Now join the annotation and fasta files together
echo
echo "Now join the annotation and fasta files together"
echo "cd $BDIR" 
cd $BDIR

echo "grep -h -v Boundary_ID $BDIR/temp/*.txt | sort -n > $BDIR/tmp.txt"
grep -h -v Boundary_ID $BDIR/temp/*.txt | sort -n > $BDIR/tmp.txt

B_FILE_FA=$(echo $BDIR/exonBoundaries_$SIZE\mers.fa)
B_FILE=$(echo $BDIR/exonBoundaries_$SIZE\mers.txt)
FINAL_FILE=$(echo $BDIR/exonBoundaries_$SIZE\mers_annotated.txt)
FINAL_FILE_GZ=$(echo $BDIR/exonBoundaries_$SIZE\mers_annotated.txt.gz)

echo "head -q -n 1 $BDIR/temp/*.txt | sort | uniq | cat - $BDIR/tmp.txt > $FINAL_FILE"
head -q -n 1 $BDIR/temp/*.txt | sort | uniq | cat - $BDIR/tmp.txt > $FINAL_FILE

echo "rm -f $BDIR/tmp.txt"
rm -f $BDIR/tmp.txt

echo
echo "Compress results file"
echo "gzip $FINAL_FILE"
gzip $FINAL_FILE
echo "gzip $B_FILE_FA"
gzip $B_FILE_FA

echo
echo "Clean-up temp files..."
echo "cd $BDIR"
cd $BDIR

echo "rm -f $B_FILE"
rm -f $B_FILE

echo "rm -f $BDIR/temp/*.txt"
rm -f $BDIR/temp/*.txt

echo
echo "Add the FID column"
echo "$SCRIPT_DIR/alternativeExpressionDatabase/appendFeatureId.pl  --infile=$FINAL_FILE_GZ  --outfile=$BDIR/tmp.txt  --prefix='EB'"
$SCRIPT_DIR/alternativeExpressionDatabase/appendFeatureId.pl  --infile=$FINAL_FILE_GZ  --outfile=$BDIR/tmp.txt  --prefix='EB'

