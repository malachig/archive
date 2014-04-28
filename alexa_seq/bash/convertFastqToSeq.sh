#!/bin/bash

E_BADARGS=65
if [ ! -n "$1" ]
then
  echo "Usage: `basename $0` /projects/alexa/public_wtss_data/illumina_hiseq/50bp_PE_mRNA_Seq HR0007 5 [i.e. basedir library lane]"
  exit $E_BADARGS
fi

#Commandline args
DIR=$1
LIB=$2
LANE=$3

R1_FILE_IN=$(echo "$DIR/$LIB/s_""$LANE""_1_sequence.txt.gz")
R2_FILE_IN=$(echo "$DIR/$LIB/s_""$LANE""_2_sequence.txt.gz")

R1_FILE_OUT=$(echo "$DIR/$LIB/s_""$LANE""_1_tmp.txt")
R2_FILE_OUT=$(echo "$DIR/$LIB/s_""$LANE""_2_tmp.txt")

#Parse read1 file and read name columns and first read seq
echo
echo START
echo Args: $DIR $LIB $LANE
echo Read1 file in: $R1_FILE_IN
echo Read1 file out: $R1_FILE_OUT
zcat $R1_FILE_IN | perl -ne 'chomp($_); $c++; if ($c == 1){if ($_ =~ /\:(\d+)\:(\d+)\:(\d+)\:(\d+)\#/){print "$1\t$2\t$3\t$4\t";}} if ($c == 2){print "$_\n";} if ($c == 4){$c = 0;}' > $R1_FILE_OUT

#Parse read2 file and get second read seq
echo Read2 file in: $R2_FILE_IN
echo Read2 file out: $R2_FILE_OUT
zcat $R2_FILE_IN | perl -ne 'chomp($_); $c++; if ($c == 2){print "$_\n";} if ($c == 4){$c = 0;}' > $R2_FILE_OUT

#Paste read1 and read2 temp files
SEQ_FILE_OUT=$(echo "$DIR/$LIB/s_""$LANE""_0001_seq.txt")
echo Join Read1 and Read2 files: $SEQ_FILE_OUT 
paste -d "" $R1_FILE_OUT $R2_FILE_OUT > $SEQ_FILE_OUT
echo bzip $SEQ_FILE_OUT
bzip2 -f $SEQ_FILE_OUT

#Delete the tmp files
rm -f $R1_FILE_OUT
rm -f $R2_FILE_OUT
echo DONE


