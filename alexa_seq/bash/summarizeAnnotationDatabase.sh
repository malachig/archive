#!/bin/bash

#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Summarize an annotation database by each feature type
#For Genes, Transcripts, Exon regions, Junctions (known & novel), Boundaries (known & novel), Introns, Silent intron regions, Active intron regions, Intergenics, Silent intergenic regions, Active intergenic regions
#Get the following:
#1.) Total events
#2.) Total events supported by EnsEMBL transcripts
#3.) Total events supported by ESTs/mRNAs
#4.) Total events supported by xESTs/mRNAs (i.e. conserved at expression level in at least one non-human species)
#5.) Total base count
#6.) Total and Percent of all bases that are masked
#7.) Total and Percent of all events with protein coding content

E_BADARGS=65
if [ ! -n "$1" ]
then
  echo "Usage: `basename $0` hs_55_37 55 62 /projects/malachig/sequence_databases (i.e. db, ensembl_version, junction_length, path)"
  exit $E_BADARGS
fi

#Commandline args
DB=$1
VERSION=$2
SEQ_SIZE=$3
DB_PATH=$4
OUTFILE=$(echo "$DB_PATH/$DB/Stats_$1.txt")

echo starting $DB
j_file=$(echo exonJunctions_"$SEQ_SIZE"\mers_annotated.txt.gz)
j_file_known="known_junctions.txt"
j_file_novel="novel_junctions.txt"
b_file=$(echo exonBoundaries_"$SEQ_SIZE"mers_annotated.txt.gz)
b_file_known="known_boundaries.txt"
b_file_novel="novel_boundaries.txt"

cd $DB_PATH/$DB/exonJunctions/
ENSG_SUPPORT_COL=$(zcat $j_file | head -n 1 | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "Supporting_EnsEMBL_Count"){print "$count";}}')
export ENSG_SUPPORT_COL
zcat $j_file | perl -ne '@line=split("\t", $_); if ($line[$ENV{"ENSG_SUPPORT_COL"}-1] > 0 || $line[$ENV{"ENSG_SUPPORT_COL"}-1] eq "Supporting_EnsEMBL_Count"){print "$_"}' > $j_file_known
zcat $j_file | perl -ne '@line=split("\t", $_); if ($line[$ENV{"ENSG_SUPPORT_COL"}-1] == 0 || $line[$ENV{"ENSG_SUPPORT_COL"}-1] eq "Supporting_EnsEMBL_Count"){print "$_"}' > $j_file_novel
gzip $j_file_known
gzip $j_file_novel

cd $DB_PATH/$DB/exonBoundaries/
ENSG_SUPPORT_COL=$(zcat $b_file | head -n 1 | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "Supporting_EnsEMBL_Count"){print "$count";}}')
export ENSG_SUPPORT_COL
zcat $b_file | perl -ne '@line=split("\t", $_); if ($line[$ENV{"ENSG_SUPPORT_COL"}-1] > 0 || $line[$ENV{"ENSG_SUPPORT_COL"}-1] eq "Supporting_EnsEMBL_Count"){print "$_"}' > $b_file_known
zcat $b_file | perl -ne '@line=split("\t", $_); if ($line[$ENV{"ENSG_SUPPORT_COL"}-1] == 0 || $line[$ENV{"ENSG_SUPPORT_COL"}-1] eq "Supporting_EnsEMBL_Count"){print "$_"}' > $b_file_novel
gzip $b_file_known
gzip $b_file_novel
cd $DB_PATH/$DB/

rm -f $OUTFILE
echo $OUTFILE

FILES="
genes/genes_annotated.txt.gz
transcripts/transcripts_annotated.txt.gz
exonRegions/exonRegions_annotated.txt.gz
exonJunctions/$j_file
exonJunctions/known_junctions.txt.gz
exonJunctions/novel_junctions.txt.gz
exonBoundaries/$b_file
exonBoundaries/known_boundaries.txt.gz
exonBoundaries/novel_boundaries.txt.gz
introns/introns_annotated.txt.gz
introns/silentIntronRegions.txt.gz
introns/activeIntronRegions.txt.gz
intergenics/intergenics_annotated.txt.gz
intergenics/silentIntergenicRegions.txt.gz
intergenics/activeIntergenicRegions.txt.gz
"

echo "Feature_Type	Feature_Count	EnsEMBL_Supported_Features	mRNA-EST_Supported_Features	Conserved_Features	Base_Count	UnMasked_Base_Count	Coding_Base_Count" >> $OUTFILE
for f in $FILES
do
   #Get the columns of interest
  ENSG_SUPPORT_COL=$(zcat $f | head -n 1 | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "Supporting_EnsEMBL_Count"){print "$count";}}')
  MRNA_SUPPORT_COL=$(zcat $f | head -n 1 | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "Supporting_mRNA_Count"){print "$count";}}')
  EST_SUPPORT_COL=$(zcat $f | head -n 1 | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "Supporting_EST_Count"){print "$count";}}')
  CONSERVED_COL=$(zcat $f | head -n 1 | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "Conserved_Species_Count"){print "$count";}}')
  BASE_COUNT_COL=$(zcat $f | head -n 1 | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "Base_Count"){print "$count";}}')
  UNMASKED_BASE_COUNT_COL=$(zcat $f | head -n 1 | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "UnMasked_Base_Count"){print "$count";}}')
  CODING_BASE_COUNT_COL=$(zcat $f | head -n 1 | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "Coding_Base_Count"){print "$count";}}')
  #Get the total number of features
  TOTAL_FEATURES=$(zcat $f | grep -v Base_Count | wc -l)
  #Get the number of ensembl supported sequences ('NA' for feature types that do not have this column)
  ENSG_SUPPORT="NA"
  if [ -n "$ENSG_SUPPORT_COL" ]
    then
      ENSG_SUPPORT=$(zcat $f | cut -f $ENSG_SUPPORT_COL | perl -ne 'if ($_ =~ /(\d+)/){if ($1 > 0){$c++;}}; if (eof){if ($c > 0){print "$c"}else{print "0"}}')
  fi

  #Get the number of EST or mRNA supported sequences ('NA' for feature types that do not have these columns)
  SEQ_SUPPORT="NA"
  if [ -n "$MRNA_SUPPORT_COL" ]
    then
      SEQ_SUPPORT=$(zcat $f | cut -f $MRNA_SUPPORT_COL,$EST_SUPPORT_COL | perl -ne 'if ($_ =~ /(\d+)\s+(\d+)/){if ($1 > 0 | $2 > 0){$c++;}}; if (eof){if ($c > 0){print "$c"}else{print "0"}}')
  fi
  #Get the number of conserved sequence ('NA' for feature types that do not have this column)
  CONSERVED="NA"
  if [ -n "$CONSERVED_COL" ]
    then
      CONSERVED=$(zcat $f | cut -f $CONSERVED_COL | perl -ne 'if ($_ =~ /(\d+)/){if ($1 > 0){$c++;}}; if (eof){if ($c > 0){print "$c"}else{print "0"}}')
  fi
  #Get the grand base count
  BASE_COUNT="NA"
  if [ -n "$BASE_COUNT_COL" ]
    then
      BASE_COUNT=$(zcat $f | cut -f $BASE_COUNT_COL | perl -ne 'if ($_ =~ /(\d+)/){if ($1 >= 0){$c+=$1;}}; if (eof){if ($c > 0){print "$c"}else{print "0"}}')
  fi
  #Get the grand unmasked base count
  UNMASKED_BASE_COUNT="NA"
  if [ -n "$UNMASKED_BASE_COUNT_COL" ]
    then
      UNMASKED_BASE_COUNT=$(zcat $f | cut -f $UNMASKED_BASE_COUNT_COL | perl -ne 'if ($_ =~ /(\d+)/){if ($1 >= 0){$c+=$1;}}; if (eof){if ($c > 0){print "$c"}else{print "0"}}')
  fi
  #Get the grand protein coding base count
  CODING_BASE_COUNT="NA"
  if [ -n "$CODING_BASE_COUNT_COL" ]
    then
      CODING_BASE_COUNT=$(zcat $f | cut -f $CODING_BASE_COUNT_COL | perl -ne 'if ($_ =~ /(\d+)/){if ($1 >= 0){$c+=$1;}}; if (eof){if ($c > 0){print "$c"}else{print "0"}}')
  fi
   #print out the results
  echo "$f	$TOTAL_FEATURES	$ENSG_SUPPORT	$SEQ_SUPPORT	$CONSERVED	$BASE_COUNT	$UNMASKED_BASE_COUNT	$CODING_BASE_COUNT" >> $OUTFILE
done

#Clean up the temp files
rm -f exonJunctions/known_junctions.txt.gz
rm -f exonJunctions/novel_junctions.txt.gz
rm -f exonBoundaries/known_boundaries.txt.gz
rm -f exonBoundaries/novel_boundaries.txt.gz

perl -ne '$_ =~ s/\s+/\t/g; print "$_\n"' $OUTFILE > tmp.txt
mv -f tmp.txt $OUTFILE
echo completed $DB
echo
















