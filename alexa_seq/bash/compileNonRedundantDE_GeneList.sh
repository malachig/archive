#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Generate a list of all genes that are significantly DE in any feature category (Gene, transcript, junction, etc.)
#Store the largest magnitude DE value from any feature in that gene

#Change to the root directory for DE values
cd /projects/malachig/solexa/figures_and_stats/DE/

#Define all the input files of interest
FILES="
ENST_v53/Mip101_vs_Mip5FuR_Gene_DE_Values_Significant_MTC.txt
Transcripts_v53/Mip101_vs_Mip5FuR_Transcript_DE_Values_Significant_MTC.txt
ENST_v53/Mip101_vs_Mip5FuR_ExonRegion_DE_Values_Significant_MTC.txt
Junctions_v53/Mip101_vs_Mip5FuR_KnownJunction_DE_Values_Significant_MTC.txt
Junctions_v53/Mip101_vs_Mip5FuR_NovelJunction_DE_Values_Significant_MTC.txt
Boundaries_v53/Mip101_vs_Mip5FuR_KnownBoundary_DE_Values_Significant_MTC.txt
Boundaries_v53/Mip101_vs_Mip5FuR_NovelBoundary_DE_Values_Significant_MTC.txt
Introns_v53/Mip101_vs_Mip5FuR_ActiveIntronRegion_DE_Values_Significant_MTC.txt
Introns_v53/Mip101_vs_Mip5FuR_SilentIntronRegion_DE_Values_Significant_MTC.txt
"

#Define a tmp file and results file
TMP1_OUT="temp1.txt"
TMP2_OUT="temp2.txt"
HEADER_OUT="header.txt"
RESULT_FILE="DE_List_AllFeatureTypes_HighestDEPerGene.txt"

rm -f $TMP1_OUT

#Get columns of interest from all files and join into a single file: stored in $TMP_OUT
for FILE in $FILES
do

  #Get the columns of interest
  NAME_COL=$(head -n 1 $FILE | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "Gene_Name"){print "$count";}}')
  ENSG_COL=$(head -n 1 $FILE | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "EnsEMBL_Gene_ID"){print "$count";}}')
  DE_COL=$(head -n 1 $FILE | perl -ne 'chomp($_); @line = split("\t", $_); foreach my $col (@line){$count++; if ($col eq "Fold_Change"){print "$count";}}')
  
  #Get the current data type
  DATA_TYPE=$(echo $FILE | perl -ne 'if ($_ =~ /Mip5FuR\_(\w+)\_DE\_/){print "$1"}')

  echo $NAME_COL $ENSG_COL $DE_COL $FILE $DATA_TYPE

  export DATA_TYPE
  cut -f $DE_COL,$NAME_COL,$ENSG_COL $FILE | perl -ne 'chomp($_); if ($_ =~ /Gene_Name/){next();} print "$_\t"; print $ENV{"DATA_TYPE"}; print "\n";' >> $TMP1_OUT
done


perl -ne 'chomp($_); if ($_ =~ /^ENSG/){@line = split("\t", $_); $ensg=$line[0]; $name=$line[1]; $fc=$line[2]; $type=$line[3]; $fca=abs($fc); if ($fca > $data{$name}{fca}){$data{$name}{fca}=$fca; $data{$name}{fc}=$fc; $data{$name}{ensg}=$ensg; $data{$name}{type}=$type; $data{$name}{count}++;}else{$data{$name}{count}++;}} if(eof){foreach my $name (keys %data){print "$name\t$data{$name}{ensg}\t$data{$name}{type}\t$data{$name}{count}\t$data{$name}{fc}\t$data{$name}{fca}\n"}}' $TMP1_OUT | sort -nr -k 6 | perl -ne '$r++; print "$r\t$_";' > $TMP2_OUT

echo "Rank	Gene_Name	EnsEMBL_Gene_ID	Top_Feature_Type	DE_Feature_Count	Fold_Change	Fold_Change_ABS" | cat > $HEADER_OUT
cat $HEADER_OUT $TMP2_OUT > $RESULT_FILE

#Clear the temp file
rm -f $TMP1_OUT
rm -f $TMP2_OUT
rm -f $HEADER_OUT


