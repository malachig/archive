#!/bin/bash
#isub -i 'docker(zlskidmore/samtools:1.10)'
#~/git/archive/misc/filterLargeIntronAlignments.sh /storage1/fs1/mgriffit/Active/scrna_mcb6c/Mouse_Bladder_MCB6C_Arora/scRNA/Rep1/v5/M_HC-JKB216-ICB_MCB6C_Rep_1-lib1/outs/ possorted_genome_bam.bam /scratch1/fs1/mgriffit/tmp/process/ /scratch1/fs1/mgriffit/tmp/regions.txt /storage1/fs1/mgriffit/Active/tmp/ 5000

#regions format perline: chr16:33952785-33969003

indir=$1
filename=$2
scratchdir=$3
regions=$4
finaldir=$5
maxsize=$6

echo Input dir: $indir
echo BAM file name: $filename
echo Working dir: $scratchdir
echo Regions list file: $regions
echo Output dir: $finaldir
echo Max intron size: $maxsize

export MAXSIZE=$maxsize

#create BAM covering only the regions of interest
cat $regions | xargs samtools view -b $indir/$filename > $scratchdir/$filename.subset.bam

#sort the BAM file
samtools sort -O bam -o $scratchdir/$filename.subset.sort.bam $scratchdir/$filename.subset.bam 
mv $scratchdir/$filename.subset.sort.bam $scratchdir/$filename.subset.bam

#create BAM with intron alignments greater than specified size removed
samtools view -h $scratchdir/$filename.subset.bam | perl -ne 'chomp; $maxlength=0; @l=split("\t",$_); if($l[5] =~ /(\d+)N/){@ops=split(/\d+/,$l[5]); shift @ops; @lengths=split(/\D+/, $l[5]); foreach $op(@ops){$length=shift(@lengths); if($op eq "N" && $length>$maxlength){$maxlength=$length}}; if($maxlength<$ENV{MAXSIZE}){print "$_\n"}}else{print "$_\n"}' 2>/dev/null | samtools view -bS - > $scratchdir/$filename

#add header back to BAM
samtools view -H $indir/$filename > $scratchdir/$filename.header
samtools reheader $scratchdir/$filename.header $scratchdir/$filename > $scratchdir/$filename.tmp
mv $scratchdir/$filename.tmp $scratchdir/$filename

#index the new BAM
samtools index $scratchdir/$filename

#copy files to final dir
cp $scratchdir/$filename $finaldir/$filename
cp $scratchdir/$filename.bai $finaldir/$filename.bai

#clean up temp files
rm -f $scratchdir/$filename.header
rm -f $scratchdir/$filename.subset.bam
rm -f $scratchdir/$filename
rm -f $scratchdir/$filename.bai

