#!/bin/bash
#isub -i 'docker(zlskidmore/samtools:1.10)'
#filterLargeIntronAlignments.sh /storage1/fs1/mgriffit/Active/scrna_mcb6c/Mouse_Bladder_MCB6C_Arora/scRNA/Rep3/v5/regtools/ Rep3_Control_tumor.bam /scratch1/fs1/mgriffit/filtered_scrna_bams/bams /storage1/fs1/mgriffit/Active/temp

indir=$1
filename=$2
scratchdir=$3
regions=$4
finaldir=$5

echo $indir
echo $filename
echo $scratchdir
echo $regions
echo $finaldir

#create BAM covering only the regions of interest
cat $regions | xargs samtools view -b $indir/$filename > $scratchdir/$filename.subset.bam

#sort the BAM file
samtools sort -O bam -o $scratchdir/$filename.subset.sort.bam $scratchdir/$filename.subset.bam 
mv $scratchdir/$filename.subset.sort.bam $scratchdir/$filename.subset.bam

#create BAM with intron alignments >50000 removed
samtools view -h $scratchdir/$filename.subset | perl -ne 'chomp; @l=split("\t",$_); if($l[5] =~ /(\d+)N/){if($1<50000){print "$_\n"}}else{print "$_\n"}' 2>/dev/null | samtools view -bS - > $scratchdir/$filename

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
rm -f $scratchdir/$filename.subset
rm -f $scratchdir/$filename
rm -f $scratchdir/$filename.bai

