#sub -i 'docker(zlskidmore/samtools:1.10)'
#filterLargeIntronAlignments.sh /storage1/fs1/mgriffit/Active/scrna_mcb6c/Mouse_Bladder_MCB6C_Arora/scRNA/Rep3/v5/regtools/ Rep3_Control_tumor.bam /scratch1/fs1/mgriffit/filtered_scrna_bams/

indir=$1
filename=$2
outdir=$2

#create BAM with intron alignments >50000 removed
samtools view -h $indir/$filename | perl -ne 'chomp; @l=split("\t",$_); if($l[5] =~ /(\d+)N/){if($1<50000){print "$_\n"}}else{print "$_\n"}' 2>/dev/null | samtools view -bS - > $outdir/$filename

#add header back to BAM
samtools view -H $indir/$filename > $outdir/$filename.header
samtools reheader $outdir/$filename.header $outdir/$filename

#index the new BAM
samtools index $outdir/$filename

#clean up temp files
rm -f $outdir/$filename.header

