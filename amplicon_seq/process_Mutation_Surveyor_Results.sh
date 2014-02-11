#Create commands to summarize SNPS from mutation surveyor reports
cd /home/samba/public/jpaul/Mutation_Surveyor/
ll /home/samba/public/jpaul/Mutation_Surveyor/*.txt | perl -ne 'chomp($_); if ($_ =~ /(MG\w+\.txt)/){print "/home/rmorin/genbank/snp_check.pl /home/samba/public/jpaul/UMPS_Traces/UMPS_GenomicGenBank_plus1KbFlank_MutSurv3.gbk $1 $1 > /home/samba/public/jpaul/Mutation_Surveyor/Summary/$1.out\n"}' > tmp.sh

bash tmp.sh

#Join the output into one big file and sort
cd /home/samba/public/jpaul/Mutation_Surveyor/Summary

#Move OTB and StPauls results into separate folders

#OTB_BCCA
cd /home/samba/public/jpaul/Mutation_Surveyor/Summary/OTB_BCCA/
cat *.out | sort > Combined.out

#Identify those position in this file that correspond to known SNPs
cut -f 1 Combined.out | sort | join - /projects/rmorin/SNPs/all_known_snps.txt.current > Observed_KnownSNPs.txt

#Strip out the unique position that correspond to observed known SNPS and use these to produce a list of only the NOVEL records found
cut -f 1 -d " " Observed_KnownSNPs.txt | sort | uniq | grep -v -f - Combined.out > Combined_NOVEL_OTB_BCCA.out
cut -f 1 -d " " Observed_KnownSNPs.txt | sort | uniq | grep -f - Combined.out > Combined_KNOWN_OTB_BCCA.out

#StPauls
cd /home/samba/public/jpaul/Mutation_Surveyor/Summary/StPauls
cat *.out | sort > Combined.out
cut -f 1 Combined.out | sort | join - /projects/rmorin/SNPs/all_known_snps.txt.current > Observed_KnownSNPs.txt
cut -f 1 -d " " Observed_KnownSNPs.txt | sort | uniq | grep -v -f - Combined.out > Combined_NOVEL_StPauls.out
cut -f 1 -d " " Observed_KnownSNPs.txt | sort | uniq | grep -f - Combined.out > Combined_KNOWN_StPauls.out






