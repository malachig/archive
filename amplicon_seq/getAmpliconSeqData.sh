
#Store plate to amplicon info in
/home/malachig/svn/amplicon_seq_analysis/plate_to_amplicon.txt

#Create directories - OTB_BCCA
perl -ne 'chomp($_); @line=split(" ", $_); if($_ =~ /^\#/){next();} unless($_ =~ /\w+/){next();} print "mkdir /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/$line[2]/\n";' /home/malachig/svn/amplicon_seq_analysis/plate_to_amplicon.txt

#Create directories - StPauls
perl -ne 'chomp($_); @line=split(" ", $_); if($_ =~ /^\#/){next();} unless($_ =~ /\w+/){next();} print "mkdir /home/samba/public/jpaul/UMPS_Traces/StPauls/$line[2]/\n";' /home/malachig/svn/amplicon_seq_analysis/plate_to_amplicon.txt


#Copy Forward and Reverse reads - OTB_BCCA
perl -ne 'chomp($_); @line=split(" ", $_); if($_ =~ /^\#/){next();} unless($_ =~ /\w+/){next();} print "cp /home/aldente/private/Projects/Marco_Research/$line[0]/AnalyzedData/$line[0]$line[1].CR/chromat_dir/*ab1 /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/$line[2]/\n";' /home/malachig/svn/amplicon_seq_analysis/plate_to_amplicon.txt
perl -ne 'chomp($_); @line=split(" ", $_); if($_ =~ /^\#/){next();} unless($_ =~ /\w+/){next();} print "cp /home/aldente/private/Projects/Marco_Research/$line[0]/AnalyzedData/$line[0]$line[1].C21/chromat_dir/*ab1 /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/$line[2]/\n";' /home/malachig/svn/amplicon_seq_analysis/plate_to_amplicon.txt

#Copy Forward and Reverse reads - St. Pauls
perl -ne 'chomp($_); @line=split(" ", $_); if($_ =~ /^\#/){next();} unless($_ =~ /\w+/){next();} @rows=split(",", $line[3]); foreach my $row(@rows){print "cp /home/aldente/private/Projects/Marco_Research/$line[0]/AnalyzedData/$line[0]$line[1].CR/chromat_dir/*_$row*ab1 /home/samba/public/jpaul/UMPS_Traces/StPauls/$line[2]/\n";}' /home/malachig/svn/amplicon_seq_analysis/plate_to_amplicon.txt
perl -ne 'chomp($_); @line=split(" ", $_); if($_ =~ /^\#/){next();} unless($_ =~ /\w+/){next();} @rows=split(",", $line[3]); foreach my $row(@rows){print "cp /home/aldente/private/Projects/Marco_Research/$line[0]/AnalyzedData/$line[0]$line[1].C21/chromat_dir/*_$row*ab1 /home/samba/public/jpaul/UMPS_Traces/StPauls/$line[2]/\n";}' /home/malachig/svn/amplicon_seq_analysis/plate_to_amplicon.txt

#Check file counts - St. Pauls
ls */* | perl -ne 'if ($_ =~ /(F\d+\_R\d+)/){print "\n$1";}' | sort | uniq -c


#Rename reads - OTB_BCCA
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F41_R41/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F41_R41/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F41_R41/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F41_R41/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F42_R42/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F42_R42/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F42_R42/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F42_R42/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F17_R17/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F17_R17/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F17_R17/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F17_R17/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F16_R15/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F16_R15/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F16_R15/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F16_R15/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F38_R38/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F38_R38/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F38_R38/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F38_R38/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F39_R39/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F39_R39/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F39_R39/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F39_R39/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F43_R43/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F43_R43/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F43_R43/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F43_R43/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F44_R44/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F44_R44/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F44_R44/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F44_R44/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F45_R45/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F45_R45/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F45_R45/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F45_R45/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F46_R46/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F46_R46/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F46_R46/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F46_R46/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F47_R47/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F47_R47/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F47_R47/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F47_R47/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F48_R48/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F48_R48/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F48_R48/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F48_R48/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F50_R50/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F50_R50/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F50_R50/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F50_R50/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F51_R51/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F51_R51/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F51_R51/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F51_R51/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F52_R52/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F52_R52/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F52_R52/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F52_R52/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F53_R53/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F53_R53/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F53_R53/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F53_R53/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F17_R17b/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F17_R17b/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F17_R17b/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F17_R17b/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F16_R15b/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F16_R15b/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F16_R15b/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F16_R15b/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F49_R49/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F49_R49/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F49_R49/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F49_R49/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F54_R54/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F54_R54/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F54_R54/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F54_R54/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F55_R55/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F55_R55/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F55_R55/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F55_R55/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F56_R56/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F56_R56/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F56_R56/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F56_R56/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F57_R57/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F57_R57/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F57_R57/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F57_R57/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F59_R59/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F59_R59/tmp.sh
ll -l /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F59_R59/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/OTB_BCCA/F59_R59/tmp.sh


#Rename reads - OTB_BCCA
ll -l /home/samba/public/jpaul/UMPS_Traces/StPauls/*/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.CR)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_R$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/StPauls/tmp.sh

ll -l /home/samba/public/jpaul/UMPS_Traces/StPauls/*/*.ab1 | perl -ne 'chomp($_); if ($_ =~ /(\/home.*)(\/MG.*)(\.C21)(.*ab1)/){print "\nmv $1$2$3$4 $1$2_F$4"}' > tmp.sh; bash tmp.sh; rm -f /home/samba/public/jpaul/UMPS_Traces/StPauls/tmp.sh

#Check file counts - St. Pauls
ls */* | perl -ne 'if ($_ =~ /(F\d+\_R\d+)/){print "\n$1";}' | sort | uniq -c











