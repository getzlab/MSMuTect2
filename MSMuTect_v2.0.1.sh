#!/bin/bash

# define bams
TBAM="~/git_repos/MSMuTect2/104R16869.bam"
TBAI="~/git_repos/MSMuTect2/104R16869.bai"

NBAM="~/git_repos/MSMuTect2/104R16905.bam"
NBAI="~/git_repos/MSMuTect2/104R16905.bai"


PAIRID="test"

######
# Input
# $1 Pair_ID
# $2 Tumor  BAM. Full PATH
# $3 Normal BAM. Full PATH
# $4 Loci file list

#Calling Tumor alleles
python3 bin/msmutect.py -I $2 -l $4 -O $1.Tumor.hist
sh bin/hist2py.sh $1.Tumor.hist
python bin/get_all.py  $1.Tumor.hist.mot bin/P_A.csv > $1.Tumor.hist.mot.all



#Calling Normal alleles
python3 bin/msmutect.py -I $3 -l $4 -O $1.Normal.hist
sh bin/hist2py.sh $1.Normal.hist
python bin/get_all.py  $1.Normal.hist.mot bin/P_A.csv > $1.Normal.hist.mot.all

#Taking the shared loci
sh bin/Shared_loci_v3.sh $1.Tumor.hist.mot.all  $1.Normal.hist.mot.all A


#Calling mutations
python bin/Find_mutations2.py  $1.Tumor.hist.mot.all.tmp.par.reg  $1.Normal.hist.mot.all.tmp.par.reg bin/P_A.csv 8 0.3 0.031    > $1.mut

#Simple output format 

tr '\n' ' ' < $1.mut | awk '{gsub("@","\n");print $0}' > $1.mut.cln
        

awk 'BEGIN{print ("Locus\tDecision\tNornal_histogram\tNormal_alleles\tNormal_frequencies\tTumor_histogram\tTumor_alleles\tTumor_frequencies\tmotif\tmotif_size\tref_length\tNorm_num_alleles\tTum_num_alleles\tNorm_allele1\tNorm_allele2\tTum_allele1\tTum_allele2\tTum_allele3\tTum_allele4")}{gsub("]","");gsub("@","") split($0,a,"[");split($2,b,":");split(b[4],c,"");split(a[3],d," ");split(a[6],e," ");printf $2"\t"$1"\t"a[2]"\t"a[3]"\t"a[4]"\t"a[5]"\t"a[6]"\t"a[7]"\t"b[4]"\t"length(c)"\t"b[5]"\t"length(d)"\t"length(e)"\t"; if(length(d)==1){printf d[1]"\t-9\t"}else {printf d[1]"\t"d[2]"\t"}; if(length(e)==1){printf e[1]"\t-9\t-9\t-9"}if(length(e)==2){printf e[1]"\t"e[2]"\t-9\t-9"} if(length(e)==3){printf e[1]"\t"e[2]"\t"e[3]"\t-9"} if(length(e)==4){printf e[1]"\t"e[2]"\t"e[3]"\t"e[4]} ;printf "\n"}' $1.mut.cln | grep -v All | awk 'BEGIN{FS="\t"}{if (NF==19){print $0}}' > $1.mut.maf_like

sh bin/add_altered_base.sh $1.mut.maf_like $4

awk 'BEGIN{FS="\t"}{n=n+1;if(n==1){print $0};if($2==1){print $0}}' $1.mut.maf_like > $1.mut.maf_like.dec 

