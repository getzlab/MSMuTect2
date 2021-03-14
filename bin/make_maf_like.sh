#####
# Making MSMuTect output a maf_like_file

# Inpuet
# $1 = *.cln file.


#awk 'BEGIN{print ("Locus\tDecision\tNornal_histogram\tNormal_alleles\tNormal_frequencies\tTumor_histogram\tTumor_alleles\tTumor_frequencies")}{gsub("]",""); split($0,a,"[");print $2"\t"$1"\t"a[2]"\t"a[3]"\t"a[4]"\t"a[5]"\t"a[6]"\t"a[7] }' $1 | awk '{split($3,a," ");for(i=1;i<length(a)+1;i=i+2){his=his a[i]"-"a[i+1]","};$3=his;print $1"\t"$2"\t"$3"\t"$4"\t"$5"\"$6"\t"$7}' > $1.maf_like
#awk 'BEGIN{print ("Locus\tDecision\tNornal_histogram\tNormal_alleles\tNormal_frequencies\tTumor_histogram\tTumor_alleles\tTumor_frequencies")}{gsub("]",""); split($0,a,"[");print $2"\t"$1"\t"a[2]"\t"a[3]"\t"a[4]"\t"a[5]"\t"a[6]"\t"a[7] }' $1 | awk '{split($3,a," ");printf $1"\t"$2"\t";for(i=1;i<length(a)+1;i=i+2){printf a[i]"-"a[i+1]","};printf "\t"$4"\t"$5"\t"$6"\t"$7"\n}' > $1.maf_like
awk 'BEGIN{print ("Locus\tDecision\tNornal_histogram\tNormal_alleles\tNormal_frequencies\tTumor_histogram\tTumor_alleles\tTumor_frequencies")}{gsub("]",""); split($0,a,"[");print $2"\t"$1"\t"a[2]"\t"a[3]"\t"a[4]"\t"a[5]"\t"a[6]"\t"a[7] }' $1 | awk 'BEGIN{FS="\t"}{n=n+1;if(n==1){print $0}else{split($3,a," ");split($6,b," ");printf $1"\t"$2"\t";for(i=1;i<length(a)-1;i=i+2){printf a[i]"-"a[i+1]","};printf a[length(a)-1]"-"a[length(a)]"\t"$4"\t"$5"\t"; for(i=1;i<length(b)-1;i=i+2){printf b[i]"-"b[i+1]","}; printf b[length(b)-1]"-"b[length(b)]"\t"$7"\t"$8"\n"}}' | awk 'BEGIN{while((getline a< "bin/Loci_list_gene_type.txt")>0){split(a,b," ");data[b[1]]=b[2]"\t"b[3]}}{split($1,c,":");locus=c[1]":"c[2]":"c[3];n=n+1;if(n==1){print $0"\tGene\tType"}else{print $0"\t"data[locus] } }'> $1.maf_like
