# Adding the exact base change 
#
# Input
# $1 *maf_like file. Use a full PATH do not use the ~ symbol 
# $2 Refernece file in a phobos format. Must use the same reference file that was used to call the MS-indels


head -n 1 $1 | awk '{print $0"\tNormal_bases_allele1\tNormal_bases_allele2\tTumor_bases_allele1\tTumor_bases_allele2\tTumor_bases_allele3\tTumor_bases_allele4"}' > $1.base_change

awk -v mut_file=$1 'BEGIN{while((getline a < mut_file)>0){gsub("chr","",a);split(a,b,"\t");split(b[1],q,":");mut[q[1]":"q[2]":"q[3]]=a};}{locus=$1":"$4":"$5;if(mut[locus]!=""){split($13,b,"");split($14,c,"");motif="";for(k=1;k<length(b)+1;k=k+1){motif=motif c[k]};split(mut[locus],d,"\t");printf mut[locus]"\t"; for(k=1;k<d[14]+1;k=k+1){printf motif};printf "\t";for(k=1;k<d[15]+1;k=k+1){printf motif}; printf "\t";for(k=1;k<d[16]+1;k=k+1){printf motif}; printf "\t";for(k=1;k<d[17]+1;k=k+1){printf motif}; printf "\t";for(k=1;k<d[18]+1;k=k+1){printf motif}; printf "\t";for(k=1;k<d[19]+1;k=k+1){printf motif};  printf "\n"}}' $2 >> $1.base_change
