
####
# input $1
# $1= file1 .e.g REBC-ACAI-NB1-A-1-0-D-A49V-36/A_impure_5_100_fa.aln_NoDup.hist.py.all
# $2= file2 e.g. REBC-ACAI-NT1-A-1-1-D-A49V-36/A_impure_5_100_fa.aln_NoDup.hist.py.all
# $3 motif. e.g. C
#awk -v motif=$3 '{data[$1" "FILENAME]=$0;count[$1]=count[$1]+1;files[FILENAME]=1}END{for(i in count){if(count[i]==2){for(j in files){n=n+1;qw[n]=j;out[n]=j".reg2"}; split(out[1],a,"/"); split(out[2],b,"/"); split(b[2],c,".reg");d=motif"_impure_5_100_fa.aln_NoDup.hist.py.all.bed.par.reg";file[1]= a[1]"/"a[1]"_"b[1]"_"d;file[2]=b[1]"/"b[1]"_"a[1]"_"d; print data[i" "qw[1]]> file[1]; print data[i" "qw[2]]>file[2]}}}' $1 $2



 awk '{gsub(" -999 "," -99999999 "); print $0}' $1 | awk '{split($0,a,"-99999999");split(a[1],b,":");split(b[5],c," ");printf $1" ";for(i=2;i<length(c);i=i+2){su[c[i]]=su[c[i]]+c[i+1]};for(i in su){printf i" "su[i]" "};printf "-99999999"a[2]"-99999999"a[3]"-99999999"a[4];printf "\n";delete su}' > $1.tmp


 awk '{gsub(" -999 "," 99999999 "); print $0}' $2 | awk '{split($0,a,"99999999");split(a[1],b,":");split(b[5],c," ");printf $1" ";for(i=2;i<length(c);i=i+2){su[c[i]]=su[c[i]]+c[i+1]};for(i in su){printf i" "su[i]" "};printf "-99999999"a[2]"-99999999"a[3]"-99999999"a[4];printf "\n";delete su}'   > $2.tmp


awk -v motif=$3 '{data[$1" "FILENAME]=$0;count[$1]=count[$1]+1;files[FILENAME]=1}END{for (i in count){if(count[i]==2){for(j in files){file=j".par.reg";print data[i" "j]> file}}}}' $1.tmp $2.tmp

awk '{gsub(/\[/,"");gsub(/\]/,"");print $0}' $1.tmp.par.reg > $1.tmp.par.reg.fx
cp $1.tmp.par.reg.fx $1.tmp.par.reg

awk '{gsub(/\[/,"");gsub(/\]/,"");print $0}' $2.tmp.par.reg > $2.tmp.par.reg.fx
cp $2.tmp.par.reg.fx $2.tmp.par.reg

