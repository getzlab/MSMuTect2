

awk 'BEGIN{FS="\t"}{print $1}' $1 | awk '{split($1,a,":");if(a[4]!="yosi"){gsub("_",", ");a[5]=int(a[5]);$1=a[1]":"a[2]":"a[3]":"a[4]":"a[5]",";if(NF>1&&a[5]<30){print $0}}}'  | awk 'BEGIN{FS=","}{for(i=2;i<NF+1;i=i+1){$i=int($i)}; for(i=1;i<NF;i=i+1){printf $i", "};printf $NF"\n"}'  > $1.mot
