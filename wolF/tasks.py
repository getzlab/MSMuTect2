import wolf

def msmutect_sample(bam, bai, loci, sample_type, output_prefix,):
    return wolf.Task(
        name = 'msmutect_' + sample_type,
        inputs= locals(),
        outputs = {'hist': '*.hist.mot.all'},
        script = [
            """
python3 /app/msmutect.py -I ${bam} -l ${loci} -O ${output_prefix}.${sample_type}.hist
bash /app/hist2py.sh ${output_prefix}.${sample_type}.hist
python /app/get_all.py  ${output_prefix}.${sample_type}.hist.mot /app/src/P_A.csv > ${output_prefix}.${sample_type}.hist.mot.all
            """
        ],
        docker = "gcr.io/broad-getzlab-workflows/msmutect2_wolf:v2",
        resources={"mem": '8G', 'cpus-per-task': '4'}
    )

def postprecess_msindel(tumor_hist, normal_hist, output_prefix, loci_file):
    return wolf.Task(
        name = "msmutect_postprocess",
        inputs = locals(),
        script=[
            r"""
#Taking the shared loci
bash /app/Shared_loci_v3.sh ${tumor_hist}  ${normal_hist} A

#Calling mutations
python3 /app/Find_mutations2.py  ${output_prefix}.Tumor.hist.mot.all.tmp.par.reg  ${output_prefix}.Normal.hist.mot.all.tmp.par.reg /app/P_A.csv 8 0.3 0.031  > ${output_prefix}.mut

#Simple output format 
tr '\n' ' ' < ${output_prefix}.mut | awk '{gsub("@","\n");print $0}' > ${output_prefix}.mut.cln
awk 'BEGIN{print ("Locus\tDecision\tNornal_histogram\tNormal_alleles\tNormal_frequencies\tTumor_histogram\tTumor_alleles\tTumor_frequencies\tmotif\tmotif_size\tref_length\tNorm_num_alleles\tTum_num_alleles\tNorm_allele1\tNorm_allele2\tTum_allele1\tTum_allele2\tTum_allele3\tTum_allele4")}{gsub("]","");gsub("@","") split($0,a,"[");split($2,b,":");split(b[4],c,"");split(a[3],d," ");split(a[6],e," ");printf $2"\t"$1"\t"a[2]"\t"a[3]"\t"a[4]"\t"a[5]"\t"a[6]"\t"a[7]"\t"b[4]"\t"length(c)"\t"b[5]"\t"length(d)"\t"length(e)"\t"; if(length(d)==1){printf d[1]"\t-9\t"}else {printf d[1]"\t"d[2]"\t"}; if(length(e)==1){printf e[1]"\t-9\t-9\t-9"}if(length(e)==2){printf e[1]"\t"e[2]"\t-9\t-9"} if(length(e)==3){printf e[1]"\t"e[2]"\t"e[3]"\t-9"} if(length(e)==4){printf e[1]"\t"e[2]"\t"e[3]"\t"e[4]} ;printf "\n"}' $1.mut.cln | grep -v All | awk 'BEGIN{FS="\t"}{if (NF==19){print $0}}' > ${output_prefix}.mut.maf_like
bash /app/add_altered_base.sh ${output_prefix}.mut.maf_like ${loci_file}
awk 'BEGIN{FS="\t"}{n=n+1;if(n==1){print $0};if($2==1){print $0}}' ${output_prefix}.mut.maf_like > ${output_prefix}.mut.maf_like.dec 
            """
        ],
        outputs={
            "MS_indels_file" : "*.mut.cln",
            "maf_life": "*.mut.maf_like",
            "maf_like_dec" : "*.mut.maf_like.dec" 
        },
        docker = "gcr.io/broad-getzlab-workflows/msmutect2_wolf:v2",
        resources={"mem": '8G'}
    )



def msmutect_workflow(
    tumor_bam, 
    tumor_bai,
    normal_bam,
    normal_bai, 
    pair_name,
    loci_file,
    ):
    
    # call tumor
    msT = msmutect_sample(
        bam = tumor_bam,
        bai = tumor_bai,
        loci = loci_file,
        sample_type = "tumor",
        output_prefix= pair_name
    )

    # call normal
    msN = msmutect_sample(
        bam = normal_bam,
        bai = normal_bai,
        loci = loci_file,
        sample_type = "normal",
        output_prefix= pair_name
    )

    # merge
    postprecess_msindel(
        tumor_hist=msT["hist"],
        normal_hist=msN["hist"],
        output_prefix = pair_name,
        loci_file=loci_file
    )

if __name__ == "__main__":
    with wolf.Workflow(workflow = msmutect_workflow) as e:
        e.run(
            tumor_bam = "/mnt/nfs/wgs_ref/hg38_bam/1205d516-e6e7-49b4-9681-8b706a7b211d_wgs_gdc_realn.bam", 
            tumor_bai = "/mnt/nfs/wgs_ref/hg38_bam/1205d516-e6e7-49b4-9681-8b706a7b211d_wgs_gdc_realn.bai",
            normal_bam = "/mnt/nfs/wgs_ref/hg38_bam/d5011161-3d6b-4802-8a64-614f6d2298f9_wgs_gdc_realn.bam",
            normal_bai = "/mnt/nfs/wgs_ref/hg38_bam/d5011161-3d6b-4802-8a64-614f6d2298f9_wgs_gdc_realn.bai", 
            pair_name = "testrun",
            loci_file = "/mnt/nfs/wgs_ref/hg38_1_to_15_loci.phobos",
            )