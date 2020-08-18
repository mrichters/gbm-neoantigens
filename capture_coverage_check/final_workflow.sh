# 08.17.20
# variant coverage check for custom capture sequencing
#!/bin/bash

name=$1
cram=$2
ref_fasta='/gscmnt/gc2560/core/model_data/00b5d92847b94b06ac93e56ff551aae8/builde171e4340cb1492ebca9b5e3fc669b05/all_sequences.fa'
# cram paths: /gscmnt/gc2547/griffithlab/mrichters/gbm/bam_directories/ts_cram_paths.tsv

# step 1 sort vcf
# /gscuser/cmiller/usr/bin/isub
# vcf-sort -c gbm_unique_variants_b38.vcf > gbm_unique_variants_b38.sorted.vcf

# step 2 create a master vcf for each sample
#cat /gscmnt/gc2547/griffithlab/mrichters/gbm/bam_directories/ts_cram_paths.tsv | awk '{print "cp gbm_unique_variants_b38.sorted.vcf "$1".sorted.vcf && sed -i '\''s/TUMOR/"$1"/'\'' sample_vcfs/"$1".sorted.vcf"}' | /bin/bash

# step 3 cram to bam
LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -J $name.cram-to-bam -q research-hpc -a 'docker(zlskidmore/samtools)' -M 15000000 -R 'select[mem>=15000] span[hosts=1] rusage[mem=15000]' -oo logs/$name.cram-to-bam.log "samtools view -T $ref_fasta -b $cram -o bam_files/$name.bam"

/gscmnt/gc2547/griffithlab/mrichters/bwait -n $name.cram-to-bam -s 30

LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -J $name.bam-index -q research-hpc -a 'docker(zlskidmore/samtools)' -M 15000000 -R 'select[mem>=15000] span[hosts=1] rusage[mem=15000]' -oo logs/$name.bam-index.log "samtools index bam_files/$name.bam"

/gscmnt/gc2547/griffithlab/mrichters/bwait -n $name.bam-index -s 30

# step 4 bam-readcount
LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -J $name.bam-readcount -q research-hpc -a 'docker(mgibio/cle)' -M 10000000 -R 'select[mem>10000] span[hosts=1] rusage[mem=10000]' -oo logs/$name.bam-readcount.log python -u /usr/bin/bam_readcount_helper.py sample_vcfs/$name.sorted.vcf $name $ref_fasta bam_files/$name.bam bam-readcount_files/ 20 30

/gscmnt/gc2547/griffithlab/mrichters/bwait -n $name.bam-readcount -s 30

echo "$name complete"
 
# remove bam file
rm -rf bam_files/$name.bam
rm -rf bam_files/$name.bam.bai

echo "$name bam files removed"

# step 5 annotate vcfs with bam-readcount results
LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -J $name.vcf-annotator -q research-hpc -a 'docker(griffithlab/vatools:4.1.0)' -M 10000000 -R 'select[mem>10000] span[hosts=1] rusage[mem=10000]' -oo logs/$name.vcf-annotator.log "vcf-readcount-annotator sample_vcfs/$name.sorted.vcf bam-readcount_files/${name}_bam_readcount_snv.tsv DNA -t snv -o annotated_vcfs/$name.snv.vcf"

/gscmnt/gc2547/griffithlab/mrichters/bwait -n vcf-annotator -s 10

LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -J $name.vcf-annotator2 -q research-hpc -a 'docker(griffithlab/vatools:4.1.0)' -M 10000000 -R 'select[mem>10000] span[hosts=1] rusage[mem=10000]' -oo logs/$name.vcf-annotator2.log "vcf-readcount-annotator annotated_vcfs/$name.snv.vcf bam-readcount_files/${name}_bam_readcount_indel.tsv DNA -t indel -o annotated_vcfs/$name.annotated.vcf"

/gscmnt/gc2547/griffithlab/mrichters/bwait -n $name.vcf-annotator2 -s 10

# remove intermediate vcf
rm -rf annotated_vcfs/$name.snv.vcf

# step 6 convert annotated vcf to tsv file
LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -J $name.vtt -q research-hpc -a 'docker(mgibio/gatk-cwl:3.6.0)' -M 10000000 -R 'select[mem>10000] span[hosts=1] rusage[mem=10000]' -oo logs/$name.vtt.log /usr/bin/java -Xmx8g -jar /opt/GenomeAnalysisTK.jar -T VariantsToTable --variant annotated_vcfs/$name.annotated.vcf -F CHROM -F POS -F REF -F ALT -GF GT -GF DP -GF AD --out tsv_files/$name.annotated.tsv --reference_sequence $ref_fasta

echo "$name final tsv complete"
