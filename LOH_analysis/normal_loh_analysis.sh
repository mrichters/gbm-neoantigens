# 08.19.20
#!/bin/bash

# germline paths batch3 /gscmnt/gc2547/griffithlab/mrichters/gbm/bam_directories/batch3_germline_paths.tsv
name=$1
input_dir=$2
ref_fasta="/gscmnt/gc2560/core/model_data/2887491634/build21f22873ebe0486c8e6f69c15435aa96/all_sequences.fa"
working_dir="/gscmnt/gc2547/griffithlab/mrichters/gbm/gbm_antigens_twck/loh_analysis/results/$name"
mkdir -p $working_dir

#GATK SelectVariants w/ germline vcf - selecting only SNPs
LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -J $name.selectvar -q research-hpc -a 'docker(broadinstitute/gatk)' -M 10000000 -R 'select[mem>10000] span[hosts=1] rusage[mem=10000]' -oo logs/$name.selectvariants.log gatk SelectVariants -R $ref_fasta -V $input_dir/annotated.filtered.final.vcf.gz -select-type SNP -O $working_dir/snps.vcf.gz

/gscmnt/gc2547/griffithlab/mrichters/bwait -n $name.selectvar -s 30

#GATK VariantsToTable w/ germline vcf snps
LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -J $name.vtt -q research-hpc -a 'docker(broadinstitute/gatk)' -M 10000000 -R 'select[mem>10000] span[hosts=1] rusage[mem=10000]' -oo logs/$name.variantstt.log gatk VariantsToTable -V $working_dir/snps.vcf.gz -F CHROM -F POS -GF GT -GF AD -GF DP -O $working_dir/snps.table

/gscmnt/gc2547/griffithlab/mrichters/bwait -n $name.vtt -s 30

#Calculate normal VAFs/find good het positions (40-60%)
LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -J $name.vafs_r -q research-hpc -a 'docker(rocker/tidyverse)' -M 10000000 -R 'select[mem>10000] span[hosts=1] rusage[mem=10000]' -oo logs/$name.vafs_r.log Rscript find_normal_vafs.R $working_dir snps.table

/gscmnt/gc2547/griffithlab/mrichters/bwait -n $name.vafs_r -s 30

# parse heterzygous.positions.txt and heterozygous.snps.table by chromosome
chromosomes=$(cat $working_dir/het.positions.txt | cut -f1 | sort -V | uniq)
for chrom in $chromosomes; do grep -w $chrom $working_dir/het.positions.txt > $working_dir/$chrom.het.positions.txt; done
for chrom in $chromosomes; do grep -w $chrom $working_dir/het.snps.table > $working_dir/$chrom.het.snps.table; done

