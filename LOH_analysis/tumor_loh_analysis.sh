# 08.20.20
#!/bin/bash

# cram paths /gscmnt/gc2547/griffithlab/mrichters/gbm/bam_directories/final_names_somatic_exome_dirs.tsv

name=$1
cram=$(echo $2 | sed 's/$/\/tumor.cram/')
ref_fasta="/gscmnt/gc2560/core/model_data/2887491634/build21f22873ebe0486c8e6f69c15435aa96/all_sequences.fa"

gl_name=$(grep $(echo $name | cut -d'-' -f1) ../../bam_directories/batch3_germline_paths.tsv | cut -f1)
working_dir="/gscmnt/gc2547/griffithlab/mrichters/gbm/gbm_antigens_twck/loh_analysis/results/$gl_name"

mkdir -p $working_dir
mkdir -p logs
mkdir -p bam_files

# Cram to bam
LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -J $name.cram-to-bam -q research-hpc -a 'docker(zlskidmore/samtools)' -M 15000000 -R 'select[mem>=15000] span[hosts=1] rusage[mem=15000]' -oo logs/$name.cram-to-bam.log "samtools view -T $ref_fasta -b $cram -o bam_files/$name.bam"

/gscmnt/gc2547/griffithlab/mrichters/bwait -n $name.cram-to-bam -s 30

# Index bam
LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -J $name.bam-index -q research-hpc -a 'docker(zlskidmore/samtools)' -M 15000000 -R 'select[mem>=15000] span[hosts=1] rusage[mem=15000]' -oo logs/$name.bam-index.log "samtools index bam_files/$name.bam"

/gscmnt/gc2547/griffithlab/mrichters/bwait -n $name.bam-index -s 30

# for each chrom in sample
for chrom in $working_dir/c*het.positions.txt; do

	chr=$(echo $(basename $chrom) | cut -d'.' -f1)
	#x_axis=$(echo "${chr^}" | sed 's/$/ Position/')

	# Run bam-readcount
	LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -J $name.$chr.bamrc -q research-hpc -a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' -M 5000000 -R 'select[mem>5000] span[hosts=1] rusage[mem=5000]' -e logs/$name.bam-readcount.err -oo logs/$name.$chr.bam-readcount.out "bam-readcount -w 1 -b 20 -q 20 -f $ref_fasta -l $chrom bam_files/$name.bam $chr"

	/gscmnt/gc2547/griffithlab/mrichters/bwait -n $name.$chr.bamrc -s 30
	
	# move bam-reacount output from stdout to tsv in results dir
	grep "^$chr" logs/$name.$chr.bam-readcount.out > $working_dir/$name.$chr.readcounts.tsv
	
	# Calculate tumor VAFs
	LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -J $name.loh_plot -q research-hpc -a 'docker(rocker/tidyverse)' -M 5000000 -R 'select[mem>5000] span[hosts=1] rusage[mem=5000]' -oo logs/$name.loh_plot.log Rscript find_tumor_vafs.R $working_dir $name.$chr.readcounts.tsv $name.$chr.tumor.VAFs.table #$gl_name.$chr.het.snps.table $name.$chr.loh.pdf $x_axis $name

	#Segmentation
	#download scp_seg_plots folders (containing tumor.chr.tumor.VAFs.table and chr.het.snps.table per cohort) and run ~/Desktop/GBM_Project/LOH_plots/run_loh.sh to get final plots with segmentation

done

echo "tumor $name complete"

# remove temp files
rm -rf bam_files/$name.bam
rm -rf bam_files/$name.bam.bai

echo "$name bam files removed"


