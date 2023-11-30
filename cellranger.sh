#!/usr/bin/bash

#export NXF_VER=22.
export NXF_SINGULARITY_CACHEDIR=/home/cpc/jiakuan.zhao/sing_img
export TMPDIR=/localscratch
export R_TempDir=/localscratch


nextflow run nf-core/scrnaseq -r 2.1.0  \
    -w /lustre/scratch/users/jiakuan.zhao/nextflow/ \
	--input ./samplesheet.csv \
	--outdir /home/cpc/jiakuan.zhao/datasets/Zindel_EGFR_Adhension_NC_2021/raw_counts/ \
	--genome GRCm38 \
	--aligner cellranger \
	--custom_config_base /home/cpc/jiakuan.zhao/nextFlow/nfcore_config \
	-profile hmgu,cpu_normal,singularity \
	-resume \
	-bg >> ./cellranger.log.txt
