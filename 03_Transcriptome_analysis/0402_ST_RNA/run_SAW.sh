#!/usr/bash
export PATH=/share/software/singularity/3.5.3/bin/:${PATH}

dir="00.data/HHXX/"
h5="00.data/HHXX/HHXX.barcodeToPos.h5"
ipr="00.data/HHXX/HHXX.ipr"
icf="00.data/HHXX/HHXX.tar.gz"
SName="Chick"
id="CH29"
index="00.data/refGenome/Chick/index/"
gtf="00.data/refGenome/Chick/genes.gtf"
path="work/"

fq1=$(find ${dir} -name "*.fq.gz"|paste -sd,)
outDir="${path}/${id}"
mkdir -p ${outDir} && \
bash pipeline/SAW/Scripts/stereoPipeline_v7.1.sh \
    -sif pipeline/SAW_7.1.sif \
    -splitCount 16 \
    -maskFile ${h5} \
    -fq1 ${fq1}  \
    -speciesName ${SName} \
    -tissueType ${id} \
    -refIndex ${index} \
    -annotationFile ${gtf} \
    -rRNAremove Y \
    -threads 16 \
    -outDir ${outDir} \
    -imageRecordFile ${ipr} \
    -imageCompressedFile ${icf} \
    -doCellBin Y