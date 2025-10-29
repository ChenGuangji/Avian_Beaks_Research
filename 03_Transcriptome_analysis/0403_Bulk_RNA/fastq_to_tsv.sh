#### fastp ####
#!/bin/bash
set -xveo pipefail
hostname
id=$1
in1=$2
in2=$3
outpath=$4
out1="${outpath}/${id}/${id}_R1.fq.gz"
out2="${outpath}/${id}/${id}_R2.fq.gz"

mkdir -p "${outpath}/${id}" && \
        fastp --thread 10 -i ${in1} -I ${in2} -o ${out1} -O ${out2} && \
        echo "${id} done at `date`"

#### hisat2 ####
#!/usr/bash
set -xveo pipefail

if [ $# -lt 5 ]
then
    echo "sh $0 <fa> <r1> <r2> <id> <path>"
        exit 0
fi

fa=$1
r1=$2
r2=$3
id=$4
path=$5

#mkdir -p ${path}/${id}/ && cd ${path}/${id}/ && \
mkdir -p ${path}/${id}/ && cd ${path}/${id}/ && gzip -dc ${r1} > ${path}/${id}/R1.fq && gzip -dc ${r2} > ${path}/${id}/R2.fq && \
        hisat2 --dta -x ${fa} -1 ${path}/${id}/R1.fq -2 ${path}/${id}/R2.fq -q --phred33 -p 10 | samtools view -b -@ 10 - -o ${path}/${id}/${id}.hisat2.bam && \
        samtools sort -n -@ 10 -o ${path}/${id}/${id}.hisat2.sortn.bam ${path}/${id}/${id}.hisat2.bam && echo "alignment sort n done" && \
        samtools sort -@ 10 -o ${path}/${id}/${id}.hisat2.sort.bam ${path}/${id}/${id}.hisat2.bam && echo "alignment sort done" && \
        rm -rf ${path}/${id}/R1.fq ${path}/${id}/R2.fq

#### htseq-count ####
htseq-count -f bam -n 40 \
 ../02.hisat2/*/*.hisat2.sort.bam \
 ../00.ref/gene.gtf --type gene --idattr gene_name 1>gene.count.tsv 2>gene.count.log