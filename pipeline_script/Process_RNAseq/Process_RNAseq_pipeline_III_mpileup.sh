#!/bin/bash

sample="${1}"
refGenome="${2}"
mismatchN="${5}"
hetsMetaPositions="${3}"
bam="${4}"
pipelineDir="${6:-./Process_RNAseq}"

outDir="${pipelineDir}/${sample}/star2pss_EndtoEnd_wasp/mismatch${mismatchN}/mpileup"
mkdir -p ${outDir}

pileup_out="${outDir}/Allchr_hets_all_transcript.pileup"
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>> We start with sample: ${sample}"

samtools mpileup -d 0 -B -s -f ${refGenome} -l ${hetsMetaPositions} ${bam} > ${pileup_out}

