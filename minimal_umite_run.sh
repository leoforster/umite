#!/bin/bash

STAR_GENOME_DIR="" # path to STAR genome 
PATH_TO_GTF="" # path to GTF file e.g. Mus_musculus.GRCm38.100.gtf from ensembl
FASTQ_DIR="" # path to fastq directory
OUTDIR="" # path to output directory

PARALLEL_JOBS=8
FASTQ_FILE_SUFFIX=".fastq.gz"
FASTQ_R1_PREFIX="_1"
FASTQ_R2_PREFIX="_2"

echo "Starting umite run at $(date)"

mkdir -p $OUTDIR
FASTQ_R1_SUFFIX="${FASTQ_R1_PREFIX}${FASTQ_FILE_SUFFIX}"

# List cell IDs
samples=$(ls ${FASTQ_DIR}/*${FASTQ_R1_SUFFIX} | sed "s/${FASTQ_R1_SUFFIX}//" | xargs -n1 basename | sort -u)
echo $samples > samples_used.txt

# umiextract with flexible TSO
umiextract \
	-1 $(echo "$samples" | xargs -I{} echo -n "${FASTQ_DIR}/{}${FASTQ_R1_PREFIX}${FASTQ_FILE_SUFFIX} ") \
	-2 $(echo "$samples" | xargs -I{} echo -n "${FASTQ_DIR}/{}${FASTQ_R2_PREFIX}${FASTQ_FILE_SUFFIX} ") \
	-d "$OUTDIR" \
	-c "$PARALLEL_JOBS" \
	--umilen 8 \
	--anchor "ATTGCGCAATG" \
	--trailing "GGG" \
	--search_region 30 \
	--fuzzy_umi 

# STAR alignments
run_star() {
    sample="$1"
    r1="${OUTDIR}/${sample}${FASTQ_R1_PREFIX}_umiextract${FASTQ_FILE_SUFFIX}" # _1_umiextract.fastq.gz
    r2="${OUTDIR}/${sample}${FASTQ_R2_PREFIX}_umiextract${FASTQ_FILE_SUFFIX}" # _2_umiextract.fastq.gz
    outprefix="${OUTDIR}/${sample}_"

    STAR \
	--genomeDir "$GENOME_DIR" \
        --readFilesIn "$r1" "$r2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$outprefix" \
        --outSAMtype BAM Unsorted \
        --clip3pAdapterSeq CTGTCTCTTATACACATCT \
        --runThreadN 1 \
        --genomeLoad LoadAndKeep
}

export -f run_star
export sorted_samples
export GENOME_DIR FASTQ_DIR OUTDIR PARALLEL_JOBS

# load STAR genome
STAR \
	--genomeDir "$GENOME_DIR" \
	--genomeLoad LoadAndExit \
	--runThreadN "$PARALLEL_JOBS"

# run STAR in parallel
echo "$samples" | parallel -j "$PARALLEL_JOBS" run_star {}

# remove STAR genomes
STAR \
    --genomeDir "$GENOME_DIR" \
    --genomeLoad Remove

# Sort BAMs by readname
echo "$samples" | parallel -j "$PARALLEL_JOBS" "samtools sort -n $OUTDIR/{}_Aligned.out.bam > $OUTDIR/{}_Aligned.sort.bam && rm $OUTDIR/{}_Aligned.out.bam"

# prepare umicount by parsing annotation from GTF
umicount \
	-gtf "${PATH_TO_GTF}" \
	--GTF_dump "${OUTDIR}/umite_GTF_dump.pkl"

# count UMI reads with umicount
umicount \
	--bams $(echo "$samples" | xargs -I{} echo -n "${OUTDIR}/{}_Aligned.sort.bam ") \
	-d "$OUTDIR" \
	-c "$PARALLEL_JOBS" \
	--GTF_skip_parse "${OUTDIR}/umite_GTF_dump.pkl" \
	--mm_count_primary \
	--UMI_correct

echo "Finished umite run at $(date)"
