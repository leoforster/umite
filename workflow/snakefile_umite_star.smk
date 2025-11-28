from os.path import exists, basename, dirname, join
from os import rename
import pandas as pd

######################### umite SNAKEMAKE WORKFLOW #########################
# 
# this workflow outlines a umite run on a set of fastq files. inputs and 
# outputs are specified in the config file snakeconfig.yaml , along with 
# individual program args and other parameters. please first ensure the 
# required data are written to snakeconfig.yaml before running this work-
# flow via snakemake --jobs 99 --sdm conda -s snakefile_umite_star.smk

# read from the config file
configfile: 'snakeconfig.yaml'

# check required inputs
for i in [config['samples_file'], config['fastq_dir'], config['output_dir'], config['log_dir'],
          config['reference']['genome'], config['reference']['annotation']]:

    if not exists(i):
        raise FileNotFoundError(f"invalid path: {i}")

# populate run data
runID = config['runID']
cores = config['threads']

# determine output tables to generate
basecols = []
if config['umicount']['combine_unspliced']:
    basecols += ['U', 'R']
else:
    basecols += ['UE', 'UI', 'RE', 'RI']
if config['umicount']['dedupe_umis']:
    basecols += ['D']

file_ext = [f".{i}.tsv" for i in basecols] # .D.tsv, .UE.tsv, ...

# read samples list from file
with open(config['samples_file'], 'r') as f:
    samples = f.read().splitlines()

print(f"running umite snakemake workflow {runID} on {len(samples)} samples, with fastqs in {config['fastq_dir']}")

############# START WORKFLOW 

# define outputs
rule all:
    input:
        expand(join(config['output_dir'], f"{runID}_umite{{ext}}"), ext=file_ext)

rule prepare_star_indices:
    input:
        ref_genome = ancient(config['reference']['genome']),
        gene_annotation = ancient(config['reference']['annotation'])
    output:
        directory(join(dirname(config['reference']['genome']), 'star_genome_index'))
    log: join(config['log_dir'], f"{runID}.star_index.log")
    conda: 'umite_conda.yaml'
    threads: cores
    params:
        extra=config['extra_args']['star_index']
    shell:
        '''
        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.ref_genome} \
            --sjdbGTFfile {input.gene_annotation} \
            {params.extra}
        '''

def format_umiextract_output_name(suffix):
    readpart, filepart = suffix.split('.', 1) # e.g. ['_R1','fastq.gz']
    return f"{readpart}_umiextract.{filepart}" # _R1_umiextract.fastq.gz

rule umiextract:
    input:
        R1s = expand(join(config['fastq_dir'], f"{{sample}}{config['R1_suffix']}"), sample=samples),
        R2s = expand(join(config['fastq_dir'], f"{{sample}}{config['R2_suffix']}"), sample=samples)
    output:
        # these have form sample_R1_umiextract.fastq.gz
        temp(expand(join(config['output_dir'], f"{{sample}}{format_umiextract_output_name(config['R1_suffix'])}"), sample=samples)),
        temp(expand(join(config['output_dir'], f"{{sample}}{format_umiextract_output_name(config['R2_suffix'])}"), sample=samples))
    log: join(config['log_dir'], f"{runID}.umiextract.log")
    conda: 'umite_conda.yaml'
    threads: cores
    params:
        outdir=config['output_dir'],
        anchor_seq=config['umiextract']['anchor_seq'],
        trailing_seq=config['umiextract']['trailing_seq'],
        umilen=config['umiextract']['umilen'],
        search_region=config['umiextract']['search_region'],
        min_seqlen=config['umiextract']['min_seqlen'],
        only_umi=('--only_umi' if config['umiextract']['only_umi'] else ''),
        fuzzy_umi=('--fuzzy_umi' if config['umiextract']['fuzzy_umi'] else ''),
        anchor_mismatches=(f"--anchor_mismatches {config['umiextract']['anchor_mismatches']}" if config['umiextract']['fuzzy_umi'] else ''),
        anchor_indels=(f"--anchor_indels {config['umiextract']['anchor_indels']}" if config['umiextract']['fuzzy_umi'] else ''),
        trailing_hamming=(f"--trailing_hamming {config['umiextract']['trailing_hamming']}" if config['umiextract']['fuzzy_umi'] else ''),
        extra=config['extra_args']['umiextract']
    shell:
        '''
        umiextract \
            -1 {input.R1s} \
            -2 {input.R2s} \
            -d {params.outdir} \
            -c {threads} \
            --umilen {params.umilen} \
            --anchor {params.anchor_seq} \
            --trailing {params.trailing_seq} \
            --search_region {params.search_region} \
            --min_seqlen {params.min_seqlen} \
            {params.only_umi} {params.fuzzy_umi} \
            {params.anchor_mismatches} {params.anchor_indels} {params.trailing_hamming} \
            {params.extra}
        '''

rule star_genome_load:
    input:
        join(config['reference']['genome'], 'star_genome_index', 'genomeParameters.txt') # created with index
    output:
        temp(touch(join(config['log_dir'], 'star_genome_load.completed')))
    shell:
        '''
        STAR \
            --genomeLoad LoadAndExit \
            --genomeDir {{os.dirname(input)}}
        '''

rule star_alignment:
    input:
        R1_extract = join(config['output_dir'], f"{{sample}}{format_umiextract_output_name(config['R1_suffix'])}"),
        R2_extract = join(config['output_dir'], f"{{sample}}{format_umiextract_output_name(config['R2_suffix'])}"),
        indices = rules.prepare_star_indices.output
    output:
        temp(join(config['output_dir'], '{sample}_Aligned.out.bam'))
    log: join(config['log_dir'], f"{runID}_{{sample}}.star_align.log")
    conda: 'umite_conda.yaml'
    threads: cores
    params:
        outprefix=join(config['output_dir'], '{sample}_'),
        extra=config['extra_args']['star_align']
    shell:
        '''
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.indices} \
            --genomeLoad LoadAndKeep \
            --readFilesIn {input.R1_extract} {input.R2_extract} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.outprefix} \
            --outSAMtype BAM Unsorted \
            {params.extra}
        '''

rule star_genome_unload:
    input:
        bams=expand(join(config['output_dir'], '{sample}_Aligned.out.bam'), sample=samples),
        genome_load_state=join(config['log_dir'], 'star_genome_load.completed')
    output:
        touch(join(config['log_dir'], f"{runID}_{{sample}}.star_genome_unload.complete"))
    params:
        genome_dir=directory(join(dirname(config['reference']['genome']), 'star_genome_index'))
    shell:
        '''
        STAR \
            --genomeLoad Remove \
            --genomeDir {params.genome_dir}
        '''

rule sort_aligned_reads:
    input:
        rules.star_alignment.output
    output:
        join(config['output_dir'], '{sample}.namesort.bam')
    log: join(config['log_dir'], f"{runID}_{{sample}}.sort_reads.log")
    conda: 'umite_conda.yaml'
    threads: cores
    shell:
        'samtools sort -@ {threads} -n {input} > {output}'

rule parse_dump_GTF:
    input:
        ancient(config['reference']['annotation'])
    output:
        join(dirname(config['reference']['annotation']), 'umicount_GTF_dump.pkl')
    log: join(config['log_dir'], f"{runID}.gtf_dump.log")
    conda: 'umite_conda.yaml'
    shell:
        'umicount -g {input} --GTF_dump {output}'

rule umicount:
    input:
        bams = expand(join(config['output_dir'], '{sample}.namesort.bam'), sample=samples),
        gtf_dump = ancient(rules.parse_dump_GTF.output)
    output:
        expand(join(config['output_dir'], f"{runID}_umite{{ext}}"), ext=file_ext)
    log: join(config['output_dir'], f"{runID}_umicount.log")
    conda: 'umite_conda.yaml'
    threads: cores
    params:
        outdir=config['output_dir'],
        min_read_mapQ=config['umicount']['min_read_mapQ'],
        dedupe_umis=('' if config['umicount']['dedupe_umis'] else '--no_dedup'),
        combine_unspliced=('--combine_unspliced' if config['umicount']['combine_unspliced'] else ''),
        count_multimappers=('--mm_count_primary' if config['umicount']['count_multimappers'] else ''),
        multiple_primary_action=config['umicount']['multiple_primary_action'],
        correct_umis=('--UMI_correct' if config['umicount']['correct_umis'] else ''),
        hamming_threshold=(f"--hamming_threshold {config['umicount']['hamming_threshold']}" if config['umicount']['correct_umis'] else ''),
        count_ratio_threshold=(f"--count_ratio_threshold {config['umicount']['count_ratio_threshold']}" if config['umicount']['correct_umis'] else ''),
        extra=config['extra_args']['umicount']
    run:
        shell('''
        umicount \
            --bams {input.bams} \
            -d {params.outdir} \
            -c {threads} \
            --GTF_skip_parse {input.gtf_dump} \
            --min_read_mapQ {params.min_read_mapQ} \
            --multiple_primary_action {params.multiple_primary_action} \
            {params.dedupe_umis} {params.combine_unspliced} {params.count_multimappers} \
            {params.correct_umis} {params.hamming_threshold} {params.count_ratio_threshold} \
            {params.extra}
        ''')

        for ext in file_ext:
            os.rename(join(params.outdir, f"umite{ext}"),
                      join(params.outdir, f"{runID}_umite{ext}"))
