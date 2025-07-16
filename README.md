# umite

**UMI extraction & counting for Smart‑seq3 scRNA‑seq libraries**

umite *unites* two tools that turn raw Smart‑seq3 FASTQ files into gene‑by‑cell count tables for downstream analysis:

![umite_scheme](images/umite_scheme.png)

| Step | Tool              | Purpose |
|------|------------------|---------|
| ①    | **umiextract** | Detect and label UMI-containing reads in a FASTQ file. Use optional fuzzy UMI matching to incrase the number of detected UMIs. |
| ②    | *(external)*     | Align reads using e.g. the splice-aware aligner, **STAR**. Then sort reads by read-name using **samtools** |
| ③    | **umicount**   | Parse aligned reads and assign reads/UMIs to genes & exons, while deduplicate and (optionally) error‑correcting UMIs.  |

Additional information is available in [our paper](https://XXX)

---

## Usage

### Installation

```bash
git clone https://github.com/leoforster/umite.git
cd umite
pip install -e .
```

if this has issuses you can try using `python -m pip install -e .` to install via the local python interpreter instead. You can also install the `pytest` suite by calling `pip install -e .[dev]`.

umite requires **Python≥3.7** and the packages *HTSeq*, *regex* and *RapidFuzz*, which are automatically installed in the example above.
For the alignment step, you will need **samtools** and ideally the **STAR** aligner.

### Pipeline

The repository ships with [`minimal_umite_run.sh`](minimal_umite_run.sh) , a minimal script for running the umite pipeline.

#### ① Detect UMI-containing reads with umiextract

in this example we have paired-end Smart-seq3 libraries from two cells: cellA and cellB. 
We enable error-tolerant detection of UMI-containing reads with `--fuzzy_umi` , and the script takes paired-end FASTQs as input. Note that parallel processing of libraries is possible using the  `--cores` argument.

```bash
umiextract \
    -1 cellA_R1.fastq.gz cellB_R1.fastq.gz \
    -2 cellA_R2.fastq.gz cellB_R2.fastq.gz \
    --umilen 8 \
    --fuzzy_umi 	# enable error-tolerant UMI detection
```

This will output modified FASTQs for each cell, e.g. `cellA_R1_umiextract.fastq.gz` . These files are essentially idential to the input where detected UMIs are trimmed from read sequences and appended to the readnames.

#### ② Read alignment

Here using `STAR` , note that `--genomeDir` requires a pre-existing genome index.

```bash
STAR \
  --genomeDir /path/to/STAR/index \	
  --readFilesIn cellA_R1_umiextract.fastq.gz \
                 cellA_R2_umiextract.fastq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM Unsorted
```

then sort BAM files by read-name, such that paired-reads are adjacent in the file:

```bash
samtools sort -n -o cellA_sorted.bam Aligned.out.bam
```

#### ③  Quantify counts with umicount 

Finally, quantify counts per gene per cell for UMI-containing (`U`) and internal-reads (`R`) using `umicount`. This step can process multiple BAMs in parallel by setting `--cores` , and requires a GTF file containing genome annotations (e.g. from Ensembl). In the example, setting `--mm_count_primary` causes the primary alignment to be counted for multimapping reads, and `--UMI_correct` enables gene-wise UMI correction by collapsing sequencing errors using directional Hamming-distances. Other options are detailed below.

```bash
umicount \
    --bams cellA_sorted.bam cellB_sorted.bam \
    --gtf Mus_musculus.GRCm39.102.gtf \
    --mm_count_primary \
    --UMI_correct
```

This will produce the following output files:

```
umite.U.tsv   # unique UMI counts (per gene × cells)
umite.R.tsv   # internal-read (i.e. non-UMI) counts
umite.D.tsv   # counts for UMI-duplicates (for QC)
log.txt           # summary of processing & statistics (optional, enabled with --logfile)
```

The output counts matrices contain samples (cells) in rows, with columns denoted by genes parsed from the GTF. These begin with read categories (e.g. `_unmapped`, `_multimapping`, `_ambiguous`) that report the fate of every read from the `umiextract` FASTQ according to the following schema:

![umite_read_categories](images/umite_readcat_scheme.png)

---

## Command‑line reference

Run `umiextract -h` or `umicount -h` for the full list of options.

### umiextract

| Flag | Description | Default |
|------|-------------|---------|
| `-1 / --read1` | *R1* FASTQ files (space‑separated) | *required* |
| `-2 / --read2` | *R2* FASTQ files (same order as R1) | – |
| `-d / --output_dir` | Where to write processed FASTQs | `.` |
| `-c / --cores` | Parallel workers (one sample per core) | `4` |
| `-l / --logfile` | Path to log file | `sys.stdout` |
| `--umilen` | UMI length in bp | `8` |
| `--anchor` | Pre‑UMI anchor (TSO) sequence | `ATTGCGCAATG` |
| `--trailing` | Post‑UMI trailing sequence | `GGG` |
| `--search_region` | Sequence cutoff to search for UMI | `-1` (whole read) |
| `--fuzzy_umi` | Enable mismatch/indel‑tolerant search | off |
| `--anchor_mismatches` | Max mismatches in anchor | `2` |
| `--anchor_indels` | Max indels in anchor | `1` |
| `--trailing_hamming_threshold` | Max Hamming distance in trailing | `2` |
| `--min_seqlen` | Minimum remaining sequence after trimming UMI | `-1` |
| `--only_umi` | Drop reads that lack a detectable UMI | off |

### umicount

| Flag | Description | Default |
|------|-------------|---------|
| `-f / --bams` | Read‑name–sorted BAM files | *required* |
| `-d / --output_dir` | Output directory | `.` |
| `-c / --cores` | Parallel workers (one BAM per core) | all cores |
| `-l / --logfile` | Path to log file | `sys.stdout` |
| `-g / --gtf` | Ensembl‑style GTF annotation | *required* (see below) |
| `--tmp_dir` | Directory to save temporary files | `--output_dir` |
| `--no_dedup` | Skip deduplication and report all UMI-reads | off |
| `--mm_count_primary` | Count primary alignment for multimapping reads | off |
| `--multiple_primary_action` | When a read has mutliple primary alignments: `warn`, `raise`, or `skip` | `warn` |
| `--min_read_mapQ` | Min mapQ to keep read | `0` |
| `--combine_unspliced` | If set, don't distinguish intronic and exonic reads | off |
| `--UMI_correct` | Enable gene‑wise UMI collapse by Hamming distance | off |
| `--hamming_threshold` | Hamming threshold for merging UMIs | `1` |
| `--count_ratio_threshold` | Only merge UMIs if one has 2*threhsold as many counts | `2` |

Of note, as GTF parsing can take several minutes, umicount implements the option to parse from the GTF file once and dump the contents to a `pickle` file. Using `--gtf` with `--GTF_dump` will enable dumping parsed GTF data to a `pickle` file which can be used as input for `umicount` with `--GTF_skip_parse` instead of `--gtf`. This functionality is useful when running multiple repeat quantifications, however generally `--gtf` is the better option. Here a minimal example of this functionality:

```bash
umicount \
	-gtf examplefile.gtf \
	--GTF_dump umite_GTF_dump.pkl

umicount \
	--bams example.bam
	--GTF_skip_parse umite_GTF_dump.pkl
```
