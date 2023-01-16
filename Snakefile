# pipeline for data predprocessing and normalization

# command: 'snakemake --cores N --config input="/path/to/input" \
# output="/path/to/output" \
# transcriptome_annotation="/path/to/transcriptome_annotation.gff" \
# genome="/path/to/genome_fasta" \
# index="/path/to/genome_idx" \
# ref_structures="/path/to/ref_structures" \
# reads_length_threshold=N'

workdir: config['output']

reads_length_threshold = str(config['reads_length_threshold'])
samples, = glob_wildcards(os.path.join(config['input'], "{samples}_R1_001.fastq.gz"))


rule all:
    input:
        qc_check_1 = os.path.join(config['output'],"qc/00_raw/multiqc_report.html"),
        qc_check_2 = os.path.join(config['output'],"qc/01_umi_removed/multiqc_report.html"),
        read_1_umiremoved_both_gz = expand(os.path.join(config['output'],"{index}_R1.umiremoved_both"  +
                                                           ".clipped.fastq.gz"),
            index=samples),
        read_2_umiremoved_both_gz = expand(os.path.join(config['output'],"{index}_R2.umiremoved_both"  +
                                                           ".clipped.fastq.gz"),
            index=samples),
        no_overlaps = os.path.join(config['output'],'transcripts.no_overlaps.gff'),
        logs_before_filtering = expand(os.path.join(config['output'],"{index}" +
                                                                        ".sorted.dedup.stats"),
            index=samples),
        sorted_dedup_alignment = expand(os.path.join(config['output'],"{index}" +
                                                        ".sorted.dedup.bam"),
            index=samples),
        logs_after_filtering = expand(os.path.join(config['output'],"{index}" +
                                                                    ".sorted.dedup.unique.stats"),
            index=samples),
        filtered_normalized_counts = expand(os.path.join(config['output'],"{index}" +
                                                            ".normalized_counts.filtered"),
            index=samples)


# raw reads predprocessing
# preprocessing raw reads steps and generate:
# - reads for alignment step (in gz archive)
# - FastQC reports
# - logs and trim stats
include: "01_raw_reads/01_raw_reads_predprocessing/Snakefile"


# generate the alignment from umiremoved both sides reads with segemehl
include: "01_raw_reads/02_alignment_and_deduplication/Snakefile"


# generate files with raw counts for positive and negative strands separately from the alignment:

# OH-library counts preparation
include: "02_alignment/01a_OH_library_transcripts/Snakefile"

# cP-library counts preparation
include: "02_alignment/01b_cP_library_transcripts/Snakefile"


# generate files with normalized counts for high covered transcripts
include: "03_normalization/Snakefile"
