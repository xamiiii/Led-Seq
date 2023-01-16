# Led-Seq - <ins>l</ins>igation-<ins>e</ins>nhanced <ins>d</ins>ouble-end <ins>seq</ins>uence-based structure analysis of RNA
This snakemake pipeline automates the data analysis from raw sequencing reads to a normalized probing signal from Led-Seq experiments.
It includes the processing of two types of libraries that each capture one of the two RNA fragments resulting from lead-cleavage: 2'3'-cP and 5'-OH.


**Input data**

- Input directory with raw data: .fastq.gz files with forward/reverse reads for both libraries (e.g. LIBRARYNAME_R1_001.fastq.gz, LIBRARYNAME_R2_001.fastq.gz)
- genome (.fasta) and genome index (.idx) (prepare index using mapping tool segemehl)
- transcriptome annotation (.gff)
- structure (dot-bracket notation) and sequence for benchmark transcripts (see e.g. /inputfiles/benchmark_RNAs_structures.fa , can be an empty file)
- setting a minimum length for reads after adapter and UMI removal (default = 12 nt)

  > Important note: illumina adapter sequences must be individualized in the pipeline code.


**Output data**

Pipeline results include:

- adapter and UMI trimming statistics
- adapter-trimmed and umi-removed reads (in fastq.gz format)
- fastqc and multiqc reports
- alignments of UMI-deduplicated reads (.dedup.bam) and aligning statistics
- raw and normalized probing signal (in tsv format)


**Workflow description**

- Raw read processing (sub-pipeline 1-1)
  - Raw reads quality check
  - Adapter trimming with statistics collection
  - UMI extraction
  - Processed reads quality check
  - Quality trimming with statistics collecting


- Read alignment (sub-pipeline 1-2)
  - Mapping of reads to the genome
  - Sort and index alignment
  - UMI deduplication


- Raw count extraction (separate for 2'3'-cP and 5'-OH libraries) (sub-pipelines 2-1a, 2-1b)
  - Alignment filtering (keep only properly mapped pairs and only unique hits)
  - Intersect alignment with transcriptome annotations
  - Calculate the coverage of individual transcripts
  - Collect raw counts for all positions in transcriptome


- Signal normalization (sub-pipeline 3)
  - Normalize raw counts within transcripts
  - Extract transcripts that fulfill coverage criteria 
  - Filter probing signal (keep high covered transcripts) and add structure information (if reference available)


***
## Requirements

> setting up a conda env with all tools is highly recommend. 
> 
> Instructions can be found [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

python ≥ 3.6

Other requirements can be installed by

`pip install -r requirements`


**Optional:**

PyCharm plugin: SnakeCharm
***

## Usage

To run the pipeline it is required to specify all inputs
and run the following command:

_Snakemake can automatically determine which parts of the workflow can be run in parallel_

_By specifying more than one available core, i.e. `--cores 10` we can optimize the scheduling of jobs_

Usage sample:

```commandline
snakemake --cores N --config \
input=/path/to/folder_with_raw_reads \
output=path/to/output/folder \
reads_length_threshold=12 \
genome=/path/to/genome.fa \
index=/path/to/genome_index.fa.idx \
transcriptome_annotation=/path/to/transcriptome_annotation.gff \
ref_structures=/path/to/reference_RNA_structures.fa
```

 > Important note: all paths should be absolute paths
***


## Authors

Sarah von Löhneysen [lsarah@bioinf.uni-leipzig.de](lsarah@bioinf.uni-leipzig.de)

Julie Ozerova [julie@bioinf.uni-leipzig.de](julie@bioinf.uni-leipzig.de)

<!--If you use the pipeline, you may want to cite the following publication:
Tim Kolberg, Sarah von Löhneysen, Iuliia Ozerova, Karolin Wellner, Roland K. Hartmann, Peter F. Stadler and Mario Mörl (2023), "Led-Seq- ligation-enhanced double-end sequence-based structure analysis of RNA", Journal XX(XX), XX-XX.-->
