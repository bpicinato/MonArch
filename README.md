# MonArch: an automated pipeline to annotate circRNAs in prokaryotic RNA-Seq data

MonArch is an automated pipeline to annotate circRNAs in RNA-Seq reads. It searches for circRNA signatures (the chiastic alignment of the two halves of the read) in the reads and then groups close junctions into one circRNA.

If you use the contents of this repository, please cite:

Picinato, B. A., Franceschini-Santos, V. H., Zaramela, L., Vêncio, R. Z. N., Koide, T. (2025) **Archaea express circular isoforms of IS200/IS605-associated ωRNAs**. In preparation.

# Installation and requisites

MonArch was developed and tested on Ubuntu.

To use MonArch, download the git repository on your machine.

We recommend using [Conda](https://docs.conda.io/en/latest/) to install all MonArch dependencies. We provided a `.yml` file to create the environment with the dependencies in the `misc` directory. Navigate to the MonArch directory and create the environment with:
```
conda env create -f misc/monarch.yml
```
Before running MonArch, activate the environment with:
```
conda activate monarch
```

If you do not want to use Conda to install all the requisites (**not recommended**), make sure your system has Python (v=3.10.8, with the sys, os, BioPython, and datetime packages), Blast (v=2.10.1), and GNU Parallel (v=20160822) installed. Be aware that you must change how they are called in the main `monarch.sh` script since it assumes that the dependencies were installed with Conda and are found in the Conda directory. The versions listed are the ones used by MonArch; other versions of the dependencies might work, but they were not tested.

# Usage

MonArch requires three parameters: a **reference genome** file in FASTA format; a comma-separated list of FASTA files with the **reads** from which circRNAs should be detected; and the **path** to an output directory (it will be created if it does not exist).

Below is a list of all possible arguments; please take notice of the default behavior for the optional parameters and verify that they fit your analysis. The `-h` or `--help` arguments show this help message on the command-line interface.
```
                                   MonArch
         Automated pipeline to annotate circular RNAs in RNA-Seq data
------------------------------------------------------------------------------

Usage:

 monarch.sh --ref_genome <reference.fasta> --reads <reads1.fasta,reads2.fasta>
                --output <path/to/output> [--threads INT]
               [--prefix STR] [--max_circle_size INT]
               [--min_read_cov FLOAT] [--max_mismatch INT]
               [--dont_invert_strand]

Required arguments:

  --ref_genome         The reference genome to be used as subject for blast

  --reads              Reads to be used as query for blast. Or a comma separated
                       list of reads, if applicable

  --output             Path to output dir. Will be created if it does not exist

Optional arguments:

 --threads             Number of threads for the analysis (Default = 1)

 --prefix              Prefix for files created by MonArch (Default = circles)

 --max_circle_size     Maximum circle size (Default = 3500)

 --min_read_cov        Minimum read coverage by first and second hits together.
                       Must be a number between 0 and 1 (0 - 100%) (Default = 0.9)

 --max_mismatch        Maximum number of mismatches allowed for each alignment
                       (Default = 0)

 --dont_invert_strand  Do not invert input strand. Default behavior of MonArch
                       is to turn 'plus' strand turns into 'minus' and vice versa.
                       Depending on the sequencing protocol one might not want this.
```

We have provided a small test dataset in the `misc` folder. It will create a `test_output` folder in the `misc` directory. You can run it with:
```
monarch --ref_genome misc/Hsalinarum.fa --reads misc/test_input.fa --output misc/test_output
```

# Outputs

MonArch outputs can be divided into three categories: 1) alignment results (BLASTn results), 2) junction results (regarding circularization junctions before joining into _ensembles_), and 3) _ensemble_ results (final results after _ensemble_ formation - see [The Pipeline](#the-pipeline)). 

#### 1. blastn.tbl

BLASTn standard tab-delimited output compressed by gzip. Columns are:
```
sseqid  qseqid  qlen  length  qstart  qend  sstart  send  sstrand  mismatch
```
You can learn more about them with the [BLAST documentation](https://www.ncbi.nlm.nih.gov/books/NBK279684/).

#### 2. junctions.info

Tab-delimited file with information on every junction found. Each line contains information on one read that contained a circularization junction. Columns are: 
```
chromosomeID read_name	strand	anchor_half	scond_hit_count	overlap_value	read_start1	read_end1*	read_start2*	read_end2	genome_start1	genome_end1*	genome_start2*	genome_end2
```
"1" start and end coordinates are of the anchor of the alignment pair, and "2" are from the second alignment. The highlighted coordinate columns with an "*" are the ones used to determine circRNA coordinates in the subsequent processing steps. 	

#### 3. junctions.bed

Standard [BED](https://en.wikipedia.org/wiki/BED_(file_format)) file in which the score column contains information on the overlap value. When there is no overlap or gap, score = 100; the overlap value is subtracted from 100, and the gap value is added, generating flag values from 97 to 103.

#### 4. ensembles.simple

Tab-delimited file with a simplified list of the circRNA _ensembles_ annotated; each line represents an individual circRNA _ensemble_. Columns are:
```
ensemble_ID	chr	start	end	strand	flag97	flag98	flag99	flag100	flag101	flag102	flag103
```
The last 6 columns contain how many junctions of each type are present in the ensemble (see the "flag" explanation in the [junctions.bed](#3-junctionsbed) section).

#### 5. ensembles.full

Tab-delimited file with detailed information on each read/junction that comprises the circRNA _ensembles_. Each contains information about a single circRNA in the ensemble, including its junction sequence. Columns are:
```
ensemble_ID  chr  start   end  read_ID  overlap_flag  strand  read_sequence  circRNA_length  halfA_sequence  halfB_sequence  empirical_junction_sequence  real_junction_sequence  circRNA_sequence 
```
Read name is prefixed with the ensemble ID to which it belongs, followed by two underscores: `CircEnse_XXXX__`.

Double "vertical bars" (||) show where the junction occurs in junction sequence columns. Lowercase letters in the junction represent a base that was part of an overlap or a gap (see [The Pipeline](#the-pipeline) below). When there are no vertical bars with a lowercase letter in the junction (as the "t" in `CACGGGTCGtAATACCCAA`), the base is from an overlap; when it is between the vertical bars, it is a gap (as the "c" in `CTCCCGCTAC|c|AGCGGGATGG`).

"Empirical junction" is the junction sequence present in the read. "Real junction" is the junction as it is in the cell. These sequences will be different if the sequencing protocol uses dUTP (MonArch default behavior), or will be the same if  the `--dont_invert_strand` option is selected.

#### 6. ensembles.bed

Standard [BED](https://en.wikipedia.org/wiki/BED_(file_format)) file for the final _ensemble_ annotation. The score column contains information on how many reads are in that circRNA *ensemble*. 

# The Pipeline

MonArch can be divided into two main parts: (1) identification of individual circularization junctions in the reads and (2) grouping similar junctions into circRNA _ensembles_ (Figure 1).

![pipeline](https://github.com/user-attachments/assets/880b9e75-3bca-44fe-a3ec-958ec0e8383f)
***Figure 1**: Schematic of the MonArch pipeline. First, it aligns the input RNA-Seq reads in the provided reference genome using BLASTn and searches for reads with a chiastic alignment that came from a circularization junction. MonArch then groups close annotated circRNA junctions into one entity, a circRNA _ensemble_ (CircEnse).*

In the first step, the reads are aligned to the reference genome with BLASTn. Then, `MonArch-GetJunctions.py` searches for a pair of alignments from the same read that could represent a circularization junction to annotate it. Many of the input options affect this step (Figure 2A). The alignments must be uniquely aligned in the genome, have at most `--max_mismatch` mismatches to the reference genome, be no further than `--max_circle_size` bases from each other, and together cover at least `--min_read_cov` of the read. Moreover, the best of the two BLASTn alignments must cover at least half of the read, and the other one should be at least 8 bases long. We allow alignments to have at most a 3nt overlap or gap between them (Figure 2B). An overlap occurs when a base in the circularization junction can be aligned to the reference genome by either alignment of the pair, while a gap occurs when a base in the circularization junction does not align with the reference genome. The coordinates of the circularization junction are adjusted accordingly.

![FigureS1-1](https://github.com/user-attachments/assets/032154cf-1bd4-4ed7-8d0e-b6b09840118a)
***Figure 2**: Details on the MonArch pipeline. (A) Schematic of some MonArch parameters to identify circularization junctions in RNA-Seq reads. Default values are represented. (B) Reads that contain a circularization junction are allowed to have at most a 3nt "overlap" or "gap" between the two halves of the alignment. An overlap (left) occurs when a base (N) can be aligned to either side of the transcript. The final circRNA coordinate always considers that the base came from the 5' portion of the circularized transcript. A gap (right) occurs when there is a base (N) between the two portions of the alignment that does not align to the reference genome.*

The post-processing steps merge circularization junctions in _ensembles_ by `MonArch-PostProcessing.py` and get their sequence with `MonArch-GetJunctionsSeq.py`. Circularization junctions are merged into one ensemble if their start and end coordinates are no further than 3nt from the start and end coordinates of junctions that are already part of that _ensemble_.

# Known issues

## Large genomes/datasets

MonArch was developed to analyze prokaryotic RNA-Seq data and was successfully tested with multiple archaeal datasets and genomes. Although its premise of finding circularization junctions is not specific to prokaryotes, MonArch might be too slow to handle large eukaryotic datasets and genomes. This is because 1) it uses BLAST to align the RNA-Seq reads to the reference genome, and 2) the post-processing scripts are not optimized for very large datasets. This could be partly mitigated by first aligning the raw RNA-Seq reads to the reference genome with another aligner tool (such as HISAT, Bowtie, TopHat) and then using the **non-aligned reads** (where the circularization junction reads should be) as input for MonArch to reduce the input size. This approach was tested with _H. salinarum_ RNA-Seq data and yielded similar results as when the raw RNA-Seq reads were used. Both approaches resulted in the same annotated circRNAs, with the raw reads resulting in more reads supporting each circRNA.

## Paired-end sequencing

MonArch exclusively searches for the circularization junction sequence on the input reads. It does not benefit from paired-end sequencing information when the R1 and R2 reads are each on a different side of the junction and can indicate circularization. When using paired-end sequencing data with MonArch, you can group the R1 and R2 reads in one FASTA file and use it as input. If you wish to _quantify_ circRNAs with paired-end sequencing data, be aware that, in a small circRNA, both R1 and R2 reads could contain the circularization junction sequence and the circRNA will be counted twice in that _ensemble_ by MonArch. The `junctions.bed` and `ensembles.full` output files have information on read names from which the circRNA was annotated that can help mitigate this problem.

# License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)

```
   __o
 _ \<_
(_)/(_)
```
