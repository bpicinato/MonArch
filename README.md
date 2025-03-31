# MonArch: an automated pipeline to annotate circRNAs in prokaryotic RNA-Seq data

MonArch is an automated pipeline to annotate circRNAs in RNA-Seq reads. It searches for circRNA signatures (the chiastic alignment of the two halves of the read) in the reads and then group close junctions into one circRNA.

If you use the contents of this repository, please cite:

Picinato, B. A., Franceschini-Santos, V. H., Zaramela, L., Vêncio, R. Z. N., Koide, T. (2025) **Archaea express circular isoforms of IS200/IS605-associated ωRNAs**. In press.

# Instalation and requisites

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
MonArch requires three parameters: a **reference genome** file in **FASTA** format; a comma-separated list of **FASTA** files with the **reads** from which circRNAs should be detected; and the path to an output directory (it will be created if it does not exist).

An example of a minimal command to call MonArch is:
```
bash monarch.sh --ref_genome Hsalinarum.fa --reads trat1.fa --output output_dir
```

Below is a list of all possible arguments; please take notice of the default behavior for the optional parameters and check if they fit your analysis. The `-h` or `--help` arguments show this help message on the command line interface.
```
                                   MonArch
         Automated pipeline to annotate circular RNAs in RNA-Seq data
------------------------------------------------------------------------------

Usage:

 monarch.sh --ref_genome <reference.fasta> --reads <reads1.fasta,reads2.fasta>
                --output <path/to/output> [--threads INT]
               [--prefix STR] [--max_circle_size INT]
               [--min_sec_hit_size INT] [--min_read_cov FLOAT]
               [--max_mismatch INT] [--no_strandness] [--dont_invert_strand]

Required arguments:

  --ref_genome         The reference genome to be used as the subject for blast.

  --reads              Reads to be used as the query for blast. Or a comma-separated
                       list of reads, if applicable.

  --output             Path to output directory. It will be created if it does not exist.

Optional arguments:

 --threads             Number of threads for the analysis (Default = 10).

 --prefix              Prefix for files created by MonArch (Default = circles).

 --max_circle_size     Maximum circRNA size allowed (Default = 3500).

 --min_sec_hit_size    Minimum size of the second hit (Default = 8).

 --min_read_cov        Minimum read coverage by first and second hits together.
                       It must be between 0 and 1 (0 - 100%) (Default = 0.9).

 --max_mismatch        Maximum number of mismatches allowed for each alignment
                       (Default = 0).

 --no_strandness       Do not consider the sequencing protocol stranded. MonArch
                       default behavior is to consider the sequencing protocol stranded.
                       This parameter overrides --dont_invert_strand.

 --dont_invert_strand  Do not invert input strand. Default behavior of MonArch
                       is to turn 'plus' strand turns into 'minus' and vice versa.
                       Depending on the sequencing protocol one might not want this.
```

# Outputs
MonArch outputs can be divided into three categories: alignment results (Blastn results), junction results (regarding circularization junctions before joining into ensembles), and ensemble results (final results after ensemble formation). See [The Pipeline](#the-pipeline) section. The most useful files for downstream analysis are the .bed.

- **blastn.tbl**: Blastn standard tab-delimited output compressed by gzip. Columns are: sseqid, qseqid, qlen, length, qstart, qend, sstart, send, sstrand, mismatch. You can learn more about it with the [Blastn documentation](https://www.ncbi.nlm.nih.gov/books/NBK279684/), section "Options for the command-line applications".

- **junctions.info**: Tab-delimited file with information on every junction found. Each line contains information on one read that contained a circularization junction. Columns are: chromosome ID, read name, strand, anchor alignment half (A/B, see [The Pipeline](#the-pipeline)), how many "second hits" were found, overlap value (values -3 to 3), start of the anchor alignment on the read, _end of the anchor alignment on the read_, _start of the second alingment on the read_, end of the second alingment on the read, start of the anchor alignment on the genome, _end of the anchor alignment on the genome_, _start of the second alingment on the genome_, end of the second alingment on the genome. The highlighted coordinate columns are the ones used to determine circRNA coordinates. 	

- **junctions.bed**: Standard [BED](https://en.wikipedia.org/wiki/BED_(file_format)) file in which the score column contains information on the overlap value (VALUE MUST BE 0-100 - see what to do, if we want to explain our 100 flag system or not).

- **ensembles.simple**: 

- **ensembles.full**: 

- **ensembles.bed**: 

# The Pipeline
- figure

# Known issues
## Large genomes/datasets
MonArch was developed to analyze prokaryotic RNA-Seq data and was successfully tested with multiple archaeal datasets and genomes. Even though its premise of finding circularization junctions is not specific to prokaryotes, MonArch might be too slow to deal with large eukaryotic datasets and genomes. This is because 1) it uses Blast to align the RNA-Seq reads to the reference genome; and 2) the post-processing scripts are not optimized for very large datasets. This could be partly mitigated by first aligning the raw RNA-Seq reads to the reference genome with another aligner tool (such as HISAT, Bowtie, TopHat) and then using the **non-aligned reads** (where the circularization junction reads should be) as input for MonArch to reduce the input size. This approach was tested with _H. salinarum_ RNA-Seq data and yielded similar results as when the raw RNA-Seq reads were used. Both approaches resulted in the same annotated circRNAs, with the raw reads resulting in more reads supporting each circRNA.

## Paired-end sequencing
MonArch exclusively searches for the circularization junction sequence on the input reads. It does not benefit from paired-end sequencing information when the R1 and R2 reads are each on a different side of the junction and can indicate circularization. When using paired-end sequencing data with MonArch, you can group the R1 and R2 reads in one FASTA file and use it as input. If you wish to _quantify_ circRNAs with paired-end sequencing data, be aware that, in a small circRNA, both R1 and R2 reads could contain the circularization junction sequence and the circRNA will be counted twice. The X and Y outputs possuem informação de quais reads x junções vieram. **should R2 reads be inverted?**

# License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)

```
   __o
 _ \<_
(_)/(_)
```
