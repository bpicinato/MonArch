# MonArch: an automated pipeline to annotate circRNAs in prokaryotic RNA-Seq data

# Purpose and general information
MonArch was developed to annotate circRNAs from RNA-Seq (.fa files) data using the premise of 

If you use the contents of this repository, please cite:
```
citation
```

# Instalation and requisites

MonArch was developed and tested on Ubuntu distributions.

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
MonArch requires three parameters: a reference genome file in FASTA format; a comma-separated list of FASTA files with the reads; and the path to an output directory (it will be created if it does not exist).

An example of a minimal command to call MonArch is:
```
bash monarch.sh --ref_genome Hsalinarum.fa --reads trat1.fa --output output_dir
```

Below is a list of all possible arguments. This help message can be called in the command line with the `-h` or `--help` arguments.
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

  --ref_genome         The reference genome to be used as subject for blast

  --reads              Reads to be used as query for blast. Or a comma separated
                       list of reads, if applicable

  --output             Path to output dir. Will be created if does not exist

Optional arguments:

 --threads             Number of threads for the analysis (Default = 10)

 --prefix              Prefix for files created by MonArch (Default = circles)

 --max_circle_size     Maximum circle size (Default = 3500)

 --min_sec_hit_size    Minimum size of the second hit (Default = 8)

 --min_read_cov        Minimum read coverage by first and second hits together.
                       Must be  number between 0 and 1 (0 - 100%) (Default = 0.9)

 --max_mismatch        Maximum number of mismatches allowed for each alignment
                       (Default = 0)

 --no_strandness        Do not consider sequencing protocol stranded. MonArch
                       default behavior is to consider sequencing protocol stranded.
                       This paramenter overrides --dont_invert_strand.

 --dont_invert_strand  Do not invert input strand. Default behavior of MonArch
                       is to turn 'plus' strand turns into 'minus' and vice versa.
                       Depending on the sequencing protocol one might not want this.
```

# Outputs
- details of the outputs

# The Pipeline
- figure

# Known issues
MonArch was developed to analyze prokaryotic RNA-Seq data and was tested with multiple archaeal datasets and genomes. Even though its premise of finding circularization junctions is not specific to prokaryotes, MonArch might be too slow to deal with large eukaryotic datasets and genomes. This is becasue 1) it uses Blast to align the RNA-Seq reads to the reference genome; and 2) the post-processing scripts are not optimized to deal with very large datasets. This could be partly mitigated by first alinging the raw RNA-Seq reads to the reference genome with another aligner tool (such as HISAT, Bowtie, TopHat) and then using the **non-aligned reads** (where the circularization junction reads are) as input for MonArch to reduce the input size. This approach was tested with _H. salinarum_ RNA-Seq data and yielded similar results as when the raw RNA-Seq reads were used. Both approaches resulted in the same annotated circRNAs, with the raw reads resulting in more reads supporting each circRNA.

When using 


- cannot effectvely use paired-end sequencing info

# License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)

# Contact information
If you find any bug or error when using MonArch, please use the repository Issues tab to report it.
If you need other information, you can email:
- beatriz [dot] picinato [at] usp [dot] br
- tkoide
- rvencio
- vinicius

```
   __o
 _ \<_
(_)/(_)
```
