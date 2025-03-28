# MonArch: an automated pipeline to annotate circRNAs in prokaryotic RNA-Seq data

# Purpose and general information
MonArch was developed to annotate circRNAs from RNA-Seq (.fa files) data using the premise of 

If you use the contents of this repository, please cite:
```
citation
```

# Instalation and requisites

MonArch was developed and tested on Ubuntu distributions.

We recommend using [Conda](https://docs.conda.io/en/latest/) for installing all MonArch dependencies. We provided a `.yml` file in the `misc` directory with all the dependencies. Navigate to the MonArch directory and create the enviroment with the dependencies with:
```
conda env create -f misc/monarch.yml
```
Before running MonArch, activate the enviroment with:
```
conda activate monarch
```

If you do not want to use Conda to install all the requisites (**not recommended**), below is a list of all the dependencies of MonArch. Be aware that you will need to change how they are called in the main `monarch.sh` script since it assumes that the dependencies were installed with Conda and can be found in the Conda directory. The versions listed are the ones tested and used in MonArch; other versions of the dependencies might work, but they were not tested.

- Python v=3.10.8
- Python packages
  - sys 
  - os
  - BioPython
  - datetime 
- Blast v=2.10.1
- GNU Parallel v=20160822

# Usage
- detais on the command and input formats

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
