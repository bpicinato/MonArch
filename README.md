# Monarch-Seq

# Table of contents

* [Installation](#installation)
* [Input data](#input-data)
* [The pipeline](#the-pipeline)
* [Output files](#output-files)
* [Downstream analysis](#downstream-analysis)
* [Utils](#utils)
* [Known issues](#known-issues)
* [Acknowledgements](#acknowledgements)
* [Citation](#citation)
* [License](#license)

# Installation

We recommend using [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) to install dependencies of **Monarch-Seq**. 
After installing Conda, you can create the `monarch-seq` environment with:

```sh
# create the env
conda env create -f Monarch-Seq/misc/monarch-seq.yml -n monarch-seq

# activate the env to run Monarch-Seq
conda activate monarch-seq

```

# Input data

`Monarch-Seq.sh` needs the three following inputs:

* **`--ref_genome`**: this parameter should contain the path to a **reference genome** in ```.fasta``` format, which will be used as **subject for blastn**;

* **`--reads`**: this variable should contain the path to a file composed of the **reads** in ```.fasta``` format that will be used as **query for blastn**, among which the **circular RNAs (circRNAs) junctions will be found**. The present pipeline is able to find junctions by aligning reads to a genome because RNA circularization provides permutation of genomic information. In other words, the linear molecules derived from circRNAs do not align to the reference genome, as shown in Figure 1. Therefore, a more effective way to use this variable is to choose a file with **reads that did not align to the genome in past alignments**. It is also recommended that reads generated from **samples treated with RNase R** are used as ```--reads```, because this enzyme degrades many linear molecules. In this case, the query is composed of RNAs derived from circRNAs and other RNAs resistant to digestion [^1].

<p align="center">
   Figure 1
</p>

<p align="center">
  <img width="460" height="300" src="https://user-images.githubusercontent.com/82397631/126039185-92d582bf-bf54-4bf1-8c22-f14f8b991ea8.png">
</p>

<p align="center">
   Modified from [^2]
</p>

* **```--output```**: this parameter should contain the **path to the output directory**. If it does not exist, it will be created.



[1]: [Dobdele et al, 2021](https://www.embopress.org/doi/full/10.15252/embr.202052072?campaign=woletoc)
   
[2]: [Danan et al, 2012](https://pubmed.ncbi.nlm.nih.gov/22140119/)


# The Pipeline

The **pipeline** lies divided in 5 parts,

* In the first stage a reference genome BLAST/n **database is created** through ```makeblastdb -in``` .

* In the second stage for ease of processing and process optimization, **fragments the file**, Splitting input reads\n  ```pyfasta split -n``` .

* Then, on the third stage, through the command ```blastn``` is **implemented the BLAST** in all files generated in the previous step. BLAST is an algorithm used to compare biological sequences, in this case, of halobacterium salinarum.

* The fourth step consists of **Parsing Results\n** , where files are passed to python scripts ```GetCircles.py``` , which results in output files that enable the analysis of best hits, alignments, etc, in addition to providing the coordinates of the circles that have been searched.

* Finally, in the fifth step, occurs the **removal of temporary files\n**, reuniting the output files into a single file again, referring to the one that had been fragmented.


# Output files

* **```circles.bed```** is a file with BED extension (https://en.wikipedia.org/wiki/BED_(file_format)). The **first column** contains the name of the chromosomes or plasmids. The **second and the third columns** contain, respectively, the start and the end coordinates of the circRNAs found by the pipeline. The **fourth column** contains the circRNAs names, that are given based on the order of the found junctions (circrna_0001, circrna_0002 and so forth). The **fifth column** contains the score that was used to represent if the sequences that compose the junctions are from the same or from different DNA strands. If both sequences are from the same strand, score = 100. If they are from different strands, score = 50. Moreover, if there is an overlap of n nucleotides between them, the number n will be subtracted from the previous score value (50 or 100). If there is a gap of m nucleotides between them, the number m will be added to the previous score value (50 or 100). For example, if a junction is composed of two sequences from the same strand and there is an overlap of 2 nucleotides, score = 98. The **sixth column** contains “+” or “-” in each line, which represents, respectively, if the best hit (explicar em algum momento o que é best hit, porque ainda não foi falado) is from the forward or reverse strand.

(Figura para representar de forma mais esquemática a questão de overlap e gap)

# Downstream analysis

Assuming you have RNA-Seq data for RNase R-treated and no-treated libraries, you should run **Monarch-Seq** for all of them. 
However, only the ensembles found in RNase R-treated data should be considered, and the control library should be used only to increase the read count of each ensemble.
We provide the `03-MergeLibraries.sh` script to automatically make it.

Say you have run two libraries of RNase R-treated data and two control libraries with the prefix `treat1`, `treat2`, `ctrl1`, `ctrl2`, respectivelly.
You should call `03-MergeLibraries.sh` like the following:

```sh
./03-MergeLibraries.sh --treated treat1.ensembles.full.readnames,treat2.ensembles.full.readnames \
                       --control ctrl1.ensembles.full.readnames,ctrl2.ensembles.full.readnames \
                       --output final_data
                       
```


RNA-Seq libraries, you might want to only keep ensembles found in the RNase R treated libraries, using the control libraries only to count

# Utils

We provided some scripts that could be useful for make it easier to analyze the output data from Monarch-Seq.

# Known issues

Monarch-Seq was made to deal with reads from single-end sequencing!

# Acknowledgements

The project Monarch-seq  would not have been possible without **labisismi** and all laboratory collaborators, in addition to other assistance involving the Faculty of Medicine of the University of São Paulo.
   
This project use Basic Local Alignment Search Tool - BLATS compare genomic sequences, so  the authors of this algorithm have been very helpful during Monarch-seq's development. It also uses Gaggle Genome Browser - GGB that was developed by the teacher involved in the project, Prof. Dr. Tie Koide to view high-density data plotted against genome coordinates.


# License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
