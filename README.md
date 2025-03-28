# MonArch: an automated pipeline to annotate circRNAs in prokaryotic RNA-Seq data

If you use the contents of this repository, please cite:
```
citation
```

# Usage
- detais on the command and input formats

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

If you do not want to use Conda to install all the requisites (**not recommended**), below is a list of all the dependencies of MonArch. Be aware that you might need to change how they are called in the main `monarch.sh` script since it assumes that the dependencies were installed with Conda. The versions listed are the ones used in MonArch; other versions of the dependencies might work, but they were not tested.

- Python v=3.10
- Python packages
  - sys 
  - os
  - BioPython
  - datetime 
- Blast v=2.10.1
- GNU Parallel v=20160822

# Outputs
- details of the outputs

# The Pipeline
- figure

# Known issues
- small genomes/data
- cannot effectvely use paired-end sequencing info

# License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
