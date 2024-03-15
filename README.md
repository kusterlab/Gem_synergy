# GEM_synergy
## Overview
This repository contains code to analyze and plot the data from *S. Hoefer et al. 2024*. The required input files can be downloaded from Zenodo [(10.5281/zenodo.10792252)](https://zenodo.org/records/10792252). More information on required input files and details on data processing steps are provided in the code. 
## SIMSI-Transfer Input and Output
The MaxQuant results needed as input to SIMSI-Transfer [(kusterlab/SIMSI-Transfer)](https://github.com/kusterlab/SIMSI-Transfer), and the SIMSI-Transfer output tables can be downloaded here: (link will follow). 
## Annotation of Phosphorylation Sites
Before starting data processing of decryptM experiments for CurveCurator [(kusterlab/curve_curator)](https://github.com/kusterlab/curve_curator), phosphorylation site annotations were added to the SIMSI-Transfer output (p10_evidence.txt). To reproduce this, inside python_environment, run the following commands:
```shell
conda env create -f gemsynergy_env.yaml -n gemsynergy
conda activate gemsynergy
python gemsynergy_pSiteAnnotation.py
```
The required files `gemsynergy_env.yaml` and `gemsynergy_pSiteAnnotation.py` are provided in this repository. The file `Phosphosite_seq.fasta` can be downloaded here: (link will follow).
