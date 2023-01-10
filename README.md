# moorea_ati
Around-the-Island 16S analyses for Moorea Coral Reef Network (LTER)

This repository holds the workflow for the Around-the-Island (ATI) microbial sample processing and analyses for 2021 collections ONLY

This workflow is currently in progress. 
Dataframes in this repository are **preliminary** and may be **inaccurate**

File structure is consistent with "input", "output" and "procedure":

"procedure" hosts the annotated scripts for running pipelines and analyses 
"input" hosts any input files necessary (for ex. metadata files)
"output" hosts any output files from the script

Any folder prefixed with "old_" are deprecated and from previous analyses that are not up to date.


The order of operations is as follows:

1. bioinformatics
	- Follow the script in the "procedures" directory. 
	Note: rawdata are too large to include in the github repository. These are currently privately held but will be publicly accessible upon publication
	Note: primary outputs of the QIIME2 pipeline - "ati-2021-demux-paired-end.qza" and "ati-2021-demux-paired-end-trim.qza" - are too large to be held on the github repository. If necessary, please ask repository owner for access.

2. metadata
	- Follow the R script in the "procedures" directory
	Note: metadata from the wider ATI community that is merged with sequencing metadata is subject to change and may update often.

3. core_analyses
	- This is split into the core analyses that have been run according to: 1) pre-processing (including input into R and decontamination), 2) alpha diversity, 3) beta diversity, and any other additional analyses (e.g., map interpolations)
	Note: the pre-processing script MUST be run first, ahead of any other script in this folder. 

