BINF*6110 - Bioinformatics for Genomics - Assignment 2 

This analysis utilizes splice-aware sequence alignment provided by the STAR tool and the edgeR package in R to identify the differentially expressing genes across different stages of biofilm formation of yeast strains used for wine aging 

The data used for this analysis was derived from the study conducted by Marodanov et al (2020), which includes three triplicate FASTQ samples respective to the three stages of the yeast biofilm stages which inclued:

--Early
--Thin
--Mature 

For the RNA-seq analysis, two DGE lists are generated comparing Early against Mature transcriptomes and Early agaisnt Thin + Mature transcriptomes.
The methods for this analysis were primarily derived from lectures provided by Dr. Lewis Lukens at the University of Guelph.

References: 
1. Mardanov, A. V., Eldarov, M. A., Beletsky, A. V., Tanashchuk, T. N., Kishkovskaya, S. A., & Ravin, N. V. (2020). Transcriptome profile of yeast strain used for biological wine aging revealed dynamic changes of gene expression in course of flor development. Frontiers in Microbiology, 11. https://doi.org/10.3389/fmicb.2020.00538
