# Arctic_eDNA_2021

This repository contains scripts used in the analysis of Wang et al. 2021 Late Quaternary dynamics of Arctic biota from ancient environmental genomics. https://doi.org/10.1038/s41586-021-04016-x


Script | used in which section of the paper | function
--- | --- | --- 
bam_header_subset.R | Supplementary Information section 9.2.4 | to subset headers of bam file to include only reference headers appeared in the alignments
Plant_normalization.R | Supplementary Information section 9.5.2 | plant reads normalisation 
paleoDistributionModel.R | Supplementary Information section 12.2 | ancient human distribution niche modelling
ExtractAndMapModels.R | Supplementary Information section 12.2 | mapping human niche modelling results to each eDNA site
get_consensus.py | Supplementary Information section 14.1.2 | to creat a consensus sequence from the multiple sequence alignment
Animal_spatiotenporal_modelling.R | Supplementary Information section 13.3 | to model covariates of animal eDNA presence/absence 
INLA_ST_functions.R | Supplementary Information section 13.3 | functions used in Animal_spatiotenporal_modelling.R 
spde-book-functions.R  | Supplementary Information section 13.3 | functions used in Animal_spatiotenporal_modelling.R 
