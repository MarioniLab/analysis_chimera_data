This repository contains the code used to produce the results and figures for the following publication:


Strauss, Magdalena E., Mai-Linh Nu Ton, Samantha Mason, ..., Berthold Gottgens, John C. Marioni, and Carolina Guibentif. 2025. 
A single cell and tissue-scale analysis suite resolves Mixl1's role in heart development. iScience (accepted).


The scripts are organised into folders as follows:


additional_processing_analysis: additional work not covered by the other folders

AML: code to reproduce results on AML patient data 

core_functions_chimera_analysis: R scripts containing core functions to run perturbSuite for the chimera data sets, 

for mapping them to the reference data set and auxiliary functions

DA_analysis:  differential abundance analysis for cell types and lineage trajectories

DE: differential gene expression analysis

dynamic_analysis_lineage_trajectories: COSICC_kinetics analysis to test for delays along lineage trajectories

mapping_chimeras_to_atlas: scripts for mapping the chimera data to the reference atlas

Mesenchyme_subtyping: subtyping Mesenchyme cell type, detecting subtypes with signature of Juxta-cardiac-field

Mixl1_data_processing: cell calling and normalisation for Mixl1 chimeras
