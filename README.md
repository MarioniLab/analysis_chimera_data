This repository contains the code used to produce the results and figures for the following preprint:


Strauss, Magdalena E., Mai-Linh Nu Ton, Samantha Mason, ..., Berthold Gottgens, John C. Marioni, and Carolina Guibentif. 2023. 
Bespoke single cell molecular and tissue-scale analysis reveals mechanisms underpinning development and disease in complex developing cell populations. 
bioRxiv. https://www.biorxiv.org/content/10.1101/2023.10.11.561904v1


The scripts are organised into folders as follows:


additional_processing_analysis: additional work not covered by the other folders

AML: code to reproduce results on AML patient data (Fig. 5 in preprint)

core_functions_chimera_analysis: R scripts containing core functions to run perturbSuite for the chimera data sets, 

for mapping them to the reference data set and auxiliary functions

mapping_chimeras_to_atlas: scripts for mapping the chimera data to the reference atlas

Mesenchyme_subtyping: subtyping Mesenchyme cell type, detecting subtypes with signature of Juxta-cardiac-field

Mixl1_data_processing: cell calling and normalisation for Mixl1 chimeras

perturbSuite_DA_chimeras: perturbSuite_DA applied to Mixl1 and T chimeras

perturbSuite_DE_dynamic: general application of perturbSuite_DE to chimeras

perturbSuite_general: core functions to apply perturbSuite in general (not chimera specific)

perturbSuite_kinetics_chimeras: application of perturbSuite_kinetics to the chimera data sets

T_chimera_analysis: application of perturbSuite_DE to Limb mesoderm for the T chimeras
