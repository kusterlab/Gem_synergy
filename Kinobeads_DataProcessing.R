#############################################
# 1. Before CurveCurator: Processing MaxQuant output and calculating correction factors
#############################################

# To run this step, please download the file "01_Kinobeads_proteingroups.txt" from Zenodo (see README). 
# This file should be imported as a dataframes called "Kinobeads_proteingroups".

#############################################

# Filter out contaminants, reverse hits and proteins only identified by site from proteingroups of all four experiments
Kinobeads_proteingroups <- Kinobeads_proteingroups %>% dplyr::filter(Potential.contaminant != "+" & Reverse != "+" & Only.identified.by.site != "+")

# Fill up the (few) empty Gene Names with values from the Majority Protein IDs-column
Kinobeads_proteingroups$Gene.names <- ifelse(Kinobeads_proteingroups$Gene.names == "", Kinobeads_proteingroups$Majority.protein.IDs, Kinobeads_proteingroups$Gene.names)

# Rename relevant columns to fit the .toml-file for curve fitting in CurveCurator
# e.g. LFQ intensity columns are renamed to "Raw"
colnames(Kinobeads_proteingroups) <- gsub("LFQ.intensity.","Raw ", colnames(Kinobeads_proteingroups))
colnames(Kinobeads_proteingroups) <- gsub("DMSO","0", colnames(Kinobeads_proteingroups))
colnames(Kinobeads_proteingroups) <- gsub("Protein.IDs","Proteins", colnames(Kinobeads_proteingroups))
colnames(Kinobeads_proteingroups) <- gsub("Gene.names","Genes", colnames(Kinobeads_proteingroups))

# Remove proteins that were not quantified in the vehicle experiment
Kinobeads_proteingroups <- Kinobeads_proteingroups %>% dplyr::filter(`Raw 0` > 0) 
experiments <- unique(Kinobeads_proteingroups$Experiment)

# Export this and perform curve fitting in CurveCurator (https://github.com/kusterlab/curve_curator)
# Afterwards, read in output tables of all four experiments, merge into one big table
# Alternatively, the output of CurveCurator can be found as "02_Kinobeads_curves.txt" in Zenodo (see README).

# Select relevant columns from "Kinobeads_proteingroups" and calculate correction factors for each experiment, for Kd,app calculation (in step 2)
correctionfactor_df <- Kinobeads_proteingroups  %>% dplyr::select(Experiment, Genes, "Raw 0", "Raw PDPD", "Unique.peptides.0") 
correctionfactor_df$correction <- correctionfactor_df$`Raw PDPD`/correctionfactor_df$`Raw 0`
correctionfactor_df$correction <- ifelse(correctionfactor_df$correction == 0, NA, correctionfactor_df$correction)

# Per protein, calculate the mean correction factor from all four experiments. Add columns to original table.
corr_mean_df <- correctionfactor_df %>% dplyr::group_by(Genes) %>% dplyr::mutate(correction_mean = mean(correction, na.rm = T)) %>% dplyr::select(Genes, correction_mean) %>% as.data.frame
correctionfactor_df <- correctionfactor_df %>% left_join(corr_mean_df, by = 'Genes') %>% distinct

# Export Kinobeads_proteingroups for each experiment and perform curve fitting in CurveCurator (https://github.com/kusterlab/curve_curator)
# Afterwards, read in output tables of all four experiments, and merge into one big table.
# Merge these correction factors to the CurveCurator output (already done for the next step).





#############################################
# 2. After CurveCurator: Applying filter criteria
#############################################

# To run this step, please download the file "02_Kinobeads_curves.txt" from Zenodo (see README). 
# This file should be imported as a dataframe called "Kinobeads_curves".
# It already contains the correction factors calculated in Step 1.

#############################################

# Calculate apparent dissociation constants (Kd,app) by multiplying median correction factors with EC50s
Kinobeads_curves$EC50 <- 10^(-Kinobeads_curves$pEC50)*1000000000
Kinobeads_curves$Kdapp <- Kinobeads_curves$EC50 * Kinobeads_curves$correction_mean
Kinobeads_curves$pKdapp <- -log10(Kinobeads_curves$Kdapp/1000000000)

# Apply filter criteria
Kinobeads_curves$Target_CurveFit <- ifelse(Kinobeads_curves$Curve.R2 >= 0.7 & 
                                                          Kinobeads_curves$pEC50 > 6 & 
                                                          Kinobeads_curves$Unique.peptides.0 > 3 & 
                                                          Kinobeads_curves$Curve.Fold.Change <(-1) &
                                                          Kinobeads_curves$Curve.Slope >0.2, T,F) 

# After exporting Kinobeads_curves, this table was further annotated manually. 
# Targets were curated based on inspection of spectra and unique peptide counts.
# The resulting file can be found as in Zenodo as Kinobeads_plotting.txt (see README).
# This will be the input file for Kinobeads_plotting.R