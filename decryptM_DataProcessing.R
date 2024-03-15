#############################################
# 1. Before CurveCurator: Filtering for unambigiously assigned p-sites with localization probability >0.75
#############################################

# To run this step, please download the files "01_p10_evidence_annotated.txt" and "01_p10_msms.txt" from Zenodo (see README). 
# These files should be imported as dataframes called "evid_SIMSI" and "msms_SIMSI", respectively.
#
# Please note: This step includes long run time (>1h). You can also start from later steps using the intermediate input files provided at Zenodo (see README).
# 
# Before this step, SIMSI-Transfer was run (https://github.com/kusterlab/SIMSI-Transfer), using the evidence.txt output tables of the MaxQuant search.
# The data was further processed by adding P-Site annotations using the Python Package.
# The python code for P-site annotation can be found in this repository.

#############################################

# Annotate and filter evid_SIMSI for phosphopeptides and remove potential contaminants and reverse hits
evid_SIMSI$'Phospho (STY)' <- str_count(evid_SIMSI$`Modified sequence`,"Phospho \\(STY\\)")
evid_SIMSI <- evid_SIMSI %>% dplyr::filter(`Phospho (STY)` > 0 & `Potential contaminant` != '+' & Reverse != '+')

# Change column names to make it compatible with CurveCurator, and fix error in drug name
colnames(evid_SIMSI)[colnames(evid_SIMSI) == "Gene Names"] <- "Gene names"
colnames(evid_SIMSI)[colnames(evid_SIMSI) == "Protein Names"] <- "Protein names"
evid_SIMSI$Experiment[which(evid_SIMSI$Experiment == 'AsPC1_GEM_Ceralaserib_ddPTM')] <- 'AsPC1_GEM_Ceralasertib_ddPTM'

# Annotate and filter evid_SIMSI for phosphopeptides and remove reverse hits
msms_SIMSI$'Phospho (STY)' <- str_count(msms_SIMSI$`Modified sequence`,"Phospho \\(STY\\)")
msms_SIMSI <- msms_SIMSI %>% dplyr::filter(`Phospho (STY)` > 0 & Reverse != '+')

# Extract number of channels used in the analysis
channel_numbers <- 1:sum(grepl("corrected", colnames(evid_SIMSI)))

# Get unique peptide sequences with site localization probabilities from the following table
# Therefore, extract all reported site localization probabilities (values in brackets from column `Phospho (STY) Probabilities`)
# Put into list, then remove brackets, make numeric and sort from high to low
STY_Probabilities <- msms_SIMSI %>% dplyr::select(`Phospho (STY) Probabilities`, `Modified sequence`, Experiment, Charge, `Phospho (STY)`) %>% distinct
vector_siteprob <- STY_Probabilities$`Phospho (STY) Probabilities`
list_siteprob <- regmatches(vector_siteprob, gregexpr("\\(.*?\\)", vector_siteprob)) # extract the probabilities of each site per peptide
list_siteprob <- lapply(X = list_siteprob, FUN = function(t) gsub("[\\(\\)]", "", x = t)) # remove brackets
list_siteprob <- lapply(X = list_siteprob, FUN = function(t) as.numeric(t)) # make numeric
list_siteprob <- lapply(X = list_siteprob, FUN = function(t) sort(t, decreasing = T)) # sort probabilities per peptide from high to low

# Create vector of same length with site count
vector_sitecount <- STY_Probabilities$`Phospho (STY)`

# Loop over list: For each phosphorylated peptide,look if any of the 'assigned' sites have probability greater 0.25, 0.5, or 0.75
# Also, check if any of the 'assigned' sites have identical prob. with 'non-assigned' ones (i.e. are ambiguously localized) 
# (only when multiple possible sites were quantified)
for(i in 1:length(vector_sitecount)){
  
  n_reported_sites <- vector_sitecount[i] # count of assigned sites
  all_sites_prob <- list_siteprob[[i]] # count of assigned and non-assigned sites (i.e. all sites with prob. values)
  reported_sites_prob <- all_sites_prob[1:n_reported_sites] # only the prob. of the assigned sites (works because values are sorted)
  
  if(vector_sitecount[i] > 0){ # Check if the probability of any of the reported sites is below the threshold (will filter for >0.75 later)
    STY_Probabilities$Prob_Lower_0.25[i] <- any(reported_sites_prob < 0.25)
    STY_Probabilities$Prob_Lower_0.5[i] <- any(reported_sites_prob < 0.5)
    STY_Probabilities$Prob_Lower_0.75[i] <- any(reported_sites_prob < 0.75)
    
  } else {
    STY_Probabilities$Prob_Lower_0.25[i] <- NA
    STY_Probabilities$Prob_Lower_0.5[i] <- NA
    STY_Probabilities$Prob_Lower_0.75[i]<- NA
  }
  
  if(length(all_sites_prob) > n_reported_sites){
    other_sites_prob <- all_sites_prob[(n_reported_sites+1):length(all_sites_prob)] # extract other sites that are reported but not assigned
    STY_Probabilities$Ambig_Localization[i] <- any(other_sites_prob %in% reported_sites_prob) # Check if any of them has same probability as the assigned ones
  } else {
    STY_Probabilities$Ambig_Localization[i] <- NA
  }
}

# merge the new Localization Filter columns to msms_SIMSI
STY_Probabilities <- merge(msms_SIMSI, STY_Probabilities %>% 
                             dplyr::select(`Phospho (STY) Probabilities`,"Prob_Lower_0.25", "Prob_Lower_0.5", "Prob_Lower_0.75","Ambig_Localization") %>% 
                             distinct, by = 'Phospho (STY) Probabilities', all.x = T)

# Extract which peptides to exclude from analysis (i.e. localization probability < 0.75 & ambiguously localized)
peptides_to_exclude <- STY_Probabilities %>% dplyr::filter(Prob_Lower_0.75 == T | Ambig_Localization == T) 

# Now filter out these peptides from p10_evidence.txt using the 'summary_ID' from p10_msms.txt
# Create function to filter dataframe based on pattern
filter_dataframe <- function(data, pattern) {
  data %>%
    filter(!str_detect(summary_ID, pattern))
}

# Initialize the filtered dataframe
evid_SIMSI_filtered <- evid_SIMSI

# Make sure the summary_IDs will be correctly matched (before, between or after semicolon, and not within another number)
peptides_to_exclude <- paste0("(^|[^0-9])", peptides_to_exclude, "([^0-9]|$)")

# Loop through each pattern and filter the dataframe
for (pattern in peptides_to_exclude) {
  evid_SIMSI_filtered <- filter_dataframe(evid_SIMSI_filtered, pattern)
}





#############################################
# 2. Before CurveCurator: Summing up reporter intensities of p-peptides containing the same site(s)
#############################################

# To run this step, please download the file "02_p10_evidence_filtered.txt" from Zenodo (see README). 
# This file should be imported as a as a dataframe called "evid_SIMSI_filtered".
#
# Please note: This step includes long run time (>1h). You can also start from later steps using the intermediate input files provided at Zenodo (see README).

#############################################

# Replace Phospho-annotation with "(ph)", remove all other modifications from Modified Sequence column
evid_SIMSI_filtered$`Modified sequence` <- gsub("\\(Phospho \\(STY\\)\\)","\\(ph\\)", evid_SIMSI_filtered$`Modified sequence`)
evid_SIMSI_filtered$`Modified sequence` <- gsub("\\(Oxidation \\(M\\)\\)","", evid_SIMSI_filtered$`Modified sequence`)
evid_SIMSI_filtered$`Modified sequence` <- gsub("\\(Acetyl \\(Protein N-term\\)\\)","", evid_SIMSI_filtered$`Modified sequence`)
evid_SIMSI_filtered$`Modified sequence` <- gsub("_","", evid_SIMSI_filtered$`Modified sequence`)

# Fill empty rows if site positions column with modified sequence information
evid_SIMSI_filtered$`Site positions` <- ifelse(evid_SIMSI_filtered$`Site positions` == "", evid_SIMSI_filtered$`Modified sequence`, evid_SIMSI_filtered$`Site positions`)

# Create a column in evid_SIMSI_filtered called "canonical" in which in the end, only one uniprot ID is listed (ideally canonical)
# These IDs will be used to extract site position information for only one isoform (non-redundant info)
canonical <- strsplit(evid_SIMSI_filtered$Proteins, ";") 
canonical <- lapply(canonical, function(x){
  if(length(x) >1){
    x <- (x[!grepl("-", x)])[1]
  } else {
    x <- x
  }})
evid_SIMSI_filtered$canonical <- unlist(canonical)

# Create temporary Uniprot ID column and Gene name column which has as little missing values as possible
# use Proteins_temp to fill empty information in the Gene names column
evid_SIMSI_filtered$Proteins_temp <- sapply(strsplit(evid_SIMSI_filtered$`Site positions`, "_"), `[`, 1)
evid_SIMSI_filtered$Proteins_temp <- ifelse(is.na(evid_SIMSI_filtered$Proteins_temp), sapply(strsplit(evid_SIMSI_filtered$`Leading proteins`, ";"), `[`, 1), evid_SIMSI_filtered$Proteins_temp)
evid_SIMSI_filtered$Proteins_temp <- ifelse(is.na(evid_SIMSI_filtered$Proteins_temp), sapply(strsplit(evid_SIMSI_filtered$`Proteins`, ";"), `[`, 1), evid_SIMSI_filtered$Proteins_temp)
evid_SIMSI_filtered$Genes_temp <- evid_SIMSI_filtered$`Gene names`
evid_SIMSI_filtered$Genes_temp <- ifelse(evid_SIMSI_filtered$Genes_temp == "", evid_SIMSI_filtered$Proteins_temp , evid_SIMSI_filtered$Genes_temp)

evid_SIMSI_filtered$canonical[which(is.na(evid_SIMSI_filtered$canonical))] <- evid_SIMSI_filtered$Proteins_temp[which(is.na(evid_SIMSI_filtered$canonical))]

# Make column with short position names. First, split short position names into list
evid_SIMSI_filtered$SitePositionsShort <- evid_SIMSI_filtered$`Site positions`
sitepositionsshort <- strsplit(evid_SIMSI_filtered$SitePositionsShort, ";") # Strplit the individual Site positions (creates a list)
names(sitepositionsshort) <- evid_SIMSI_filtered$canonical # name with one particular isoform of protein, needed for site info extraction

# Loop through list to get site information only
for(k in 1:length(sitepositionsshort)){
  if(grepl("_", sitepositionsshort[[k]][1])){
    sitepositionsshort[[k]] <- sitepositionsshort[[k]][grepl(paste0(names(sitepositionsshort)[k], "_"), sitepositionsshort[[k]], fixed = T)]   # get site information only for one isoform of the protein, not multiple (works since only one is given in name of list elements)
    sitepositionsshort[[k]] <- gsub(".*_", "", sitepositionsshort[[k]]) # Remove protein accession (everything before "_") and keep only the site information
    sitepositionsshort[[k]] <- paste(sitepositionsshort[[k]], collapse = ";") # Concatenate multiple sites with ";"
  }
  else{
    sitepositionsshort[[k]] <- sitepositionsshort[[k]]
  }
}

# Make some additional columns for site annotation
evid_SIMSI_filtered$SitePositionsShort <- unlist(sitepositionsshort)
evid_SIMSI_filtered$ModSeqShort <- evid_SIMSI_filtered$`Modified sequence`
evid_SIMSI_filtered$SitePos_ModSeq <- paste(evid_SIMSI_filtered$SitePositionsShort, evid_SIMSI_filtered$ModSeqShort, sep = " ")

# Make additional column with site information: Gene name + reported sites. Will be used for concatenating/summing up peptide intensities
evid_SIMSI_filtered$Gene_SitePos <- paste(evid_SIMSI_filtered$Genes_temp, evid_SIMSI_filtered$SitePositionsShort, sep= '_') 

# Before summing up intensities, concatenate summary_ID (to later find probabilities in msms.txt), mod. seq. and site pos
evid_SIMSI_filtered_temp <- evid_SIMSI_filtered %>% 
  dplyr::select(`Site positions`, `Modified sequence`, summary_ID, Experiment, Gene_SitePos) %>% 
  group_by(`Gene_SitePos`, Experiment) %>%
  mutate(Prob_concat = paste0(summary_ID, collapse = ";")) %>% 
  mutate(ModSeq_concat = paste0(`Modified sequence`, collapse = ";")) %>%
  mutate(SitePos_concat = paste0(`Site positions`, collapse = ";")) %>%
  dplyr::select(-`Modified sequence`, -summary_ID, -`Site positions`) %>% distinct
evid_SIMSI_filtered <- merge(evid_SIMSI_filtered, evid_SIMSI_filtered_temp, by= c('Gene_SitePos', "Experiment"))

# Sum up intensities of same "Gene_SitePos" per Experiment, merge back to previous evidence
evid_SIMSI_filtered_mergedInt <- evid_SIMSI_filtered %>% group_by(`Gene_SitePos`, Experiment) %>%
  dplyr::summarise_at(vars(paste0("Reporter intensity corrected ", min(channel_numbers)):paste0("Reporter intensity corrected ", max(channel_numbers))), sum)
colnames(evid_SIMSI_filtered_mergedInt)[-c(1,2)] <- paste0("Merged ", colnames(evid_SIMSI_filtered_mergedInt)[-c(1,2)])
evid_SIMSI_filtered <- merge(evid_SIMSI_filtered, evid_SIMSI_filtered_mergedInt, by = c('Gene_SitePos', "Experiment"))





#############################################
# 3. Before CurveCurator: Normalization across TMT channels and export of CurveCurator Input files
#############################################

# To run this step, please download the files "03_p10_evidence_filtered_sum.txt" from Zenodo (see README). 
# This file should be imported as a as a dataframe called "evid_SIMSI_filtered".

#############################################

# Rename Merged columns to original columns, and use unique rows only
evid_fornorm <- evid_SIMSI_filtered %>% dplyr::select(-starts_with("Reporter intensity corrected"), -`Modified sequence`, -`Site positions`) %>%
  rename('Modified sequence'= "ModSeq_concat",  'Phospho (STY) Probabilities'="Prob_concat", 'Site positions'= "SitePos_concat") %>% distinct
colnames(evid_fornorm) <- gsub("Merged ", "", colnames(evid_fornorm))

# Extract channel numbers from table
channel_numbers <- 1:sum(grepl("corrected", colnames(evid_fornorm)))

# For each experiment, normalize across TMT Channels. Write tables into list.
drugs <- c("Elimusertib", "Ceralasertib", "Berzosertib", "Gartisertib")
evid_norm_list <- list()
CurveCurator_input_list <- list()

for(i in drugs){
  
  # Create matrix with the merged intensity columns for each TMT Channel
  # Save site information in rownames to merge back to original table later 
  evid_fornorm_temp <- evid_fornorm %>% dplyr::filter(Experiment == paste0("AsPC1_GEM_",i,"_ddPTM")) %>%  
    dplyr::distinct(`Gene_SitePos`, Experiment, .keep_all= TRUE) %>% 
    dplyr::select(`Gene_SitePos`, starts_with("Reporter intensity corrected"))
  rownames(evid_fornorm_temp) <- evid_fornorm_temp$`Gene_SitePos`
  evid_fornorm_temp <- evid_fornorm_temp %>% dplyr::select(-`Gene_SitePos`)
  colnames(evid_fornorm_temp) <- as.character(channel_numbers)
  
  # Remove rows containing any missing value
  evid_fornorm_temp$zerocount <- rowSums(evid_fornorm_temp == 0) # Add column with number of missing values per row
  evid_fornorm_temp <- evid_fornorm_temp %>% dplyr::filter(zerocount == 0) %>% dplyr::select(-zerocount) # filter and remove column zerocount

  # Calculate correction factors based merged intensities
  colMed <- evid_fornorm_temp %>% dplyr::summarise_all(median) %>% as.numeric
  colMedMedian <- median(colMed)
  colMedCFs <- colMedMedian/colMed %>% as.numeric
  
  # Apply on data matrix
  evid_norm_temp <- sweep(evid_fornorm_temp, 2, colMedCFs, FUN = "*")
  colnames(evid_norm_temp) <- paste0("Reporter intensity corrected ", names(evid_norm_temp))

  # Add back the row names as Site annotation column, merge back to original table based on this column
  # Thereby, replace the "Reporter intensity corrected" column, as this naming is needed for CurveCurator later
  evid_norm_temp$`Gene_SitePos` <- rownames(evid_norm_temp)
  evid_norm <- merge(evid_fornorm %>% dplyr::filter(Experiment == paste0("AsPC1_GEM_",i,"_ddPTM")) %>% 
                            dplyr::select(-starts_with("Reporter intensity corrected")), 
                     evid_norm_temp, by="Gene_SitePos", all =T )
  
  # Calculate median scores and PEPs per annotated p-site, and filter for missing values in vehicle
  evid_norm <- evid_norm %>% group_by(`Gene_SitePos`) %>% 
    mutate(Score = median(Score, na.rm =T), `Delta score` = median(`Delta score`, na.rm =T), PEP = median(PEP, na.rm =T)) %>%
    dplyr::filter(!is.na(`Reporter intensity corrected 1`))
  
  # Remove duplicated Site position rows, since we are only interested in merged intensities now
  evid_norm <- evid_norm[!duplicated(evid_norm$`Gene_SitePos`), ] 
  
  # These tables are the input for CurveCurator
  CurveCurator_input_list[[which(drugs==i)]] <- evid_norm %>% dplyr::select(c(Experiment, `Modified sequence`, `Gene names`, `Protein names`,
                                                                              Score, Proteins, starts_with("Reporter intensity corrected"), `Potential contaminant`))
  # These tables will be used to merge back to CurveCurator-Output
  evid_norm_list[[which(drugs==i)]] <- evid_norm

}

# The dataframes in CurveCurator_input_list are the evidence tables used as input for CurveCurator (https://github.com/kusterlab/curve_curator/).
# All input files for CurveCurator (processed evidence tables and .toml files) can be downloaded in Zenodo (see README).
# Please note: The name of the exported dataframes must match the names indicated in the .toml when using CurveCurator.





#############################################
# 4. After CurveCurator: Add back columns from evidence file & prepare for t-test (GEM vs. vehicle)
#############################################

# To run this step, please download the files "04_curves.txt" and "04_p10_evidence_filtered_sum_norm.txt"from Zenodo (see README).
# These files should be imported as dataframes called "curves" and "evid_norm", respectively.

#############################################

drugs <- c("Elimusertib", "Ceralasertib", "Berzosertib", "Gartisertib")

curves_list <- list()
for(i in drugs){
  evid_norm_temp <- evid_norm %>% dplyr::filter(Experiment==i)
  # Channel 11 is the Vehicle which was not used in curve fitting
  # The normalized channel is needed for t-test of GEM-only later.
  colnames(evid_norm_temp) <- gsub("Reporter intensity corrected 11", "Normalized 11", colnames(evid_norm_temp))
  evid_norm_temp <- evid_norm_temp %>%
    dplyr::select(c("Gene_SitePos", "Modified sequence", starts_with("Normalized"), "PEP", "Score", "Delta score", "Gene names", "Protein names", "Proteins", #starts_with("Raw"),
                    "Site positions","Start positions","End positions","Site sequence context", "SitePositionsShort",
                    "ModSeqShort", "Transferred spectra count", "Proteins_temp", "Genes_temp")) %>% distinct
  
  curves_temp <- curves %>% dplyr::filter(Experiment==i)
  # The "Raw" intensities in the CurveCurator were normalized beforehand (see step 3). Rename to "Normalized"
  colnames(curves_temp) <- gsub("Raw", "Normalized", colnames(curves_temp)) 
  # Select columns, merge with evidence file from before CurveCurator, write into list
  curves_temp <- curves_temp %>% dplyr::select(!c(Genes, Proteins, Name, Score)) # Will use the average Score columns from evid_norm
  curves_temp <- curves_temp  %>% full_join(evid_norm_temp, by=c("Modified sequence")) 
  curves_list[[which(drugs==i)]] <- curves_temp
}

curves <- bind_rows(curves_list)

# Add site sequence window
curves$`Site sequence context_15AA`<- sapply(strsplit( curves$`Site sequence context`, ";"), "[", 1)
curves$`Site sequence context_15AA` <- substring(curves$`Site sequence context_15AA`, 9, 23)

# In column ModSeqShort, just write the first sequence
curves$ModSeqShort <- sapply(strsplit( curves$`Modified sequence`, ";"), `[`, 1)
# Add Phospho Count
curves$'Phospho (STY)' <- str_count(sapply(strsplit(curves$`Modified sequence`, ";"), `[`, 1),"(ph)")
# Add SQTQ motif filter column
curves$SQTQ <- grepl(pattern = "S(ph)Q", x = curves$ModSeqShort, fixed = T)|grepl(pattern = "T(ph)Q", x = curves$ModSeqShort, fixed = T)

# Make columns with gene name, site name and modified sequence
curves$SitePos_ModSeq <- paste(curves$SitePositionsShort, curves$ModSeqShort, sep = " ")
curves$Gene_SitePos <- paste(curves$Genes_temp, curves$SitePositionsShort, sep= '_')
curves$Gene_SitePos_ModSeq <- paste(curves$Gene_SitePos, curves$ModSeqShort, sep = " ")





# The following code generates the input table for t-test analysis in Perseus of 1uM GEM vs. vehicle
# The t-test was two-sided and used BH correction. 
# It was performed twice on the same data: first with at least 3 valid values in each group, second with 1 valid value in each group
# The second t-test was used to obtain log2 fold changes upon 1uM GEM for each processed peptide, even for 2 or more missing values.

# Channel 11 is the vehicle not containing any drug; Channel 10 is 1uM GEM (without ATRi)
# The median of medians (between GEM and control in each experiment) was used to normalize across the four TMT batches 

GEM_ctrl_df <- curves %>% dplyr::select(`Gene_SitePos`, `Normalized 11`, `Normalized 10`, Experiment)
colnames(GEM_ctrl_df)[c(2,3)] <- c("ctrl", "GEM")
GEM_ctrl_df <- GEM_ctrl_df %>% rowwise %>% mutate(median = median(c(ctrl, GEM)))
GEM_ctrl_df_temp <- GEM_ctrl_df %>% group_by(`Gene_SitePos`) %>% summarise(medianofmedian = median(median))

GEM_ctrl_df <- merge(GEM_ctrl_df, GEM_ctrl_df_temp, by = "Gene_SitePos")
GEM_ctrl_df <- GEM_ctrl_df %>% mutate(corrfactor = median/medianofmedian)
GEM_ctrl_df <- GEM_ctrl_df %>% mutate(ctrl_corr = ctrl/corrfactor, GEM_corr = GEM/corrfactor)  %>% rowwise %>%  mutate(median_corr = median(c(ctrl_corr, GEM_corr)))

GEM_ctrl_df <- GEM_ctrl_df[!duplicated(paste0(GEM_ctrl_df$`Gene_SitePos`, GEM_ctrl_df$Experiment)),]

GEM_ctrl_df <- GEM_ctrl_df %>% dplyr::select(`Gene_SitePos`, Experiment, GEM_corr, ctrl_corr) %>%
  tidyr::pivot_wider(names_from = Experiment, values_from = c(GEM_corr, ctrl_corr))

# After Perseus, the data was exported and the results were merged back to the curves table based on Gene_SitePos column
# The newly added columns are called c("GEM_pvalue", "GEM_qvalue", "GEM_log2FC", "GEM_statistics")





#############################################
# 5. After CurveCurator and GEM t-test: Filter for regulation, add annotations
#############################################

# To run this step, please download the files "05_curves_GEM.txt" and "04_p10_evidence_filtered_sum_norm.txt"from Zenodo (see README).
# This file should be imported as as a dataframe called "curves".

#############################################

# Annotate significantly regulated peptides upon GEM-only, and add direction of regulation (up/down)
curves$GEM_sign <- curves$GEM_qvalue < 0.01 & abs(curves$`GEM_log2FC`) > 1
curves$GEM_updown <- ""
curves$GEM_updown[curves$GEM_qvalue < 0.01 & (curves$`GEM_log2FC`) > 1] <- 'up'
curves$GEM_updown[curves$GEM_qvalue < 0.01 & (curves$`GEM_log2FC`) < (-1)] <- 'down'

# Save the original curve classification result by CurveCurator in separate column
# Use "Curve Regulation" column, apply special filter criteria for curves: pEC50 >5
curves$`Curve Regulation_Curator` <- curves$`Curve Regulation`
curves$`Curve Regulation`[curves$pEC50 < 5] <- ""
curves$`Curve Regulation`[curves$`Curve Regulation` == "not"] <- "" # ignore curves classified as "not-regulated" for this study
curves$EC50 <- 10^(-curves$pEC50)*1000000000 # Add EC50

# Add categorical column describing if curve regulation is true
curves$ATRi_regulated <- curves$`Curve Regulation` == 'up' |  curves$`Curve Regulation` == 'down'

# Count how many ATRi a peptide was regulated by adding column "ATRi_count_regulated"
curves <- curves  %>% dplyr::group_by(`Gene_SitePos`) %>% dplyr::add_tally(`ATRi_regulated`, name = 'ATRi_count_regulated') %>% ungroup 
# Add categorical column that tells if peptide was regulated in one direction by at least 3 ATRi ("ATRi_updown_atleast3")
curves <- curves %>% group_by(Gene_SitePos) %>% dplyr::mutate(ATRi_updown_atleast3 = if_else(any(`Curve Regulation` == 'up') & any(ATRi_count_regulated >2), "up", 
                                                                                          if_else(any(`Curve Regulation` == 'down') & any(ATRi_count_regulated >2), "down", ""), ""))


# Add column "combi_sign" which describes the regulation of both GEM and at least 3 ATRi (without direction)
curves$combi_sign <- ifelse(curves$GEM_sign == T & curves$ATRi_count_regulated >2, "GEMreg_ATRireg",
                            ifelse(curves$GEM_sign == F & curves$ATRi_count_regulated >2, "GEMnot_ATRireg",
                                   ifelse( is.na(curves$GEM_sign) & curves$ATRi_count_regulated >2, "GEMNA_ATRireg",
                                           ifelse( is.na(curves$GEM_sign) & curves$ATRi_count_regulated <3, "GEMNA_ATRinot",
                                                   ifelse(curves$GEM_sign == T & curves$ATRi_count_regulated < 3, "GEMreg_ATRinot","")))))


# Add column "combi_regulation" which describes the regulation of both GEM and at least 3 ATRi (with direction)
curves$combi_regulation <- ifelse(curves$GEM_updown == 'up' & curves$ATRi_updown_atleast3 == 'down', "GEMup_ATRidown",
                                  ifelse(curves$GEM_updown == 'down' & curves$ATRi_updown_atleast3 == 'up', "GEMdown_ATRiup",
                                         ifelse(curves$GEM_updown == 'up' & curves$ATRi_updown_atleast3 == 'up', "GEMup_ATRiup",
                                                ifelse(curves$GEM_updown == 'down' & curves$ATRi_updown_atleast3 == 'down', "GEMdown_ATRidown",
                                                       ifelse(curves$GEM_updown == 'down' & curves$ATRi_updown_atleast3 == '', "GEMdown_ATRinot",
                                                              ifelse(curves$GEM_updown == 'up' & curves$ATRi_updown_atleast3 == '', "GEMup_ATRinot",
                                                                     ifelse(curves$GEM_updown == '' & curves$ATRi_updown_atleast3 == 'down', "GEMnot_ATRidown",
                                                                            ifelse(curves$GEM_updown == '' & curves$ATRi_updown_atleast3 == 'up', "GEMnot_ATRiup",""))))))))


# Annotate kinase-substrate-relationship from PSP. Dataset was downloaded in April 2023 from https://www.phosphosite.org/staticDownloads
# Matching was performed using the column 'Site sequence context_15AA', and only human kinases were considered.
# The new column was named "Kinase"

# Subsets of this column are provided as Supplementary Table in the Paper.
# The following subset of the table is the input for plots of this study (")

curves_forplotting <- curves %>% dplyr::select(c("Site sequence context_15AA",Gene_SitePos,Experiment, "Gene names" , ModSeqShort, starts_with("Ratio"), starts_with("Curve"),
                                                 EC50, pEC50, SQTQ, GEM_log2FC, GEM_qvalue, GEM_sign, GEM_updown, ATRi_regulated, ATRi_count_regulated, 
                                                 ATRi_updown_atleast3, combi_regulation, combi_sign))





################# The following code generates the input table for PTMNavigator, generating the plot for Fig. 5c

# Filter for GEM-only and combi-regulation. Select columns and remove missing values in GEM_pvalue
PTMnavigator <- curves %>% ungroup %>% dplyr::filter(combi_sign %in% c("GEMreg_ATRireg","GEMreg_ATRinot")) %>% 
  dplyr::select(`Modified sequence`, `Proteins_temp`,Experiment, GEM_log2FC, GEM_pvalue, GEM_qvalue,
                combi_sign, `Gene_SitePos`, SQTQ, combi_regulation)  %>% dplyr::filter(!is.na(`GEM_pvalue`))

# Call the experiment "GEM", name all other columns as needed for PTMNavigator
PTMnavigator$Experiment <- "GEM"
colnames(PTMnavigator) <- c('Modified sequence', 'Protein IDs', 'Experiment',	'Fold change'	,'P Value'	,'adjusted pvalue',	'Regulation', 'Gene_SitePos', 'SQTQ', 'combi_regulation')
PTMnavigator <- PTMnavigator[!duplicated(PTMnavigator[,c('Gene_SitePos',	'Fold change'	,'P Value'	,'adjusted pvalue',	'Regulation')]),] # remove duplicates
PTMnavigator$`Modified sequence` <- as.character(map(strsplit(PTMnavigator$`Modified sequence`, split = ";"), 1)) # if not done already, use only first sequence input
PTMnavigator$`Protein IDs` <- as.character(map(strsplit(PTMnavigator$`Protein IDs`, split = ";"), 1)) # if not done already, use only first protein name

# Here, sites regulated by both GEM and ATRi are mocked as "down" regulated (blue color in PTMnavigator)
# And sites regulated only by GEM but not ATRi are mocked as "up" regulated (red color in PTMnavigator)
# The colors were changed afterwards manually to dark red/grey and light red/grey (also based on the SQTQ motif)
# Information used for coloring in the final figure are indicated in columns "SQTQ" (red/grey) and "combi_regulation" (light/dark)
categories <- unique(PTMnavigator$Regulation)
PTMnavigator$Regulation <- ifelse(PTMnavigator$Regulation == 'GEMreg_ATRireg', "down", ifelse(PTMnavigator$Regulation == 'GEMreg_ATRinot', "up",  ""))

# Modify columns to ensure compatibility with required input format
PTMnavigator$`Gene names` <- NULL
PTMnavigator$`adjusted pvalue` <- as.numeric(PTMnavigator$`adjusted pvalue`)
PTMnavigator$`P Value` <- 10^(-as.numeric(PTMnavigator$`P Value`))

# The created file can be found in Zenodo (decryptM/PTMnavigator_input; see README)
# RUN PTMNAVIGATOR https://www.proteomicsdb.org/analytics/ptmNavigator
