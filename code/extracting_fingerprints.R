library(tidyverse)
library(stringr)
source('fingerprints_hex_to_dataframe.R')
source('calculating bond properties from graph.R')

setwd("C:/Users/annel/OneDrive - Kruvelab/Kruvelab/R codes for variouse things/Computational resources/PubChem")

#----Combining all neccessary PubChemCIDs----
compounds_SU <- read_delim("SU_adducts_CID.csv",
                               delim = ",",
                               col_names = TRUE)

compounds_UT <- read_delim("UT_adducts_CID.csv",
                           delim = ",",
                           col_names = TRUE)

compounds_Corey <- read_delim("Corey_adducts_CID.csv",
                           delim = ",",
                           col_names = TRUE)

#collect all CID values that need to be looked up and 
compounds <- compounds_Lorean %>% select(PubChemCID) %>%
  bind_rows(compounds_SU %>% select(PubChemCID)) %>%
  bind_rows(compounds_UT %>% select(PubChemCID)) %>%
  bind_rows(compounds_Corey %>% select(PubChemCID)) %>%
  na.omit() %>%
  unique()

data_path <- "C:/Users/annel/OneDrive - Kruvelab/Kruvelab/R codes for variouse things/Computational resources/PubChem/unread files"

files <- dir(data_path, pattern = "*.csv") # get file names
remove <- c("fingerprint_Compound_", ".sdf.csv") #idetify the part of filename that will be removed
Fingerprints_for_adducts_summary <- tibble() #here we will collect all matchies 

setwd(data_path)
for (filename in files) {
  if (length(Fingerprints_for_adducts_summary$PubChemCID) < length(compounds$PubChemCID)) {
    datafile <- read_delim(filename,
                           delim = "\t",
                           col_names = TRUE) %>%
      rename(PubChemCID = V2,
             SMILES = V3,
             Fingerprint = V4)
    print(filename)
    PubCID_range <- str_remove_all(filename, paste(remove, collapse = "|"))
    PubCID_range <- str_split(PubCID_range, "_")
    first_CID <- as.double(PubCID_range[[1]][1])
    last_CID <- as.double(PubCID_range[[1]][2])
    compounds_small <- compounds %>%
      filter(PubChemCID > first_CID & PubChemCID < last_CID) %>%
      left_join(datafile)
    Fingerprints_for_adducts_summary <- Fingerprints_for_adducts_summary %>%
      bind_rows(compounds_small)
  } else {
    break
  }
}
setwd('..')

write_delim(Fingerprints_for_adducts_summary,
            "Fingerprints_all.csv",
            delim = ",")
Fingerprints_for_adducts_summary <- Fingerprints_for_adducts_summary %>%
  na.omit()

Decoded_figenrprints_for_adduct_summary <- tibble()
for (fingerprint in Fingerprints_for_adducts_summary$Fingerprint) {
  current_fingerprint <- fingerprint_hex_to_dataframe(fingerprint)
  print(current_fingerprint)
  Decoded_figenrprints_for_adduct_summary <- Decoded_figenrprints_for_adduct_summary %>%
    bind_rows(current_fingerprint)
}

adducts_summary <- as_tibble(cbind(Fingerprints_for_adducts_summary, Decoded_figenrprints_for_adduct_summary))

write_delim(adducts_summary,
            "Fingerprints_all.csv",
            delim = ",")

adducts_summary <- read_delim("Fingerprints_all.csv",
                             delim = ",",
                             col_names = TRUE)

Bond_descriptors_summary <- tibble()
for (smiles in adducts_summary$SMILES) {
  print(smiles)
  bond_descriptors_SMILES <- as_tibble(t(as.integer(bond_descriptors(smiles))))
  bond_descriptors_SMILES <- bond_descriptors_SMILES %>%
    mutate(SMILES = smiles)
  Bond_descriptors_summary <- Bond_descriptors_summary %>%
    bind_rows(bond_descriptors_SMILES)
}
Bond_descriptors_summary <- Bond_descriptors_summary %>%
  rename(min_dist_2O = V1,
         max_dist_2O = V2,
         min_dist_carbonyl = V3,
         max_dist_carbonyl = V4,
         max_path_non_rotatabale_bonds_carbonyl = V5,
         min_path_non_rotatabale_bonds_carbonyl = V6,
         max_path_non_rotatabale_bonds_2O = V7,
         min_path_non_rotatabale_bonds_2O = V8)

adducts_summary <- adducts_summary %>%
  left_join(Bond_descriptors_summary)

write_delim(adducts_summary,
            "PubChem_fingerprints_bond_properties.csv",
            delim = ",")
