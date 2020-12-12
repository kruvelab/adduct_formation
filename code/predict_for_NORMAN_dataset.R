library(plotly)

NORMAN = read_delim("susdat_2019-12-06-112040.csv",
                    delim = ",",
                    col_names = TRUE) %>%
  rename(PubChemCID = PubChem_CID,
         CAS = `CAS_RN Dashboard`)

NORMAN = NORMAN %>%
  group_by(Monoiso_Mass) %>%
  mutate(count = n()) %>%
  ungroup()

ggplot(data = NORMAN) +
  geom_histogram(mapping = aes(x = count))

NORMAN = NORMAN %>%
  filter(count > 10)

NORMAN = NORMAN %>%
  select(PubChemCID, Monoiso_Mass) %>%
  unique() %>%
  mutate(PubChemCID = as.numeric(PubChemCID))

data_path <- "C:/Users/annel/OneDrive - Kruvelab/Kruvelab/computational/PubChem/fingerprints"

files <- dir(data_path, pattern = "*.csv") # get file names
remove <- c("fingerprint_Compound_", ".sdf.csv") #idetify the part of filename that will be removed
Fingerprints_for_NORMAN <- tibble() #here we will collect all matchies 

setwd(data_path)
for (filename in files) {
  if (length(Fingerprints_for_NORMAN$PubChemCID) < length(NORMAN$PubChemCID)) {
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
    compounds_small <- NORMAN %>%
      filter(PubChemCID > first_CID & PubChemCID < last_CID) %>%
      left_join(datafile)
    Fingerprints_for_NORMAN <- Fingerprints_for_NORMAN %>%
      bind_rows(compounds_small)
  } else {
    break
  }
}
setwd('..')

write_delim(Fingerprints_for_NORMAN,
            "Fingerprints_NORMAN.csv",
            delim = ",")
Fingerprints_for_NORMAN <- Fingerprints_for_NORMAN %>%
  na.omit()

Decoded_figenrprints_for_NORMAN <- tibble()
for (fingerprint in Fingerprints_for_NORMAN$Fingerprint) {
  current_fingerprint <- fingerprint_hex_to_dataframe(fingerprint)
  print(current_fingerprint)
  Decoded_figenrprints_for_NORMAN <- Decoded_figenrprints_for_NORMAN %>%
    bind_rows(current_fingerprint)
}

adducts_summary_NORMAN <- as_tibble(cbind(Fingerprints_for_NORMAN, Decoded_figenrprints_for_NORMAN))

write_delim(adducts_summary_NORMAN,
            "Fingerprints_NORMAN.csv",
            delim = ",")


adducts_summary_NORMAN = read_delim("Fingerprints_NORMAN.csv",
                              delim = ",",
                              col_names = TRUE)

Bond_descriptors_summary <- tibble()
for (smiles in adducts_summary_NORMAN$SMILES) {
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

adducts_summary_NORMAN <- adducts_summary_NORMAN %>%
  left_join(Bond_descriptors_summary)

write_delim(adducts_summary_NORMAN,
            "PubChem_fingerprints_bond_properties_NORMAN.csv",
            delim = ",")

NORMAN_pred = adducts_summary_NORMAN %>%
  na.omit() %>%
  mutate(M_Na_pred = predict(classifier_ada, newdata = adducts_summary_NORMAN, type = "prob")[,2])

NORMAN_adduct_plot = 
ggplot(data = NORMAN_short %>%
         filter(Monoiso_Mass < 300)) +
  geom_point(mapping = aes(x = Monoiso_Mass,
                           y = M_Na_pred,
                           text = Name),
             alpha = 1/5) +
  labs(x = "monoisotopic mass (Da)", y = "probability of adduct formation") +
  my_theme

ggplotly(NORMAN_adduct_plot)

ggsave("adduct_formation_NORMAN.svg",
       width = 8,
       height = 8,
       units = "cm")

NORMAN_short = NORMAN_pred %>%
  select(PubChemCID, M_Na_pred, Monoiso_Mass) %>%
  left_join(NORMAN %>%
              mutate(PubChemCID = as.numeric(PubChemCID)) %>%
              select(PubChemCID, Name, CAS))
NORMAN_short = NORMAN_short %>%
  na.omit()


write_delim(NORMAN_short,
            "adduct_formation_NORMAN_short.csv",
            delim = ",")
NORMAN_short = read_delim("adduct_formation_NORMAN_short.csv",
                          delim = ",",
                          col_names = TRUE)

available_compounds = read_delim("MMK_chemicals.csv",
                                 delim = ",",
                                 col_names = TRUE)

NORMAN_available_at_MMK = NORMAN_short %>%
  left_join(available_compounds %>%
              select(CAS, Rum))

NORMAN_available_at_MMK = NORMAN_available_at_MMK %>%
  na.omit()

NORMAN_large_dif = NORMAN_available_at_MMK %>%
  group_by(Monoiso_Mass) %>%
  mutate(prob_dif = max(M_Na_pred) - min(M_Na_pred)) %>%
  ungroup() %>%
  arrange(desc(prob_dif))

write_delim(NORMAN_large_dif,
            "NORMA_dif_available_MMK.csv",
            delim = ",")

NORMAN_adduct_plot = 
  ggplot(data = NORMAN_available_at_MMK ) +
  geom_point(mapping = aes(x = Monoiso_Mass,
                           y = M_Na_pred,
                           text = Name),
             alpha = 1/5) +
  labs(x = "monoisotopic mass (Da)", y = "probability of adduct formation") +
  my_theme

ggplotly(NORMAN_adduct_plot)
