source("code/PubChem_fingerprint_calculation.R")
source("code/calculating_bond_properties_from_graph.R")

predict_adduct_formation = function(SMILES) {
  classifier_svmPoly = readRDS("models/classifier_svmPoly.rds")
  bond_descriptors_this_SMILES = bond_descriptors(SMILES)
  pubchem_descriptors_this_SMILES = pubchem_fingerprint(SMILES)
  if (length(bond_descriptors_this_SMILES) != 0 & length(pubchem_descriptors_this_SMILES) != 0) {
    Bond_descriptors_summary <- bond_descriptors_this_SMILES %>% 
      left_join(pubchem_descriptors_this_SMILES)
  }
  Bond_descriptors_summary = Bond_descriptors_summary %>%
    mutate(M_Na_pred = predict(classifier_svmPoly, newdata = Bond_descriptors_summary, type = "prob")[,2])
  return(Bond_descriptors_summary$M_Na_pred)
}

Vectorize(predict_adduct_formation)
