library(tidyverse)
library(caret)
library(gmodels)
library(caTools)
source('./code/my_theme.R')

setwd("./data/training")

#Adduct data ----
#UT data-----
UT_data = read_delim("UT_adducts_CID.csv",
                      delim = ",",
                      col_names = TRUE) %>%
  na.omit()

UT_data = UT_data %>%
  mutate(M_H = case_when(
    Class == 0 ~ 1,
    Class == 2 ~ 1,
    TRUE ~ 0
  ),
  M_Na = case_when(
    Class == 1 ~ 1,
    Class == 2 ~ 1,
    TRUE ~ 0
  ),
  Lab = "UT")

UT_data  %>% 
  group_by(M_Na) %>%
  summarise(n())

set.seed(123)

#there are about 2x so (224 vs 127) many compounds not forming an adduct as there are ones forming an adduct
UT_data_adduct = sample_n(UT_data %>%
                            filter(M_Na == 1),
                          size = 100,
                          replace = FALSE)

UT_data_No_adduct = sample_n(UT_data %>%
                             filter(M_Na == 0),
                           size = 100,
                           replace = FALSE)
UT_data_balanced = UT_data_adduct %>%
  bind_rows(UT_data_No_adduct)

#Corey data -----

Corey_data = read_delim("Corey_adducts_CID.csv",
                      delim = ",",
                      col_names = TRUE)

Corey_data = Corey_data %>%
  mutate(M_Na = case_when(
                    Class == 1 ~ 1,
                    Class == 2 ~ 1,
                    TRUE ~ 0),
         Lab = "Corey")

Corey_data  %>% 
  group_by(M_Na) %>%
  summarise(n())

#there are about 6x so (72 vs 525) many compounds forming an adduct as there are ones forming an adduct
Corey_data_No_adduct = sample_n(Corey_data %>%
                             filter(M_Na == 0),
                           size = 100,
                           replace = TRUE)

Corey_data_adduct = sample_n(Corey_data %>%
                                   filter(M_Na == 1),
                                 size = 100,
                                 replace = FALSE)

Corey_data_balanced = Corey_data_adduct %>%
  bind_rows(Corey_data_No_adduct)

#SU data -----

SU_data = read_delim("SU_adducts_CID.csv",
                      delim = ",",
                      col_names = TRUE)

SU_data = SU_data %>%
  mutate(M_Na = case_when(
                  SlopeNa == 0 ~ 0,
                  TRUE ~ 1),
        M_H = case_when(
                  SlopeH == 0 ~ 0,
                  TRUE ~ 1),
        Lab = "SU")

SU_data  %>% 
  group_by(M_Na, M_H) %>%
  summarise(n())

#this is an almost even dataset (53 vs 41)

SU_data_adduct = sample_n(SU_data %>%
                                filter(M_Na == 1),
                              size = 50,
                              replace = TRUE)

SU_data_No_adduct = sample_n(SU_data %>%
                                filter(M_Na == 0),
                              size = 50,
                              replace = TRUE)
SU_data_balanced = SU_data_adduct %>%
  bind_rows(SU_data_No_adduct)

#Celma data----

Celma_data = read_delim("Celma_adduct_CID.csv",
                              delim = ",",
                              col_names = TRUE)

Celma_data  %>% 
  group_by(M_Na) %>%
  summarise(n())
#this is an almost even dataset (274 vs 247)

Celma_data_adduct = sample_n(Celma_data %>%
                             filter(M_Na == 1),
                           size = 200,
                           replace = FALSE)

Celma_data_No_adduct = sample_n(Celma_data %>%
                                filter(M_Na == 0),
                              size = 200,
                              replace = FALSE)
Celma_data_balanced = Celma_data_adduct %>%
  bind_rows(Celma_data_No_adduct)

#Picache data----

Picache_data = read_delim("Picache_adducts_CID.csv",
                          delim = ",",
                          col_names = TRUE)

Picache_data = Picache_data %>%
  filter(Adduct == "[M+Na]" | Adduct == "[M+H]") %>%
  select(PubChemCID, Adduct, CCS) %>%
  na.omit() %>%
  group_by(PubChemCID, Adduct) %>%
  summarise(CCS = mean(CCS)) %>%
  ungroup()

Picache_data = Picache_data %>%
  spread(key = Adduct, value = CCS)

Picache_data = Picache_data %>%
  mutate(M_Na = case_when(
    is.na(`[M+Na]`) ~ 0,
    TRUE ~ 1
  ))

Picache_data = Picache_data %>%
  select(PubChemCID, M_Na) %>%
  mutate(Lab = "Picache")


Picache_data  %>% 
  group_by(M_Na) %>%
  summarise(n())
#this is unbalance agin. 151 does not form adducts while 469 does.

Picache_data_adduct = sample_n(Picache_data %>%
                                filter(M_Na == 1),
                              size = 100,
                              replace = FALSE)

Picache_data_No_adduct = sample_n(Picache_data %>%
                                   filter(M_Na == 0),
                                 size = 100,
                                 replace = FALSE)

Picache_data_balanced = Picache_data_adduct %>%
  bind_rows(Picache_data_No_adduct)

#binding all together----

data_all = UT_data %>%
  bind_rows(Corey_data) %>%
  bind_rows(SU_data) %>%
  bind_rows(Celma_data) %>%
  bind_rows(Picache_data)

data = UT_data_balanced %>%
  bind_rows(Corey_data_balanced) %>%
  bind_rows(SU_data_balanced) %>%
  bind_rows(Celma_data_balanced) %>%
  bind_rows(Picache_data_balanced) %>%
  select(-Name, -Class, -M_H, -SlopeH, -SlopeNa)


#Fingerprints----

fingerprints <- read_delim("PubChem_fingerprints_Bond_properties.csv",
                           delim = ",",
                           col_names = TRUE) %>%
  select(-SMILES, -X1, -SDFlineStart, -SDFlineEnd, -SDFID, -Fingerprint) %>%
  select(everything(), PubChemCID)

fingerprints <- fingerprints %>%
  select(-caret::nearZeroVar(fingerprints)) #removing fingerprints that do not change significantely between samples

fingerprints_with_NA <- fingerprints %>%
  select(everything()) %>%
  summarise_all(funs(sum(is.na(.))))
#only few columns have na-s and only for 5 rows, so I decide to drop the rows.

fingerprints <- fingerprints %>%
  na.omit() %>%
  unique()

#Putting adduct data and fingerprints together----
data <- data %>%
  na.omit() %>%
  left_join(fingerprints)%>%
  na.omit()

write_delim(data,
            "balanced_adduct_data_training.csv",
            delim = ",")

data_all = data_all %>%
  left_join(fingerprints)

data_test = data_all %>%
  anti_join(data)

write_delim(data_test,
            "adduct_data_test.csv",
            delim = ",")

