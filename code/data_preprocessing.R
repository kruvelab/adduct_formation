library(tidyverse)
library(caret)
library(gmodels)
library(caTools)

setwd("C:/Users/annel/OneDrive - Kruvelab/Kruvelab/R codes for variouse things/Computational resources/PubChem")


#Adduct data ----
UT_data <- read_delim("UT_adducts_CID.csv",
                      delim = ",",
                      col_names = TRUE)
UT_data <- UT_data %>%
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

#there are about 2x so (224 vs 127) many compounds not forming an adduct as there are ones forming an adduct
UT_data_adduct <- sample_n(UT_data %>%
                            filter(M_Na == 1),
                          size = 100,
                          replace = TRUE)

UT_data_No_adduct <- sample_n(UT_data %>%
                             filter(M_Na == 0),
                           size = 100,
                           replace = TRUE)
UT_data <- UT_data_adduct %>%
  bind_rows(UT_data_No_adduct)

#Corey data -----

Corey_data <- read_delim("Corey_adducts_CID.csv",
                      delim = ",",
                      col_names = TRUE)

Corey_data <- Corey_data %>%
  mutate(M_Na = case_when(
                    Class == 1 ~ 1,
                    Class == 2 ~ 1,
                    TRUE ~ 0),
         Lab = "Corey")

Corey_data  %>% 
  group_by(M_Na) %>%
  summarise(n())

#there are about 6x so (72 vs 525) many compounds forming an adduct as there are ones forming an adduct
Corey_data_No_adduct <- sample_n(Corey_data %>%
                             filter(M_Na == 0),
                           size = 100,
                           replace = TRUE)

Corey_data_adduct <- sample_n(Corey_data %>%
                                   filter(M_Na == 1),
                                 size = 100,
                                 replace = TRUE)

Corey_data <- Corey_data_adduct %>%
  bind_rows(Corey_data_No_adduct)

#SU data -----

SU_data <- read_delim("SU_adducts_CID.csv",
                      delim = ",",
                      col_names = TRUE)

SU_data <- SU_data %>%
  mutate(M_Na = case_when(
                  SlopeNa == 0 ~ 0,
                  TRUE ~ 1),
        M_H = case_when(
                  SlopeH == 0 ~ 0,
                  TRUE ~ 1),
        Lab = "SU")

SU_data  %>% 
  group_by(M_Na) %>%
  summarise(n())

#this is an almost even dataset (53 vs 41)

SU_data_adduct <- sample_n(SU_data %>%
                                filter(M_Na == 1),
                              size = 100,
                              replace = TRUE)

SU_data_No_adduct <- sample_n(SU_data %>%
                                filter(M_Na == 0),
                              size = 100,
                              replace = TRUE)
SU_data <- SU_data_adduct %>%
  bind_rows(SU_data_No_adduct)

#binding all together

data <- UT_data %>%
  bind_rows(Corey_data) %>%
  bind_rows(SU_data) 

# write_delim(data,
#             "all_adduct_data_training.csv",
#             delim = ",")

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
  na.omit()

fingerprints_copy <- fingerprints 

fingerprints_copy <- fingerprints_copy %>%
  select(-findCorrelation(fingerprints_copy,
                cutoff = 0.8)) %>%
  bind_cols(fingerprints %>% select(PubChemCID))

#Putting adduct data and fingerprints together----
data <- data %>%
  select(-Name, -Class, -M_H, -SlopeH, -SlopeNa) %>%
  na.omit() %>%
  left_join(fingerprints_copy) %>%
  na.omit()

split = sample.split(data$M_Na, SplitRatio = 0.8)

data_train = data %>%
  filter(split == TRUE) %>%
  mutate(M_Na = as.factor(M_Na))

folds <- groupKFold(data_train$Lab, k = 3) 

fitControl <- trainControl(method = "repeatedcv",
                           number = 3,
                           repeats = 1,
                           index = folds)

classifier <- train(M_Na ~ .,
                    data = data_train %>%
                      select(-PubChemCID, -Lab),
                    method = "svmPoly",
                    trCOntrol = fitControl)

data_pred <- data %>%
  mutate(M_Na_pred = predict(classifier, newdata = data),
         split = split)

data_test_pred<- data_pred %>%
  filter(split == FALSE)

#RRF::importance(classifier$finalModel)

ggplot(data = data_pred) +
  geom_point(mapping = aes(x = M_Na,
                           y = M_Na_pred),
             size = 5,
             alpha = 1/50) +
  facet_wrap(~Lab) +
  theme_bw()

for (lab in levels(factor(data_pred$Lab))) {
  data_this <- data_test_pred %>%
    filter(Lab == lab)
  print(lab)
  print(
    CrossTable(data_this$M_Na, data_this$M_Na_pred)
  )
}
CrossTable(data_test_pred$M_Na, data_test_pred$M_Na_pred)


#read in the validation set and calculate for this