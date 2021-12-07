library(tidyverse)
library(caret)
library(gmodels)
library(caTools)
source('code/my_theme.R')

setwd("./data/training")

data = read_delim("data/balanced_adduct_data_training.csv",
                  delim = ",",
                  col_names = TRUE)

compounds = data %>%
  select(PubChemCID) %>%
  unique()

set.seed(987123)

data_train = data %>%
  mutate(M_Na = case_when(
    M_Na == 1 ~ "yes",
    TRUE ~ "no"
  )) %>%
  mutate(M_Na = factor(M_Na,
                          levels = c("no", "yes")))

#it turns out that for some labs there are still features with zero variance in and these need to be removed
data_train = data_train %>%
  select(-PubChemCID, -Lab, -M_Na, everything())
for (lab in levels(as.factor(data_train$Lab))) {
  data_train <- data_train %>%
    select(-caret::nearZeroVar(data_train %>%
                                 filter(Lab == lab) %>%
                                 select(-PubChemCID, -Lab, -M_Na)))
}

folds <- groupKFold(data_train$Lab, k = 5) 

fitControl <- trainControl(method = "repeatedcv",
                           number = 5,
                           repeats = 1,
                           index = folds,
                           classProbs = TRUE)

set.seed(987123)

classifier_svmpoly <- train(M_Na ~ .,
                    data = data_train %>%
                      select(-PubChemCID, -Lab),
                    method = "svmPoly",
                    trControl = fitControl,
                    metric = "ROC")

saveRDS(classifier_svmpoly,
        file = "classifier_svmpoly.rds")

set.seed(987123)

classifier_svmlinear <- train(M_Na ~ .,
                            data = data_train %>%
                              select(-PubChemCID, -Lab),
                            method = "svmLinear",
                            trControl = fitControl)

saveRDS(classifier_svmlinear,
        file = "classifier_svmlinear.rds")


#data format for RRF and others-----
data_train = data %>%
  mutate(M_Na = factor(M_Na)) 

#it turns out that for some labs there are still features with zero variance in and these need to be removed
data_train = data_train %>%
  select(-PubChemCID, -Lab, -M_Na, everything())
for (lab in levels(as.factor(data_train$Lab))) {
  data_train <- data_train %>%
    select(-caret::nearZeroVar(data_train %>%
                                 filter(Lab == lab) %>%
                                 select(-PubChemCID, -Lab, -M_Na)))
}

folds <- groupKFold(data_train$Lab, k = 5) 

fitControl <- trainControl(method = "repeatedcv",
                           number = 5,
                           repeats = 1,
                           index = folds)

set.seed(987123)

classifier_RRF <- train(M_Na ~ .,
                            data = data_train %>%
                              select(-PubChemCID, -Lab),
                            method = "RRF",
                            trControl = fitControl)
saveRDS(classifier_RRF,
        file = "classifier_RRF.rds")

set.seed(987123)

classifier_ada <- train(M_Na ~ .,
                        data = data_train %>%
                          select(-PubChemCID, -Lab),
                        method = "ada",
                        trControl = fitControl)

saveRDS(classifier_ada,
        file = "classifier_ada.rds")

set.seed(987123)

classifier_DT <- train(M_Na ~ .,
                        data = data_train %>%
                          select(-PubChemCID, -Lab),
                        method = "rpart",
                        trControl = fitControl)

saveRDS(classifier_DT,
        file = "classifier_DT.rds")

set.seed(987123)

classifier_knn <- train(M_Na ~ .,
                        data = data_train %>%
                          select(-PubChemCID, -Lab),
                        method = "knn",
                        trControl = fitControl)

saveRDS(classifier_knn,
        file = "classifier_knn.rds")

set.seed(987123)

classifier_naive_bayes <- train(M_Na ~ .,
                        data = data_train %>%
                          select(-PubChemCID, -Lab),
                        method = "naive_bayes",
                        trControl = fitControl)

saveRDS(classifier_naive_bayes,
        file = "classifier_naive_bayes.rds")

data_train_pred <- data %>%
  mutate(M_Na_pred_svmLinear = predict(classifier_svmlinear, newdata = data),
         M_Na_pred_svmPoly = predict(classifier_svmpoly, newdata = data),
         M_Na_pred_RRF = predict(classifier_RRF, newdata = data),
         M_Na_pred_ada = predict(classifier_ada, newdata = data),
         M_Na_pred_DT = predict(classifier_DT, newdata = data),
         M_Na_pred_knn = predict(classifier_knn, newdata = data),
         M_Na_pred_naive_bayes = predict(classifier_naive_bayes, newdata = data))

data_train_pred = data_train_pred %>%
  mutate(M_Na_pred_svmLinear = case_when(
    M_Na_pred_svmLinear == "yes" ~ 1,
    TRUE ~ 0
  ),
  M_Na_pred_svmPoly = case_when(
    M_Na_pred_svmPoly == "yes" ~ 1,
    TRUE ~ 0
  ))

write_delim(data_train_pred,
            "data/Predicted_adduct_formation_training_set.csv",
            delim = ",")

