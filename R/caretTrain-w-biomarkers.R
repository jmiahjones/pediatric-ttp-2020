## Loading the Datasets

trained.model.file <- "./results/trained-models-w-biomarkers.RData"
raw.data.file <- "./data/Time to Positivity - Merged Philly and URMC. FINAL. 12.8.20.xlsx"

main_vars <- c('Sex', 'Age', 'GAge', 'Appearance', 
               'Temperature', 'DaysIll', 'Cough', 
               'UTI', 'TTP')
sens_vars <- c(main_vars, "WBCC", "Bands")

library(readxl)
library(dplyr)
library(caret)
library(pROC)

fullData <- read_xlsx(raw.data.file) %>%
  rename(WBCC=`White Blood Cell Count`,
         Bands=`Total Bands`,
         Appearance=`Ill Appearance`,
         GAge=`Gestational Age`,
         DaysIll=`Number of Days of Illness`,
         Cough=`Cough Present`,
         TTP=`Blood Culture Time to Positivity (Hours)`,
         UTI=`Urinary Tract Inflammation`,
         Bacteremia=`Bacteremia #1 (Outcome #1)`,
         BacteremiaSens=`Bacteremia #2 (Outcome #2-Sensitivity analysis)`
  ) %>%
  mutate(TTP = ifelse(!is.na(TTP), TTP,
                      130)) %>% 
  mutate_at(vars(Sex, Appearance,
                 Cough,
                 UTI,
                 Bacteremia),
            ~factor(., exclude=c("2",NA))) %>%
  filter(TTP < 130)


##### Missing Data
mice::md.pattern(fullData)
## There are ~7 missing UTI  +/- records. Let's keep a separate NA level.
## We set these as a separate factor level.

fullData <- fullData %>%
  mutate_at(vars(UTI),
            ~addNA(.))


N.Y.U.lvls <- c("No", "Yes", "Unknown")
levels(fullData$Sex) <- c("M", "F")
levels(fullData$Appearance) <- c("No", "Yes")
levels(fullData$Cough) <- c("Yes", "No")
levels(fullData$UTI) <- N.Y.U.lvls # set NA level
levels(fullData$Bacteremia) <- c("No", "Yes")

mice::md.pattern(fullData)

# trainingData <- fullData # %>% select(-site)

train_props <- fullData %>% pull(Bacteremia) %>% 
  table %>% prop.table

prop_wt <- train_props[1]/train_props[2]

record_ids <- fullData$`Record ID`
trainingData <- dplyr::select(fullData, all_of(sens_vars), Bacteremia)


############# Begin Training ################

library(parallel)
library(doParallel)
cl <- parallel::makeCluster(4, "FORK")
registerDoParallel(cl)
num_folds=5
num_repeats=10
set.seed(321597)
seeds=sample(100000)[1:(num_folds*num_repeats + 1)]
seeds=lapply(seeds, function(seed){
  set.seed(seed)
  sample(100000, size=5000)
})

set.seed(321597)
folds <- caret::createMultiFolds(trainingData$Bacteremia,
                                 k=num_folds,
                                 times=num_repeats)

trControl <- caret::trainControl(method="repeatedcv",
                                 search = "grid",
                                 number=num_folds, repeats=num_repeats,
                                 index=folds,
                                 classProbs=TRUE,
                                 summaryFunction=twoClassSummary,
                                 savePredictions = TRUE,
                                 seeds=seeds,
                                 returnResamp = "final",
                                 trim=TRUE)

trainingX <- dplyr::select(trainingData, -Bacteremia)

message("Beginning training.")

#--------------- RandomForest ---------------#
rf_grid <- expand.grid(mtry=1:11)
rfFit <- train(Bacteremia ~ ., data = trainingData,
               method = "rf",
               trControl = trControl,
               # verbose = FALSE,
               tuneGrid = rf_grid,
               ## Specify which metric to optimize
               metric = "ROC")

message("Trained RF")

#--------------- Glmnet ---------------#
glmnetFit <- train(Bacteremia ~ ., data = trainingData,
                   ## specify binomial likelihood
                   family="binomial",
                   method = "glmnet",
                   trControl = trControl,
                   # verbose = FALSE,
                   tuneLength=50,
                   ## Specify which metric to optimize
                   metric = "ROC")


message("Trained GLMNET")

#--------------- SVM ---------------#
sigmas <- kernlab::sigest(model.matrix(~.-1, trainingX),
                          na.action = na.omit, scaled = TRUE)
svmWeightGrid <- expand.grid(
  Weight=c(prop_wt,1),
  sigma=seq(sigmas[1], sigmas[2], length.out = 100),
  C=2^seq(-2,4, length.out = 10)
)

svmWeightFit <- train(Bacteremia ~ ., data = trainingData,
                      method = "svmRadialWeights",
                      trControl = trControl,
                      # verbose = FALSE,
                      tuneGrid=svmWeightGrid,
                      ## Specify which metric to optimize
                      metric = "ROC")

message("Trained SVM")

stopImplicitCluster()

#-------------- xgbTree --------------#
xgbtrControl <- caret::trainControl(method="repeatedcv",
                                    search = "grid",
                                    number=num_folds, repeats=num_repeats,
                                    index=folds,
                                    classProbs=TRUE,
                                    summaryFunction=twoClassSummary,
                                    savePredictions = TRUE,
                                    allowParallel = FALSE,
                                    seeds=seeds,
                                    returnResamp = "final",
                                    trim=TRUE)

xgbGrid = expand.grid(
  nrounds = c(100, 500, 1000, 2000, 3000),
  max_depth = c(1, 3, 5),
  eta = c(0.01,0.001),
  gamma = c(1, 2, 3),
  colsample_bytree = c(0.5, 1.0),
  min_child_weight = c(0.5, 1),
  subsample = c(0.8, 1.0)
)
xgbFit <- train(x=model.matrix(~.-1, data=trainingX),
                y=trainingData$Bacteremia,
                method = "xgbTree",
                trControl = xgbtrControl,
                scale_pos_weight = prop_wt,
                tuneGrid=xgbGrid,
                metric = "ROC")

message("Trained xgboost")

save(list=ls(), file=trained.model.file)
