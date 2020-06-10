####################################################
########### Bootstrap Results ######################
####################################################
library(xgboost)
library(glmnet)
library(kernlab)
library(randomForest)
library(caret)
library(dplyr)

cluster.file <- "./cluster-no-biomarkers.log"
trained.model.file <- "./results/trained-models-no-biomarkers.RData"
oos.prediction.file <- "./results/trained-model-predictions-no-biomarkers.RData"
boot.res.file <- "./results/results-bbccv-no-biomarkers.RData"
boot.ci.csv <- "./results/ci-no-biomarkers.csv"

load(file=trained.model.file)


#######################
####### 24-Hour Rule
#######################
library(pROC)
naive_roc = roc(as.numeric(trainingData$Bacteremia) - 1, 1*(trainingData$TTP <= 24))
naive_roc
ci.auc(naive_roc)

naive_char=coords(naive_roc, 
                  0.5,
                  ret=c("sensitivity", "specificity", "npv", "ppv"),
                  transpose = T)

print(naive_char)


#######################
####### BBC-CV Results
#######################

trained_models <- list(
  rf=rfFit,
  glmnet=glmnetFit,
  # gam=gamFit,
  svm=svmWeightFit,
  xgb=xgbFit
)
resamples(trained_models) %>% summary(metric="ROC")

xgboost::xgb.importance(model=xgbFit$finalModel) %>% select(Feature, Gain) %>% 
  mutate(Gain=round(100*Gain, 1)) %>% 
  write.csv(file="./results/feature_importance_no_biomarkers.csv", row.names = F)


### Bootstrap Approach
set.seed(530981)
B <- 2000
rowIdxs <- 1:nrow(trainingData)
bootIdxs <- replicate(B, sample(rowIdxs, size=nrow(trainingData), replace=T), simplify = F)
fold.names <- names(folds)

extract.repeat <- function(fold.name){
  strsplit(fold.name, ".", fixed=T) %>% sapply(function(x) x[2])
}

# ######################################
library(doParallel)
library(abind)
library(foreach)

if(file.exists(oos.prediction.file)){
  load(file=oos.prediction.file)
} else {
  
  cl <- makeCluster(4, "FORK", outfile=cluster.file)
  registerDoParallel(cl)
  
  pred.arr.list <- foreach::foreach(mod.idx=seq_along(trained_models)) %dopar% {
    # lapply(seq_along(trained_models), function(mod.idx){
    mod <- trained_models[[mod.idx]]
    pred.df <- mod$pred
    
    configs.df <- pred.df %>% 
      select(-pred, -obs, -No, -Yes, -rowIndex, -Resample) %>% unique %>% 
      as_tibble
    
    config.cols <- colnames(configs.df)
    
    pred.df <- pred.df %>% mutate(Rep = extract.repeat(Resample)) %>% 
      select(-Resample)
    
    cv.reps <- unique(pred.df$Rep)
    
    pred.arr <- sapply(cv.reps, function(rep.name){
      this.rep.pred.df <- pred.df %>% filter(Rep == rep.name)
      sapply(1:nrow(configs.df), function(config.idx){
        this.config <- configs.df[config.idx,]
        this.rep.pred.df %>% semi_join(this.config, by=config.cols) %>% 
          arrange(rowIndex) %>% 
          pull(Yes)
      })
    }, simplify="array")
  }
  stopImplicitCluster()
  
  pred.arr <- abind::abind(pred.arr.list, along=2)
  # pred.arr[row.idx, config.idx, cv.rep.idx]
  rm(pred.arr.list)
  
  save(pred.arr, file=oos.prediction.file)
}

dim.pred.arr <- dim(pred.arr)
nconfig <- dim.pred.arr[2]
nrepeat <- dim.pred.arr[3]

config.selection <- function(boot.mat, max=TRUE){
  # accepts matrix of aucs, which we maximize
  summ.mat <- apply(boot.mat, 2, function(col){
    c(median(col, na.rm=T), min(col, na.rm=T))
  })
  
  if(max){
    best.med <- max(summ.mat[1,])
  } else {
    best.med <- min(summ.mat[1,])
  }
  best.idxs <- which(summ.mat[1,] == best.med)
  
  if(length(best.idxs) > 1){
    if(max){
      final.best.idx <- best.idxs[which.max(summ.mat[2, best.idxs])]
    } else {
      final.best.idx <- best.idxs[which.min(summ.mat[2, best.idxs])]
    }
  } else {
    final.best.idx <- best.idxs
  }
  return(final.best.idx)
}

cut.pts <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
library(Rcpp)
library(aucC)
library(doParallel)
library(abind)
library(foreach)

message("Beginning Bootstrap...")
cl <- makeCluster(5, "FORK", outfile=cluster.file)
registerDoParallel(cl)
boot_stats <- foreach(b=1:B, 
                      .noexport = "aucC"
                      ) %dopar% {
  res <- try({
    print(paste0("Beginning Boot: ", b))
    boot.idx <- bootIdxs[[b]]
    out.boot.idx <- setdiff(rowIdxs, boot.idx)
    boot.pred.arr <- pred.arr[boot.idx, ,]
    outcome.vec <- trainingData$Bacteremia
    boot.outcome.vec <- outcome.vec[boot.idx]
    out.boot.outcome.vec <- outcome.vec[out.boot.idx]
    
    # resampled aucs
    # boot.char <- matrix(data=NA, dim=c(nconfig, nrepeat))
    boot.char <- foreach(repeat.idx=1:nrepeat, .combine="cbind") %:%
      foreach(config.idx=1:nconfig, .combine="c") %do% {
        aucC_wrap(label=boot.outcome.vec,
                  scores=boot.pred.arr[, config.idx, repeat.idx])
      }
    
    print(paste0("Selecting config for boot ", b))
    best.config.idx <- config.selection(boot.char, max=TRUE)
    set.seed(530981 + b)
    k <- sample(nrepeat, size=1)
    out.boot.pred <- pred.arr[out.boot.idx, best.config.idx, k]
    auc <- mean(boot.char[best.config.idx, ])
    
    sapply(cut.pts, function(cut.pt){
      confusion_mat <- factor(out.boot.pred < cut.pt, levels=c(T, F), labels=c("No", "Yes")) %>% 
        table(out.boot.outcome.vec, exclude = F) %>%
        caret::confusionMatrix(positive="Yes")
      
      c(auc=auc, confusion_mat$byClass[1:4])
    })
  })
  if(inherits(res, "try-error")){
    res <- NA
  }
  return(res)
}
stopCluster(cl)

save(boot_stats, file=boot.res.file)

print(paste0("# Errors: ", sum(is.na(boot_stats))))
boot_stats <- abind(boot_stats[!is.na(boot_stats)], along=3)

apply(boot_stats[,1,], 1, function(x) sum(is.na(x)))
apply(boot_stats[,7,], 1, function(x) sum(is.na(x)))

# boot_stats <- boot_stats[,-1,]
cis <- apply(boot_stats, c(1,2), quantile, probs=c(0.025, 0.975), na.rm=T)

pretty_cis <- sapply(1:(dim(boot_stats)[2]), function(i){
  sapply(seq_along(cis[1,,i]), function(j){
    estim <- mean(boot_stats[j,i,], na.rm=T)
    conf.inter <- cis[,j,i]
    paste0(round(estim, 3), " (", paste(round(conf.inter,3), collapse=", "), ")")
  })
}) %>% t
rownames(pretty_cis) <- cut.pts#[-1]
colnames(pretty_cis) <- rownames(boot_stats[,,1])

write.csv((pretty_cis), file=boot.ci.csv)

naive_est <- (factor(1*(trainingData$TTP <= 24), levels=c(0, 1), labels=c("No", "Yes")) %>% 
  table(trainingData$Bacteremia, exclude = F) %>%
  caret::confusionMatrix(positive="Yes"))$byClass[1:4]

naive_boot <- foreach(b=1:B, .combine=rbind) %do% {
  boot.idx <- bootIdxs[[b]]
  confusion_mat <- factor(1*(trainingData$TTP[boot.idx] <= 24), levels=c(0, 1), labels=c("No", "Yes")) %>% 
    table(trainingData$Bacteremia[boot.idx], exclude = F) %>%
    caret::confusionMatrix(positive="Yes")
  
  c(confusion_mat$byClass[1:4])
}
naive_cis <- apply(naive_boot, 2, quantile, probs=c(0.025, 0.975), na.rm=T)

pretty_naive_cis <- sapply(seq_along(naive_est), function(i){
  paste0(round(naive_est[i], 3), " (", paste(round(naive_cis[,i],3), collapse=", "), ")")
})

print("Naive Results:")
print(pretty_naive_cis)
