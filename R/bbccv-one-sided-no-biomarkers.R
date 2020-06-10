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
boot.res.file <- "./results/one-sided-results-bbccv-no-biomarkers.RData"
boot.ci.csv <- "./results/one-sided-ci-no-biomarkers.csv"

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
####### BBC-CV One-sided test
#######################


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

load(file=oos.prediction.file)

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
cut.pt <- 0.1
library(Rcpp)
library(aucC)
library(doParallel)
library(abind)
library(foreach)

message("Beginning Bootstrap...")
cl <- makeCluster(5, "FORK", outfile=cluster.file)
registerDoParallel(cl)
boot_diffs <- foreach(b=1:B, 
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
    
    ml_confusion_mat <- factor(out.boot.pred < cut.pt, levels=c(T, F), labels=c("No", "Yes")) %>% 
      table(out.boot.outcome.vec, exclude = F) %>%
      caret::confusionMatrix(positive="Yes")
    
    naive_confusion_mat <- 
      factor((trainingData$TTP[out.boot.idx] <= 24), levels=c(F, T), labels=c("No", "Yes")) %>% 
      table(out.boot.outcome.vec, exclude = F) %>%
      caret::confusionMatrix(positive="Yes")
    
    c(ml_confusion_mat$byClass[1:4] - naive_confusion_mat$byClass[1:4])
  })
  if(inherits(res, "try-error")){
    res <- NA
  }
  return(res)
}
stopCluster(cl)

save(boot_diffs, file=boot.res.file)

print(paste0("# Errors: ", sum(is.na(boot_diffs))))
boot_diffs <- do.call(rbind, boot_diffs[!is.na(boot_diffs)])

apply(boot_diffs, 2, function(x) sum(is.na(x)))

# two_sided_cis <- apply(boot_diffs, 2, quantile, probs=c(0.025, 0.975), na.rm=T)
one_sided_cis <- apply(boot_diffs, 2, quantile, probs=c(0.05,0.95), na.rm=T)

pretty_one_sided_cis <- sapply(1:(ncol(boot_diffs)), function(i){
  estim <- mean(boot_diffs[,i], na.rm=T)
  # conf.inter <- one_sided_cis[i]
  # paste0(round(estim, 3), " (", round(conf.inter,3), ")")
  conf.inter <- one_sided_cis[,i]
  paste0(round(estim, 3), " (", paste(round(conf.inter,3), collapse=", "), ")")
}) 
names(pretty_one_sided_cis) <- colnames(boot_diffs)

write.csv(pretty_one_sided_cis, file=boot.ci.csv)

one_sided_greater_p_values <- apply(boot_diffs, 2, function(x) mean(x <= 0)); print(one_sided_greater_p_values)
one_sided_lesser_p_values <- apply(boot_diffs, 2, function(x) mean(x >= 0)); print(one_sided_lesser_p_values)
