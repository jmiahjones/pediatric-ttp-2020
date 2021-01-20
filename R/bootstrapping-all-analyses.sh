#!/bin/bash
# Bootstrap trained model to get analysis results
nohup Rscript ./R/bbccv-w-biomarkers.R >> ./output-w-biomarkers.log 2>&1 &&
  echo "Completed bootstrap with Biomarkers!" &&
  nohup Rscript ./R/bbccv-sens-no-biomarkers.R >> ./output-sens-no-biomarkers.log 2>&1 &&
  echo "Completed sensitivity bootstrap without Biomarkers!" &&
  nohup Rscript ./R/bbccv-no-biomarkers.R >> ./output-no-biomarkers.log 2>&1 &&
  nohup Rscript ./R/bbccv-one-sided-no-biomarkers.R >> ./output-no-biomarkers.log 2>&1 &&
  echo "Completed bootstrap without Biomarkers!"

echo "Completed all bootstraps. Check in case of errors!"