#!/bin/bash
# Simple script to run the whole analysis pipeline,
# with and without the biomarkers.

# Note: Run from the main project directory -- e.g. ./R/pipeline.sh

#nohup Rscript ./R/caretTrain-w-biomarkers.R > ./output-w-biomarkers.log 2>&1 &&
#  echo "Completed training with Biomarkers!" &&
#  nohup Rscript ./R/caretTrain-no-biomarkers.R > ./output-no-biomarkers.log 2>&1 &&
#  echo "Completed training without Biomarkers!" &&
  nohup Rscript ./R/caretTrain-sens-no-biomarkers.R > ./output-sens-no-biomarkers.log 2>&1 &&
  echo "Completed training Sensitivity without Biomarkers!" &&
  nohup Rscript ./R/bbccv-sens-no-biomarkers.R >> ./output-sens-no-biomarkers.log 2>&1 &&
  echo "Completed sensitivity bootstrap without Biomarkers!"

echo "Finished training pipeline. Check in case of errors."
