#!/bin/bash
# Simple script to run the whole analysis pipeline,
# with and without the biomarkers.

# Note: Run from the main project directory -- e.g. ./R/pipeline.sh

nohup Rscript ./R/caretTrain-w-biomarkers.R > ./output-w-biomarkers.log 2>&1 &&
  nohup Rscript ./R/bbccv-w-biomarkers.R >> ./output-w-biomarkers.log 2>&1 &&
  echo "Completed training with Biomarkers!" &&
  nohup Rscript ./R/caretTrain-no-biomarkers.R > ./output-no-biomarkers.log 2>&1 &&
  nohup Rscript ./R/bbccv-no-biomarkers.R >> ./output-no-biomarkers.log 2>&1 &&
  nohup Rscript ./R/bbccv-one-sided-no-biomarkers.R >> ./output-no-biomarkers.log 2>&1

echo "Finished pipeline. Check in case of errors."
