#!/bin/bash


set -e
# Setting arguments
DRYRUN="--dry-run" # "--sub-sample", "--fake-data"
SUBSIDY=0.2

if [ $DRYRUN == "--dry-run" ]
then
    printf "\n\n\n Dry Run \n\n\n"
fi

printf "\n\n\n Running 0 Subsidy \n\n\n"
# Run with 0 subsidy 
Rscript ./optim/optimal_allocation.R  ${DRYRUN} \
    --min-cost  \
    --target-constraint=0.85 \
    --output-path=optim/data \
    --output-filename=dry-run \
    --dry-run-subsidy=0 # i.e. Social planner doesn't know
printf "\n\n\n 0 Subsidy Finished \n\n\n"

printf "\n\n\n Running ${SUBSIDY} Subsidy \n\n\n"
# Run with actual subsidy
Rscript ./optim/optimal_allocation.R  ${DRYRUN} \
    --min-cost  \
    --target-constraint=0.85 \
    --output-path=optim/data \
    --output-filename=dry-run \
    --dry-run-subsidy=0.2 # i.e. Social planner doesn't know
printf "\n\n\n ${SUBSIDY} Subsidy Finished \n\n\n"

# Now postprocessing

printf "\n\n\n Running postprocessing 1 \n\n\n"
# Social planner doesn't know 
Rscript ./optim/postprocess_allocation.R  \
    --min-cost  \
    --target-constraint=0.85  \
    --input-path=optim/data \
    --optim-input-filename=dry-run-subsidy-0-optimal-allocation.csv \
    --village-input-filename=dry-run-subsidy-0-village-locations.csv \
    --pot-input-filename=dry-run-subsidy-0-pot-locations.csv \
    --demand-input-filename=dry-run-subsidy-0-demand-data.csv \
    --true-demand-input-filename=dry-run-subsidy-${SUBSIDY}-demand-data.csv \
    --output-path=optim \
    --output-filename=problem-a-optimal-pot-plot.png \
    --comp-demand
printf "\n\n\n Postprocessing 1 Finished \n\n\n"

printf "\n\n\n Running postprocessing 2 \n\n\n"
# Social planner does know 
Rscript ./optim/postprocess_allocation.R  \
    --min-cost  \
    --target-constraint=0.85 \
    --input-path=optim/data \
    --optim-input-filename=dry-run-subsidy-${SUBSIDY}-optimal-allocation.csv \
    --village-input-filename=dry-run-subsidy-${SUBSIDY}-village-locations.csv \
    --pot-input-filename=dry-run-subsidy-${SUBSIDY}-pot-locations.csv \
    --demand-input-filename=dry-run-subsidy-${SUBSIDY}-demand-data.csv \
    --true-demand-input-filename=dry-run-subsidy-${SUBSIDY}-demand-data.csv \
    --misspecified-optim-input-filename=dry-run-subsidy-0-optimal-allocation.csv \
    --output-path=optim \
    --output-filename=problem-c-optimal-pot-plot.png 
printf "\n\n\n Postprocessing 2 Finished \n\n\n"
