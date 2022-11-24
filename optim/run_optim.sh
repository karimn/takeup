#!/bin/bash


set -e
# Setting arguments
DRYRUN="" # "--sub-sample", "--fake-data"
SKIPOPTIMISATION=""
# SKIPOPTIMISATION="dont-skip-pls"
SUBSIDY=0.25
TARGETCONSTRAINT=0.32
VERSION=66



Rscript ./optim/predict-takeup-for-optim.R 66 \
    --num-post-draws=400 \
    --from-csv \
    bracelet \
    control

Rscript ./optim/optimal_allocation.R  \
                             --posterior-median \
                             --num-cores=12 \
                             --min-cost  \
                             --target-constraint=${TARGETCONSTRAINT} \
                             --output-path=optim/data \
                             --output-filename=structural \
                             --input-path=optim/data  \
                             --village-input-filename=village-df.csv \
                             --pot-input-filename=pot-df.csv \
                             --time-limit=10000 \
                             --demand-input-filename=pred_demand_dist_fit${VERSION}.csv

Rscript ./optim/optimal_allocation.R  \
                             --num-cores=12 \
                             --min-cost  \
                             --target-constraint=${TARGETCONSTRAINT} \
                             --output-path=optim/data \
                             --output-filename=structural \
                             --input-path=optim/data  \
                             --village-input-filename=village-df.csv \
                             --pot-input-filename=pot-df.csv \
                             --time-limit=10000 \
                             --demand-input-filename=pred_demand_dist_fit${VERSION}.csv

Rscript ./optim/postprocess_allocation.R  \
                             --min-cost \
                             --posterior-median \
                             --target-constraint=0.33 \
                             --input-path=optim/data \
                             --optim-input-a-filename=structural-median-optimal-allocation.rds \
                             --village-input-filename=village-df.csv \
                             --pot-input-filename=pot-df.csv \
                             --demand-input-a-filename=pred_demand_dist_fit66.csv \
                             --output-path=optim \
                             --output-basename=structural-median-test \
                             --map-plot

Rscript ./optim/postprocess_allocation.R  \
                             --min-cost \
                             --target-constraint=0.33 \
                             --input-path=optim/data \
                             --optim-input-a-filename=structural-post-draws-optimal-allocation.rds \
                             --village-input-filename=village-df.csv \
                             --pot-input-filename=pot-df.csv \
                             --demand-input-a-filename=pred_demand_dist_fit66.csv \
                             --output-path=optim \
                             --output-basename=structural-post-draws-test \
                             --map-plot

# Rscript ./optim/postprocess_allocation.R --min-cost \
#                              --target-constraint=0.33 \
#                              --input-path=optim/data \
#                              --optim-input-a-filename=init-control-optimal-allocation.csv \
#                              --village-input-filename=village-df.csv \
#                              --pot-input-filename=pot-df.csv \
#                              --demand-input-a-filename=approx-control-demand.csv \
#                              --output-path=optim \
#                              --output-basename=init-control \
#                              --map-plot 

# Rscript ./optim/postprocess_allocation.R --min-cost \
#                              --target-constraint=0.33 \
#                              --input-path=optim/data \
#                              --optim-input-a-filename=init-bracelet-optimal-allocation.csv \
#                              --village-input-filename=village-df.csv \
#                              --pot-input-filename=pot-df.csv \
#                              --demand-input-a-filename=approx-bracelet-demand.csv \
#                              --output-path=optim \
#                              --output-basename=init-bracelet 


# Rscript ./optim/postprocess_allocation.R --min-cost \
#                              --target-constraint=0.33 \
#                              --input-path=optim/data \
#                              --optim-input-a-filename=init-reduced-optimal-allocation.csv \
#                              --village-input-filename=village-df.csv \
#                              --pot-input-filename=pot-df.csv \
#                              --demand-input-a-filename=approx-reduced-bracelet-demand.csv \
#                              --output-path=optim \
#                              --output-basename=init-reduced-bracelet 

# if [ $DRYRUN == "--dry-run" ]
# then
#     printf "\n\n\n Dry Run \n\n\n"
# fi


# if [ $SKIPOPTIMISATION != "--skip" ] 
# then
#     printf "\n\n\n Running 0 Subsidy \n\n\n"
#     # Run with 0 subsidy 
#     Rscript ./optim/optimal_allocation.R  ${DRYRUN} \
#         --min-cost  \
#         --target-constraint=${TARGETCONSTRAINT} \
#         --output-path=optim/data \
#         --output-filename=dry-run \
#         --dry-run-subsidy=0 # i.e. Social planner doesn't know
#     printf "\n\n\n 0 Subsidy Finished \n\n\n"

#     printf "\n\n\n Running ${SUBSIDY} Subsidy \n\n\n"
#     # Run with actual subsidy
#     Rscript ./optim/optimal_allocation.R  ${DRYRUN} \
#         --min-cost  \
#         --target-constraint=${TARGETCONSTRAINT} \
#         --output-path=optim/data \
#         --output-filename=dry-run \
#         --dry-run-subsidy=${SUBSIDY} # i.e. Social planner doesn't know
#     printf "\n\n\n ${SUBSIDY} Subsidy Finished \n\n\n"
# fi

# # Now postprocessing

# printf "\n\n\n Running postprocessing 1 \n\n\n"
# # Social planner doesn't know 
# Rscript ./optim/postprocess_allocation.R  \
#     --min-cost  \
#     --target-constraint=${TARGETCONSTRAINT}  \
#     --input-path=optim/data \
#     --optim-input-a-filename=dry-run-subsidy-0-optimal-allocation.csv \
#     --village-input-filename=dry-run-subsidy-0-village-locations.csv \
#     --pot-input-filename=dry-run-subsidy-0-pot-locations.csv \
#     --demand-input-a-filename=dry-run-subsidy-0-demand-data.csv \
#     --demand-input-b-filename=dry-run-subsidy-${SUBSIDY}-demand-data.csv \
#     --output-path=optim \
#     --output-basename=sim-optim-example-ms \
#     --comp-demand \
#     --map-plot
# printf "\n\n\n Postprocessing 1 Finished \n\n\n"

# printf "\n\n\n Running postprocessing 2 \n\n\n"
# # Social planner does know 
# Rscript ./optim/postprocess_allocation.R  \
#     --min-cost  \
#     --target-constraint=${TARGETCONSTRAINT} \
#     --input-path=optim/data \
#     --optim-input-a-filename=dry-run-subsidy-${SUBSIDY}-optimal-allocation.csv \
#     --optim-input-b-filename=dry-run-subsidy-0-optimal-allocation.csv \
#     --village-input-filename=dry-run-subsidy-${SUBSIDY}-village-locations.csv \
#     --pot-input-filename=dry-run-subsidy-${SUBSIDY}-pot-locations.csv \
#     --demand-input-a-filename=dry-run-subsidy-${SUBSIDY}-demand-data.csv \
#     --output-path=optim \
#     --output-basename=sim-optim-example
# printf "\n\n\n Postprocessing 2 Finished \n\n\n"
