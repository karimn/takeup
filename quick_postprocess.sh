#!/usr/bin/env bash


postprocess_struct_models () {
    echo "RUNNING MODEL: $1"
    echo "RUNNING VERSION: $VERSION"

    Rscript --no-save \
            --no-restore \
            --verbose \
            quick_ate_postprocess.R \
            ${VERSION} \
            ${POSTPROCESS_INOUT_ARGS} \
            --model=$1 \
            1 2 3 4 > temp/log/struct-postprocess_$1_${VERSION}.txt 2>&1 

    Rscript --no-save \
            --no-restore \
            --verbose \
            quick_ate_postprocess.R \
            ${VERSION} \
            ${POSTPROCESS_INOUT_ARGS} \
            --model=$1 \
            --prior \
            1 2 3 4 > temp/log/struct-postprocess_$1_${VERSION}.txt 2>&1 

    Rscript --no-save \
            --no-restore \
            --verbose \
            quick_roc_postprocess.R \
            ${VERSION} \
            ${POSTPROCESS_INOUT_ARGS} \
            --model=$1 \
            1 2 3 4 > temp/log/struct-postprocess_$1_${VERSION}.txt 2>&1 

    Rscript --no-save \
            --no-restore \
            --verbose \
            quick_submodel_postprocess.R \
            ${VERSION} \
            ${POSTPROCESS_INOUT_ARGS} \
            --model=$1 \
            1 2 3 4 > temp/log/struct-postprocess_$1_${VERSION}.txt 2>&1 
}

postprocess_rf_models () {
    echo "RUNNING: $1"
    echo "RUNNING VERSION: $VERSION"
    Rscript --no-save \
            --no-restore \
            --verbose \
            quick_ate_postprocess.R \
            ${VERSION} \
            ${POSTPROCESS_INOUT_ARGS} \
            --model=$1 \
            1 2 3 4 > temp/log/rf-postprocess_$1_${VERSION}.txt 2>&1 

    Rscript --no-save \
            --no-restore \
            --verbose \
            quick_ate_postprocess.R \
            ${VERSION} \
            ${POSTPROCESS_INOUT_ARGS} \
            --model=$1 \
            --prior \
            1 2 3 4 > temp/log/rf-postprocess_$1_${VERSION}.txt 2>&1 



}

