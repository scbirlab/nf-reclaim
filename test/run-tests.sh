#!/usr/bin/env bash

set -exuo pipefail

GITHUB=${1:-no}
script_dir="$(dirname $0)"

if [ "$GITHUB" == "gh" ]
then
    export NXF_CONTAINER_ENGINE=docker
    export NXF_DOCKER_OPTS="-u $(id -u):$(id -g)"
    docker_flag='-profile gh'
else
    export SINGULARITY_FAKEROOT=1
    docker_flag=''
fi

# Examples without sample sheet
nextflow run "$script_dir"/.. \
    -resume $docker_flag \
    --test \
    --outputs output-83332 \
    --organism_id 83332

# Examples with sample sheet
nextflow run "$script_dir"/.. \
    -resume $docker_flag \
    --test \
    -c nextflow.config \
    --sample_sheet inputs/sample-sheet.csv \
    --inputs inputs \
    --outputs output-sample-sheet