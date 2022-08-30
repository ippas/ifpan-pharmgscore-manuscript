#!/usr/bin/env bash

#SBATCH -p plgrid
#SBATCH -t 3-0:0:0
#SBATCH --cpus-per-task=2
#SBATCH --mem=20GB
#SBATCH --output=/net/scratch/people/plgjacekh/slurm-log/%j.out
#SBATCH --error=/net/scratch/people/plgjacekh/slurm-log/%j.err

set -ex

PROJ=/net/people/plgjacekh/projects/ifpan-gosborcz-ukb/
TOOLS_DIR=/net/archive/groups/plggneuromol/tools/


# sites
$TOOLS_DIR/gsutil/gsutil -m cp -r \
    "gs://gcp-public-data--gnomad/release/3.1.1/ht/genomes/gnomad.genomes.v3.1.1.sites.ht/" \
    $PROJ/raw/gnomad/

# coverage
$TOOLS_DIR/gsutil/gsutil -m cp -r \
    "gs://gcp-public-data--gnomad/release/3.0.1/coverage/genomes/gnomad.genomes.r3.0.1.coverage.ht/" \
    $PROJ/raw/gnomad/
