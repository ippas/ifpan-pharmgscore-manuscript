#!/usr/bin/env bash

#SBATCH -p plgrid
#SBATCH -t 4:0:0
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --output=/net/archive/groups/plggneuromol/jhajto/slurm-std/%j.out
#SBATCH --error=/net/archive/groups/plggneuromol/jhajto/slurm-std/%j.err

set -ex

PROJ=/net/people/plgjacekh/projects/ifpan-gosborcz-ukb/
WD=$PROJ/raw/adr/plink/

cd $WD

# Population-level FE variants, PLINK format - initial 50k release
FIELD=23160

# Downloading only for chr 1, but it seems to contain all data
# *.bed
srun $PROJ/ukb-utils/gfetch $FIELD -a../../dataset/k62979r47888.key -c1
# *.fam
srun $PROJ/ukb-utils/gfetch $FIELD -a../../dataset/k62979r47888.key -c1 -m
# *.bim
wget -nv -nd biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/ukb_fe_exm_chrall_v1.bim
