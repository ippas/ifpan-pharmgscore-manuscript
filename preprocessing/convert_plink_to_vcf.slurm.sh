#!/usr/bin/env bash

#SBATCH -p plgrid
#SBATCH -t 4:0:0
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --output=/net/archive/groups/plggneuromol/jhajto/slurm-std/%j.out
#SBATCH --error=/net/archive/groups/plggneuromol/jhajto/slurm-std/%j.err

set -e

PROJ=/net/people/plgjacekh/projects/ifpan-gosborcz-ukb/
WD=$PROJ/raw/adr/plink/

cd $WD

# Unify file names
# ukb23160_c1_b0_v1.bed
mv -v ukb_fe_exm_chrall_v1.bim ukb23160_c1_b0_v1.bim
mv -v ukb23160_c1_b0_v1_s49945.fam ukb23160_c1_b0_v1.fam

TOOLS_DIR=/net/archive/groups/plggneuromol/jhajto/tools/plink1.9
$TOOLS_DIR/plink --bfile ukb23160_c1_b0_v1 --recode vcf --out ukb23160
