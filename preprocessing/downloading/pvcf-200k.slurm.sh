#!/usr/bin/env bash

#SBATCH -p plgrid
#SBATCH -t 3:0:0
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB
#SBATCH --array=1-977%3
#SBATCH --output=/net/scratch/people/plgjacekh/slurm-log/%A_%a.out
#SBATCH --error=/net/scratch/people/plgjacekh/slurm-log/%A_%a.err

set -x

export PROJ=/net/people/plgjacekh/projects/ifpan-gosborcz-ukb/
# use SCRATCH as a download folder
# export WD=$PROJ/raw/pvcf-200k/
export WD=$SCRATCH/pvcf-200k/
TOOLS_DIR=/net/archive/groups/plggneuromol/tools/


if [ ! -f $WD/pvcf_blocks.txt ]
then
    flock --nonblock \
        $WD/wget.lock \
        bash -c '
            wget -nd -P $WD biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/pvcf_blocks.txt
            cp $PROJ/raw/dataset/k62979r48815.key $WD/
        '
    sleep 10
fi


cd $WD


# Population level exome OQFE variants, pVCF format - interim 200k release
FIELD=23156

LINE=$(grep -P "^$SLURM_ARRAY_TASK_ID\t" $WD/pvcf_blocks.txt)

chromosome=$(echo "$LINE" | cut -f 2)
block=$(echo "$LINE" | cut -f 3)
echo $chromosome $block
# use SCRATCH as a download folder
# $PROJ/ukb-utils/gfetch $FIELD -a../dataset/k62979r48815.key -c$chromosome -b$block
$PROJ/ukb-utils/gfetch $FIELD -ak62979r48815.key -c$chromosome -b$block
