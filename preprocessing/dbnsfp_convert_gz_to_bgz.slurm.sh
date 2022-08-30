#!/usr/bin/env bash

#SBATCH -p plgrid
#SBATCH -t 8:0:0
#SBATCH --cpus-per-task=1
#SBATCH --output=/net/scratch/people/plgjacekh/slurm-log/%j.out
#SBATCH --error=/net/scratch/people/plgjacekh/slurm-log/%j.err


DB_PATH=/net/people/plgjacekh/projects/ifpan-gosborcz-ukb/raw/dbNSFP4.3a/

for FILE in "$DB_PATH/*chr*.gz"
do
    BGZ_FILE="${FILE%gz}bgz"

    if [ ! -s "$BGZ_FILE" ]
    then
        echo "Converting $(basename -- $FILE)..."
        zcat "$FILE" | bgzip -c > "$BGZ_FILE"
    fi
done
