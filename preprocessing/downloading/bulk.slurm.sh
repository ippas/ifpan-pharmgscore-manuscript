#!/usr/bin/env bash

#SBATCH -p plgrid
#SBATCH -t 1-0:0:0
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --output=/net/archive/groups/plggneuromol/jhajto/slurm-std/%j.out
#SBATCH --error=/net/archive/groups/plggneuromol/jhajto/slurm-std/%j.err

set -e

if [ -z "$1" ] ; then
	echo "Must specify batch file"
	exit 1
elif ! [ -e "$1" ] ; then
	echo "Invalid batch file"
	exit 1
else
	FETCH_CAP=50000
	if [ $(wc -l <$1) -ge $FETCH_CAP ]; then
		echo "WARNING: Batch file too long (>= $FETCH_CAP)"
		echo "Trim the file or specify -s and -m flag in ukbfetch"
	fi
fi


PROJ=/net/people/plgjacekh/projects/ifpan-gosborcz-ukb/
WD=$PROJ/raw/adr/gvcf/

cp $1 $WD
cd $WD
BULK_FILE=$(basename "$1")


srun $PROJ/ukb-utils/ukbfetch -a../../dataset/k62979r47888.key -b$BULK_FILE
rm $BULK_FILE
find . -maxdepth 1 -type f -regextype posix-egrep -regex ".*.gvcf.gz(.tbi)?" \
	-exec bash -c 'mv -v $0 ${0/.gvcf.gz/.g.vcf.gz}' {} \;
