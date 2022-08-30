#!/usr/bin/env bash

#SBATCH -A plgdepresja
#SBATCH --partition plgrid
#SBATCH --mem=20GB
#SBATCH --time 1-00:00:00
#SBATCH --job-name ukb-r
#SBATCH --output=/net/scratch/people/plgjacekh/slurm-log/%j.out
#SBATCH --error=/net/scratch/people/plgjacekh/slurm-log/%j.err


ID="49742"


echo "MD5..."
ukb-utils/ukbmd5 raw/dataset/ukb$ID.enc

echo "Decrypt..."
ukb-utils/ukbunpack raw/dataset/ukb$ID.enc raw/dataset/k62979r$ID.key

echo "Convert..."
ukb-utils/ukbconv raw/dataset/ukb$ID.enc_ukb docs -eukb-utils/encoding.ukb
ukb-utils/ukbconv raw/dataset/ukb$ID.enc_ukb r -eukb-utils/encoding.ukb
