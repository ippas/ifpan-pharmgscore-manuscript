#!/usr/bin/env bash


DIR=data/pgxpop/

for GENE in CYP2C9 CYP2D6 CYP4F2 CYP2C19 CYP3A5 DPYD SLCO1B1 CYP2B6 NUDT15
do
    for PART in {0..1}
    do
        echo "[$(date +"%H:%M:%S")] Doing $GENE, part $PART..."
        python PGxPOP/bin/PGxPOP.py \
            --vcf "$DIR/pgxpop-part-$PART.vcf.bgz" \
            -g "$GENE" \
            --build grch38 \
            > "pgxpop-part-$PART-$GENE.log" 2>&1
        mv pgxpop_results.txt "$DIR/pgxpop-part-$PART-$GENE.txt"
    done
done
