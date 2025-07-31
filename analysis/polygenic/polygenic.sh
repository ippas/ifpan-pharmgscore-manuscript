#!/bin/bash

export TOOLS_DIR="tools/"
export PROJ_DIR="data/"

wget "https://workflows.intelliseq.com/intelliseq/workflows/-/raw/dev/src/main/wdl/modules/gvcf-pgx-genotype/gvcf-pgx-genotype.wdl" \
    -O tmp/gvcf-pgx-genotype.wdl_wget


SAMPLE_ID="NA00001"
X_SAMPLE_ID="${SAMPLE_ID//-/X}"

# because of mv function at the end two array jobs can't write to the same directory
export OUT_DIR=data/polygenic/cromwell-output/sample
INPUTS_JSON="{
    \"gvcf_pgx_genotype.sample_id\": \"$X_SAMPLE_ID\",
    \"gvcf_pgx_genotype.vcf_gz\": \"$PROJ_DIR/sample.vcf.gz\",
    \"gvcf_pgx_genotype.vcf_gz_tbi\": \"$PROJ_DIR/sample.vcf.gz.tbi\",
    \"gvcf_pgx_genotype.selected_genes\": [
        \"CYP2B6\", \"CYP2C19\", \"CYP2C9\", \"CYP2D6\", \"CYP3A5\",
        \"CYP4F2\", \"DPYD\", \"NUDT15\", \"SLCO1B1\"
    ]
}"
OPTIONS_JSON="{
    \"final_workflow_outputs_dir\": \"$OUT_DIR\",
    \"use_relative_output_paths\": true
}"

# Run the Cromwell job for each sample
bash -c 'java \
    -Dconfig.file=$PROJ_DIR/analysis/polygenic/cromwell-wd/config.conf \
    -jar $TOOLS_DIR/cromwell run \
        $PROJ_DIR/analysis/polygenic/cromwell-wd/gvcf-pgx-genotype.wdl \
        --inputs /dev/fd/3 \
        --options /dev/fd/4' 3< <(echo $INPUTS_JSON) 4< <(echo $OPTIONS_JSON)

# Move output files and rename them based on the sample ID
mv -- $OUT_DIR/openpgx.json $OUT_DIR/${SAMPLE_ID}-openpgx.json
mv -- $OUT_DIR/pharmgkb.json $OUT_DIR/${SAMPLE_ID}-pharmgkb.json
mv -- $OUT_DIR/${X_SAMPLE_ID}-merged_models.json $OUT_DIR/${SAMPLE_ID}-merged_models.json
mv -- $OUT_DIR/${X_SAMPLE_ID}-phased.vcf.gz $OUT_DIR/${SAMPLE_ID}-phased.vcf.gz
mv -- $OUT_DIR/${X_SAMPLE_ID}-phased.vcf.gz.tbi $OUT_DIR/${SAMPLE_ID}-phased.vcf.gz.tbi

# Get diplotype phenotype
OPENPGX_PATH_STEM="$OUT_DIR/$SAMPLE_ID-openpgx"
$TOOLS_DIR/jq -r 'to_entries[] | [.key, .value] | @tsv' \
    "$OPENPGX_PATH_STEM.json" > "$OPENPGX_PATH_STEM.tsv"
docker run --rm \
  -v "$OUT_DIR":"$OUT_DIR" \
  -v "$(dirname "$OPENPGX_PATH_STEM.tsv")":"$(dirname "$OPENPGX_PATH_STEM.tsv")" \
  pgkb/pharmcat:2.13.0 \
  pharmcat --phenotyper \
    --phenotyper-outside-call-file "$OPENPGX_PATH_STEM.tsv" \
    --output-dir "$OUT_DIR"
