#!/usr/bin/env bash
# Reproduce the annotation/ analysis end to end (Phases 1-5).
# Requires: python3 + annotation/requirements.txt, curl, and network access to
# Ensembl Plants BioMart, Ensembl Compara REST, and NCBI E-utilities.
#   python3 -m pip install -r annotation/requirements.txt
#   bash annotation/reproduce.sh
set -euo pipefail
cd "$(dirname "$0")"                       # -> annotation/
BM="https://plants.ensembl.org/biomart/martservice"

echo "[1/6] Pull GO + Pfam for B. rapa (Chiifu, dataset brapa_eg_gene) from Ensembl Plants BioMart"
GO_XML='<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="plants_mart" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6"><Dataset name="brapa_eg_gene" interface="default"><Attribute name="ensembl_gene_id"/><Attribute name="go_id"/><Attribute name="name_1006"/><Attribute name="namespace_1003"/></Dataset></Query>'
PFAM_XML='<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="plants_mart" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6"><Dataset name="brapa_eg_gene" interface="default"><Attribute name="ensembl_gene_id"/><Attribute name="pfam"/></Dataset></Query>'
curl -s -m 300 -G "$BM" --data-urlencode "query=$GO_XML"   -o brapa_go.tsv
curl -s -m 300 -G "$BM" --data-urlencode "query=$PFAM_XML" -o brapa_pfam.tsv
echo "    brapa_go.tsv=$(wc -l < brapa_go.tsv) rows  brapa_pfam.tsv=$(wc -l < brapa_pfam.tsv) rows"

echo "[2/6] Build ID crosswalk + GO/Pfam annotation layer"
python3 build_annotation.py

echo "[3/6] Curate the tiered floral-scent gene set"
python3 curate_scent_geneset.py

echo "[4/6] Radiation gene-set test (permutation + BH) and power/sensitivity"
python3 phase3_robust.py
python3 phase3_power.py

echo "[5/6] Cross-species orthology (Ensembl Compara REST) + model-species records (NCBI)"
python3 scent_orthology.py
python3 phase4c_link.py

echo "[6/6] Tables + figures"
python3 make_tables.py
python3 make_figures.py
echo "Done. See annotation/tables/, annotation/figures/, and docs/results.html."
