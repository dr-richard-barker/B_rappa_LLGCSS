# Phase 1 — gene-ID reconciliation & annotation layer

This directory builds the join key and functional annotation that the
scent-under-radiation analysis depends on. See `../README_PLAN.md` for the full plan.

## The key finding that unblocked everything

The two experiments looked un-joinable because they use different-looking gene IDs:

- **scent** dataset → `Bra######` (Chiifu / Ensembl Plants `Brapa_1.0`)
- **radiation** dataset → `BRA######`

These are the **same loci differing only in case**. After capitalising, **all 31,756
scent genes map 1:1 onto radiation genes** (100%). So the datasets join directly on
gene ID — no cross-assembly mapping is needed. (`Brara.A00001` FPsc IDs and the
`BRADO####` *Bradyrhizobium* KEGG file are unrelated third namespaces — do not use them.)

## Scent gene set — defined in *B. rapa* space

Because this Ensembl Plants release exposes no Arabidopsis orthologs or descriptions for
*B. rapa*, the floral-volatile gene set is defined directly from **GO terms + Pfam
domains** (terpenoid, benzenoid/phenylpropanoid, fatty-acid/green-leaf-volatile,
apocarotenoid routes, plus SABATH volatile esterases). Curation is in
`build_annotation.py` (`SCENT_GO`, `SCENT_PFAM`); non-scent methyltransferases
(histone/RNA/protein) are deliberately excluded.

Result: **413 scent-candidate genes**, 322 present in both experiments.

## Files

Committed (small, the deliverables):
- `build_annotation.py` — builds the crosswalk + annotation + scent set.
- `scent_radiation_test.py` — preliminary enrichment test (Phase 3 preview).
- `scent_geneset.tsv` — the 413 curated scent-candidate genes.

Regenerable / git-ignored (large; on disk, ship via deposit):
- `gene_annotation.tsv` — per-gene GO/Pfam + scent flag (41,018 genes).
- `id_crosswalk.tsv` — `Bra` ↔ `BRA` + presence flags.
- `brapa_go.tsv`, `brapa_pfam.tsv` — raw Ensembl Plants BioMart pulls.

## Reproduce

```bash
# 1. Pull GO + Pfam from Ensembl Plants BioMart (dataset brapa_eg_gene = Brapa_1.0):
#    (see the curl commands in build_annotation.py header / git history)
# 2. Build the annotation layer + scent set:
python3 annotation/build_annotation.py
# 3. Preliminary radiation × scent test:
python3 annotation/scent_radiation_test.py
```

## Preliminary answer to "does GCR radiation influence floral scent?"

*Exploratory — our radiation counts are an in-house re-analysis (the primary radiation
study is collaborator-led); treat as a signal to follow up, not a final result.*

- The **radiation effect is very weak** at 40 cGy: only **13 genes** are DE in WT
  (40 vs 0 cGy, adjP<0.1) and **0** in the *anthocyaninless* mutant. The dominant signal
  in the radiation experiment is genotype, not dose.
- The scent gene set **validates** on the scent axis: 29/310 candidates are DE between
  High- and Low-scent lines.
- **Scent × radiation:** only **1** of the 13 radiation DEGs is a scent-candidate
  (`Bra029041`, an O-methyltransferase / SABATH-adjacent volatile-ester enzyme;
  radiation log2FC +1.8). Fold-enrichment is 8.3× but **Fisher p = 0.11 — not
  significant**. The test is under-powered by the tiny radiation DEG count, so this is a
  *"no significant evidence, but under-powered"* result, not a clean negative.

**Defensible statement:** this repository now establishes the join and the scent gene set
needed to ask the question; on the current in-house radiation counts there is no
significant preferential perturbation of scent genes, with `Bra029041` as the single
follow-up candidate. Re-running against the collaborators' radiation counts (and higher
doses, if available) is the logical next step (Phase 3/4).
