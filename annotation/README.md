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

## Scent gene set — defined in *B. rapa* space, tiered (Phase 2)

Because this Ensembl Plants release exposes no Arabidopsis orthologs or descriptions for
*B. rapa*, the floral-volatile gene set is defined directly from **GO terms + Pfam
domains**. Phase 2 (`curate_scent_geneset.py`) structures it by biosynthetic **route**
and **confidence tier**, and drops over-broad terms (e.g. GO:0008299 "isoprenoid
biosynthetic process", which also captures sterols / photosynthetic carotenoids).

- **Tier 1 (core volatile-forming enzymes, 108 genes):** terpene synthases
  (Pfam PF01397/PF03936), SABATH methyltransferases (PF03492), lipoxygenase (PF00305),
  carotenoid cleavage dioxygenases (PF03055).
- **Tier 2 (supporting / route-level, 255 genes):** BAHD acyltransferases (PF02458),
  O-methyltransferases (PF00891 / GO:0008171), PAL (PF00221), terpenoid-route GO,
  jasmonate metabolism.
- **Total: 363 genes** (269 present in both experiments). Family sources cited in the
  script header (Dudareva 2013; Chen 2011; D'Auria 2006; Effmert 2005).

`scent_geneset.tsv` columns: `bra_id, BRA_id, in_scent_data, in_radiation_data, tier,
routes, evidence`.

## Files

Committed (small, the deliverables):
- `build_annotation.py` — Phase 1: builds the crosswalk + GO/Pfam annotation.
- `curate_scent_geneset.py` — Phase 2: tiered, route-classified scent gene set.
- `scent_radiation_test.py` — tier-aware enrichment test (Phase 3 preview).
- `scent_orthology.py` — Phase 4: cross-species conservation via Ensembl Compara.
- `scent_geneset.tsv` — the 363 curated scent genes (tier + route).
- `scent_orthology_matrix.tsv` — scent gene × species ortholog-count matrix.
- `scent_reference_species.tsv` — Phase 4b: literature-curated scent genes in
  Petunia / snapdragon / rose (species absent from Ensembl Plants).

Regenerable / git-ignored (large; on disk, ship via deposit):
- `gene_annotation.tsv` — per-gene GO/Pfam + scent flag (41,018 genes).
- `id_crosswalk.tsv` — `Bra` ↔ `BRA` + presence flags.
- `brapa_go.tsv`, `brapa_pfam.tsv` — raw Ensembl Plants BioMart pulls.

## Phase 4 — cross-species conservation of the scent module (Ensembl orthology)

Ensembl Plants Compara does not include our Chiifu `Bra######` assembly (only
`brassica_rapa_ro18`), so conservation is **anchored on Arabidopsis** — a curated panel
of 12 canonical floral-scent genes, one+ per route — with orthologs pulled across all
plant genomes via the Compara REST homology endpoint (`scent_orthology.py` →
`scent_orthology_matrix.tsv`; raw JSON cached in `ortho_cache/`).

**Result — the scent module is an ancient toolkit with lineage-specific tailoring:**
- **All four biosynthetic routes are deeply conserved** across flowering plants — every
  panel gene has orthologs from Brassicaceae through Solanaceae (tomato/potato) to
  monocots (rice). The core scent machinery predates the crucifers.
- **The volatile-*tailoring* enzymes expand lineage-specifically** — this is the known
  engine of scent diversity:
  - benzenoid/ester O-methyltransferases: **COMT1 = 29 Brassica copies** (vs 1 in
    tomato/potato), **BSMT1/SABATH = 13** (many-to-many) — strong Brassicaceae expansion.
  - terpene synthases expand elsewhere (e.g. **TPS10 = 27 orthologs in grape**, and in
    monocots) — a different lineage's elaboration.
  - `brassica_rapa_ro18` has **lost the JMT ortholog** (0 copies) — a contraction.
- Part of the Brassica copy inflation reflects the **Brassica whole-genome triplication**,
  not necessarily scent-specific selection — stated as a caveat, not a result.

**Tie back to the radiation result:** `Bra029041` — the single gene that moved under
radiation — is an **O-methyltransferase**, i.e. it belongs to the *most Brassica-expanded,
most regulatorily flexible* scent-tailoring family (COMT/SABATH), not the conserved core.
So the "floral scent may be altered in space" hypothesis is best framed as: *if* radiation
perturbs scent, the expanded ester/benzenoid-methyltransferase families are the likely
point of action — a specific, testable prediction for Phase 3 on higher-dose data.

**Scope caveats:** the classic scent-model plants (Petunia, snapdragon, rose) are **not in
Ensembl Plants** — "other scent producers" here means Brassicas + Solanaceae + the broad
plant tree. A Petunia/snapdragon comparison would need Sol Genomics / dedicated orthology
(Phase 4b). Compara's B. rapa is R-o-18, so its counts are a proxy for our Chiifu data.

## Phase 4b — the model scent species Ensembl can't reach (curated)

The classic floral-scent models — **Petunia, snapdragon (*Antirrhinum*), rose (*Rosa*)** —
are **not in Ensembl Plants Compara**, and snapdragon is **absent from OrthoDB** too, so
their scent orthologs cannot be pulled programmatically. `scent_reference_species.tsv` is a
**literature-curated** table of the *characterized* scent genes in these species, cross-
referenced to our B. rapa routes. (Anchor genes web-verified this session: PhODO1 /
Verdonk 2005, AmNES/LIS / Nagegowda 2008, RhNUDX1 / Magnard 2015; supporting citations are
indicative — verify before manuscript use.)

**Cross-lineage synthesis — the same "tailoring" node keeps defining scent:**
- The **methyltransferase / acyltransferase tailoring step is the scent-defining node in
  every lineage**: Petunia **PhBSMT** (methylbenzoate), snapdragon **AmBAMT** (methylbenzoate),
  rose **RhOMT** (3,5-dimethoxytoluene) — all SABATH/OMT or BAHD enzymes decorating a
  conserved precursor. This is the **same family** that is most Brassica-expanded in Phase 4
  (COMT1 29 copies, BSMT1/SABATH 13) **and** that contains `Bra029041`, the one gene that
  moved under radiation. Across rosids (Brassica, rose) and asterids (Petunia, snapdragon),
  scent identity is set at this labile decorating step — a coherent, cross-clade thread.
- **The same volatile can be made by non-orthologous enzymes** — rose makes geraniol via a
  Nudix hydrolase (**RhNUDX1**), not the usual terpene synthase. So **orthology alone will
  miss functionally convergent scent genes**; a scent comparison must combine orthology
  (Phase 4) with pathway/enzyme-function curation (this table).
- Phylogenetic spread is good: rosids (B. rapa, rose) + asterids (Petunia, snapdragon)
  bracket the eudicots, so conserved-vs-labile calls are not an artefact of one clade.

**To make this programmatic (future Phase 4c):** pull proteomes from Sol Genomics
(Petunia axillaris/inflata), the *Antirrhinum majus* genome (Li 2019), and *Rosa chinensis*
(Raymond 2018), and run reciprocal-best-hit / OrthoFinder against the B. rapa scent set.

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
- The scent gene set **validates** on the scent axis, across all four routes: DE between
  High/Low lines = terpenoid 5/34, ester 13/117, apocarotenoid 4/15, benzenoid 4/71,
  GLV 4/42 (Tier1+2 overall 28/258; Tier1 core 11/72). The set captures real scent
  variation.
- **Scent × radiation (tiered):** the **Tier-1 core volatile enzymes** (terpene
  synthases, SABATH, LOX, CCD) are **0/55 among radiation DEGs** — no radiation response.
  In the broader Tier1+2 set, only **1/232** is radiation-DE: `Bra029041`, an
  O-methyltransferase (radiation log2FC +1.8), fold ≈10× but **Fisher p = 0.094 — not
  significant**.

**Defensible statement:** the repository now establishes the join and a validated,
tiered scent gene set. On the current in-house radiation counts there is **no significant
preferential perturbation of floral-scent biosynthesis genes** at 40 cGy — and notably
the *core* volatile-forming enzymes are entirely unresponsive — with the broad-tier
O-methyltransferase `Bra029041` as the single, non-significant follow-up candidate. The
test is under-powered by the weak overall radiation response, so this is a *"no
significant effect, under-powered"* result. Re-running against the collaborators'
radiation counts (and higher doses, if available) is the logical next step (Phase 3/4).
