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
- `scent_radiation_test.py` — tier-aware Fisher enrichment (Phase 1/2 preview).
- `phase3_geneset_test.py` — Phase 3: powered rank-based gene-set test + confounds.
- `scent_orthology.py` — Phase 4: cross-species conservation via Ensembl Compara.
- `scent_geneset.tsv` — the 363 curated scent genes (tier + route).
- `scent_orthology_matrix.tsv` — scent gene × species ortholog-count matrix.
- `scent_reference_species.tsv` — Phase 4b: literature-curated scent genes in
  Petunia / snapdragon / rose (species absent from Ensembl Plants).
- `phase4c_link.py` + `scent_query_proteins.faa` + `scent_query_accessions.tsv` —
  Phase 4c: real NCBI protein records for the landmark genes + family bridge.
- `phase4c_protocol.md` — Phase 4c: reproducible OrthoFinder/DIAMOND recipe.

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

## Phase 4c — making 4b programmatic (executed part + protocol)

No alignment toolchain (BLAST/DIAMOND/OrthoFinder) exists in this sandbox, so Phase 4c is
split into the part that is **executable and grounded here** and a **reproducible protocol**
for the user's compute (`phase4c_protocol.md`).

**Executed (`phase4c_link.py`):** fetched **real NCBI/UniProt protein records** for the
seven landmark model-species scent genes (no invented accessions) → `scent_query_proteins.faa`
+ `scent_query_accessions.tsv`:

| gene | route | accession |
|---|---|---|
| PhBSMT (Petunia) | ester | `A4ZDG8` (BSMT3_PETHY) |
| PhBPBT (Petunia) | ester | `Q6E593` (BEBT1_PETHY) |
| PhODO1 (Petunia) | regulator | `Q50EX6` (ODO1_PETHY) |
| AmBAMT (snapdragon) | ester | `Q9FYZ9` (BAMT_ANTMA) |
| AmNES/LIS (snapdragon) | terpenoid | `ABR24418` |
| RhNUDX1 (rose) | terpenoid | `M4I1C6` (NUDT1_ROSHC) |
| RhOMT (rose) | benzenoid | `AAM23005` |

The script then bridges each to its **candidate B. rapa orthologs** = the scent-set members
of the same route (ester → 174 B. rapa genes / 39 Tier-1; terpenoid → 52 / 39; benzenoid →
97). This is the honest link at family resolution; single-gene orthology across ~120 My
needs the alignment step below.

**Protocol (`phase4c_protocol.md`):** grounded proteome sources (Rosa `GCF_002994745.2`
RefSeq; Petunia via Sol Genomics; Antirrhinum via Li 2019 portal; Chiifu B. rapa via Ensembl
Plants) + OrthoFinder / DIAMOND-RBH commands + the intersection step. **Prediction to test:**
the ester-MT queries (PhBSMT/AmBAMT/RhOMT) should land in a B. rapa SABATH/O-MT orthogroup
containing `Bra029041`; RhNUDX1 (a Nudix hydrolase) should have **no** terpene-synthase
ortholog — a built-in positive control for pathway convergence.

## Reproduce

```bash
# 1. Pull GO + Pfam from Ensembl Plants BioMart (dataset brapa_eg_gene = Brapa_1.0):
#    (see the curl commands in build_annotation.py header / git history)
# 2. Build the annotation layer + scent set:
python3 annotation/build_annotation.py
# 3. Preliminary radiation × scent test:
python3 annotation/scent_radiation_test.py
```

## Does GCR radiation influence floral scent? — Phase 3 powered test

*Exploratory — our radiation counts are an in-house re-analysis (the primary radiation
study is collaborator-led); treat as hypothesis-generating, not a final result.*

The Phase-1 preview used a Fisher test on the 13 hard DEGs — badly under-powered.
`phase3_geneset_test.py` instead uses a **rank-based competitive gene-set test**
(Mann-Whitney U of the scent set's radiation dose effect vs the genome-wide background,
iDEP interaction model, main 40-vs-0 cGy term, 31,008 genes) — every gene's effect size
counts, no DEG cutoff.

- **Core volatile enzymes (Tier-1: TPS, SABATH, LOX, CCD) show no shift** (signed-log2FC
  MWU p = 0.40). The enzymes that directly *make* the scent are unresponsive to 40 cGy.
- **The broad scent set (Tier-1+2, n=287) shows a weak but nominally significant
  down-shift** relative to background (signed-log2FC MWU **p = 0.032**, z = −2.14; scent
  median ≈ 0 vs background +0.13 — i.e. scent genes do *not* join the mild genome-wide
  up-shift at 40 cGy). Magnitude (|log2FC|) is unchanged (p = 0.95), so scent genes are
  not *more* variable, just directionally lower.
- **The signal is concentrated in the ester / methyltransferase route** (SABATH/BAHD/OMT,
  n=142): median −0.11, **MWU p = 0.012** — the same *tailoring* family that the
  comparative genomics (Phase 4/4b) flagged as the cross-lineage scent node. Terpenoid,
  benzenoid, GLV and apocarotenoid routes show nothing (p ≥ 0.58).
- **Confounds are clean:** 0 scent genes have a significant dose×genotype or
  dose×preservative interaction — the (weak) scent response does not depend on the
  *anthocyaninless* background or the preservative.
- **Dual-hit candidates** (radiation-responsive *and* scent-associated High/Low):
  `Bra013161` (GLV; dose +2.5, scent adjP 0.01), `Bra028224` and `Bra039555` (ester).
  `Bra029041` is dose-**up** (+1.6) — an outlier *against* its route's overall down-trend.

**Defensible statement:** at 40 cGy the core scent-forming enzymes are unaffected, but the
broad scent set — specifically the **ester/methyltransferase tailoring route** — shows a
weak, directional, nominally significant (uncorrected; borderline after multiple-testing)
**down-shift** relative to the genomic radiation response, landing exactly where the
cross-species conservation analysis predicted the labile scent node to be. This is a
small, hypothesis-generating effect, not a headline result. The decisive test is to re-run
`phase3_geneset_test.py` on the **collaborators' radiation counts** (higher power, and
ideally higher doses), for which the ID join and gene set are already prepared.
