# Legacy / archived files

These files were moved out of the active repository during the 2026-07 tidy-up.
They are **preserved on disk under `docs/legacy/`** but are **not tracked in git**
(so they stay out of the clean deposit). If you want any of them back in the
repository, run `git add -f docs/legacy/<file>`.

| Archived file | Why archived | Canonical replacement |
|---|---|---|
| `01-RNAseq_processing_Brappa_test-checkpoint.ipynb` | Stale Jupyter autosave; mixed mouse/Brapa reference chain (RSeQC step still points at `Mus_musculus.GRCm39`) | `01-RNAseq_processing_Brapa_FPsc_v1_3.ipynb` |
| `02-RNAseq_analysis_UPDATE_needed.ipynb` | Un-ported predecessor; still the mouse OSD-104 template (`org.Mm.eg.db`, `ENSMUSG…`, FLT vs GC) | `02-RNAseq_analysis_Brapa.ipynb` |
| `02-RNAseq_analysis_Brapa_executed.ipynb` | Near-identical duplicate of the canonical notebook (no saved outputs either) | `02-RNAseq_analysis_Brapa.ipynb` |
| `02-RNAseq_analysis_Brapa.ipynb.bak` | Editor backup | `02-RNAseq_analysis_Brapa.ipynb` |
| `brapa_all_genes.tsv` | **Wrong organism** — contains *Bradyrhizobium* genes (`bra:BRADO…`, dnaA/gyrB/recF), not *Brassica rapa*. Anything derived from it is suspect. | Rebuild from KEGG `bra` / Ensembl Plants (Phase 1 of `README_PLAN.md`) |
| `test_fixes.py`, `validate_fixes.md` | Jules/AI bug-fix scratch notes from an earlier session | — |

See `README_PLAN.md` for the consolidation plan these moves are part of.
