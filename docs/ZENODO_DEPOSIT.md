# Zenodo deposit guide

Prepared for R. Barker. **Nothing here mints a DOI automatically** — the deposit is
triggered by *you* publishing a GitHub release (or uploading to Zenodo) under your own
account. This file is the checklist + the exact steps.

## Recommended route: GitHub → Zenodo release archiving (code + derived tables)

This captures the repository (code, small tables, docs) and mints a versioned DOI.

1. **One-time:** log in to <https://zenodo.org> with GitHub, go to **Account → GitHub**,
   and flip the switch **ON** for `dr-richard-barker/B_rappa_LLGCSS`.
2. **Confirm metadata is present** (already in the repo, so Zenodo reads it):
   - `.zenodo.json` — title, description, creators (ORCID), keywords, license CC0-1.0.
   - `CITATION.cff` — citation metadata (GitHub "Cite this repository" button).
3. **Publish a release:**
   ```bash
   git tag -a v1.0.0 -m "Scent-under-radiation framework: phases 0-4c"
   git push origin v1.0.0
   # then on GitHub: Releases → Draft a new release → choose tag v1.0.0 → Publish
   ```
   Zenodo auto-archives the release and issues a DOI (plus a version-independent "concept" DOI).
4. **Back-fill the DOI** into `README.md`, `CITATION.cff` (`doi:`), and `.zenodo.json`
   (`related_identifiers`), then commit. (Optional: add the Zenodo DOI badge to the README.)

## What is / isn't in this deposit

- **In (git-tracked):** all analysis code (`annotation/*.py`, R scripts, notebooks), the
  curated small tables (`scent_geneset.tsv`, `scent_orthology_matrix.tsv`,
  `scent_reference_species.tsv`, `scent_query_*`), docs, figures.
- **Not in git (see `.gitignore` / `DATA.md`):** reference genomes, STAR/RSEM indexes, and
  the large derived DE / count / annotation tables. For a **data** deposit of those, either:
  - point to the primary sources (SRA `SRR4417237–244`; radiation OSDR accession), **or**
  - make a **separate Zenodo data upload** of the large tables and cross-link the two DOIs
    (`related_identifiers`: software `isSupplementedBy` data).

## Pre-deposit checklist (Zenodo readiness §6, updated)

- [x] `README.md` rewritten to match real contents
- [x] `CITATION.cff` present with ORCID `0000-0001-5681-9857`
- [x] `.zenodo.json` present (authors, keywords, license)
- [x] `LICENSE` present (CC0-1.0)
- [x] Large binaries handled (`.gitignore` + `DATA.md`); history purged (241→98 MB)
- [x] `environment.yml` present
- [x] Data provenance documented (`DATA.md`; scent SRA accessions)
- [x] One canonical notebook/script per step; legacy archived (`docs/LEGACY.md`)
- [ ] **Radiation OSDR accession/DOI** filled into `DATA.md`, `.zenodo.json`, README — *needs you*
- [ ] **Confirm collaborator authorship / consent** to include the radiation re-analysis publicly
- [ ] Strip leftover mouse/OSD-104 text from the kept notebooks (cosmetic; pre-submission)
- [ ] Decide: combined repo, or split scent (deposit) vs radiation (collaborator) — see plan §5.6
- [ ] Tag `v1.0.0` and publish the release → DOI; then back-fill the DOI

## Before you publish — two things only you can decide
1. **Radiation data sharing.** The radiation experiment is collaborator-led. Either include the
   in-house counts/results (with their sign-off) or make the deposit scent-focused and cite their
   OSDR DOI. This choice drives what goes in the archive.
2. **Authorship.** Confirm the author list and contributions before minting a citable DOI.
