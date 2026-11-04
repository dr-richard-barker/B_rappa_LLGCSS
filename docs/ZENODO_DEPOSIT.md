# Zenodo Deposit Guide: Brassica rapa Lunar LEAF Radiation & Redox Decoder

Prepared for R. Barker. **Nothing here mints a DOI automatically** — the deposit is triggered by *you* publishing a GitHub release (or uploading to Zenodo) under your own account. This file is the checklist + the exact steps.

## Recommended Route: GitHub → Zenodo Release Archiving (Code + Derived Tables)

This captures the repository (code, small tables, docs) and mints a versioned DOI.

1. **One-time:** log in to <https://zenodo.org> with GitHub, go to **Account → GitHub**, and flip the switch **ON** for `dr-richard-barker/B_rappa_LLGCSS`.
2. **Confirm metadata is present** (already in the repo, so Zenodo reads it):
   - `.zenodo.json` — title, description, creators (ORCID), keywords, license CC0-1.0.
   - `CITATION.cff` — citation metadata (GitHub "Cite this repository" button).
3. **Publish a release:**
   ```bash
   git tag -a v1.0.0 -m "Lunar LEAF B. rapa cosmic radiation and redox decoding release"
   git push origin v1.0.0
   # then on GitHub: Releases → Draft a new release → choose tag v1.0.0 → Publish
   ```
   Zenodo auto-archives the release and issues a DOI.
4. **Back-fill the DOI** into `README.md`, `CITATION.cff` (`doi:`), and `.zenodo.json` (`related_identifiers`), then commit.

## What is in this Deposit

- **In (git-tracked):** all analysis code, orthology tables (18,770 RBH pairs), derived statistics, docs, figures (Fig 1–6).
- **Primary Data:** NASA Open Science Data Repository (OSDR) Lunar LEAF candidate investigation in *Brassica rapa*.
