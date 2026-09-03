#!/usr/bin/env python3
"""Rebuild this repo's ortholog and KEGG annotation from source.

Two faults are fixed here, both of which silently produced plausible-looking output.

**1. The KEGG annotation was for the wrong organism.** The KEGG organism code used
throughout was `bra`, which is *Bradyrhizobium* sp. ORS 278 -- a soil bacterium -- not
*Brassica rapa*. KEGG's own record settles it:

    $ curl https://rest.kegg.jp/info/bra
    bra   Bradyrhizobium sp. ORS 278

The correct code is `brp`. The error is visible in the artefacts once you look: the gene
symbols in `brapa_symbol_to_kegg_id_map.csv` were `dnaA`, `dnaN`, `recF` (bacterial
replication genes), the ids were `BRADO#####`, and the pathway list included "Microbial
metabolism in diverse environments". Every downstream pathway result built on those files
was describing a bacterium.

**2. The ortholog file was an error message.** `brapa_to_arabidopsis_orthologs.tsv` was
100 bytes containing

    Query ERROR: caught BioMart::Exception::Usage:
    Attribute athaliana_eg_homolog_ensembl_gene NOT FOUND

Ensembl Plants exposes no Arabidopsis-homolog attribute for *B. rapa*, so that query
cannot be made to work. The replacement comes from reciprocal best hit (DIAMOND, TAIR10
proteome), computed in the sibling arabidopsis-drem-osdr repository and copied in here.

**Correcting the organism is necessary but not sufficient.** KEGG `brp` is keyed by NCBI
GeneID (`brp:103844700`) while this repo's RSEM matrices are indexed by Ensembl Plants
`Bra######`; the two namespaces do not intersect at all, and a join test against the
repo's own counts returns 0. `brp` is also almost entirely symbol-free -- 495 of 44,411
genes carry one -- so a `symbol` column built from it degenerates to LOC numbers.

KEGG `ath` has neither problem: it is keyed on AGI codes, which is exactly what the
ortholog map emits, and its genes carry real symbols. So `Bra###### -> AT#G##### ->
pathway` joins with no bridging step and every link is checkable.

The files the analysis scripts read are therefore built on the Arabidopsis route, because
that is the one that reaches the data. The `brp` files are still written, under
`annotation/kegg/`, so the repo holds correct-organism reference data and nothing anywhere
still points at a bacterium.

  python3 annotation/fix_kegg_and_orthologs.py

(Distinct from `build_annotation.py`, which is Phase 1: the Bra/BRA crosswalk,
GO/Pfam pull and scent gene set. This script only rebuilds the ortholog map and
the KEGG annotation.)
"""
from __future__ import annotations

import csv
import datetime as dt
import json
import re
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ANN = ROOT / "annotation"
KEGG_DIR = ANN / "kegg"
META = ROOT / "Brapa_analysis" / "Metadata"
RBH_SOURCE = (Path.home() / "Documents" / "arabidopsis-drem-osdr" / "data" / "brapa"
              / "brapa_arabidopsis_rbh.tsv")
REST = "https://rest.kegg.jp"
AGI = re.compile(r"^AT[1-5CM]G\d{5}$", re.I)


def get(path: str) -> str:
    url = f"{REST}/{path}"
    with urllib.request.urlopen(url, timeout=180) as r:
        return r.read().decode("utf-8", "replace")


def rows(text: str) -> list[list[str]]:
    return [ln.split("\t") for ln in text.strip("\n").split("\n") if ln]


def write_tsv(path: Path, header: list[str] | None, data) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        if header:
            w.writerow(header)
        for r in data:
            w.writerow(r)
            n += 1
    return n


def main() -> int:
    qc: dict = {"generated": dt.date.today().isoformat(), "source": REST}

    # ---------------------------------------------------------- 1. orthologs
    if not RBH_SOURCE.exists():
        sys.exit(f"ortholog source not found: {RBH_SOURCE}\n"
                 "Run scripts/26_brapa_orthologs.py in arabidopsis-drem-osdr first.")
    pairs = []
    with RBH_SOURCE.open() as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5:
                # Restore the repo's own capitalisation: counts use `Bra000001`.
                pairs.append([f[0].capitalize(), f[1].upper(), f[2], f[3], f[4]])
    pairs.sort()
    n_orth = write_tsv(ANN / "brapa_to_arabidopsis_orthologs.tsv",
                       ["brapa_gene", "arabidopsis_agi", "identity_pct", "evalue",
                        "bitscore"], pairs)
    bra2at = {p[0].upper(): p[1] for p in pairs}
    print(f"orthologs                     {n_orth:>7} pairs")
    qc["orthologs"] = {"n_pairs": n_orth,
                       "method": "DIAMOND reciprocal best hit vs TAIR10, identity >=50%",
                       "replaces": "a 100-byte file containing a BioMart error"}

    # ---------------------------------------------------------- 2. KEGG brp
    brp_genes = rows(get("list/brp"))
    brp_paths = rows(get("list/pathway/brp"))
    brp_link = rows(get("link/pathway/brp"))

    # Same four-column shape the old (wrong-organism) brapa_all_genes.tsv had, so any
    # existing reader keeps working -- process_genes.py included.
    write_tsv(ROOT / "brapa_all_genes.tsv", None, brp_genes)
    write_tsv(KEGG_DIR / "brp_pathway_names.tsv", ["pathway_id", "pathway_name"],
              [[p[0].replace("path:", ""), p[1].split(" - ")[0]] for p in brp_paths])
    write_tsv(KEGG_DIR / "brp_gene_to_pathway.tsv", ["pathway_id", "ncbi_gene_id"],
              [[l[1].replace("path:", ""), l[0].split(":")[1]] for l in brp_link])
    n_sym = sum(1 for g in brp_genes
                if len(g) > 3 and ";" in g[3] and " " not in g[3].split(";")[0])
    print(f"brapa_all_genes.tsv           {len(brp_genes):>7} genes  (was Bradyrhizobium)")
    print(f"annotation/kegg/brp_*.tsv     {len(brp_paths):>7} pathways, "
          f"{len(brp_link)} links, {n_sym} genes with a symbol")
    qc["kegg_brp"] = {"organism": "brp (Brassica rapa, field mustard)",
                      "replaced": "bra (Bradyrhizobium sp. ORS 278)",
                      "n_genes": len(brp_genes), "n_pathways": len(brp_paths),
                      "n_links": len(brp_link),
                      "id_space": "NCBI GeneID; joins 0 of this repo's Bra###### counts",
        "n_genes_with_symbol": n_sym,
        "kept_as": "annotation/kegg/brp_*.tsv (reference only)"}

    # ------------------------------------------- 3. the join that works: via Arabidopsis
    ath_paths = rows(get("list/pathway/ath"))
    ath_link = rows(get("link/pathway/ath"))
    ath_genes = rows(get("list/ath"))

    at_symbol = {}
    for g in ath_genes:
        agi = g[0].split(":")[1].upper()
        desc = g[3] if len(g) > 3 else ""
        first = desc.split(";")[0].strip() if ";" in desc else ""
        at_symbol[agi] = first.split(",")[0].strip() if first else agi

    by_at: dict[str, list[str]] = {}
    for l in ath_link:
        by_at.setdefault(l[0].split(":")[1].upper(), []).append(l[1].replace("path:", ""))

    joined, covered = [], set()
    for bra, at in sorted(bra2at.items()):
        for pw in by_at.get(at, []):
            joined.append([pw, bra.capitalize(), at, at_symbol.get(at, at)])
            covered.add(bra)
    write_tsv(ANN / "brapa_gene_to_kegg_via_arabidopsis.tsv",
              ["pathway_id", "brapa_gene", "arabidopsis_agi", "symbol"], joined)

    # These two are what brapa_gsea.R reads for TERM2GENE/TERM2NAME. They must be keyed
    # by the ids the DGE tables actually carry, which are Bra######, not NCBI GeneIDs.
    path_names = [[p[0].replace("path:", ""), p[1].split(" - ")[0]] for p in ath_paths]
    write_tsv(META / "brapa_kegg_pathway_names.tsv", None, path_names)
    write_tsv(META / "brapa_kegg_link.tsv", None, [[r[0], r[1]] for r in joined])
    write_tsv(KEGG_DIR / "ath_pathway_names.tsv", ["pathway_id", "pathway_name"],
              path_names)

    # `kegg_id,symbol` schema preserved for brapa_gsea.R, which joins on $symbol and
    # reads $kegg_id. The join keys are added as extra columns; R ignores what it does
    # not ask for. Built on ath because brp carries essentially no symbols.
    seen, sym_rows = set(), []
    for bra, at in sorted(bra2at.items()):
        if at in seen:
            continue
        seen.add(at)
        sym_rows.append([f"ath:{at}", at_symbol.get(at, at), bra.capitalize(), at])
    with (ROOT / "brapa_symbol_to_kegg_id_map.csv").open("w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["kegg_id", "symbol", "brapa_gene", "arabidopsis_agi"])
        w.writerows(sym_rows)

    real_sym = sum(1 for r in sym_rows if not AGI.match(r[1]))
    print(f"via-Arabidopsis join          {len(joined):>7} links "
          f"({len(covered)} B. rapa genes, {len(ath_paths)} pathways)")
    print(f"brapa_kegg_link.tsv           {len(joined):>7} links   (Bra-keyed, usable)")
    print(f"brapa_symbol_to_kegg_id_map   {len(sym_rows):>7} rows   "
          f"({real_sym} with a real symbol)")
    qc["kegg_via_arabidopsis"] = {
        "n_links": len(joined), "n_brapa_genes": len(covered),
        "n_pathways": len(ath_paths),
        "why": "KEGG ath is keyed on AGI codes, which the ortholog map already emits, so "
               "Bra###### -> AT#G##### -> pathway joins with no bridging step. KEGG brp "
               "is keyed on NCBI GeneID and cannot join to this repo's counts.",
    }

    (ANN / "annotation_qc.json").write_text(json.dumps(qc, indent=2) + "\n")
    print(f"\nwrote {ANN / 'annotation_qc.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
