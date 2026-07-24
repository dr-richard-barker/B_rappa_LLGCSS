#!/usr/bin/env python3
"""
Phase 4 — cross-species conservation of the floral-scent gene set (Ensembl orthology).

Ensembl Plants Compara does not include the Chiifu `Bra######` assembly our expression
data uses (only `brassica_rapa_ro18`). So conservation is anchored on *Arabidopsis*
(the reference for floral-scent biology, and present in Compara): for a curated panel of
canonical Arabidopsis scent genes — one per family across the four biosynthetic routes —
we pull orthologs across all plant genomes via the Compara REST homology endpoint and
tabulate ortholog COUNT per species (conservation + lineage-specific copy-number).

The B. rapa link: `brassica_rapa_ro18` counts are the B. rapa copy numbers; our Chiifu
scent set (annotation/scent_geneset.tsv) contains the corresponding family members.

Curated Arabidopsis floral-scent panel (sources: Dudareva 2013 New Phytol; Chen 2011
Plant J [TPS]; Effmert 2005 [SABATH]; Tholl 2006; Vogel 2010 [BSMT]; Auldridge 2006 [CCD]).
Endpoint: https://rest.ensembl.org/homology/id/{sp}/{gene}?compara=plants (Accept: json)
"""
import json, os, time, urllib.request, urllib.error
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "ortho_cache"); os.makedirs(CACHE, exist_ok=True)

PANEL = [  # (AT id, gene, route)
    ("AT5G23960", "TPS21",  "terpenoid"),
    ("AT2G24210", "TPS10",  "terpenoid"),
    ("AT4G16740", "TPS03",  "terpenoid"),
    ("AT1G61120", "TPS04",  "terpenoid"),
    ("AT2G37040", "PAL1",   "benzenoid_phenylpropanoid"),
    ("AT5G54160", "COMT1",  "benzenoid_phenylpropanoid"),
    ("AT3G11480", "BSMT1",  "ester(SABATH)"),
    ("AT1G19640", "JMT",    "ester(SABATH)"),
    ("AT3G45140", "LOX2",   "fatty_acid_GLV"),
    ("AT4G15440", "HPL1",   "fatty_acid_GLV"),
    ("AT3G63520", "CCD1",   "apocarotenoid"),
    ("AT4G19170", "CCD4",   "apocarotenoid"),
]
# species of interest to highlight (exact Compara production names; anchor
# arabidopsis_thaliana is the query species so is not among its own orthologs)
FOCUS = ["arabidopsis_lyrata", "arabidopsis_halleri", "brassica_rapa_ro18",
         "brassica_oleracea", "brassica_napus", "brassica_juncea",
         "solanum_lycopersicum_gca000188115v5cm", "solanum_tuberosum",
         "vitis_vinifera", "oryza_sativa"]
TOMATO = "solanum_lycopersicum_gca000188115v5cm"

def fetch(at):
    fp = os.path.join(CACHE, at + ".json")
    if os.path.exists(fp) and os.path.getsize(fp) > 50:
        return json.load(open(fp))
    url = (f"https://rest.ensembl.org/homology/id/arabidopsis_thaliana/{at}"
           f"?compara=plants;type=orthologues;format=condensed")
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                d = json.loads(r.read().decode())
            json.dump(d, open(fp, "w")); return d
        except urllib.error.HTTPError as e:
            if e.code == 429: time.sleep(2 + attempt * 2); continue
            return {"error": f"HTTP {e.code}"}
        except Exception as e:
            time.sleep(1 + attempt)
    return {"error": "failed"}

rows = []          # (gene, route, species -> count, types)
species_all = set()
for at, gene, route in PANEL:
    d = fetch(at)
    homs = (d.get("data") or [{}])[0].get("homologies", []) if "error" not in d else []
    by_sp = Counter(h["species"] for h in homs)
    types = Counter(h.get("type", "?") for h in homs)
    species_all |= set(by_sp)
    rows.append((gene, at, route, by_sp, types, len(homs)))
    print(f"{gene:7s} {at}  {route:26s} orthologs={len(homs):3d} in {len(by_sp):2d} species "
          f"[{','.join(f'{k}:{v}' for k,v in types.items())}]")

# ---- write conservation matrix (gene x focus-species -> ortholog count) ----
with open(os.path.join(HERE, "scent_orthology_matrix.tsv"), "w") as fo:
    fo.write("gene\tAT_id\troute\ttotal_orthologs\tn_species\t" + "\t".join(FOCUS) + "\n")
    for gene, at, route, by_sp, types, tot in rows:
        cells = "\t".join(str(by_sp.get(s, 0)) for s in FOCUS)
        fo.write(f"{gene}\t{at}\t{route}\t{tot}\t{len(by_sp)}\t{cells}\n")

# ---- conservation summary --------------------------------------------------
print("\n=== conservation summary ===")
print(f"plant genomes seen across panel: {len(species_all)}")
for gene, at, route, by_sp, types, tot in rows:
    focus_present = [s.split('_')[0][:3]+s.split('_')[-1][:3] for s in FOCUS if by_sp.get(s)]
    brassica = sum(by_sp.get(s,0) for s in FOCUS if s.startswith("brassica"))
    deep = "yes" if (by_sp.get(TOMATO) or by_sp.get("solanum_tuberosum") or by_sp.get("oryza_sativa")) else "no"
    print(f"  {gene:7s} {route:26s} Brassica copies={brassica:2d}  conserved-to-tomato/rice={deep}")
print("\nwrote annotation/scent_orthology_matrix.tsv (+ raw JSON in annotation/ortho_cache/)")
