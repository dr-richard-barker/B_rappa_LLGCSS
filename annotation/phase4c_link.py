#!/usr/bin/env python3
"""
Phase 4c — programmatic linkage of model-species scent genes to the B. rapa scent set.

No alignment toolchain is available in this sandbox (no BLAST/DIAMOND/OrthoFinder), and
Petunia/snapdragon/rose are absent from Ensembl Plants. So this script does the part that
IS executable and grounded:
  1. Fetch REAL protein accessions + sequences for the landmark model-species scent genes
     from NCBI (E-utilities) — upgrades the Phase-4b table from names to verifiable records
     and produces a query FASTA for the user's OrthoFinder/DIAMOND run.
  2. Bridge each to its biosynthetic route/family and to the candidate B. rapa orthologs =
     the members of that family in our scent set (annotation/scent_geneset.tsv).
The full reciprocal-ortholog run is documented as a protocol in annotation/phase4c_protocol.md.
"""
import urllib.request, urllib.parse, json, os, time, re
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# landmark genes: (label, route, ncbi_protein_search_term)
QUERIES = [
    ("PhBSMT (Petunia)",     "ester",                     'Petunia hybrida[Organism] AND benzoic acid carboxyl methyltransferase'),
    ("PhBPBT (Petunia)",     "ester",                     'Petunia hybrida[Organism] AND benzoyl-CoA benzyl alcohol benzoyltransferase'),
    ("PhODO1 (Petunia)",     "regulator",                 'Petunia hybrida[Organism] AND ODORANT1'),
    ("AmBAMT (snapdragon)",  "ester",                     'Antirrhinum majus[Organism] AND benzoic acid carboxyl methyltransferase'),
    ("AmNES/LIS (snapdragon)","terpenoid",                'Antirrhinum majus[Organism] AND nerolidol linalool synthase'),
    ("RhNUDX1 (rose)",       "terpenoid",                 'Rosa[Organism] AND nudix hydrolase RhNUDX1'),
    ("RhOMT (rose)",         "benzenoid_phenylpropanoid", 'Rosa[Organism] AND orcinol O-methyltransferase'),
]

def get(url):
    for a in range(4):
        try:
            with urllib.request.urlopen(url, timeout=45) as r: return r.read().decode()
        except Exception: time.sleep(1 + a)
    return ""

def esearch(term):
    u = f"{EUTILS}/esearch.fcgi?db=protein&retmax=1&retmode=json&term=" + urllib.parse.quote(term)
    try: return (json.loads(get(u)).get("esearchresult", {}).get("idlist") or [None])[0]
    except Exception: return None

def efetch(uid):
    fa = get(f"{EUTILS}/efetch.fcgi?db=protein&id={uid}&rettype=fasta&retmode=text")
    if not fa.startswith(">"): return None, None, None
    hdr = fa.splitlines()[0][1:]
    seq = "".join(fa.splitlines()[1:])
    acc = hdr.split()[0]
    return acc, hdr, seq

# ---- fetch ----
recs = []
for label, route, term in QUERIES:
    uid = esearch(term); time.sleep(0.4)
    if not uid:
        print(f"{label:26s} route={route:26s} NO NCBI HIT for query"); recs.append((label, route, "", "no hit", 0, "")); continue
    acc, hdr, seq = efetch(uid); time.sleep(0.4)
    if not acc:
        print(f"{label:26s} efetch failed uid={uid}"); recs.append((label, route, "", "efetch failed", 0, "")); continue
    print(f"{label:26s} route={route:26s} -> {acc}  len={len(seq)}  {hdr[:60]}")
    recs.append((label, route, acc, hdr, len(seq), seq))

# ---- write query FASTA + accession table ----
with open(os.path.join(HERE, "scent_query_proteins.faa"), "w") as ff:
    for label, route, acc, hdr, n, seq in recs:
        if seq: ff.write(f">{acc} {label} | route={route} | {hdr}\n{seq}\n")
with open(os.path.join(HERE, "scent_query_accessions.tsv"), "w") as ft:
    ft.write("gene\troute\tncbi_accession\tprotein_length\tdefline\n")
    for label, route, acc, hdr, n, seq in recs:
        ft.write(f"{label}\t{route}\t{acc}\t{n}\t{hdr}\n")

# ---- bridge to B. rapa scent set (candidate orthologs = same-route members) ----
route_members = Counter(); route_t1 = Counter()
with open(os.path.join(HERE, "scent_geneset.tsv")) as fh:
    import csv
    for row in csv.DictReader(fh, delimiter='\t'):
        if row["in_radiation_data"] == "1" or row["in_scent_data"] == "1":
            for rt in row["routes"].split(";"):
                route_members[rt] += 1
                if row["tier"] == "1": route_t1[rt] += 1
print("\n=== gene -> route -> candidate B. rapa scent-set orthologs (same route) ===")
for label, route, acc, hdr, n, seq in recs:
    if route == "regulator":
        print(f"  {label:26s} regulator (MYB) — not in enzyme scent set; compare via TF orthology"); continue
    print(f"  {label:26s} {route:26s} candidates: {route_members.get(route,0)} B. rapa genes "
          f"({route_t1.get(route,0)} Tier-1 core)")
print("\nNote: Bra029041 (radiation-moved O-MT) is an ester/benzenoid-MT candidate ortholog of "
      "PhBSMT / AmBAMT / RhOMT — the cross-lineage tailoring node.")
print("wrote annotation/scent_query_proteins.faa + scent_query_accessions.tsv")
