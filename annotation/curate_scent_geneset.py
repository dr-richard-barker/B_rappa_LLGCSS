#!/usr/bin/env python3
"""
Phase 2 — curated, tiered floral-scent (volatile) biosynthesis gene set for B. rapa.

Refines the Phase-1 GO/Pfam flag into a structured set with a biosynthetic ROUTE and a
confidence TIER, so downstream tests can run on a specific (Tier 1) or sensitive
(Tier 1+2) definition. Over-broad terms from v1 (e.g. GO:0008299 "isoprenoid
biosynthetic process", which also captures sterols / photosynthetic carotenoids /
ubiquinone) are demoted or dropped.

Gene-family sources (floral VOC biosynthesis):
  - Dudareva et al. 2013, New Phytologist 198:16-32 (VOC biosynthesis review)
  - Pichersky & Gershenzon 2002; Muhlemann et al. 2014 PCE (benzenoids)
  - Chen et al. 2011 Plant J. (terpene synthase / TPS family)
  - D'Auria 2006 (BAHD acyltransferases); Effmert et al. 2005 (SABATH methyltransferases)

Routes: terpenoid | benzenoid_phenylpropanoid | fatty_acid_GLV | apocarotenoid | ester
Inputs (git-ignored, regenerable): annotation/brapa_pfam.tsv, annotation/brapa_go.tsv,
                                    annotation/id_crosswalk.tsv
Output: annotation/scent_geneset.tsv  (overwrites the Phase-1 version, now tiered)
"""
import csv, os, re
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)

def norm(g):
    m = re.fullmatch(r'(?i)bra(\d+)', g.strip().strip('"')); return f"Bra{m.group(1)}" if m else None

# ---- rules -----------------------------------------------------------------
# Tier 1: core volatile-forming enzyme families (Pfam = high structural specificity)
PFAM_T1 = {
    "PF01397": ("terpenoid", "terpene synthase (TPS, N-term)"),
    "PF03936": ("terpenoid", "terpene synthase (TPS, C-term)"),
    "PF03492": ("ester", "SABATH methyltransferase (methyl-benzoate/-salicylate/-jasmonate)"),
    "PF00305": ("fatty_acid_GLV", "lipoxygenase (LOX; green-leaf-volatile route)"),
    "PF03055": ("apocarotenoid", "carotenoid cleavage dioxygenase (CCD; caveat: incl. NCED/ABA)"),
}
# Tier 2: supporting / route-level (broader; medium specificity)
PFAM_T2 = {
    "PF02458": ("ester", "BAHD acyltransferase (volatile esters; broad family)"),
    "PF00891": ("benzenoid_phenylpropanoid", "O-methyltransferase (benzenoid volatiles; also lignin)"),
    "PF00221": ("benzenoid_phenylpropanoid", "phenylalanine ammonia-lyase (PAL; phenylpropanoid entry)"),
}
GO_T2 = {
    "GO:0016114": ("terpenoid", "terpenoid biosynthetic process"),
    "GO:0016102": ("terpenoid", "diterpenoid biosynthetic process"),
    "GO:0051762": ("terpenoid", "sesquiterpene biosynthetic process"),
    "GO:0016099": ("terpenoid", "monoterpenoid biosynthetic process"),
    "GO:0010333": ("terpenoid", "terpene synthase activity"),
    "GO:0010334": ("terpenoid", "sesquiterpene synthase activity"),
    "GO:0009699": ("benzenoid_phenylpropanoid", "phenylpropanoid biosynthetic process"),
    "GO:0008171": ("benzenoid_phenylpropanoid", "O-methyltransferase activity"),
    "GO:0009694": ("fatty_acid_GLV", "jasmonic acid metabolic process"),
    "GO:0009695": ("fatty_acid_GLV", "jasmonic acid biosynthetic process"),
    "GO:0010436": ("apocarotenoid", "carotenoid dioxygenase activity"),
    "GO:0080031": ("ester", "methyl salicylate esterase activity"),
    "GO:0080032": ("ester", "methyl jasmonate esterase activity"),
}

# ---- load annotation -------------------------------------------------------
pfam_by = defaultdict(set); go_by = defaultdict(set)
with open(os.path.join(HERE, "brapa_pfam.tsv")) as fh:
    r = csv.reader(fh, delimiter='\t'); next(r)
    for row in r:
        if len(row) >= 2 and row[1] and norm(row[0]): pfam_by[norm(row[0])].add(row[1])
with open(os.path.join(HERE, "brapa_go.tsv")) as fh:
    r = csv.reader(fh, delimiter='\t'); next(r)
    for row in r:
        if len(row) >= 2 and row[1] and norm(row[0]): go_by[norm(row[0])].add(row[1])
present = {}
with open(os.path.join(HERE, "id_crosswalk.tsv")) as fh:
    r = csv.DictReader(fh, delimiter='\t')
    for row in r: present[row["bra_id"]] = (row["in_scent_data"], row["in_radiation_data"])

# ---- classify --------------------------------------------------------------
rows = []
for g in sorted(present):
    routes, ev, tier = set(), [], None
    for pf in pfam_by.get(g, ()):
        if pf in PFAM_T1:
            tier = 1; routes.add(PFAM_T1[pf][0]); ev.append(f"T1 Pfam:{pf} {PFAM_T1[pf][1]}")
    for pf in pfam_by.get(g, ()):
        if pf in PFAM_T2:
            routes.add(PFAM_T2[pf][0]); ev.append(f"T2 Pfam:{pf} {PFAM_T2[pf][1]}"); tier = tier or 2
    for go in go_by.get(g, ()):
        if go in GO_T2:
            routes.add(GO_T2[go][0]); ev.append(f"T2 {go} {GO_T2[go][1]}"); tier = tier or 2
    if tier:
        ins, inr = present[g]
        rows.append([g, "BRA" + g[3:], ins, inr, tier, ";".join(sorted(routes)), " | ".join(ev)])

with open(os.path.join(HERE, "scent_geneset.tsv"), "w", newline='') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(["bra_id", "BRA_id", "in_scent_data", "in_radiation_data", "tier", "routes", "evidence"])
    w.writerows(rows)

# ---- summary ---------------------------------------------------------------
from collections import Counter
tot = len(rows); t1 = sum(1 for r in rows if r[4] == 1)
both = sum(1 for r in rows if r[2] == "1" and r[3] == "1")
print(f"curated scent genes: {tot}  (Tier1 core: {t1}, Tier2: {tot-t1})")
print(f"  present in both experiments: {both}")
route_ct = Counter();
for r in rows:
    for rt in r[5].split(";"): route_ct[rt] += 1
print("by route:")
for rt, c in route_ct.most_common(): print(f"  {rt:28s} {c}")
print("\nby route x tier1:")
rt1 = Counter()
for r in rows:
    if r[4] == 1:
        for rt in r[5].split(";"): rt1[rt] += 1
for rt, c in rt1.most_common(): print(f"  {rt:28s} {c}")
print("\nwrote annotation/scent_geneset.tsv")
