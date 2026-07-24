#!/usr/bin/env python3
"""
Phase 1 — gene-ID reconciliation + annotation layer for B_rappa_LLGCSS.

Key empirical finding this script is built on:
  The scent dataset uses  Bra######  gene IDs (Chiifu / Ensembl Plants Brapa_1.0)
  The radiation dataset uses BRA######  gene IDs (same loci, uppercased).
  -> after case-normalisation, all 31,756 scent genes map 1:1 onto radiation genes.
  So the two experiments JOIN DIRECTLY on gene ID; no cross-assembly mapping needed.

Inputs (in repo):
  NewTest/differential_expression_GLbulkRNAseq (1).csv   scent DE (Bra IDs)
  Radiation model/Results_LFC_Pval_DESeq2.csv            radiation DE (BRA IDs)
  annotation/brapa_go.tsv     Bra -> GO   (pulled from Ensembl Plants BioMart)
  annotation/brapa_pfam.tsv   Bra -> Pfam (pulled from Ensembl Plants BioMart)

Outputs (annotation/):
  id_crosswalk.tsv       canonical Bra id <-> BRA id, presence flags
  gene_annotation.tsv    per-gene GO/Pfam + scent-candidate flag
  scent_geneset.tsv      the curated floral-volatile biosynthesis candidate genes

Scent gene set is defined in B. rapa space (this species is not well served by
Arabidopsis-ortholog tables here) using curated GO terms + Pfam domains for the
four canonical floral-volatile routes: terpenoid, benzenoid/phenylpropanoid,
fatty-acid/green-leaf-volatile, and apocarotenoid, plus SABATH volatile esterases.
"""
import csv, re, os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# ---- curated floral-volatile biosynthesis gene set -------------------------
# GO accessions (specific; deliberately EXCLUDES histone/RNA/protein methyltransferases)
SCENT_GO = {
    # terpenoid / isoprenoid
    "GO:0016114": "terpenoid biosynthetic process",
    "GO:0008299": "isoprenoid biosynthetic process",
    "GO:0010333": "terpene synthase activity",
    "GO:0016102": "diterpenoid biosynthetic process",
    "GO:0016099": "monoterpenoid biosynthetic process",
    "GO:0051762": "sesquiterpene biosynthetic process",
    "GO:0010334": "sesquiterpene synthase activity",
    "GO:0016104": "triterpenoid biosynthetic process",
    "GO:0046246": "terpene biosynthetic process",
    # benzenoid / phenylpropanoid
    "GO:0009699": "phenylpropanoid biosynthetic process",
    "GO:0009698": "phenylpropanoid metabolic process",
    "GO:0045551": "cinnamyl-alcohol dehydrogenase activity",
    # fatty-acid derived / green leaf volatiles / jasmonate
    "GO:0009695": "jasmonic acid biosynthetic process",
    "GO:0009694": "jasmonic acid metabolic process",
    "GO:0016702": "oxidoreductase activity, acting on single donors (lipoxygenase)",
    # SABATH volatile esters (methyl benzoate / salicylate / jasmonate)
    "GO:0008171": "O-methyltransferase activity",
    "GO:0080031": "methyl salicylate esterase activity",
    "GO:0080032": "methyl jasmonate esterase activity",
    # apocarotenoid volatiles
    "GO:0010436": "carotenoid dioxygenase activity",
}
# Pfam domains for scent enzymes (high-confidence structural evidence)
SCENT_PFAM = {
    "PF01397": "Terpene_synth (terpene synthase N)",
    "PF03936": "Terpene_synth_C (terpene synthase C)",
    "PF00305": "Lipoxygenase",
    "PF00891": "Methyltransf_2 (O-MT, incl. SABATH-adjacent)",
    "PF02458": "Transferase (BAHD acyltransferase, e.g. BEAT/BEBT)",
}

def norm(gid):
    """Canonical gene id: capitalise to 'Bra######' regardless of source case."""
    g = gid.strip().strip('"')
    m = re.fullmatch(r'(?i)bra(\d+)', g)
    return f"Bra{m.group(1)}" if m else None

# ---- gene id sets from the two experiments ---------------------------------
def ids_from(path, col):
    out = set()
    with open(os.path.join(ROOT, path), newline='') as fh:
        r = csv.reader(fh); h = next(r); j = h.index(col)
        for row in r:
            n = norm(row[j]) if j < len(row) else None
            if n: out.add(n)
    return out

scent_ids = ids_from("NewTest/differential_expression_GLbulkRNAseq (1).csv", "gene_id")
rad_ids   = ids_from("Radiation model/Results_LFC_Pval_DESeq2.csv", "ensembl_ID")
all_ids   = sorted(scent_ids | rad_ids)

# ---- functional annotation -------------------------------------------------
go_by  = defaultdict(list)   # bra -> [(go, name, domain)]
with open(os.path.join(HERE, "brapa_go.tsv")) as fh:
    r = csv.reader(fh, delimiter='\t'); next(r)
    for row in r:
        if len(row) >= 4 and row[1]:
            n = norm(row[0])
            if n: go_by[n].append((row[1], row[2], row[3]))
pfam_by = defaultdict(set)
with open(os.path.join(HERE, "brapa_pfam.tsv")) as fh:
    r = csv.reader(fh, delimiter='\t'); next(r)
    for row in r:
        if len(row) >= 2 and row[1]:
            n = norm(row[0])
            if n: pfam_by[n].add(row[1])

# ---- write crosswalk + annotation + scent set ------------------------------
os.makedirs(HERE, exist_ok=True)
scent_rows = []
with open(os.path.join(HERE, "id_crosswalk.tsv"), "w", newline='') as fx, \
     open(os.path.join(HERE, "gene_annotation.tsv"), "w", newline='') as fa:
    wx = csv.writer(fx, delimiter='\t'); wa = csv.writer(fa, delimiter='\t')
    wx.writerow(["bra_id", "BRA_id", "in_scent_data", "in_radiation_data"])
    wa.writerow(["bra_id", "BRA_id", "in_scent_data", "in_radiation_data",
                 "n_go", "go_ids", "pfam_ids", "is_scent_candidate", "scent_evidence"])
    for g in all_ids:
        BRA = "BRA" + g[3:]
        in_s = g in scent_ids; in_r = g in rad_ids
        wx.writerow([g, BRA, int(in_s), int(in_r)])
        gos = go_by.get(g, [])
        go_ids = sorted({x[0] for x in gos})
        pfams = sorted(pfam_by.get(g, set()))
        ev = []
        hit_go = [SCENT_GO[x] for x in go_ids if x in SCENT_GO]
        hit_pf = [SCENT_PFAM[x] for x in pfams if x in SCENT_PFAM]
        if hit_go: ev += [f"GO:{n}" for n in hit_go]
        if hit_pf: ev += [f"Pfam:{n}" for n in hit_pf]
        is_scent = bool(ev)
        wa.writerow([g, BRA, int(in_s), int(in_r), len(go_ids),
                     ";".join(go_ids), ";".join(pfams),
                     int(is_scent), " | ".join(ev)])
        if is_scent:
            scent_rows.append([g, BRA, int(in_s), int(in_r), " | ".join(ev)])

with open(os.path.join(HERE, "scent_geneset.tsv"), "w", newline='') as fs:
    w = csv.writer(fs, delimiter='\t')
    w.writerow(["bra_id", "BRA_id", "in_scent_data", "in_radiation_data", "scent_evidence"])
    w.writerows(sorted(scent_rows))

# ---- summary ---------------------------------------------------------------
print(f"canonical genes (scent ∪ radiation): {len(all_ids)}")
print(f"  in scent data:      {len(scent_ids)}")
print(f"  in radiation data:  {len(rad_ids)}")
print(f"  in BOTH:            {len(scent_ids & rad_ids)}")
print(f"genes with GO:   {len(go_by)}   genes with Pfam: {len(pfam_by)}")
print(f"scent-candidate genes: {len(scent_rows)}")
print(f"  ...present in both experiments: {sum(1 for r in scent_rows if r[2] and r[3])}")
print("wrote annotation/{id_crosswalk,gene_annotation,scent_geneset}.tsv")
