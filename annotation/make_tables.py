#!/usr/bin/env python3
"""
Build manuscript Tables 1-3 from the committed data, and emit:
  docs/TABLES.md            markdown (manuscript / GitHub)
  annotation/tables/*.csv   machine-readable
  docs/results.html         a Pages results page (figure gallery + Tables 1-3)
Re-run after phase3_robust.py / curate_scent_geneset.py / scent_orthology.py.
"""
import csv, os, re, html
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
DOCS = os.path.join(ROOT, "docs"); TBL = os.path.join(HERE, "tables")
os.makedirs(TBL, exist_ok=True)
def norm(g):
    m = re.fullmatch(r'(?i)bra(\d+)', g.strip().strip('"')); return f"Bra{m.group(1)}" if m else None
def f(x, p=2):
    try: return f"{float(x):+.{p}f}"
    except: return ""

ROUTES = ["terpenoid","benzenoid_phenylpropanoid","ester","fatty_acid_GLV","apocarotenoid"]
RLAB = {"terpenoid":"Terpenoid","benzenoid_phenylpropanoid":"Benzenoid / phenylpropanoid",
        "ester":"Ester / methyltransferase","fatty_acid_GLV":"Fatty-acid / GLV",
        "apocarotenoid":"Apocarotenoid"}

# ---- load ----
sg = list(csv.DictReader(open(os.path.join(HERE,"scent_geneset.tsv")), delimiter='\t'))
tier = {r["bra_id"]:int(r["tier"]) for r in sg}
routes = {r["bra_id"]:r["routes"].split(";") for r in sg}
evid = {r["bra_id"]:r["evidence"] for r in sg}
dose = {}
if os.path.exists(os.path.join(HERE,"phase3_dose_stats.tsv")):
    for r in csv.DictReader(open(os.path.join(HERE,"phase3_dose_stats.tsv")), delimiter='\t'):
        dose[r["bra_id"]] = float(r["dose_effect"])
scentDE = {}
rows = list(csv.reader(open(os.path.join(ROOT,"NewTest/differential_expression_GLbulkRNAseq (1).csv"))))
h=rows[0]; gi=h.index("gene_id"); li=h.index("Log2fc_(High)v(Low)"); pi=h.index("Adj.p.value_(High)v(Low)")
for rr in rows[1:]:
    g=norm(rr[gi])
    if g:
        try: scentDE[g]=(float(rr[li]), float(rr[pi]))
        except: pass
robust = {r["name"]:r for r in csv.DictReader(open(os.path.join(HERE,"phase3_robust_results.tsv")), delimiter='\t')}

# ================= TABLE 1: scent gene-set composition =================
t1 = [["Biosynthetic route","Tier 1 (core enzymes)","Tier 2 (supporting)","Total","Example families"]]
EXAMP = {"terpenoid":"terpene synthases (TPS)","benzenoid_phenylpropanoid":"PAL, O-methyltransferases",
         "ester":"SABATH MTs, BAHD acyltransferases","fatty_acid_GLV":"lipoxygenase (LOX), HPL",
         "apocarotenoid":"carotenoid cleavage dioxygenases (CCD)"}
for rt in ROUTES:
    a=sum(1 for g in routes if rt in routes[g] and tier[g]==1)
    b=sum(1 for g in routes if rt in routes[g] and tier[g]==2)
    t1.append([RLAB[rt], str(a), str(b), str(a+b), EXAMP[rt]])
tot1=sum(1 for g in tier if tier[g]==1); tot2=sum(1 for g in tier if tier[g]==2)
t1.append(["**All routes**", f"**{tot1}**", f"**{tot2}**", f"**{len(tier)}**", "—"])

# ================= TABLE 2: radiation x scent candidate genes =================
# scent-set genes that are scent-associated (High/Low adjP<0.1), ranked by |radiation dose effect|
cand=[]
for g in routes:
    sd=scentDE.get(g); d=dose.get(g)
    if sd and sd[1] is not None and sd[1]<0.1 and d is not None:
        cand.append((g,d,sd))
cand.sort(key=lambda x: abs(x[1]), reverse=True)
t2=[["Gene (Bra)","Route","Tier","Radiation dose effect (40 vs 0 cGy, log2)","Scent High/Low log2FC","Scent adjP","Enzyme evidence"]]
for g,d,sd in cand[:12]:
    t2.append([g, RLAB[routes[g][0]], str(tier[g]), f(d), f(sd[0]), f"{sd[1]:.2g}", evid[g].split('|')[0].strip()[:42]])

# ================= TABLE 3: characterised scent genes in model species =================
acc={}
for r in csv.DictReader(open(os.path.join(HERE,"scent_query_accessions.tsv")), delimiter='\t'):
    key=r["gene"].split()[0]  # e.g. 'PhBSMT'
    acc[key]=r["ncbi_accession"]
t3=[["Species","Gene","Enzyme / role","Route","Volatile product","Family","NCBI/UniProt","Reference"]]
for r in csv.reader(open(os.path.join(HERE,"scent_reference_species.tsv")), delimiter='\t'):
    if not r or r[0].startswith("#") or r[0]=="species": continue
    sp,clade,gene,enz,route,vol,fam,ref = r[:8]
    a = next((acc[k] for k in acc if gene.split('/')[0].replace('1','').replace('2','') in k or k in gene), "")
    t3.append([sp,gene,enz,route.replace('_',' '),vol,fam,a,ref])

# ---- writers ----
def md_table(t):
    out=["| "+" | ".join(t[0])+" |", "|"+"|".join("---" for _ in t[0])+"|"]
    for row in t[1:]: out.append("| "+" | ".join(row)+" |")
    return "\n".join(out)
def csv_write(path,t):
    with open(path,"w",newline='') as fo:
        w=csv.writer(fo)
        for row in t: w.writerow([c.replace('**','') for c in row])
def html_table(t, cap):
    th="".join(f"<th>{html.escape(c)}</th>" for c in t[0])
    body=""
    for row in t[1:]:
        cells="".join(f"<td>{html.escape(c).replace('**','')}</td>" for c in row)
        body+=f"<tr>{cells}</tr>"
    return (f'<p style="font-weight:600;margin:1.6rem 0 .3rem">{html.escape(cap)}</p>'
            f'<table><thead><tr>{th}</tr></thead><tbody>{body}</tbody></table>')

CAP1="Table 1. Composition of the curated floral-scent gene set in B. rapa (363 genes)."
CAP2="Table 2. Candidate genes: scent-associated (High vs Low, adjP<0.1) and their radiation dose effect. Note: the scent set as a whole is NOT significantly radiation-responsive after correction (see Fig 3); these are individual candidates, not a validated hit list."
CAP3="Table 3. Characterised floral-scent genes in model scent species absent from Ensembl Plants (literature-curated; NCBI/UniProt accessions verified where fetched)."

with open(os.path.join(DOCS,"TABLES.md"),"w") as fo:
    fo.write("# Tables\n\n_Generated by `annotation/make_tables.py` from the committed data._\n\n")
    fo.write("## "+CAP1+"\n\n"+md_table(t1)+"\n\n## "+CAP2+"\n\n"+md_table(t2)+"\n\n## "+CAP3+"\n\n"+md_table(t3)+"\n")
csv_write(os.path.join(TBL,"table1_geneset_composition.csv"),t1)
csv_write(os.path.join(TBL,"table2_candidates.csv"),t2)
csv_write(os.path.join(TBL,"table3_model_species.csv"),t3)

# ---- supplementary tables: S1 full scent set, S2 full ortholog matrix ----
def tsv_to_csv(src, dst):
    with open(src) as fi, open(dst,"w",newline='') as fo:
        w=csv.writer(fo)
        for row in csv.reader(fi, delimiter='\t'): w.writerow(row)
tsv_to_csv(os.path.join(HERE,"scent_geneset.tsv"), os.path.join(TBL,"TableS1_scent_geneset_full.csv"))
tsv_to_csv(os.path.join(HERE,"scent_orthology_matrix.tsv"), os.path.join(TBL,"TableS2_ortholog_matrix.csv"))

# ---- results.html (Pages) ----
FIGS=[("figures/Fig1_design.png","Fig 1. Two-experiment design and the Bra↔BRA gene-ID join."),
      ("figures/Fig2_geneset.png","Fig 2. (A) The 363-gene scent set by route and tier. (B) Scent-axis validation: % of each route's genes differentially expressed between High- and Low-scent lines."),
      ("figures/Fig3_radiation_test.png","Fig 3. Radiation dose effect (40 vs 0 cGy) by scent route vs genomic background. No route is significant after correlation-aware permutation and BH correction; the ester/methyltransferase route gives the strongest (still non-significant) trend."),
      ("figures/Fig4_conservation.png","Fig 4. Cross-species conservation of an Arabidopsis-anchored scent panel (Ensembl Compara ortholog counts). Note the Brassicaceae-specific expansion of COMT1 and the terpene-synthase expansions in grape/tomato."),
      ("figures/FigS1_power.png","Fig S1. Sensitivity of the gene-set test: power to detect a coordinated shift δ in the scent set's dose response. The set-level test is well-powered (~3–4% detectable at 80% power); observed effects (dashed) fall below, so the null is informative, not merely under-powered.")]
gallery="".join(f'<figure><img src="{s}" alt="{html.escape(c)}"><figcaption>{html.escape(c)}</figcaption></figure>' for s,c in FIGS)
page=f"""<!doctype html><html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Results · B. rapa floral scent under radiation</title>
<meta name="description" content="Figures and tables for the B. rapa scent-under-radiation framework. Preliminary, non-significant at 40 cGy after correction.">
<link rel="stylesheet" href="assets/style.css"></head><body>
<header class="site"><nav><span class="brand">🌸 B. rapa · Scent × Radiation</span>
<a href="index.html">Home</a><a href="#figures">Figures</a><a href="#tables">Tables</a>
<a href="https://github.com/dr-richard-barker/B_rappa_LLGCSS">GitHub ↗</a></nav></header>
<main><div class="hero" style="padding-bottom:.6rem">
<h1>Results — figures &amp; tables</h1>
<p class="lede">Draft manuscript figures and tables, generated from the committed analysis
(<code>annotation/</code>). All values regenerate via the scripts.</p></div>
<p class="note"><b>Read the effect honestly.</b> At 40 cGy there is <b>no statistically significant
effect</b> of radiation on the scent gene set after correlation-aware permutation and
multiple-testing correction (Fig 3). The tables list a validated gene set and individual
candidates — not a confirmed radiation-responsive hit list.</p>
<div class="content">
<h2 id="figures">Figures</h2>
<div class="gallery">{gallery}</div>
<h2 id="tables">Tables</h2>
{html_table(t1,CAP1)}
{html_table(t2,CAP2)}
{html_table(t3,CAP3)}
<p class="note">Sources: <code>scent_geneset.tsv</code>, <code>phase3_robust_results.tsv</code>,
<code>scent_orthology_matrix.tsv</code>, <code>scent_reference_species.tsv</code>,
<code>scent_query_accessions.tsv</code>. Markdown: <a href="TABLES.md"><code>docs/TABLES.md</code></a>;
CSV: <a href="https://github.com/dr-richard-barker/B_rappa_LLGCSS/tree/main/annotation/tables"><code>annotation/tables/</code></a>.
Supplementary (Fig S1, Tables S1–S2): <a href="SUPPLEMENTARY.md"><code>docs/SUPPLEMENTARY.md</code></a>;
reproduce all: <a href="https://github.com/dr-richard-barker/B_rappa_LLGCSS/blob/main/annotation/reproduce.sh"><code>annotation/reproduce.sh</code></a>.</p>
</div></main>
<footer class="site">B. rapa — floral scent under space radiation (Lunar LEAF) · work in progress · CC0 ·
<a href="https://github.com/dr-richard-barker/B_rappa_LLGCSS">GitHub</a></footer></body></html>"""
with open(os.path.join(DOCS,"results.html"),"w") as fo: fo.write(page)

print("wrote docs/TABLES.md, docs/results.html, annotation/tables/*.csv")
print(f"Table 1 rows: {len(t1)-1}  Table 2 candidates: {len(t2)-1}  Table 3 rows: {len(t3)-1}")
