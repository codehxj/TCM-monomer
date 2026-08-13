from __future__ import annotations

import re
import textwrap
import unicodedata
import html
import time
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import requests
import seaborn as sns


ROOT = Path(__file__).resolve().parent
PAPER = ROOT / "论文"
FIG = PAPER / "figures"

COLS = {
    "compound": "单体",
    "source": "来源",
    "model": "实验模型",
    "pathway": "作用通路",
    "effect": "核心生物学效应",
    "target": "主要靶点（或备注上升下降反应）",
    "clinical": "临床试验数据",
    "toxicity": "主要副作用",
    "pmid": "参考文献（PubMed ID）",
}

CATEGORY = {
    "1,4,6-tri-o-galloyl-beta-d-glucopyranose": "Hydrolysable tannins",
    "amygdalin": "Cyanogenic glycosides",
    "astragalin": "Flavonoids",
    "astragaloside iv": "Triterpenoid saponins",
    "beta-lapachone": "Quinones",
    "beta-sitosterol": "Phytosterols",
    "caffeic acid": "Phenolic acids",
    "catechin": "Flavonoids",
    "chrysophanol": "Anthraquinones",
    "curdione": "Terpenoids",
    "dictamnine": "Alkaloids",
    "ellagic acid": "Polyphenols",
    "emodin aloe-emodin": "Anthraquinones",
    "eucalyptol": "Terpenoids",
    "formononetin": "Flavonoids",
    "glabridin": "Flavonoids",
    "hederagenin": "Triterpenoids",
    "honokiol": "Neolignans",
    "isorhamnetin": "Flavonoids",
    "isovitexin": "Flavonoids",
    "jaranol": "Flavonoids",
    "kaempferol": "Flavonoids",
    "luteolin": "Flavonoids",
    "macrocarpal i": "Meroterpenoids",
    "magnolol": "Neolignans",
    "myricetin": "Flavonoids",
    "naamidine j": "Alkaloids",
    "neothalfine": "Alkaloids",
    "obovatol": "Neolignans",
    "pachymic acid": "Triterpenoids",
    "pinocembrin": "Flavonoids",
    "rhein": "Anthraquinones",
    "sanguinarine": "Alkaloids",
    "secoemestrin c": "Fungal alkaloid-like metabolites",
    "silibinin": "Flavonoids",
    "stigmasterol": "Phytosterols",
    "syringic acid": "Phenolic acids",
    "tagitinin c": "Sesquiterpene lactones",
    "tanshinone iia": "Diterpenoid quinones",
    "wighteone": "Flavonoids",
}

PATHWAYS = {
    "PI3K/Akt/mTOR": r"PI3K|Akt|mTOR|p70S6K",
    "Wnt/beta-catenin": r"Wnt|catenin|β-catenin|LEF|TCF",
    "NF-kB/inflammation": r"NF|炎症|inflamm|TNF|IL-|COX|PTGS|NLRP",
    "STAT3": r"STAT3",
    "MAPK/JNK/ERK": r"MAPK|JNK|ERK|p38",
    "Ferroptosis/GPX4": r"ferropt|铁死亡|GPX4|Nrf2|SLC7A11|GSH|lipid peroxid",
    "Apoptosis/caspases": r"apopt|凋亡|caspase|Bax|Bcl|PARP",
    "EMT/ECM/MMP": r"EMT|ECM|MMP|侵袭|迁移|metasta",
    "ROS/oxidative stress": r"ROS|氧化|oxidative|mitochond",
    "Cell cycle": r"cell cycle|细胞周期|Cyclin|CDK|G2|G1",
    "Autophagy": r"autophagy|自噬|LC3|Beclin",
    "Immune microenvironment": r"PD-L1|macrophage|巨噬|CD8|immune|免疫|T cell",
    "Gut microbiota": r"microbiota|菌群|肠道微生物",
}

HIGHLIGHTS = {"luteolin", "kaempferol", "tanshinone iia", "astragaloside iv"}


def norm(x: object) -> str:
    s = str(x or "").strip()
    s = unicodedata.normalize("NFKD", s)
    s = s.replace("β", "beta").replace("Ⅱ", "II")
    s = re.sub(r"[\s_]+", " ", s)
    s = re.sub(r"[^0-9a-zA-Z,.\- ]+", " ", s)
    return re.sub(r"\s+", " ", s).strip().lower()


def clean(x: object) -> str:
    if pd.isna(x):
        return ""
    return re.sub(r"\s+", " ", str(x).strip())


def latex(x: object) -> str:
    s = clean(x)
    s = s.replace("β", "beta").replace("α", "alpha").replace("κ", "kappa")
    s = s.replace("↓", " down").replace("↑", " up").replace("²", "2")
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode("ascii")
    repl = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    return "".join(repl.get(ch, ch) for ch in s)


def short(s: object, n: int) -> str:
    s = clean(s)
    return s if len(s) <= n else s[: n - 3].rstrip() + "..."


def pmids(x: object) -> list[str]:
    return re.findall(r"\b\d{6,9}\b", str(x or ""))


def model_class(s: str) -> str:
    t = s.lower()
    if re.search(r"网络药理|分子对接|in silico|bioinform|comput|molecular docking", t):
        return "Network/in silico"
    if re.search(r"小鼠|大鼠|aom|dss|xenograft|nude|裸鼠|animal|mouse|mice|rat", t):
        return "In vivo"
    if re.search(r"hct|ht-29|sw480|sw620|caco|ct26|mc38|rko|lovo|cell|细胞", t):
        return "In vitro"
    return "Not specified" if not t else "Other experimental evidence"


def hits(row: pd.Series) -> list[str]:
    s = " ".join(str(row.get(c, "")) for c in ["pathway", "effect", "target"])
    return [k for k, p in PATHWAYS.items() if re.search(p, s, re.I)]


def cells(s: str) -> list[str]:
    pats = [r"HCT[- ]?116", r"HT[- ]?29", r"SW480", r"SW620", r"Caco[- ]?2", r"LoVo", r"RKO", r"CT26", r"MC38", r"T84", r"DLD[- ]?1"]
    out = []
    for p in pats:
        out += re.findall(p, s, re.I)
    return sorted(set(x.upper().replace(" ", "-") for x in out))


def doses(s: str) -> list[str]:
    pat = r"\b\d+(?:\.\d+)?\s?(?:nM|uM|µM|μM|mM|mg/kg|mg\/kg|mg mL-1|mg/mL|ug/mL|µg/mL|μg/mL|mg·kg−1|mg kg-1)\b"
    return sorted(set(re.findall(pat, s, re.I)))


def load() -> pd.DataFrame:
    df = pd.read_excel(ROOT / "单体.xlsx", sheet_name="Sheet1")
    df = df[[v for v in COLS.values() if v in df.columns]].rename(columns={v: k for k, v in COLS.items()})
    for c in df.columns:
        df[c] = df[c].map(clean)
    df = df[df["compound"] != ""].drop_duplicates()
    df["key"] = df["compound"].map(norm)
    df["category"] = df["key"].map(CATEGORY).fillna("Other natural products")
    df["model_class"] = df["model"].map(model_class)
    df["pathway_hits"] = df.apply(hits, axis=1)
    df["cell_lines"] = df.apply(lambda r: ", ".join(cells(" ".join([r["model"], r["effect"], r["pathway"], r["target"]]))), axis=1)
    df["dose"] = df.apply(lambda r: ", ".join(doses(" ".join([r["model"], r["effect"], r["pathway"], r["target"]]))), axis=1)
    return df


def table_outputs(df: pd.DataFrame) -> None:
    PAPER.mkdir(exist_ok=True)
    out = df.rename(columns={
        "category": "化学类别", "compound": "单体", "source": "来源", "model": "实验模型",
        "pathway": "作用通路", "effect": "核心效应", "pmid": "参考文献",
        "cell_lines": "细胞系", "dose": "剂量",
    })[["化学类别", "单体", "来源", "实验模型", "细胞系", "剂量", "作用通路", "核心效应", "参考文献"]]
    out.sort_values(["化学类别", "单体", "参考文献"]).to_csv(PAPER / "monomer_table.csv", index=False, encoding="utf-8-sig")
    out.sort_values(["化学类别", "单体", "参考文献"]).to_excel(PAPER / "monomer_table.xlsx", index=False)


def figures(df: pd.DataFrame) -> None:
    FIG.mkdir(exist_ok=True)
    sns.set_theme(style="whitegrid", font="DejaVu Sans")
    cats = Counter(df.drop_duplicates("key")["category"])
    items = sorted(cats.items(), key=lambda x: (-x[1], x[0]))
    plt.figure(figsize=(8.8, 5))
    ax = sns.barplot(x=[v for _, v in items], y=[k for k, _ in items], color="#4C78A8")
    ax.bar_label(ax.containers[0], padding=3)
    ax.set_xlabel("Number of distinct monomers")
    ax.set_ylabel("Chemical category")
    ax.set_title("Distribution of TCM monomers by chemical category")
    plt.tight_layout()
    plt.savefig(FIG / "category_distribution.png", dpi=300)
    plt.close()

    models = Counter(df["model_class"])
    plt.figure(figsize=(7, 5.3))
    plt.pie(list(models.values()), labels=list(models.keys()), autopct="%1.1f%%", startangle=110, colors=sns.color_palette("Set2", len(models)))
    plt.title("Experimental model types")
    plt.tight_layout()
    plt.savefig(FIG / "model_type_pie.png", dpi=300)
    plt.close()

    ph = Counter(x for hs in df["pathway_hits"] for x in hs)
    top = ph.most_common(12)
    plt.figure(figsize=(9, 5.3))
    ax = sns.barplot(x=[v for _, v in top], y=[k for k, _ in top], color="#59A14F")
    ax.bar_label(ax.containers[0], padding=3)
    ax.set_xlabel("Number of records")
    ax.set_ylabel("Pathway/process")
    ax.set_title("High-frequency signaling pathways and phenotypic processes")
    plt.tight_layout()
    plt.savefig(FIG / "pathway_frequency.png", dpi=300)
    plt.close()
    network()


def network() -> None:
    reps = [
        ("Luteolin", "Flavonoids", "PI3K/Akt/mTOR; GPX4", "apoptosis, ferroptosis", "growth and immune restraint"),
        ("Kaempferol", "Flavonoids", "Wnt/beta-catenin; ER stress", "cell-cycle arrest", "chemosensitization"),
        ("Tanshinone IIA", "Diterpenoid quinones", "survivin; Ang2/Tie2; ferroptosis", "apoptosis, ferroptosis", "tumor inhibition"),
        ("Astragaloside IV", "Triterpenoid saponins", "PI3K/Akt/mTOR; circRNA/miRNA", "migration/invasion inhibition", "tumor suppression"),
        ("Honokiol", "Neolignans", "STAT3/NF-kB", "stemness and inflammation restraint", "sensitization"),
        ("Rhein", "Anthraquinones", "ROS/ER stress", "apoptosis", "growth inhibition"),
        ("Sanguinarine", "Alkaloids", "MAPK/JNK/ERK", "ferroptosis and migration inhibition", "anti-metastatic effect"),
        ("beta-sitosterol", "Phytosterols", "Wnt/beta-catenin", "proliferation inhibition", "growth inhibition"),
        ("Caffeic acid", "Phenolic acids", "NF-kB/inflammation", "inflammatory cytokine reduction", "CAC risk modulation"),
    ]
    colors = {"Flavonoids": "#4C78A8", "Diterpenoid quinones": "#E15759", "Triterpenoid saponins": "#59A14F", "Neolignans": "#B07AA1", "Anthraquinones": "#F28E2B", "Alkaloids": "#9C755F", "Phytosterols": "#76B7B2", "Phenolic acids": "#EDC948"}
    fig, ax = plt.subplots(figsize=(14, 8), dpi=180)
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    def box(x, y, text, fc, w=.18, h=.06, fs=8.5, bold=False):
        ax.add_patch(plt.Rectangle((x-w/2, y-h/2), w, h, facecolor=fc, edgecolor="#444", lw=1.1))
        ax.text(x, y, text, ha="center", va="center", fontsize=fs, weight="bold" if bold else "normal", wrap=True)
    ys = [0.88 - i * 0.085 for i in range(len(reps))]
    for (compound, cat, pathway, phenotype, outcome), y in zip(reps, ys):
        box(.12, y, compound, colors[cat], bold=norm(compound) in HIGHLIGHTS)
        box(.40, y, pathway, "#F2F5F9", w=.24, fs=8)
        box(.67, y, phenotype, "#F6F0E8", w=.22, fs=8)
        ax.annotate("", xy=(.28, y), xytext=(.21, y), arrowprops=dict(arrowstyle="->", lw=1))
        ax.annotate("", xy=(.55, y), xytext=(.52, y), arrowprops=dict(arrowstyle="->", lw=1))
        ax.annotate("", xy=(.84, .5), xytext=(.78, y), arrowprops=dict(arrowstyle="->", lw=0.9, color="#666"))
    box(.91, .5, "Anti-CRC effect\nreduced growth, invasion,\ninflammation and resistance", "#D9EAD3", w=.18, h=.13, fs=9, bold=True)
    for x, title in [(.12, "Monomers"), (.40, "Pathways"), (.67, "Cell phenotypes"), (.91, "Anti-CRC")]:
        ax.text(x, .965, title, ha="center", va="center", fontsize=11, weight="bold")
    for i, (cat, c) in enumerate(colors.items()):
        x = .05 + (i % 4) * .235
        y = .065 - (i // 4) * .035
        ax.add_patch(plt.Rectangle((x, y), .018, .018, facecolor=c))
        ax.text(x + .024, y + .009, cat, va="center", fontsize=8)
    ax.set_title("Monomer -> pathway -> cell phenotype -> anti-CRC effect", fontsize=14, weight="bold", pad=20)
    plt.tight_layout()
    plt.savefig(FIG / "mechanism_network.png", dpi=300, bbox_inches="tight")
    plt.close()


def cite_list(series: pd.Series, n=3) -> str:
    seen = []
    for v in series:
        for p in pmids(v):
            if p not in seen:
                seen.append(p)
    return ", ".join([fr"\cite{{PMID{x}}}" for x in seen[:n]]) or "Not indexed"


def evidence_sentence(sub: pd.DataFrame, max_rows=4) -> str:
    rows = []
    for _, r in sub.head(max_rows).iterrows():
        bits = []
        if r["model"]:
            bits.append("model: " + latex(short(r["model"], 110)))
        if r["cell_lines"]:
            bits.append("cell line(s): " + latex(r["cell_lines"]))
        if r["dose"]:
            bits.append("dose: " + latex(r["dose"]))
        else:
            bits.append("dose: not reported in the curated record")
        if r["pathway"]:
            bits.append("pathway: " + latex(short(r["pathway"], 110)))
        if r["effect"]:
            bits.append("effect: " + latex(short(r["effect"], 120)))
        cs = cite_list(pd.Series([r["pmid"]]), 1)
        rows.append("; ".join(bits) + f" {cs}.")
    return " ".join(rows)


def category_sections(df: pd.DataFrame) -> str:
    parts = [r"\section{Evidence by Chemical Category}"]
    order = [k for k, _ in Counter(df.drop_duplicates("key")["category"]).most_common()]
    for cat in order:
        sub = df[df["category"] == cat].copy()
        compounds = ", ".join([x for x, _ in Counter(sub["compound"]).most_common(8)])
        models = ", ".join([f"{k} ({v})" for k, v in Counter(sub["model_class"]).most_common()])
        paths = ", ".join([k for k, _ in Counter(x for hs in sub["pathway_hits"] for x in hs).most_common(5)]) or "pathways not consistently annotated"
        parts.append(fr"\subsection{{{latex(cat)}}}")
        parts.append(
            f"This category contained {sub['key'].nunique()} curated monomer(s), including {latex(compounds)}. "
            f"The extracted experimental context was {latex(models)}. Frequently annotated mechanisms were {latex(paths)}. "
            "The following statements are direct evidence summaries from the curated model, pathway, target and effect fields: "
            + evidence_sentence(sub.sort_values(["compound", "pmid"]), 5)
        )
    return "\n\n".join(parts)


def highlight_sections(df: pd.DataFrame) -> str:
    names = ["Luteolin", "Kaempferol", "Tanshinone IIA", "Astragaloside IV"]
    parts = [r"\section{Priority Monomer Evidence}"]
    for name in names:
        sub = df[df["key"] == norm(name)].copy()
        paths = ", ".join([k for k, _ in Counter(x for hs in sub["pathway_hits"] for x in hs).most_common(6)])
        cell = ", ".join([x for x, _ in Counter(", ".join(sub["cell_lines"]).split(", ")).most_common(6) if x])
        dose_vals = sorted({d for x in sub["dose"] for d in x.split(", ") if d})
        dose_text = ", ".join(dose_vals[:8]) if dose_vals else "not reported in the curated records"
        parts.append(fr"\subsection{{{latex(name)}}}")
        parts.append(
            f"{latex(name)} had {len(sub)} curated record(s). Extracted cell/model signals included {latex(cell or 'cell lines not consistently specified')}; "
            f"recorded dose information was {latex(dose_text)}. Recurrent mechanisms were {latex(paths or 'not consistently annotated')}. "
            + evidence_sentence(sub.sort_values("pmid"), 6)
        )
    return "\n\n".join(parts)


def monomer_table_tex(df: pd.DataFrame) -> str:
    grouped = []
    for key, sub in df.groupby("key"):
        grouped.append((
            sub["category"].iloc[0], sub["compound"].iloc[0],
            "; ".join([x for x, _ in Counter(sub["source"]).most_common(2) if x]),
            "; ".join([x for x, _ in Counter(sub["model_class"]).most_common(2) if x]),
            ", ".join([k for k, _ in Counter(x for hs in sub["pathway_hits"] for x in hs).most_common(3)]) or short("; ".join(sub["pathway"].head(2)), 95),
            short("; ".join([x for x, _ in Counter(sub["effect"]).most_common(2) if x]), 125),
            cite_list(sub["pmid"], 3),
        ))
    grouped.sort(key=lambda x: (x[0], x[1].lower()))
    out = [r"\begin{scriptsize}", r"\begin{longtable}{L{0.13\textwidth}L{0.16\textwidth}L{0.14\textwidth}L{0.20\textwidth}L{0.20\textwidth}L{0.12\textwidth}}", r"\caption{Monomers organized by chemical category. Highlighted rows denote priority monomers.}\label{tab:monomers}\\", r"\toprule", r"Monomer & Source & Experimental model & Pathway & Core effect & References \\", r"\midrule", r"\endfirsthead", r"\toprule", r"Monomer & Source & Experimental model & Pathway & Core effect & References \\", r"\midrule", r"\endhead"]
    cur = None
    for cat, comp, src, mod, path, eff, cites in grouped:
        if cat != cur:
            out.append(fr"\multicolumn{{6}}{{l}}{{\textbf{{{latex(cat)}}}}}\\")
            cur = cat
        prefix = r"\rowcolor{yellow!16} " if norm(comp) in HIGHLIGHTS else ""
        name = fr"\textbf{{{latex(comp)}}}" if norm(comp) in HIGHLIGHTS else latex(comp)
        out.append(f"{prefix}{name} & {latex(short(src or 'Not specified', 75))} & {latex(short(mod, 60))} & {latex(short(path, 85))} & {latex(short(eff, 110))} & {cites} \\\\")
    out += [r"\bottomrule", r"\end{longtable}", r"\end{scriptsize}"]
    return "\n".join(out)


def mechanism_table(df: pd.DataFrame) -> str:
    out = [r"\begin{scriptsize}", r"\begin{longtable}{L{0.18\textwidth}L{0.22\textwidth}L{0.38\textwidth}L{0.15\textwidth}}", r"\caption{Mechanism-oriented evidence map.}\label{tab:mechanisms}\\", r"\toprule", r"Mechanism & Monomers & Extracted model/effect evidence & References \\", r"\midrule", r"\endfirsthead", r"\toprule", r"Mechanism & Monomers & Extracted model/effect evidence & References \\", r"\midrule", r"\endhead"]
    for mech, _ in Counter(x for hs in df["pathway_hits"] for x in hs).most_common(12):
        sub = df[df["pathway_hits"].map(lambda hs: mech in hs)]
        mons = ", ".join([x for x, _ in Counter(sub["compound"]).most_common(6)])
        evidence = "; ".join([short(x, 80) for x, _ in Counter((sub["model_class"] + ": " + sub["effect"]).dropna()).most_common(3)])
        out.append(f"{latex(mech)} & {latex(short(mons, 105))} & {latex(short(evidence, 190))} & {cite_list(sub['pmid'], 4)} \\\\")
    out += [r"\bottomrule", r"\end{longtable}", r"\end{scriptsize}"]
    return "\n".join(out)


def write_tex(df: pd.DataFrame) -> None:
    cats = Counter(df.drop_duplicates("key")["category"])
    models = Counter(df["model_class"])
    ph = Counter(x for hs in df["pathway_hits"] for x in hs)
    intro_cites = cite_list(df["pmid"], 7)
    tex = rf"""
\documentclass[pdflatex]{{sn-jnl}}
\usepackage{{graphicx}}
\usepackage{{amsmath,amssymb}}
\usepackage{{booktabs}}
\usepackage{{array}}
\usepackage{{tabularx}}
\usepackage{{longtable}}
\usepackage{{float}}
\usepackage{{placeins}}
\usepackage[table]{{xcolor}}
\usepackage[numbers,sort&compress]{{natbib}}
\hypersetup{{hypertexnames=false}}
\newcolumntype{{L}}[1]{{>{{\raggedright\arraybackslash}}p{{#1}}}}
\raggedbottom

\begin{{document}}
\title[TCM monomers and colorectal cancer]{{Evidence-based review of traditional Chinese medicine monomers in colorectal cancer: extracted models, pathways, phenotypes and translational priorities}}
\author[1]{{Xianjun Han}}\email{{hxj@ahu.edu.cn}}
\author[1]{{Jincheng Fang}}\email{{E24301275@stu.ahu.edu.cn}}
\equalcont{{These authors contributed equally to this work.}}
\author[2]{{Zijian Wu}}\email{{wuzijian@ahtcm.edu.cn}}
\author*[2,3]{{Can Bai}}\email{{baican@ahtcm.edu.cn}}
\author[3]{{Renbao Huang}}\email{{hhmoxue@163.com}}
\affil[1]{{\orgdiv{{School of Computer Science and Technology}}, \orgname{{Anhui University}}, \orgaddress{{\street{{111 Jiulong Road}}, \city{{Hefei}}, \postcode{{230601}}, \state{{Anhui}}, \country{{China}}}}}}
\affil[2]{{\orgdiv{{School of Acupuncture and Tuina College}}, \orgname{{Anhui University of Chinese Medicine}}, \orgaddress{{\street{{350 Longzihu Road}}, \city{{Hefei}}, \postcode{{230012}}, \state{{Anhui}}, \country{{China}}}}}}
\affil[3]{{\orgdiv{{First Clinical Medical College}}, \orgname{{Anhui University of Chinese Medicine}}, \orgaddress{{\street{{350 Longzihu Road}}, \city{{Hefei}}, \postcode{{230012}}, \state{{Anhui}}, \country{{China}}}}}}

\abstract{{We curated {len(df)} records covering {df['key'].nunique()} natural monomers related to colorectal cancer (CRC). Evidence was extracted directly from the structured spreadsheet and the monomer-folder PMID organization, including experimental model, cell-line mentions, dose mentions when present, pathways, targets, biological effects and PubMed identifiers. The dataset spans {len(cats)} chemical categories. The most frequent pathway/effect annotations were {latex(', '.join(f'{k} ({v})' for k, v in ph.most_common(7)))}. Clinical trial information was absent from the curated clinical column; therefore, the review focuses on preclinical and computational evidence and marks missing dose information explicitly.}}
\keywords{{colorectal cancer, traditional Chinese medicine, natural monomer, PubMed, mechanism network, flavonoids, ferroptosis}}
\maketitle

\section{{Introduction}}
CRC therapy still faces metastasis, recurrence, acquired resistance and toxicity. TCM-derived monomers are attractive because they frequently modulate connected cancer phenotypes, including proliferation, apoptosis, ferroptosis, EMT, inflammation, oxidative stress and immune context. This manuscript does not infer unreported experiments. Instead, it extracts model descriptions, cell-line mentions, dose strings, pathways, core effects and PubMed identifiers from the curated data table and links them to the prepared monomer literature folders. Representative evidence is cited with PMID-based BibTeX keys such as {intro_cites}.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.98\textwidth]{{figures/mechanism_network.png}}
\caption{{Mechanism network summarizing extracted monomer-to-pathway-to-phenotype-to-anti-CRC relationships. Monomer colors indicate chemical categories.}}
\label{{fig:network}}
\end{{figure}}

\section{{Dataset Overview and Statistics}}
The curated monomer-level chemical distribution was {latex(', '.join(f'{k}: {v}' for k, v in cats.most_common()))}. Experimental model categories were {latex(', '.join(f'{k}: {v}' for k, v in models.most_common()))}. Dose strings were extracted by regular expression from model, pathway, effect and target fields; absence of a dose string is reported as not reported rather than estimated.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.88\textwidth]{{figures/category_distribution.png}}
\caption{{Number of distinct monomers by chemical category.}}
\label{{fig:category}}
\end{{figure}}

{monomer_table_tex(df)}

{category_sections(df)}

{highlight_sections(df)}

\section{{Mechanism-Level Synthesis}}
Across categories, the extracted annotations converged on apoptosis/caspases, NF-kB/inflammation, PI3K/Akt/mTOR, EMT/ECM/MMP, cell-cycle regulation, ROS/oxidative stress, MAPK/JNK/ERK and immune-microenvironment modules. Fig.~\ref{{fig:pathways}} shows the high-frequency signals, and Table~\ref{{tab:mechanisms}} links mechanisms to monomers, model/effect phrases and references.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.88\textwidth]{{figures/pathway_frequency.png}}
\caption{{High-frequency signaling pathways and phenotypic processes extracted from pathway, target and effect fields.}}
\label{{fig:pathways}}
\end{{figure}}

{mechanism_table(df)}

\section{{Experimental Model Context}}
Model selection shapes the meaning of each observation. In vitro cell-line evidence is useful for cytotoxicity, apoptosis, migration, invasion and pathway perturbation. In vivo AOM/DSS, xenograft and syngeneic contexts are more informative for inflammation, immune activity, pharmacokinetics and tissue toxicity. Network pharmacology and molecular docking identify candidate mechanisms but require orthogonal validation.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.72\textwidth]{{figures/model_type_pie.png}}
\caption{{Share of experimental model types in the extracted records.}}
\label{{fig:models}}
\end{{figure}}

\section{{Clinical Translation, Safety and Future Work}}
The clinical trial column was empty in the source table, so clinical efficacy should not be claimed from this dataset. Future work should prioritize: standardized compound identity and purity; concentration-response testing in CRC cells, patient-derived organoids and matched normal intestinal organoids; rescue experiments for ferroptosis, apoptosis and pathway causality; immune-competent AOM/DSS or syngeneic models for inflammation and microenvironment mechanisms; and pharmacokinetic/toxicity packages for each candidate. Combination studies should be mechanism-matched: luteolin and kaempferol with fluoropyrimidine or oxaliplatin sensitization assays, tanshinone IIA with ferroptosis/redox rescue experiments, and astragaloside IV with inflammation-associated CRC and invasion models.

\section{{Data Availability}}
The extracted table is available as \texttt{{monomer\_table.csv}} and \texttt{{monomer\_table.xlsx}} in the manuscript directory. Figures are stored in \texttt{{figures/}}. Bibliographic records are stored in \texttt{{sn-bibliography.bib}}.

\section{{Acknowledgements}}
This work was supported by the National Natural Science Foundation of China (grant number 62106005).

\bibliographystyle{{plainnat}}
\bibliography{{sn-bibliography}}
\end{{document}}
"""
    (PAPER / "sn-article.tex").write_text(textwrap.dedent(tex).strip() + "\n", encoding="utf-8")


def ensure_bib() -> None:
    root_bib = ROOT / "sn-bibliography.bib"
    paper_bib = PAPER / "sn-bibliography.bib"
    if (not paper_bib.exists() or paper_bib.stat().st_size < 1000) and root_bib.exists():
        paper_bib.write_text(root_bib.read_text(encoding="utf-8"), encoding="utf-8")


def bib_clean(s: object) -> str:
    s = html.unescape(clean(s))
    s = s.replace("β", "beta").replace("α", "alpha").replace("κ", "kappa").replace("γ", "gamma")
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode("ascii")
    repl = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_", "{": "", "}": ""}
    return "".join(repl.get(ch, ch) for ch in s)


def fetch_pubmed_records(ids: list[str]) -> dict[str, dict[str, object]]:
    if not ids:
        return {}
    out: dict[str, dict[str, object]] = {}
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    for i in range(0, len(ids), 80):
        batch = ids[i : i + 80]
        r = requests.get(url, params={"db": "pubmed", "id": ",".join(batch), "retmode": "xml"}, timeout=30)
        r.raise_for_status()
        root = ET.fromstring(r.text)
        for article in root.findall(".//PubmedArticle"):
            pmid = article.findtext(".//MedlineCitation/PMID") or ""
            art = article.find(".//Article")
            if not pmid or art is None:
                continue
            authors = []
            for a in art.findall("AuthorList/Author"):
                last = a.findtext("LastName") or ""
                ini = a.findtext("Initials") or ""
                coll = a.findtext("CollectiveName") or ""
                if last:
                    authors.append(f"{last}, {ini}".strip())
                elif coll:
                    authors.append(coll)
            year = art.findtext("Journal/JournalIssue/PubDate/Year") or art.findtext("Journal/JournalIssue/PubDate/MedlineDate") or ""
            m = re.search(r"\d{4}", year)
            doi = ""
            for aid in article.findall(".//ArticleIdList/ArticleId"):
                if aid.attrib.get("IdType") == "doi":
                    doi = aid.text or ""
            out[pmid] = {
                "authors": authors or ["PubMed indexed record"],
                "title": art.findtext("ArticleTitle") or f"PubMed record {pmid}",
                "journal": art.findtext("Journal/Title") or art.findtext("Journal/ISOAbbreviation") or "PubMed",
                "year": m.group(0) if m else "n.d.",
                "doi": doi,
            }
        time.sleep(0.35)
    return out


def ensure_cited_bib_entries() -> None:
    tex = (PAPER / "sn-article.tex").read_text(encoding="utf-8")
    cited = sorted(set(re.findall(r"PMID(\d{6,9})", tex)))
    bib_path = PAPER / "sn-bibliography.bib"
    bib = bib_path.read_text(encoding="utf-8") if bib_path.exists() else ""
    present = set(re.findall(r"@article\{PMID(\d{6,9})", bib))
    missing = [x for x in cited if x not in present]
    if not missing:
        return
    records = fetch_pubmed_records(missing)
    lines = [bib.rstrip(), ""]
    for pmid in missing:
        rec = records.get(pmid)
        if not rec:
            continue
        lines.append(f"@article{{PMID{pmid},")
        lines.append(f"  author = {{{' and '.join(bib_clean(a) for a in rec['authors'])}}},")
        lines.append(f"  title = {{{bib_clean(rec['title'])}}},")
        lines.append(f"  journal = {{{bib_clean(rec['journal'])}}},")
        lines.append(f"  year = {{{bib_clean(rec['year'])}}},")
        lines.append(f"  pmid = {{{pmid}}},")
        if rec.get("doi"):
            lines.append(f"  doi = {{{bib_clean(rec['doi'])}}},")
        lines.append("}")
        lines.append("")
    bib_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    df = load()
    table_outputs(df)
    figures(df)
    ensure_bib()
    write_tex(df)
    ensure_cited_bib_entries()
    print(f"wrote {PAPER / 'sn-article.tex'}")
    print(f"wrote {PAPER / 'monomer_table.csv'}")
    print(f"wrote figures to {FIG}")


if __name__ == "__main__":
    main()
