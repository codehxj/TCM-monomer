from __future__ import annotations

import html
import re
import shutil
import textwrap
import time
import unicodedata
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import requests
import seaborn as sns


ROOT = Path(__file__).resolve().parent
PAPER_DIR = ROOT / "论文"
FIG_DIR = PAPER_DIR / "figures"
HIGHLIGHT_MONOMERS = {"luteolin", "kaempferol", "tanshinone iia", "astragaloside iv"}

REQUIRED_COLUMNS = {
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

CATEGORY_MAP = {
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

PATHWAY_PATTERNS = {
    "PI3K/Akt/mTOR": r"PI3K|Akt|mTOR|p70S6K",
    "Wnt/beta-catenin": r"Wnt|catenin|β-catenin|LEF",
    "NF-kB/inflammation": r"NF|炎症|inflamm|TNF|IL-|COX|PTGS",
    "STAT3": r"STAT3",
    "MAPK/JNK/ERK": r"MAPK|JNK|ERK|p38",
    "Ferroptosis/GPX4": r"ferropt|铁死亡|GPX4|Nrf2|SLC7A11|GSH",
    "Apoptosis/caspases": r"apopt|凋亡|caspase|Bax|Bcl|PARP",
    "EMT/ECM/MMP": r"EMT|ECM|MMP|侵袭|迁移|metasta",
    "ROS/oxidative stress": r"ROS|氧化|oxidative|mitochond",
    "Cell cycle": r"cell cycle|细胞周期|Cyclin|CDK|G2|G1",
    "Autophagy": r"autophagy|自噬|LC3|Beclin",
    "Immune microenvironment": r"PD-L1|macrophage|巨噬|CD8|immune|免疫|T cell",
    "Gut microbiota": r"microbiota|菌群|肠道微生物",
}


def normalize_name(value: object) -> str:
    text = str(value or "").strip()
    text = unicodedata.normalize("NFKD", text)
    text = text.replace("β", "beta").replace("Ⅱ", "II").replace("IIA", "iia")
    text = re.sub(r"[\s_]+", " ", text)
    text = re.sub(r"[^0-9a-zA-Z,.\- ]+", " ", text)
    text = re.sub(r"\s+", " ", text).strip().lower()
    return text


def clean_text(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    text = re.sub(r"\s+", " ", text)
    return text


def ascii_latex(text: object) -> str:
    text = clean_text(text)
    text = text.replace("β", "beta").replace("α", "alpha").replace("κ", "kappa")
    text = text.replace("↓", " down").replace("↑", " up").replace("²", "2")
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii")
    repl = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(repl.get(ch, ch) for ch in text)


def short_text(text: str, limit: int = 120) -> str:
    text = clean_text(text)
    if len(text) <= limit:
        return text
    return text[: limit - 3].rstrip() + "..."


def read_data() -> pd.DataFrame:
    xlsx = None
    for candidate in ROOT.glob("*.xlsx"):
        try:
            preview = pd.read_excel(candidate, sheet_name="Sheet1", nrows=1)
        except Exception:
            continue
        if REQUIRED_COLUMNS["compound"] in preview.columns and REQUIRED_COLUMNS["pmid"] in preview.columns:
            xlsx = candidate
            break
    if xlsx is None:
        raise FileNotFoundError("Could not find the source workbook with Sheet1 and the expected monomer columns.")
    df = pd.read_excel(xlsx, sheet_name="Sheet1")
    df = df.dropna(how="all")
    cols = list(REQUIRED_COLUMNS.values())
    df = df[[c for c in cols if c in df.columns]].copy()
    df = df.rename(columns={v: k for k, v in REQUIRED_COLUMNS.items()})
    for col in df.columns:
        df[col] = df[col].map(clean_text)
    df = df[df["compound"] != ""]
    df = df.drop_duplicates()
    folder_names = {
        normalize_name(p.name): p.name
        for p in ROOT.iterdir()
        if p.is_dir() and p.name not in {"figures", "论文", "__pycache__"}
    }
    df["compound_key"] = df["compound"].map(normalize_name)
    df["folder_matched"] = df["compound_key"].isin(folder_names)
    df["category"] = df["compound_key"].map(CATEGORY_MAP).fillna("Other natural products")
    return df


def classify_model(text: str) -> str:
    t = text.lower()
    if re.search(r"网络药理|分子对接|in silico|bioinform|comput|molecular docking", t):
        return "Network/in silico"
    if re.search(r"小鼠|大鼠|aom|dss|xenograft|nude|裸鼠|mc38|animal|mouse|mice|rat", t):
        return "In vivo"
    if re.search(r"hct|ht-29|sw480|sw620|cell|细胞|caco|ct26|rko|loVo".lower(), t):
        return "In vitro"
    if t:
        return "Other experimental evidence"
    return "Not specified"


def pathway_hits(row: pd.Series) -> list[str]:
    text = " ".join(str(row.get(c, "")) for c in ["pathway", "effect", "target"])
    hits = []
    for label, pattern in PATHWAY_PATTERNS.items():
        if re.search(pattern, text, flags=re.I):
            hits.append(label)
    return hits


def extract_pmids(series: pd.Series) -> list[str]:
    pmids: list[str] = []
    for value in series:
        for match in re.findall(r"\b\d{6,9}\b", str(value)):
            if match not in pmids:
                pmids.append(match)
    return pmids


def fetch_pubmed(pmids: list[str]) -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    for i in range(0, len(pmids), 80):
        batch = pmids[i : i + 80]
        params = {"db": "pubmed", "id": ",".join(batch), "retmode": "xml"}
        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            root = ET.fromstring(response.text)
        except Exception as exc:
            print(f"PubMed fetch failed for batch {i // 80 + 1}: {exc}")
            continue
        for article in root.findall(".//PubmedArticle"):
            pmid = article.findtext(".//MedlineCitation/PMID") or ""
            art = article.find(".//Article")
            if not pmid or art is None:
                continue
            title = "".join(art.findtext("ArticleTitle") or "").strip()
            journal = art.findtext("Journal/Title") or art.findtext("Journal/ISOAbbreviation") or ""
            year = (
                art.findtext("Journal/JournalIssue/PubDate/Year")
                or art.findtext("Journal/JournalIssue/PubDate/MedlineDate")
                or ""
            )
            year_match = re.search(r"\d{4}", year)
            year = year_match.group(0) if year_match else ""
            authors = []
            for author in art.findall("AuthorList/Author"):
                last = author.findtext("LastName") or ""
                initials = author.findtext("Initials") or ""
                collective = author.findtext("CollectiveName") or ""
                if last:
                    authors.append(f"{last}, {initials}".strip())
                elif collective:
                    authors.append(collective)
            doi = ""
            for aid in article.findall(".//ArticleIdList/ArticleId"):
                if aid.attrib.get("IdType") == "doi":
                    doi = aid.text or ""
                    break
            records[pmid] = {
                "title": html.unescape(title),
                "journal": html.unescape(journal),
                "year": year,
                "authors": authors,
                "doi": doi,
            }
        time.sleep(0.35)
    return records


def bibtex_escape(text: object) -> str:
    text = html.unescape(clean_text(text))
    text = text.replace("β", "beta").replace("α", "alpha").replace("κ", "kappa")
    text = text.replace("γ", "gamma").replace("δ", "delta").replace("≤", "<=")
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii")
    repl = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": "",
        "}": "",
    }
    text = "".join(repl.get(ch, ch) for ch in text)
    return text


def write_bib(pmids: list[str], records: dict[str, dict[str, object]]) -> None:
    lines = []
    for pmid in pmids:
        rec = records.get(pmid, {})
        authors = rec.get("authors") or ["PubMed indexed record"]
        title = rec.get("title") or f"PubMed record {pmid}"
        journal = rec.get("journal") or "PubMed"
        year = rec.get("year") or "n.d."
        doi = rec.get("doi") or ""
        lines.append(f"@article{{PMID{pmid},")
        lines.append(f"  author = {{{' and '.join(map(bibtex_escape, authors))}}},")
        lines.append(f"  title = {{{bibtex_escape(title)}}},")
        lines.append(f"  journal = {{{bibtex_escape(journal)}}},")
        lines.append(f"  year = {{{bibtex_escape(year)}}},")
        lines.append(f"  pmid = {{{pmid}}},")
        if doi:
            lines.append(f"  doi = {{{bibtex_escape(doi)}}},")
        lines.append("}")
        lines.append("")
    text = "\n".join(lines)
    (ROOT / "sn-bibliography.bib").write_text(text, encoding="utf-8")
    (PAPER_DIR / "sn-bibliography.bib").write_text(text, encoding="utf-8")


def save_clean_tables(df: pd.DataFrame) -> None:
    tidy = df[
        [
            "compound",
            "category",
            "source",
            "model",
            "pathway",
            "effect",
            "target",
            "clinical",
            "toxicity",
            "pmid",
            "folder_matched",
        ]
    ].copy()
    tidy.to_csv(ROOT / "monomer_crc_tidy.csv", index=False, encoding="utf-8-sig")
    tidy.to_excel(ROOT / "monomer_crc_tidy.xlsx", index=False)
    paper_table = tidy.rename(
        columns={
            "category": "化学类别",
            "compound": "单体",
            "source": "来源",
            "model": "实验模型",
            "pathway": "作用通路",
            "effect": "核心效应",
            "pmid": "参考文献",
        }
    )[["化学类别", "单体", "来源", "实验模型", "作用通路", "核心效应", "参考文献"]]
    paper_table = paper_table.sort_values(["化学类别", "单体", "参考文献"])
    paper_table.to_csv(PAPER_DIR / "monomer_table.csv", index=False, encoding="utf-8-sig")
    paper_table.to_excel(PAPER_DIR / "monomer_table.xlsx", index=False)


def make_figures(df: pd.DataFrame) -> dict[str, Counter]:
    FIG_DIR.mkdir(exist_ok=True)
    sns.set_theme(style="whitegrid", font="DejaVu Sans")
    category_counts = Counter(df.drop_duplicates("compound_key")["category"])
    model_counts = Counter(df["model_class"])
    pathway_counts: Counter[str] = Counter()
    for hits in df["pathway_hits"]:
        pathway_counts.update(hits)

    plt.figure(figsize=(8.6, 4.8))
    cat_items = sorted(category_counts.items(), key=lambda x: (-x[1], x[0]))
    ax = sns.barplot(x=[v for _, v in cat_items], y=[k for k, _ in cat_items], color="#4C78A8")
    ax.set_xlabel("Number of distinct monomers")
    ax.set_ylabel("Chemical category")
    ax.set_title("Distribution of natural monomers by chemical category")
    for container in ax.containers:
        ax.bar_label(container, padding=3)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "category_distribution.png", dpi=300)
    plt.close()

    plt.figure(figsize=(6.8, 5.2))
    labels = list(model_counts.keys())
    sizes = list(model_counts.values())
    colors = sns.color_palette("Set2", len(labels))
    plt.pie(sizes, labels=labels, autopct="%1.1f%%", startangle=110, colors=colors)
    plt.title("Experimental model categories")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "model_type_pie.png", dpi=300)
    plt.close()

    top_pathways = pathway_counts.most_common(12)
    plt.figure(figsize=(8.8, 5.2))
    ax = sns.barplot(x=[v for _, v in top_pathways], y=[k for k, _ in top_pathways], color="#59A14F")
    ax.set_xlabel("Number of records mentioning the pathway/process")
    ax.set_ylabel("Pathway/process")
    ax.set_title("Frequently reported mechanisms")
    for container in ax.containers:
        ax.bar_label(container, padding=3)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "pathway_frequency.png", dpi=300)
    plt.close()
    make_mechanism_network(df)
    return {"category": category_counts, "model": model_counts, "pathway": pathway_counts}


def make_mechanism_network(df: pd.DataFrame) -> None:
    """Draw a compact monomer -> pathway -> effect -> anti-CRC mechanism map."""
    FIG_DIR.mkdir(exist_ok=True)
    reps = [
        ("Luteolin", "Flavonoids", "PI3K/Akt/mTOR", "Apoptosis and ferroptosis"),
        ("Kaempferol", "Flavonoids", "Wnt/beta-catenin", "Chemosensitization"),
        ("Tanshinone IIA", "Diterpenoid quinones", "Ferroptosis/GPX4", "Cell-cycle arrest"),
        ("Astragaloside IV", "Triterpenoid saponins", "PI3K/Akt/mTOR", "Invasion suppression"),
        ("Honokiol", "Neolignans", "STAT3/NF-kB", "Stemness reduction"),
        ("Rhein", "Anthraquinones", "ROS/ER stress", "Apoptosis"),
        ("Sanguinarine", "Alkaloids", "MAPK/JNK/ERK", "Migration inhibition"),
        ("beta-sitosterol", "Phytosterols", "Wnt/beta-catenin", "Proliferation inhibition"),
        ("Caffeic acid", "Phenolic acids", "NF-kB/inflammation", "Inflammation restraint"),
    ]
    palette = {
        "Flavonoids": "#4C78A8",
        "Diterpenoid quinones": "#E15759",
        "Triterpenoid saponins": "#59A14F",
        "Neolignans": "#B07AA1",
        "Anthraquinones": "#F28E2B",
        "Alkaloids": "#9C755F",
        "Phytosterols": "#76B7B2",
        "Phenolic acids": "#EDC948",
    }
    pathway_colors = "#F2F5F9"
    effect_colors = "#F6F0E8"
    anti_crc_color = "#D9EAD3"

    fig, ax = plt.subplots(figsize=(14, 8), dpi=180)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    def node(x, y, text, fc, ec="#4D4D4D", w=0.18, h=0.055, fs=9, weight="normal"):
        box = plt.Rectangle((x - w / 2, y - h / 2), w, h, facecolor=fc, edgecolor=ec, linewidth=1.2)
        ax.add_patch(box)
        ax.text(x, y, text, ha="center", va="center", fontsize=fs, weight=weight, wrap=True)

    ys = [0.88 - i * 0.085 for i in range(len(reps))]
    pathway_y = {}
    effect_y = {}
    for i, (compound, category, pathway, effect) in enumerate(reps):
        y = ys[i]
        node(0.12, y, compound, palette.get(category, "#BAB0AC"), w=0.18, fs=8.5, weight="bold" if normalize_name(compound) in HIGHLIGHT_MONOMERS else "normal")
        pathway_y.setdefault(pathway, y)
        effect_y.setdefault(effect, y)
        ax.annotate("", xy=(0.31, y), xytext=(0.21, y), arrowprops=dict(arrowstyle="->", lw=1.1, color="#555555"))

    p_items = list(pathway_y.items())
    for j, (pathway, y) in enumerate(p_items):
        node(0.42, y, pathway, pathway_colors, w=0.22, fs=8.3)
        ax.annotate("", xy=(0.58, y), xytext=(0.53, y), arrowprops=dict(arrowstyle="->", lw=1.1, color="#555555"))

    e_items = list(effect_y.items())
    for effect, y in e_items:
        node(0.69, y, effect, effect_colors, w=0.22, fs=8.3)
        ax.annotate("", xy=(0.84, 0.50), xytext=(0.80, y), arrowprops=dict(arrowstyle="->", lw=1.0, color="#777777", alpha=0.75))

    node(0.91, 0.50, "Anti-CRC effects\nreduced growth, invasion,\ninflammation and resistance", anti_crc_color, w=0.18, h=0.12, fs=9, weight="bold")
    ax.text(0.12, 0.96, "Representative monomers", ha="center", va="center", fontsize=11, weight="bold")
    ax.text(0.42, 0.96, "Signal pathways", ha="center", va="center", fontsize=11, weight="bold")
    ax.text(0.69, 0.96, "Cellular effects", ha="center", va="center", fontsize=11, weight="bold")
    ax.text(0.91, 0.66, "Integrated outcome", ha="center", va="center", fontsize=11, weight="bold")

    legend_y = 0.06
    x0 = 0.05
    for idx, (category, color) in enumerate(palette.items()):
        x = x0 + (idx % 4) * 0.235
        y = legend_y - (idx // 4) * 0.035
        ax.add_patch(plt.Rectangle((x, y), 0.018, 0.018, facecolor=color, edgecolor="none"))
        ax.text(x + 0.024, y + 0.009, category, va="center", fontsize=8.2)
    ax.set_title("TCM monomers regulate CRC through convergent pathway-effect modules", fontsize=14, weight="bold", pad=18)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "mechanism_network.png", dpi=300, bbox_inches="tight")
    plt.close()


def representative_cites(df: pd.DataFrame) -> dict[str, str]:
    cites = {}
    for name in ["Luteolin", "Kaempferol", "Honokiol", "Tanshinone IIA", "Astragaloside IV", "rhein", "beta-sitosterol"]:
        rows = df[df["compound"].str.lower() == name.lower()]
        pmids = extract_pmids(rows["pmid"]) if not rows.empty else []
        if pmids:
            cites[name] = f"PMID{pmids[0]}"
    all_pmids = extract_pmids(df["pmid"])
    cites["general"] = f"PMID{all_pmids[0]}" if all_pmids else ""
    return cites


def build_summary_table(df: pd.DataFrame) -> str:
    grouped = []
    for key, sub in df.groupby("compound_key", sort=True):
        compound = sub["compound"].iloc[0]
        category = sub["category"].iloc[0]
        source = "; ".join([x for x, _ in Counter(sub["source"]).most_common(2) if x]) or "Not specified"
        model = "; ".join([x for x, _ in Counter(sub["model_class"]).most_common(2) if x]) or "Not specified"
        hits = Counter(x for hits in sub["pathway_hits"] for x in hits)
        mechanisms = ", ".join([k for k, _ in hits.most_common(3)]) or short_text("; ".join(sub["pathway"].dropna().astype(str).head(2)), 95) or "Mechanism not specified"
        effects = short_text("; ".join([x for x, _ in Counter(sub["effect"]).most_common(2) if x]), 115) or "Not specified"
        pmids = extract_pmids(sub["pmid"])[:3]
        grouped.append((compound, category, source, model, mechanisms, effects, ", ".join([f"\\cite{{PMID{x}}}" for x in pmids])))
    grouped.sort(key=lambda x: (x[1], x[0].lower()))
    rows = [
        r"\begin{scriptsize}",
        r"\begin{longtable}{L{0.13\textwidth}L{0.16\textwidth}L{0.14\textwidth}L{0.20\textwidth}L{0.20\textwidth}L{0.12\textwidth}}",
        r"\caption{Chemical-category-organized monomers, sources, models, mechanisms and representative references. Highlighted rows denote priority monomers for follow-up validation.}\label{tab:monomers}\\",
        r"\toprule",
        r"Monomer & Source & Experimental model & Pathway & Core effect & References \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"Monomer & Source & Experimental model & Pathway & Core effect & References \\",
        r"\midrule",
        r"\endhead",
    ]
    current_category = None
    for compound, category, source, model, mechanisms, effects, cites in grouped:
        if category != current_category:
            rows.append(rf"\multicolumn{{6}}{{l}}{{\textbf{{{ascii_latex(category)}}}}}\\")
            current_category = category
        prefix = r"\rowcolor{yellow!16} " if normalize_name(compound) in HIGHLIGHT_MONOMERS else ""
        monomer = rf"\textbf{{{ascii_latex(compound)}}}" if normalize_name(compound) in HIGHLIGHT_MONOMERS else ascii_latex(compound)
        rows.append(
            f"{prefix}{monomer} & {ascii_latex(short_text(source, 90))} & {ascii_latex(short_text(model, 80))} & "
            f"{ascii_latex(short_text(mechanisms, 95))} & {ascii_latex(effects)} & {cites} \\\\"
        )
    rows.extend([r"\bottomrule", r"\end{longtable}", r"\end{scriptsize}"])
    return "\n".join(rows)


def build_mechanism_reference_table(df: pd.DataFrame) -> str:
    rows = [
        r"\begin{scriptsize}",
        r"\begin{longtable}{L{0.19\textwidth}L{0.25\textwidth}L{0.34\textwidth}L{0.15\textwidth}}",
        r"\caption{Mechanism-oriented reference map linking recurrent pathways to representative monomers and cellular readouts.}\label{tab:mechanism_refs}\\",
        r"\toprule",
        r"Mechanism & Representative monomers & Typical cellular readouts & References \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"Mechanism & Representative monomers & Typical cellular readouts & References \\",
        r"\midrule",
        r"\endhead",
    ]
    for mechanism, _count in Counter(x for hits in df["pathway_hits"] for x in hits).most_common(12):
        sub = df[df["pathway_hits"].map(lambda hits: mechanism in hits)]
        monomers = ", ".join([x for x, _ in Counter(sub["compound"]).most_common(5)])
        effects = "; ".join([x for x, _ in Counter(sub["effect"]).most_common(3) if x])
        pmids = extract_pmids(sub["pmid"])[:4]
        cites = ", ".join([f"\\cite{{PMID{x}}}" for x in pmids])
        rows.append(f"{ascii_latex(mechanism)} & {ascii_latex(short_text(monomers, 120))} & {ascii_latex(short_text(effects, 170))} & {cites} \\\\")
    rows.extend([r"\bottomrule", r"\end{longtable}", r"\end{scriptsize}"])
    return "\n".join(rows)


def write_tex(df: pd.DataFrame, stats: dict[str, Counter]) -> None:
    cites = representative_cites(df)
    n_records = len(df)
    n_monomers = df["compound_key"].nunique()
    n_categories = df.drop_duplicates("compound_key")["category"].nunique()
    matched = int(df.drop_duplicates("compound_key")["folder_matched"].sum())
    cat_sentence = ", ".join(f"{k} ({v})" for k, v in stats["category"].most_common(6))
    model_sentence = ", ".join(f"{k}: {v}" for k, v in stats["model"].most_common())
    pathway_sentence = ", ".join(f"{k} ({v})" for k, v in stats["pathway"].most_common(8))
    table = build_summary_table(df)
    mechanism_table = build_mechanism_reference_table(df)

    tex = rf"""
% Auto-generated from monomer_crc_tidy.xlsx by generate_crc_review.py.
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
\newcolumntype{{Y}}{{>{{\raggedright\arraybackslash}}X}}
\raggedbottom

\begin{{document}}

\title[Traditional Chinese medicine monomers and colorectal cancer]{{Natural monomers from traditional Chinese medicine in colorectal cancer: chemical classes, experimental models and mechanistic evidence}}

\author[1]{{Xianjun Han}}\email{{hxj@ahu.edu.cn}}
\author[1]{{Jincheng Fang}}\email{{E24301275@stu.ahu.edu.cn}}
\equalcont{{These authors contributed equally to this work.}}
\author[2]{{Zijian Wu}}\email{{wuzijian@ahtcm.edu.cn}}
\author*[2,3]{{Can Bai}}\email{{baican@ahtcm.edu.cn}}
\author[3]{{Renbao Huang}}\email{{hhmoxue@163.com}}

\affil[1]{{\orgdiv{{School of Computer Science and Technology}}, \orgname{{Anhui University}}, \orgaddress{{\street{{111 Jiulong Road}}, \city{{Hefei}}, \postcode{{230601}}, \state{{Anhui}}, \country{{China}}}}}}
\affil[2]{{\orgdiv{{School of Acupuncture and Tuina College}}, \orgname{{Anhui University of Chinese Medicine}}, \orgaddress{{\street{{350 Longzihu Road}}, \city{{Hefei}}, \postcode{{230012}}, \state{{Anhui}}, \country{{China}}}}}}
\affil[3]{{\orgdiv{{First Clinical Medical College}}, \orgname{{Anhui University of Chinese Medicine}}, \orgaddress{{\street{{350 Longzihu Road}}, \city{{Hefei}}, \postcode{{230012}}, \state{{Anhui}}, \country{{China}}}}}}

\abstract{{Traditional Chinese medicine (TCM) contains chemically diverse small molecules with growing experimental evidence in colorectal cancer (CRC). We curated {n_records} records covering {n_monomers} monomers from a structured spreadsheet and matched {matched} monomers to the prepared compound folders. The resulting evidence map spans {n_categories} chemical categories, cell-based assays, animal models, network pharmacology, molecular docking and formulation studies. Flavonoids, terpenoid-derived compounds, anthraquinones, alkaloids, phytosterols and phenolic acids recur as major groups. Mechanistically, the literature concentrates on apoptosis, ferroptosis, PI3K/Akt/mTOR, Wnt/beta-catenin, NF-kB-driven inflammation, STAT3, EMT/ECM remodeling, oxidative stress and immune microenvironment remodeling. This review integrates chemical classification, model context, mechanism networks, clinically relevant safety gaps and experimentally actionable development strategies for TCM monomers in CRC.}}

\keywords{{colorectal cancer, traditional Chinese medicine, natural monomers, flavonoids, ferroptosis, PI3K/Akt, Wnt/beta-catenin, PubMed}}

\maketitle

\section{{Introduction}}
Colorectal cancer (CRC) remains a major cause of cancer morbidity and mortality worldwide, and its clinical management still faces recurrence, metastatic spread, chemoresistance and treatment-limiting toxicity. Natural monomers from TCM provide a chemically rich source of candidate modulators because many compounds act on convergent cancer phenotypes rather than on a single isolated target. In the curated dataset used here, the evidence base includes direct cell viability and migration assays, AOM/DSS-associated carcinogenesis models, xenografts, network pharmacology, molecular docking and drug-delivery studies. This spectrum is valuable for hypothesis generation, but it also requires careful separation of mechanistic plausibility from clinical readiness.

The present review is based on a structured curation of \texttt{{monomer\_crc\_tidy.xlsx}} and the associated compound folders. The dataset includes repeated reports for well-studied monomers such as luteolin, kaempferol, honokiol, tanshinone IIA, astragaloside IV, rhein and beta-sitosterol \cite{{{cites.get("Luteolin", "PMID40488304")},{cites.get("Kaempferol", "PMID38493716")},{cites.get("Honokiol", "PMID34192810")},{cites.get("Tanshinone IIA", "PMID33535870")},{cites.get("Astragaloside IV", "PMID38831823")},{cites.get("rhein", "PMID33946531")},{cites.get("beta-sitosterol", "PMID36603684")}}}. Rather than treating these monomers as interchangeable anti-cancer agents, we organize them by chemical class, model system and recurrent mechanism.

Conceptually, the curated evidence can be read as a four-layer pharmacological network: monomers, signaling pathways, cellular phenotypes and integrated anti-CRC outcomes. Fig.~\ref{{fig:network}} summarizes this logic. The diagram emphasizes that chemically distinct monomers can converge on a limited set of pathway-effect modules, while structurally related compounds can still diverge in their dominant readouts depending on dose, model and delivery system.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.98\textwidth]{{figures/mechanism_network.png}}
\caption{{Network schematic linking representative TCM monomers to signaling pathways, cellular effects and anti-CRC outcomes. Node colors in the monomer column denote chemical categories, and highlighted monomers correspond to priority candidates discussed in the text.}}
\label{{fig:network}}
\end{{figure}}

\section{{Natural Monomer Classification}}
Chemical classification helps clarify why structurally related monomers often converge on similar biological effects. The curated monomer-level distribution is shown in Fig.~\ref{{fig:category}}. The largest classes were {ascii_latex(cat_sentence)}. Flavonoids dominated the collection, reflecting the frequent study of luteolin, kaempferol, myricetin, isorhamnetin, formononetin, glabridin, pinocembrin, catechin and silibinin in CRC-related systems. These compounds commonly affect oxidative stress, inflammatory signaling, apoptosis and epithelial-mesenchymal transition.

Terpenoid and terpenoid-derived molecules formed another important axis, including astragaloside IV, curdione, eucalyptol, hederagenin, pachymic acid, tagitinin C and tanshinone IIA. Anthraquinones such as rhein, emodin, aloe-emodin and chrysophanol were repeatedly linked to apoptosis, ER stress, mitochondrial dysfunction and ferroptosis-like responses. Alkaloids and alkaloid-like products, including dictamnine, sanguinarine, neothalfine, naamidine J and secoemestrin C, added a mechanistically distinct group with evidence for cytotoxicity, kinase modulation and stress-response pathways. Phenolic acids, tannins, neolignans and phytosterols rounded out the landscape and often appeared in studies focused on inflammation, gut barrier function, Wnt signaling or formulation.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.88\textwidth]{{figures/category_distribution.png}}
\caption{{Distribution of distinct natural monomers across chemical categories after folder-matched curation.}}
\label{{fig:category}}
\end{{figure}}

{table}

\section{{Pharmacological Actions and Mechanisms}}
The curated evidence emphasizes several recurring pharmacological themes. First, many monomers reduce CRC cell proliferation and promote programmed cell death. Flavonoids such as luteolin and kaempferol are frequently associated with caspase activation, Bcl-2 family rebalancing, mitochondrial dysfunction and cell-cycle arrest. Anthraquinones and quinones broaden this theme through redox stress and ER-stress-associated apoptosis, while tanshinone IIA and tagitinin C are often discussed in relation to ferroptosis.

Second, signal transduction pathways repeatedly converge on PI3K/Akt/mTOR, Wnt/beta-catenin, STAT3, NF-kB and MAPK modules. These pathways connect growth-factor signaling, inflammatory tone, survival, invasion and stemness. The frequency map in Fig.~\ref{{fig:pathways}} highlights the most common mechanistic labels extracted from the pathway, effect and target fields: {ascii_latex(pathway_sentence)}. This does not imply that every report provides equivalent causal evidence; network pharmacology and docking studies should be interpreted as prioritization evidence, whereas perturbation experiments in cell or animal models provide stronger mechanistic support.

Third, several monomers target the tumor microenvironment. The dataset contains studies linking natural products to macrophage polarization, CD8-positive T-cell activation, PD-L1 regulation, cytokine suppression, ECM remodeling and MMP down-regulation. These effects are especially relevant for colitis-associated CRC, where inflammation, epithelial injury, dysbiosis and tumor initiation are biologically intertwined.

\subsection{{Flavonoids and polyphenolic regulators}}
Flavonoids represent the densest evidence cluster. Luteolin is repeatedly connected with PI3K/Akt/mTOR inhibition, STAT3/NF-kB attenuation, GPX4-related ferroptosis, apoptosis and immune activation \cite{{{cites.get("Luteolin", "PMID40488304")}}}. These findings suggest a dual role: direct suppression of malignant epithelial phenotypes and reshaping of inflammatory or immune context. Kaempferol has a similarly broad profile, with reports implicating Wnt/beta-catenin, ER stress, DNA damage, 5-fluorouracil sensitization and gut-microbiota-linked bile acid signaling \cite{{{cites.get("Kaempferol", "PMID38493716")}}}. For these two priority flavonoids, the next decisive experiments should quantify exposure-response relationships in matched CRC organoids, normal intestinal organoids and immune-competent animal models.

\subsection{{Terpenoids, quinones and ferroptosis-oriented mechanisms}}
Tanshinone IIA is a prominent diterpenoid quinone with evidence for apoptosis, cell-cycle arrest, survivin targeting, Ang2/Tie2 modulation and ferroptosis-related activity \cite{{{cites.get("Tanshinone IIA", "PMID33535870")}}}. Astragaloside IV, a triterpenoid saponin, appears more frequently in PI3K/Akt/mTOR, circRNA/miRNA and inflammatory-microenvironment contexts \cite{{{cites.get("Astragaloside IV", "PMID38831823")}}}. This distinction is useful translationally: tanshinone IIA may be prioritized for redox and ferroptosis-combination screens, whereas astragaloside IV is better suited to inflammation-associated CRC and invasion-suppression models.

\subsection{{Anthraquinones, alkaloids, sterols and lignans}}
Anthraquinones such as rhein, emodin and aloe-emodin are commonly linked to mitochondrial stress, ER stress, apoptosis and migration inhibition \cite{{{cites.get("rhein", "PMID33946531")}}}. Alkaloids such as sanguinarine and dictamnine add kinase, MAPK, NF-kB and ferroptosis-associated mechanisms, although their therapeutic window needs especially careful evaluation. Phytosterols, represented by beta-sitosterol and stigmasterol, recur in Wnt/beta-catenin, inflammation and network-pharmacology studies \cite{{{cites.get("beta-sitosterol", "PMID36603684")}}}. Neolignans such as honokiol and magnolol connect anti-inflammatory, anti-stemness and immune-modulating activities, and therefore may be better framed as sensitizers or microenvironment modulators than as stand-alone cytotoxins.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.88\textwidth]{{figures/pathway_frequency.png}}
\caption{{Frequently reported pathways and biological processes across curated records. A record can contribute to more than one mechanism.}}
\label{{fig:pathways}}
\end{{figure}}

{mechanism_table}

\section{{Experimental Models}}
Experimental context strongly affects how the evidence should be read. The model distribution is summarized in Fig.~\ref{{fig:models}} and was classified as {ascii_latex(model_sentence)}. In vitro CRC cell lines, including HCT116, HT-29, SW480, SW620, CT26, MC38 and related systems, remain the most common setting for estimating cytotoxicity, migration, invasion, apoptosis and pathway modulation. These studies are efficient and mechanistically useful, but they cannot fully represent metabolism, bioavailability, immune context or intestinal ecology.

In vivo models provide a more integrated view. AOM/DSS models are particularly informative for inflammation-driven tumorigenesis, whereas xenograft and syngeneic models are useful for tumor growth, immune response and delivery-system evaluation. Computational studies, including network pharmacology, molecular docking and bioinformatic target prioritization, were also common. Their main value is to connect multi-target monomers with candidate pathways, but their conclusions require experimental validation.

A practical experimental pipeline should therefore separate four questions. First, does the monomer reach effective intracellular concentrations in CRC cells or organoids? Second, does pathway modulation occur before cell death rather than as a late stress artifact? Third, is the effect maintained in co-culture or immune-competent settings? Fourth, does the compound preserve selectivity against normal epithelial, hepatic and hematopoietic cells? Applying this sequence to luteolin, kaempferol, tanshinone IIA and astragaloside IV would turn the current evidence map into testable translational programs.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.72\textwidth]{{figures/model_type_pie.png}}
\caption{{Distribution of experimental model categories among curated records.}}
\label{{fig:models}}
\end{{figure}}

\section{{Clinical Translation and Safety}}
The clinical field is still immature. In the source spreadsheet, the clinical trial column was empty across all curated records, indicating that the current dataset mainly supports preclinical and mechanistic conclusions rather than clinical efficacy claims. This gap is important because many natural monomers have limited aqueous solubility, variable oral absorption, extensive metabolism or context-dependent cytotoxicity. Nanocarriers, prodrugs, combination strategies and microbiota-aware delivery may improve exposure, but they also add regulatory and safety complexity.

Safety data were sparse, with only a small number of rows containing explicit adverse-effect notes. The absence of reported toxicity in a curation table should not be read as evidence of safety. Several monomers can interact with chemotherapy, redox balance, cytochrome metabolism or immune signaling; therefore, translational development should include pharmacokinetics, maximum tolerated dose, organ toxicity, genotoxicity where relevant, and herb-drug interaction testing. For compounds proposed as adjuvants to 5-fluorouracil, oxaliplatin, irinotecan or immunotherapy, synergy should be evaluated together with normal intestinal epithelial and hematopoietic toxicity.

\section{{Future Research Directions}}
Future work should move from descriptive pathway lists toward experimentally ranked development strategies. For flavonoids, a high-yield direction is combination testing with 5-fluorouracil, oxaliplatin or irinotecan using Bliss, Loewe or ZIP synergy models, followed by validation of apoptosis, DNA-damage response and Wnt/beta-catenin readouts. For ferroptosis-oriented compounds such as tanshinone IIA, beta-lapachone and tagitinin C, screens should combine GPX4, SLC7A11, Nrf2 and lipid-peroxidation assays with rescue experiments using ferrostatin-1 or liproxstatin-1. For inflammation-associated monomers such as astragaloside IV, honokiol, glabridin and caffeic acid, AOM/DSS models, cytokine panels, barrier-function assays and microbiota profiling should be prioritized.

Combination strategies should be mechanism-matched rather than merely empirical. PI3K/Akt/mTOR-modulating monomers can be paired with chemotherapy or radiotherapy to test survival-pathway blockade; Wnt/beta-catenin modulators can be paired with stemness and organoid-renewal assays; NF-kB/STAT3 modulators should be tested in colitis-associated CRC and macrophage-rich co-cultures; ferroptosis inducers should be evaluated with redox-active drugs and iron-metabolism biomarkers. Formulation studies should report not only tumor inhibition but also plasma exposure, tissue distribution, release kinetics and off-target toxicity.

\section{{Conclusion and Outlook}}
This curated review shows that TCM monomers relevant to CRC are chemically diverse but mechanistically convergent. Flavonoids, terpenoids, anthraquinones, alkaloids, phytosterols, neolignans and phenolic compounds repeatedly modulate proliferation, apoptosis, ferroptosis, inflammation, EMT, oxidative stress, Wnt/beta-catenin, PI3K/Akt/mTOR, STAT3 and immune-microenvironment pathways. The strongest near-term value of this evidence map is not to claim clinical efficacy, but to prioritize monomers and combinations for rigorous validation.

Future work should emphasize standardized compound identity and purity, dose-exposure relationships, orthogonal pathway validation, matched normal-cell toxicity, organoid and patient-derived models, microbiota-aware animal studies and transparent negative results. Clinically useful development will require moving beyond isolated pathway diagrams toward reproducible pharmacology, biomarker-defined indications and safety packages compatible with modern oncology standards. The most immediate candidates for structured follow-up are luteolin and kaempferol for flavonoid-based chemosensitization, tanshinone IIA for ferroptosis/redox programs, and astragaloside IV for inflammation-associated and invasion-suppression models.

\section{{Data Availability}}
The curated table generated for this review is provided as \texttt{{monomer\_table.csv}} in the manuscript directory. Figures are stored in \texttt{{figures/}}, including the mechanism network and three statistical charts. The PubMed-derived bibliography is stored as \texttt{{sn-bibliography.bib}}.

\section{{Acknowledgements}}
This work was supported by the National Natural Science Foundation of China (grant number 62106005).

\bibliographystyle{{plainnat}}
\bibliography{{sn-bibliography}}

\end{{document}}
"""
    text = textwrap.dedent(tex).strip() + "\n"
    (ROOT / "sn-article.tex").write_text(text, encoding="utf-8")
    (PAPER_DIR / "sn-article.tex").write_text(text, encoding="utf-8")


def copy_template_assets() -> None:
    for name in ["sn-jnl.cls", "sn-mathphys-num.bst"]:
        src = PAPER_DIR / name
        dst = ROOT / name
        if src.exists() and not dst.exists():
            shutil.copy2(src, dst)
    bst_src = PAPER_DIR / "bst"
    bst_dst = ROOT / "bst"
    if bst_src.exists() and not bst_dst.exists():
        shutil.copytree(bst_src, bst_dst)
    FIG_DIR.mkdir(exist_ok=True)


def main() -> None:
    copy_template_assets()
    df = read_data()
    df["model_class"] = df["model"].map(classify_model)
    df["pathway_hits"] = df.apply(pathway_hits, axis=1)
    save_clean_tables(df)
    stats = make_figures(df)
    pmids = extract_pmids(df["pmid"])
    print(f"Fetching {len(pmids)} unique PubMed records")
    records = fetch_pubmed(pmids)
    print(f"Fetched {len(records)} PubMed records")
    write_bib(pmids, records)
    write_tex(df, stats)
    print(f"Generated {ROOT / 'sn-article.tex'}")
    print(f"Generated {ROOT / 'sn-bibliography.bib'}")
    print(f"Generated figures in {FIG_DIR}")


if __name__ == "__main__":
    main()
