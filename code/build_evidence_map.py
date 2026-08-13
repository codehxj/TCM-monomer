from __future__ import annotations

import re
import shutil
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


PAPER = Path(__file__).resolve().parent
ROOT = PAPER.parent
SOURCE_XLSX = ROOT / "\u5355\u4f53.xlsx"
SOURCE_BIB = PAPER / "sn-bibliography.bib"
SOURCE_BIB_BACKUP = PAPER / "sn-bibliography.source.bib"
FIGURES = PAPER / "figures"

CATEGORY = {
    "luteolin": "Flavonoids",
    "kaempferol": "Flavonoids",
    "catechin": "Flavonoids",
    "formononetin": "Flavonoids",
    "glabridin": "Flavonoids",
    "astragalin": "Flavonoids",
    "isovitexin": "Flavonoids",
    "isorhamnetin": "Flavonoids",
    "jaranol": "Flavonoids",
    "myricetin": "Flavonoids",
    "pinocembrin": "Flavonoids",
    "silibinin": "Flavonoids",
    "wighteone": "Flavonoids",
    "tanshinone iia": "Diterpenoid quinones",
    "astragaloside iv": "Triterpenoid saponins",
    "honokiol": "Neolignans",
    "magnolol": "Neolignans",
    "obovatol": "Neolignans",
    "rhein": "Anthraquinones",
    "chrysophanol": "Anthraquinones",
    "emodin/aloe-emodin": "Anthraquinones",
    "emodin   aloe-emodin": "Anthraquinones",
    "beta-sitosterol": "Phytosterols",
    "stigmasterol": "Phytosterols",
    "beta-lapachone": "Quinones",
    "caffeic acid": "Phenolic acids",
    "syringic acid": "Phenolic acids",
    "ellagic acid": "Polyphenols",
    "curdione": "Terpenoids",
    "eucalyptol": "Terpenoids",
    "dictamnine": "Alkaloids",
    "sanguinarine": "Alkaloids",
    "neothalfine": "Alkaloids",
    "naamidine j": "Alkaloids",
    "hederagenin": "Triterpenoids",
    "pachymic acid": "Triterpenoids",
    "amygdalin": "Cyanogenic glycosides",
    "macrocarpal i": "Meroterpenoids",
    "tagitinin c": "Sesquiterpene lactones",
    "secoemestrin c": "Fungal alkaloid-like metabolites",
    "1,4,6-tri-o-galloyl-\u03b2-d-glucopyranose": "Hydrolysable tannins",
}

CRC_TERMS = (
    "colorectal", "colon cancer", "rectal cancer", "colitis-associated",
    "crc", "cac", "\u7ed3\u76f4\u80a0", "\u7ed3\u80a0\u764c", "\u7ed3\u80a0\u708e\u76f8\u5173",
    "hct116", "hct-116", "ht29", "ht-29", "sw480", "sw620", "lovo",
    "dld-1", "dld1", "caco-2", "ct26", "mc38",
)
ANIMAL_TERMS = (
    "xenograft", "syngeneic", "aom", "dss", "mouse", "mice", "rat",
    "\u5c0f\u9f20", "\u5927\u9f20", "\u88f8\u9f20", "\u79fb\u690d\u7624", "\u4f53\u5185", "\u52a8\u7269",
)
CELL_TERMS = (
    "cell", "hct116", "hct-116", "ht29", "ht-29", "sw480", "sw620",
    "lovo", "dld-1", "dld1", "caco-2", "ct26", "mc38", "\u7ec6\u80de", "\u4f53\u5916",
)
COMPUTATIONAL_TERMS = (
    "network pharmacology", "molecular docking", "in silico",
    "\u7f51\u7edc\u836f\u7406", "\u5206\u5b50\u5bf9\u63a5", "\u7eaf\u8ba1\u7b97", "\u8ba1\u7b97\u673a\u6a21\u62df",
)
REVIEW_TERMS = (
    "review", "perspective", "\u7efc\u8ff0", "advancement",
    "herbal medicine for colorectal cancer", "pharmacological potential",
)
MIXTURE_TERMS = (
    "formula", "decoction", "granules", "powder", "extract", "essential oil",
    "prescription", "mixture", "green tea catechins", "ginseng, quercetin, and tea",
    "\u590d\u65b9", "\u6c64", "\u6563", "\u9897\u7c92", "\u63d0\u53d6\u7269",
)

PATHWAY_RULES = {
    "Apoptosis/cell cycle": ("apoptosis", "caspase", "bax", "bcl-2", "\u51cb\u4ea1", "\u5468\u671f"),
    "PI3K/Akt/mTOR/metabolism": ("pi3k", "akt", "mtor", "pkm2", "glycol", "\u7cd6\u9175\u89e3"),
    "Wnt/beta-catenin/stemness": ("wnt", "beta-catenin", "\u03b2-catenin", "stem", "lef1", "\u5e72\u6027"),
    "NF-kB/STAT3/inflammation": ("nf-\u03bab", "nf-kb", "nf-kappa", "stat3", "inflamm", "\u708e\u75c7"),
    "MAPK/JNK/ERK": ("mapk", "jnk", "erk", "p38"),
    "ROS/ferroptosis": ("ferropt", "gpx4", "slc7a11", "ros", "oxidative", "\u94c1\u6b7b\u4ea1", "\u6c27\u5316"),
    "EMT/ECM/metastasis": ("emt", "mmp", "migration", "invasion", "metasta", "\u8fc1\u79fb", "\u4fb5\u88ad", "\u8f6c\u79fb"),
    "Immune microenvironment": ("macrophage", "t cell", "pd-1", "pd-l1", "immune", "tam", "\u5de8\u566c", "\u514d\u75ab"),
}

# Only values directly checked in local PDFs are added here. They supplement,
# but never replace, the original extraction workbook.
VERIFIED_DOSE = {
    "33757400": "10-20 \u00b5M luteolin (HCT116; PDF abstract)",
    "35453311": "12.5 or 25 \u00b5M luteolin with 1 \u00b5M oxaliplatin in vitro; luteolin 50 mg/kg/day with oxaliplatin 10 mg/kg three times/week in xenograft mice",
    "31180555": "100 \u00b5M kaempferol with 50 \u00b5M 5-FU for 48 h",
    "33535870": "AOM 10 mg/kg and 2.5% DSS; tanshinone IIA 200 mg/kg intraperitoneally",
    "41051621": "Tanshinone IIA 5 mg/kg every 2 days intraperitoneally in xenograft mice",
    "39276049": "Sanguinarine concentration series: SW620 0-2 \u00b5mol/L; HCT-116 0-3 \u00b5mol/L",
    "39310106": "20 \u00b5M honokiol for 48 h in the initial CRC cell screen",
    "30199885": "Aloe-emodin 10, 20, and 40 \u00b5M in SW620 and HT29 cells",
    "37735401": "Curdione 12.5, 25, and 50 \u00b5M in CRC cells",
    "40469706": "Rhein 10, 20, and 50 \u00b5M in HT-29 and SW480 cells",
}
VERIFIED_SAFETY = {
    "33946531": "No obvious toxicity reported in the HCT116 xenograft experiment (PDF abstract).",
}
VERIFIED_MODEL = {
    "35408903": "HCT8 and 5-FU-resistant HCT8-R colorectal cancer cells (verified in local PDF)",
    "40812220": "Colon cancer cells and MC38 syngeneic mouse model (verified in local PDF)",
}
VERIFIED_TITLE = {
    "39310106": "Honokiol enhances the sensitivity of cetuximab in KRASG13D mutant colorectal cancer by destroying SNX3-retromer.",
    "34192810": "Honokiol inhibits proliferation of colorectal cancer cells by targeting Anoctamin 1/TMEM16A Ca2+-activated Cl- channels.",
    "34309458": "Quercetin and Luteolin Improve the Anticancer Effects of 5-Fluorouracil in Human Colorectal Adenocarcinoma In Vitro Model: A Mechanistic Insight.",
    "31164916": "A tannin compound from Sanguisorba officinalis blocks Wnt/beta-catenin signaling pathway and induces apoptosis of colorectal cancer cells.",
}


def text(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def normalized_name(value: object) -> str:
    return text(value).strip().lower()


def find_pdf_index() -> dict[str, str]:
    result: dict[str, str] = {}
    for pdf in ROOT.rglob("*.pdf"):
        if PAPER in pdf.parents or pdf.parent.name in {"figures", "bst"}:
            continue
        match = re.search(r"(?<!\d)(\d{8})(?!\d)", pdf.name)
        if match and match.group(1) not in result:
            result[match.group(1)] = str(pdf.relative_to(ROOT))
    return result


def parse_bib(bib_text: str) -> dict[str, str]:
    entries: dict[str, str] = {}
    for block in re.findall(r"@article\{PMID\d+,.*?\n\}", bib_text, flags=re.S):
        key_match = re.search(r"@article\{(PMID\d+),", block)
        if key_match:
            entries[key_match.group(1).replace("PMID", "")] = block
    return entries


def bib_title(block: str) -> str:
    match = re.search(r"title\s*=\s*\{(.*?)\},\n", block, flags=re.S)
    return re.sub(r"\s+", " ", match.group(1)).strip() if match else ""


def complete_title(title: str) -> bool:
    low = title.lower()
    if len(title) < 25 or "pubmed record" in low:
        return False
    if not title.endswith((".", "?", "!")):
        return False
    return not any(low.endswith(suffix) for suffix in (" of", " the", " a", " and", ":"))


def captured_dose(joined: str, pmid: str) -> str:
    if pmid in VERIFIED_DOSE:
        return VERIFIED_DOSE[pmid]
    units = r"(?:\u00b5M|\u03bcM|uM|mM|\u03bcmol/L|\u00b5mol/L|mg/kg|mg/mL|\u03bcg/mL|%)"
    matches = re.findall(rf"\d+(?:\.\d+)?(?:\s*[-\u2013]\s*\d+(?:\.\d+)?)?\s*{units}", joined, flags=re.I)
    return "; ".join(dict.fromkeys(matches)) if matches else "not captured in source extraction"


def model_group(model: str) -> str:
    low = model.lower()
    animal = any(term in low for term in ANIMAL_TERMS)
    cell = any(term in low for term in CELL_TERMS)
    computational = any(term in low for term in COMPUTATIONAL_TERMS)
    if animal and cell:
        return "In vitro + in vivo"
    if animal:
        return "In vivo"
    if cell:
        return "In vitro"
    if computational:
        return "Network/in silico"
    return "Not specified"


def pathway_modules(joined: str) -> str:
    low = joined.lower()
    found = [name for name, terms in PATHWAY_RULES.items() if any(term in low for term in terms)]
    return "; ".join(found)


def strength(model: str, pathway: str, phenotype: str, title: str) -> str:
    joined = f"{model} {pathway} {phenotype} {title}".lower()
    mgroup = model_group(model)
    title_lower = title.lower()
    if any(term in title_lower for term in REVIEW_TERMS):
        return "Contextual review-only"
    if any(term in title_lower for term in MIXTURE_TERMS):
        return "Contextual mixture/formulation evidence"
    if mgroup == "Network/in silico":
        return "Level 4: in silico-only"
    informative = bool(pathway or phenotype)
    if mgroup == "In vitro + in vivo" and informative:
        return "Level 1: animal-validated mechanistic evidence"
    if mgroup == "In vivo" and informative:
        return "Level 2: in vivo evidence"
    if mgroup == "In vitro" and informative:
        return "Level 3: cell-based mechanistic evidence"
    if any(term in joined for term in COMPUTATIONAL_TERMS):
        return "Level 4: in silico/contextual evidence"
    return "Excluded: insufficient extractable evidence"


def build_table() -> tuple[pd.DataFrame, dict[str, str]]:
    if not SOURCE_BIB_BACKUP.exists():
        shutil.copy2(SOURCE_BIB, SOURCE_BIB_BACKUP)
    bib_entries = parse_bib(SOURCE_BIB_BACKUP.read_text(encoding="utf-8"))
    pdf_index = find_pdf_index()
    raw = pd.read_excel(SOURCE_XLSX, sheet_name=0)
    rows = []
    for _, row in raw.iterrows():
        monomer = text(row.iloc[1])
        source = text(row.iloc[2])
        model = text(row.iloc[3])
        pathway = text(row.iloc[4])
        phenotype = text(row.iloc[5])
        target = text(row.iloc[6])
        clinical = text(row.iloc[7])
        safety = text(row.iloc[8])
        pmid = re.sub(r"\D", "", text(row.iloc[10]))
        block = bib_entries.get(pmid, "")
        title = VERIFIED_TITLE.get(pmid, bib_title(block))
        if pmid in VERIFIED_MODEL:
            model = VERIFIED_MODEL[pmid]
        joined = " ".join([monomer, source, model, pathway, phenotype, target, title])
        relevant = any(term in joined.lower() for term in CRC_TERMS)
        pdf_path = pdf_index.get(pmid, "")
        evidence_strength = strength(model, pathway, phenotype, title)
        usable_bib = bool(block) and complete_title(title)
        include = (
            relevant and bool(pdf_path) and usable_bib
            and not evidence_strength.startswith("Excluded")
        )
        if not relevant:
            reason = "Excluded: no CRC-specific model or outcome extracted"
        elif not pdf_path:
            reason = "Excluded: no local PDF matched"
        elif not usable_bib:
            reason = "Excluded: incomplete bibliographic title requiring verification"
        elif evidence_strength.startswith("Excluded"):
            reason = evidence_strength
        else:
            reason = "Included"
        rows.append({
            "Monomer": monomer.strip(),
            "Class": CATEGORY.get(normalized_name(monomer), "Other natural products"),
            "Source": source or "not reported",
            "Model": model or "not reported",
            "Model type": model_group(model),
            "Dose/concentration/administration": captured_dose(joined, pmid),
            "Pathway/target (source extraction)": pathway or target or "not reported",
            "Phenotype/result (source extraction)": phenotype or "not reported",
            "In vivo evidence": "yes" if "In vivo" in model_group(model) else "not captured",
            "Safety/toxicity": VERIFIED_SAFETY.get(pmid, safety or "not reported"),
            "Clinical evidence": clinical or "not reported",
            "Evidence level": evidence_strength,
            "Mechanism module tags": pathway_modules(joined) or "not assigned",
            "PMID": pmid or "not reported",
            "BibTeX key": f"PMID{pmid}" if usable_bib else "not available",
            "Local PDF": pdf_path or "not matched",
            "Bibliographic title": title if usable_bib else "not retained - incomplete bibliographic metadata",
            "Inclusion status": reason,
        })
    result = pd.DataFrame(rows)
    result.to_csv(PAPER / "evidence_extraction_table.csv", index=False, encoding="utf-8-sig")
    result.to_excel(PAPER / "evidence_extraction_table.xlsx", index=False)
    return result, bib_entries


def write_bibliography(table: pd.DataFrame, bib_entries: dict[str, str]) -> int:
    included_pmids = sorted(set(table.loc[table["Inclusion status"] == "Included", "PMID"]))
    blocks = []
    for pmid in included_pmids:
        if pmid not in bib_entries:
            continue
        block = bib_entries[pmid]
        if pmid in VERIFIED_TITLE:
            block = re.sub(
                r"title\s*=\s*\{.*?\},\n",
                f"  title = {{{VERIFIED_TITLE[pmid]}}},\n",
                block,
                flags=re.S,
            )
        blocks.append(block)
    SOURCE_BIB.write_text("\n\n".join(blocks) + "\n", encoding="utf-8")
    return len(blocks)


def plot_bar(series: pd.Series, title: str, xlabel: str, filename: str, color: str, limit: int | None = None) -> None:
    if limit:
        series = series.head(limit)
    fig, ax = plt.subplots(figsize=(8.0, max(4.0, len(series) * 0.30 + 1.2)))
    ordered = series.sort_values()
    ordered.plot(kind="barh", ax=ax, color=color)
    ax.set_title(title, fontsize=11, pad=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("")
    ax.spines[["top", "right"]].set_visible(False)
    for y, val in enumerate(ordered.values):
        ax.text(val + max(series.max() * 0.01, 0.1), y, str(int(val)), va="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / filename, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_figures(table: pd.DataFrame) -> None:
    FIGURES.mkdir(exist_ok=True)
    included = table[table["Inclusion status"] == "Included"].copy()
    direct = included[~included["Evidence level"].str.contains("Contextual|Level 4")]
    plot_bar(direct["Class"].value_counts(), "Chemical classes in direct CRC evidence", "Included records", "chemical_category_distribution.png", "#376996")
    plot_bar(direct["Model type"].value_counts(), "Model types in direct CRC evidence", "Included records", "model_type_distribution.png", "#4d8a65")
    modules = Counter()
    for tags in direct["Mechanism module tags"]:
        for tag in tags.split("; "):
            if tag != "not assigned":
                modules[tag] += 1
    plot_bar(pd.Series(modules), "Mechanistic annotations in direct CRC evidence", "Records with annotation", "pathway_frequency_review.png", "#a2534c")
    plot_bar(included["Evidence level"].value_counts(), "Evidence-strength distribution", "Included records", "evidence_strength_distribution.png", "#8a6eac")
    plot_bar(direct["Monomer"].value_counts(), "Top monomers by direct CRC-related records", "Included records", "top_monomers_crc_records.png", "#c27a30", limit=15)


def main() -> None:
    table, entries = build_table()
    bib_count = write_bibliography(table, entries)
    write_figures(table)
    included = table[table["Inclusion status"] == "Included"]
    direct = included[~included["Evidence level"].str.contains("Contextual|Level 4")]
    summary = [
        f"Source workbook records: {len(table)}",
        f"Included evidence records: {len(included)}",
        f"Direct experimental records: {len(direct)}",
        f"Curated BibTeX entries: {bib_count}",
    ]
    (PAPER / "evidence_build_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
