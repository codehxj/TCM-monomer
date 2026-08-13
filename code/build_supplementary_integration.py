from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import pandas as pd


PAPER = Path(__file__).resolve().parent
FIGURES = PAPER / "figures"
INPUT = PAPER / "evidence_extraction_table.csv"

FORMULATION = (
    "nanopart", "nano ", "nanostruct", "hydrogel", "drug delivery",
    "delivery system", "encapsulat", "niosom", "liposome", "film",
    "composite", "nanoplatform", "nanosystem", "hydroxyapatite",
)
EXTRACT_FORMULA = (
    "decoction", "formula", "granule", "powder", "extract", "essential oil",
    "capsule", "churna", "green tea catechins", "juice concentrate",
    "total flavonoids", "honokiol-magnolol-baicalin", "magnolia officinalis",
    "zuojin", "shancigu", "herba patriniae", "lignans and polyphenols",
    "dietary flavonols", "phenolic compounds of manuka", "walnut constituents",
    "pomegranate juice", "arecaceae seeds", "micronutrient-enriched",
    "quercetin, luteolin, and xanthohumol", "sorafenib and plant-derived",
)
IN_SILICO = (
    "network pharmacology", "molecular docking", "in silico", "in-silico",
    "bioinformatics", "reverse pharmacology network", "targets prediction",
)
REVIEW = (
    "review", "perspective", "therapeutic potential", "herbal medicine",
    "targeting colorectal cancer using dietary", "targeting cancer stem cells and signalling",
    "interplay between traditional", "natural alkaloids in cancer therapy",
    "modulation of the canonical wnt", "plausible paradigm",
    "the wnt/beta-catenin signaling pathway", "tumor-associated macrophages",
)
COMBINATION = (
    "5-fluorouracil", "5fluorouracil", "5-fu", "oxaliplatin", "cisplatin",
    "irinotecan", "cetuximab", "anti-pd-1", "checkpoint", "erastin",
    "chemoresistance", "drug-induced", "sensit", "synerg",
)
IMMUNE = (
    "macrophage", "tam", "immune", "pd-1", "pd-l1", "cd8", "extracellular vesicle",
    "microenvironment",
)
FERROPTOSIS = (
    "ferropt", "gpx4", "slc7a11", "xct", "ferritinophagy", "lipid peroxid",
    "oxeiptosis", "ros", "oxidative",
)


def has_any(value: str, terms: tuple[str, ...]) -> bool:
    low = value.lower()
    return any(term in low for term in terms)


def yn(flag: bool) -> str:
    return "Yes" if flag else "No"


def classify(group: pd.DataFrame) -> dict[str, object]:
    first = group.iloc[0]
    title = str(first["Bibliographic title"])
    merged = " ".join(
        group[col].fillna("").astype(str).str.cat(sep=" ")
        for col in [
            "Bibliographic title", "Source", "Model", "Evidence level",
            "Mechanism module tags", "Pathway/target (source extraction)",
            "Phenotype/result (source extraction)",
        ]
    )
    formulation = has_any(title, FORMULATION)
    extract_formula = has_any(title, EXTRACT_FORMULA)
    in_silico = has_any(title, IN_SILICO) or all(
        "Level 4" in x for x in group["Evidence level"].astype(str)
    )
    review = has_any(title, REVIEW) or all(
        "review-only" in x for x in group["Evidence level"].astype(str)
    )
    contextual = formulation or extract_formula or in_silico or review
    experimental = any(
        str(level).startswith(("Level 1", "Level 2", "Level 3"))
        for level in group["Evidence level"]
    )
    direct = experimental and not contextual
    model_text = " ".join(group["Model"].fillna("").astype(str)).lower()
    explicit_animal = any(term in model_text for term in (
        "xenograft", "syngeneic", "aom/dss", "aom+dss", "apcmin",
        "mouse model", "mice", "rat model", "动物模型", "移植瘤",
        "小鼠模型", "大鼠", "裸鼠", "体内模型", "肝转移模型",
    ))
    animal = direct and explicit_animal
    cell_only = direct and not animal
    combination = has_any(merged, COMBINATION)
    immune = has_any(merged, IMMUNE)
    ferroptosis = has_any(merged, FERROPTOSIS)
    if formulation:
        primary = "Formulation/nanodelivery contextual evidence"
        use = "Context only; the delivery system cannot be attributed to isolated monomer activity."
    elif extract_formula:
        primary = "Extract/formula/decoction contextual evidence"
        use = "Context only; preparation effects are not assigned to an isolated monomer."
    elif in_silico:
        primary = "Network pharmacology/in silico evidence"
        use = "Hypothesis generation only unless supported by a separate direct record."
    elif review:
        primary = "Review/background evidence"
        use = "Background synthesis only; not counted as new efficacy evidence."
    elif animal:
        primary = "Animal-validated isolated monomer evidence"
        use = "High-information direct preclinical evidence."
    elif cell_only:
        primary = "Cell-only isolated monomer evidence"
        use = "Direct mechanistic evidence requiring in vivo validation."
    else:
        primary = "Other CRC-relevant contextual evidence"
        use = "Retained for mapping; not used for an isolated-monomer conclusion."
    return {
        "PMID": first["PMID"],
        "BibTeX key": first["BibTeX key"],
        "Monomer(s) mapped in source table": "; ".join(sorted(set(group["Monomer"].astype(str)))),
        "Class(es)": "; ".join(sorted(set(group["Class"].astype(str)))),
        "Bibliographic title": title,
        "Primary evidence classification": primary,
        "Direct isolated monomer experimental evidence": yn(direct),
        "Animal-validated monomer evidence": yn(animal),
        "Cell-only monomer evidence": yn(cell_only),
        "Combination therapy evidence": yn(combination),
        "Immune microenvironment evidence": yn(immune),
        "Ferroptosis/redox evidence": yn(ferroptosis),
        "Formulation/nanodelivery evidence": yn(formulation),
        "Extract/formula/decoction contextual evidence": yn(extract_formula),
        "Network pharmacology/in silico evidence": yn(in_silico),
        "Review/background evidence": yn(review),
        "Models recorded": " | ".join(group["Model"].dropna().astype(str).drop_duplicates()),
        "Exposure recorded": " | ".join(group["Dose/concentration/administration"].dropna().astype(str).drop_duplicates()),
        "Mechanism tags": "; ".join(sorted(set(
            tag for tags in group["Mechanism module tags"].fillna("")
            for tag in str(tags).split("; ") if tag and tag != "not assigned"
        ))),
        "Use in review": use,
    }


def make_supplement() -> pd.DataFrame:
    data = pd.read_csv(INPUT)
    retained = data[data["Inclusion status"].eq("Included")].copy()
    rows = [classify(group) for _, group in retained.groupby("PMID", sort=True)]
    result = pd.DataFrame(rows)
    result.to_csv(PAPER / "supplementary_evidence_table.csv", index=False, encoding="utf-8-sig")
    with pd.ExcelWriter(PAPER / "supplementary_evidence_table.xlsx", engine="openpyxl") as writer:
        result.to_excel(writer, sheet_name="Table S1 evidence map", index=False)
        summary = result["Primary evidence classification"].value_counts().rename_axis("Classification").reset_index(name="References")
        summary.to_excel(writer, sheet_name="Classification summary", index=False)
    return result


def box(ax, xy, width, height, text, color, fontsize=9, edge="#27435c"):
    x, y = xy
    patch = FancyBboxPatch(
        (x, y), width, height, boxstyle="round,pad=0.015,rounding_size=0.015",
        facecolor=color, edgecolor=edge, linewidth=1.15
    )
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height / 2, text, ha="center", va="center", fontsize=fontsize, color="#10202c")


def arrow(ax, start, end):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=14, lw=1.2, color="#536878"))


def make_flow(data: pd.DataFrame, supplement: pd.DataFrame) -> None:
    excluded = data[~data["Inclusion status"].eq("Included")]["Inclusion status"].value_counts()
    strict_direct = (supplement["Direct isolated monomer experimental evidence"] == "Yes").sum()
    fig, ax = plt.subplots(figsize=(12.6, 8.4))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    box(ax, (0.31, 0.88), 0.38, 0.075, "303 source records\n(workbook rows linked to local collection)", "#dcebf4", 10)
    box(ax, (0.31, 0.76), 0.38, 0.075, "Screening sequence\ndeduplication/metadata -> CRC relevance -> PDF match -> title completeness", "#eaf2f7", 9)
    arrow(ax, (0.50, 0.88), (0.50, 0.837))
    box(ax, (0.34, 0.63), 0.32, 0.075, "217 retained records\n168 unique PMID-linked references", "#cae7db", 10)
    arrow(ax, (0.50, 0.76), (0.50, 0.707))
    box(ax, (0.015, 0.62), 0.27, 0.115,
        "86 excluded records\n37 no CRC-specific outcome\n27 incomplete metadata\n16 insufficient extractable information\n6 no local PDF", "#f8dddd", 8.5)
    arrow(ax, (0.31, 0.795), (0.285, 0.68))
    box(ax, (0.06, 0.39), 0.40, 0.135,
        "141 provisionally direct experimental\nmonomer-annotated records\n51 animal-validated mechanistic\n13 additional in vivo\n77 cell-based mechanistic", "#cfe5f3", 9)
    box(ax, (0.54, 0.39), 0.40, 0.135,
        "76 contextual records\n45 mixture/formulation\n22 in silico/contextual\n9 review-only", "#f6e7c9", 9)
    arrow(ax, (0.43, 0.63), (0.27, 0.527))
    arrow(ax, (0.57, 0.63), (0.73, 0.527))
    box(ax, (0.17, 0.18), 0.66, 0.115,
        f"Supplementary Table S1 attribution audit: {strict_direct} references support\ndirect isolated-monomer interpretation after formulation/extract/review/in silico flags;\nborderline records remain visible as contextual evidence.", "#edf1f3", 9)
    arrow(ax, (0.27, 0.39), (0.38, 0.297))
    arrow(ax, (0.73, 0.39), (0.62, 0.297))
    ax.text(0.5, 0.985, "Evidence curation and attribution workflow", ha="center", va="top", fontsize=15, weight="bold")
    ax.text(0.5, 0.08, "Counts represent records or retained references, not effect size or therapeutic efficacy.", ha="center", fontsize=9, color="#536878")
    fig.tight_layout()
    fig.savefig(FIGURES / "evidence_curation_flow.png", dpi=320, bbox_inches="tight")
    plt.close(fig)


def make_mechanism_panels() -> None:
    fig, ax = plt.subplots(figsize=(13, 6.5))
    ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0.5, 0.96, "A. Signaling modules with directly discussed monomer evidence", ha="center", fontsize=14, weight="bold")
    modules = [
        ("PI3K/Akt/mTOR\nmetabolism", "kaempferol; rhein;\ntanshinone IIA; formononetin"),
        ("Wnt/beta-catenin\nstemness", "beta-sitosterol; silibinin;\nkaempferol; galloyl-glucose"),
        ("NF-kB/STAT3\ninflammation", "tanshinone IIA; rhein;\nformononetin; glabridin*"),
        ("MAPK/JNK/ERK\nstress", "tanshinone IIA; luteolin;\ndictamnine; beta-lapachone"),
        ("GPX4/SLC7A11\nredox death", "luteolin; sanguinarine;\ncurdione; tagitinin C"),
    ]
    for i, (head, names) in enumerate(modules):
        x = 0.03 + i * 0.194
        box(ax, (x, 0.56), 0.165, 0.18, head, "#dcebf4", 9)
        box(ax, (x, 0.25), 0.165, 0.18, names, "#eaf2f7", 8)
        arrow(ax, (x + 0.082, 0.56), (x + 0.082, 0.43))
    ax.text(0.03, 0.10, "* Glabridin is retained as preparation-associated contextual evidence in the corpus.", fontsize=9, color="#536878")
    fig.tight_layout()
    fig.savefig(FIGURES / "mechanism_signaling_pathways.png", dpi=320, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(13, 6.5))
    ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0.5, 0.96, "B. Phenotypic and tumor-microenvironment endpoints", ha="center", fontsize=14, weight="bold")
    outcomes = [
        ("Apoptosis/cell cycle", "luteolin; aloe-emodin;\nkaempferol; tanshinone IIA"),
        ("Chemotherapy\nsensitization", "luteolin; kaempferol;\ntanshinone IIA; beta-sitosterol"),
        ("Migration/EMT/\nmetastasis", "rhein; pinocembrin;\nastragaloside IV; kaempferol"),
        ("Macrophage/TME", "astragaloside IV; luteolin;\nmacrocarpal I; dictamnine"),
        ("Immunotherapy\ncombination", "macrocarpal I;\nH-M-B mixture*"),
    ]
    for i, (head, names) in enumerate(outcomes):
        x = 0.03 + i * 0.194
        box(ax, (x, 0.56), 0.165, 0.18, head, "#d7eadf", 9)
        box(ax, (x, 0.25), 0.165, 0.18, names, "#eef6f1", 8)
        arrow(ax, (x + 0.082, 0.56), (x + 0.082, 0.43))
    ax.text(0.03, 0.10, "* Honokiol-magnolol-baicalin combination is contextual, not isolated honokiol evidence.", fontsize=9, color="#536878")
    fig.tight_layout()
    fig.savefig(FIGURES / "mechanism_phenotypes_tme.png", dpi=320, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIGURES.mkdir(exist_ok=True)
    raw = pd.read_csv(INPUT)
    supplement = make_supplement()
    make_flow(raw, supplement)
    make_mechanism_panels()
    print("Supplementary Table S1 references:", len(supplement))
    print(supplement["Primary evidence classification"].value_counts().to_string())
    print("Direct isolated-monomer references after attribution audit:",
          (supplement["Direct isolated monomer experimental evidence"] == "Yes").sum())


if __name__ == "__main__":
    main()
