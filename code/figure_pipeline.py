
from __future__ import annotations

import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
SUPPLEMENTARY = ROOT / "supplementary"
DATA = ROOT / "data"
FIGURES = ROOT / "submission_figures"

EXTRACTION_FILE = SUPPLEMENTARY / "evidence_extraction_table.csv"
POOL_FILE = SUPPLEMENTARY / "evidence_pool.csv"

BLUE = "#0072B2"
SKY = "#56B4E9"
ORANGE = "#E69F00"
GREEN = "#009E73"
VERMILLION = "#D55E00"
PURPLE = "#CC79A7"
GRAY = "#666666"
LIGHT_GRAY = "#E6E6E6"

CORE_COMPOUNDS = ["Luteolin", "Kaempferol", "Tanshinone IIA", "Astragaloside IV", "Honokiol"]
COMPARATOR = "Macrocarpal I"
COMPOUNDS_WITH_COMPARATOR = CORE_COMPOUNDS + [COMPARATOR]

TABLE1_CLASS_ORDER = [
    "Flavonoids",
    "Neolignans",
    "Diterpenoid quinones",
    "Polyphenols",
    "Phytosterols",
    "Phenolic acids",
    "Anthraquinones",
    "Alkaloids",
    "Other classes",
]

DOMAIN_RULES = {
    "Treatment resistance\nor sensitization": (
        "5-fluorouracil", "5-fu", "fluorouracil", "oxaliplatin", "cetuximab",
        "chemoresistan", "drug resistan", "resistant", "sensiti", "synerg", "erastin",
    ),
    "Metastatic or\nimmune niche": (
        "metasta", "macrophage", "immune", "extracellular vesicle", "m2", "cd8", "emt", "liver", "hepatic",
    ),
    "Inflammation-associated\ntumorigenesis": (
        "aom", "dss", "colitis", "inflamm", "cac", "nf-kb", "nf-?b", "barrier",
    ),
}

TREATMENT_RULES = {
    "5-FU": ("5-fluorouracil", "5-fu", "fluorouracil"),
    "Oxaliplatin": ("oxaliplatin",),
    "Cetuximab": ("cetuximab",),
    "Immune-checkpoint\ntherapy": ("anti-pd-1", "anti-pd1", "pd-1", "checkpoint", "pembrolizumab"),
}

MODEL_RULES = {
    "Conventional xenograft": ("xenograft", "nude mouse", "nude mice", "??", "???"),
    "Syngeneic / immune-competent": ("syngeneic", "immune-competent", "immunocompetent", "mc38", "ct26"),
    "AOM/DSS or inherited tumor model": ("aom", "dss", "apcmin", "apc min"),
    "Organ-specific / liver-metastasis": (
        "liver metast", "hepatic metast", "spleen-liver", "spleen liver",
        "orthotopic", "\u809d\u8f6c\u79fb", "\u810f\u5668\u7279\u5f02",
        "\u813e\u810f-\u809d\u810f", "\u813e\u810f\u809d\u810f",
    ),
    "Patient-derived CRC organoid": ("patient-derived organoid", "patient derived organoid", "crc organoid", "???????"),
}

POSITIVE_DESIGN_FEATURE_RULES = {
    "Conventional xenograft": MODEL_RULES["Conventional xenograft"],
    "Immune-competent model": MODEL_RULES["Syngeneic / immune-competent"],
    "AOM/DSS or inherited tumor model": MODEL_RULES["AOM/DSS or inherited tumor model"],
    "Matched normal-intestinal comparison": (
        "normal intestinal", "normal colon", "normal colonic", "fhc", "ncm460", "ccd 841", "normal organoid",
    ),
    "Organ-specific / liver-metastasis model": MODEL_RULES["Organ-specific / liver-metastasis"],
}


def _contains(text: str, terms: tuple[str, ...]) -> bool:
    low = text.lower()
    return any(term in low for term in terms)


def _normalise_compound(value: object) -> str:
    text = str(value or "").lower().replace("?", "ii").replace("?", "beta")
    return re.sub(r"\s+", " ", text).strip()


def _compound_mask(series: pd.Series, compound: str) -> pd.Series:
    key = _normalise_compound(compound)
    return series.fillna("").map(_normalise_compound).map(
        lambda value: key in [part.strip() for part in re.split(r"[;|]", value)]
    )


def _combined_text(frame: pd.DataFrame) -> pd.Series:
    cols = [
        "Bibliographic title", "Models recorded", "Exposure recorded", "Mechanism tags",
        "Use in review", "Material actually tested", "Source medicinal material",
    ]
    present = [c for c in cols if c in frame.columns]
    return frame[present].fillna("").astype(str).agg(" ".join, axis=1)


def load_and_validate() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    extraction = pd.read_csv(EXTRACTION_FILE)
    pool = pd.read_csv(POOL_FILE)
    retained = extraction[extraction["Inclusion status"].eq("Included")].copy()
    direct = pool[pool["Direct isolated monomer experimental evidence"].eq("Yes")].copy()
    animal = int(pool["Animal-validated monomer evidence"].eq("Yes").sum())
    cell_only = int(pool["Cell-only monomer evidence"].eq("Yes").sum())
    first_pass = int(retained["Evidence level"].str.match(r"^Level [123]").sum())
    contextual_unique = len(pool) - len(direct)

    expected = {
        "initial row-level records": (len(extraction), 303),
        "excluded row-level records": (len(extraction) - len(retained), 86),
        "retained row-level records": (len(retained), 217),
        "unique PMID-linked references": (len(pool), 168),
        "first-pass row-level annotations": (first_pass, 141),
        "contextual row-level annotations": (len(retained) - first_pass, 76),
        "direct unique references": (len(direct), 77),
        "animal-validated unique references": (animal, 29),
        "cell-only unique references": (cell_only, 48),
        "contextual unique references": (contextual_unique, 91),
    }
    mismatches = [f"{name}: observed {obs}, expected {exp}" for name, (obs, exp) in expected.items() if obs != exp]
    if mismatches:
        raise ValueError("Evidence-count validation failed:\n" + "\n".join(mismatches))
    if pool["PMID"].nunique() != 168:
        raise ValueError("evidence_pool.csv must contain one row per unique PMID.")
    DATA.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    return extraction, retained, pool


def set_style() -> None:
    sns.set_theme(style="whitegrid", context="paper")
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.unicode_minus": False,
        "figure.dpi": 150,
        "savefig.dpi": 600,
    })


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURES / f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def _box(ax, xy, width, height, text, facecolor, edgecolor="#333333", fontsize=8):
    patch = FancyBboxPatch(
        xy, width, height, boxstyle="round,pad=0.012,rounding_size=0.012",
        facecolor=facecolor, edgecolor=edgecolor, linewidth=1, zorder=2,
    )
    ax.add_patch(patch)
    ax.text(xy[0] + width/2, xy[1] + height/2, text, ha="center", va="center", fontsize=fontsize, zorder=3)


def _arrow(ax, start, end, mutation_scale=9):
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=mutation_scale, linewidth=0.9,
        color="#444444", shrinkA=7, shrinkB=7, zorder=1, clip_on=False,
    ))


def _connector(ax, points):
    x, y = zip(*points)
    ax.plot(x, y, color="#555555", linewidth=0.9, solid_capstyle="round", zorder=1)


def panel_label(ax: plt.Axes, label: str, x: float = -0.08, y: float = 1.05) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left", clip_on=False)


def generate_evidence_curation_flow() -> None:
    extraction, retained, pool = load_and_validate()
    first_pass = retained[retained["Evidence level"].str.match(r"^Level [123]")]
    direct = pool[pool["Direct isolated monomer experimental evidence"].eq("Yes")]
    context = pool[pool["Direct isolated monomer experimental evidence"].ne("Yes")]
    flow_data = pd.DataFrame([
        ("Initial row-level records", len(extraction)),
        ("Excluded row-level records", len(extraction) - len(retained)),
        ("Retained row-level records", len(retained)),
        ("Unique PMID-linked references", len(pool)),
        ("First-pass isolated-compound-annotated rows", len(first_pass)),
        ("Contextual row-level annotations", len(retained) - len(first_pass)),
        ("Direct purified-compound unique references", len(direct)),
        ("Animal-validated direct unique references", int(direct["Animal-validated monomer evidence"].eq("Yes").sum())),
        ("Cell-only direct unique references", int(direct["Cell-only monomer evidence"].eq("Yes").sum())),
        ("Contextual unique references", len(context)),
    ], columns=["Stage", "Count"])
    flow_data.to_csv(DATA / "figure_1_evidence_flow.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.48, 5.05))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(
        0.5, 0.975, "Evidence curation and source-level attribution workflow",
        ha="center", va="top", fontsize=10, fontweight="bold",
    )

    def flow_arrow(start, end):
        ax.add_patch(FancyArrowPatch(
            start, end, arrowstyle="-|>", mutation_scale=9.5,
            linewidth=1.05, color="#444444", shrinkA=0, shrinkB=0,
            zorder=4, clip_on=False,
        ))

    def group_frame(x, y, width, height, title):
        frame = FancyBboxPatch(
            (x, y), width, height,
            boxstyle="round,pad=0.008,rounding_size=0.010",
            facecolor="#FBFCFD", edgecolor="#8A8A8A",
            linewidth=0.85, zorder=0,
        )
        ax.add_patch(frame)
        ax.text(
            x + width / 2, y + height - 0.025, title,
            ha="center", va="center", fontsize=7.4,
            fontweight="bold", color="#333333", zorder=3,
        )

    _box(ax, (0.33, 0.845), 0.34, 0.075,
         "303 initial row-level records", "#DCEAF4", BLUE, 8.0)
    _box(ax, (0.22, 0.695), 0.56, 0.095,
         "Metadata verification, CRC relevance,\nsource matching, and eligibility review",
         "#F3F6F8", "#333333", 7.6)
    _box(ax, (0.82, 0.705), 0.15, 0.075,
         "86 excluded\nrows", "#F7D7D4", VERMILLION, 7.3)
    _box(ax, (0.28, 0.565), 0.44, 0.080,
         "217 retained rows\n168 unique PMID-linked references",
         "#D9EEE7", GREEN, 7.7)

    group_frame(0.055, 0.335, 0.89, 0.165, "Row-level evidence annotation")
    _box(ax, (0.080, 0.355), 0.390, 0.095,
         "141 first-pass isolated-compound\nrow annotations",
         "#DCEAF4", BLUE, 7.3)
    _box(ax, (0.530, 0.355), 0.390, 0.095,
         "76 contextual row annotations",
         "#FBE8C6", ORANGE, 7.4)

    _box(ax, (0.30, 0.230), 0.40, 0.060,
         "Strict reference-level\nmaterial-attribution audit",
         "#EDF1F3", GRAY, 7.4)

    group_frame(0.055, 0.015, 0.89, 0.165, "Reference-level source attribution")
    _box(ax, (0.080, 0.035), 0.390, 0.095,
         "77 direct purified-compound references\n29 animal-validated + 48 cell-only",
         "#D9EEE7", GREEN, 7.1)
    _box(ax, (0.530, 0.035), 0.390, 0.095,
         "91 contextual unique references",
         "#F3F3F3", GRAY, 7.4)

    flow_arrow((0.50, 0.835), (0.50, 0.802))
    flow_arrow((0.50, 0.685), (0.50, 0.655))
    flow_arrow((0.790, 0.7425), (0.810, 0.7425))
    flow_arrow((0.50, 0.555), (0.50, 0.510))
    flow_arrow((0.50, 0.325), (0.50, 0.300))
    flow_arrow((0.50, 0.207), (0.50, 0.187))

    save_figure(fig, "Figure_1")


def compute_publication_trend(pool: pd.DataFrame) -> pd.DataFrame:
    direct_flag = pool["Direct isolated monomer experimental evidence"].eq("Yes")
    trend = (pool.assign(Evidence=np.where(direct_flag, "Direct purified-compound", "Contextual"))
             .groupby(["Publication year", "Evidence"]).size().unstack(fill_value=0)
             .reindex(range(2013, 2026), fill_value=0).reset_index())
    trend.to_csv(DATA / "figure_2a_publication_trend.csv", index=False)
    return trend


def compute_chemical_class(retained: pd.DataFrame, pool: pd.DataFrame) -> pd.DataFrame:
    first_pass = retained[retained["Evidence level"].str.match(r"^Level [123]")].copy()
    first_counts = first_pass["Class"].fillna("Not specified").value_counts()
    direct = pool[pool["Direct isolated monomer experimental evidence"].eq("Yes")].copy()
    direct_class = direct["Class(es)"].fillna("Not specified").str.split(";").str[0].str.strip().value_counts()
    rows = []
    for name in TABLE1_CLASS_ORDER[:-1]:
        rows.append({"Chemical class": name, "First-pass row-level annotations": int(first_counts.get(name, 0)), "Direct unique references": int(direct_class.get(name, 0))})
    named = set(TABLE1_CLASS_ORDER[:-1])
    rows.append({"Chemical class": "Other classes", "First-pass row-level annotations": int(first_counts[~first_counts.index.isin(named)].sum()), "Direct unique references": int(direct_class[~direct_class.index.isin(named)].sum())})
    df = pd.DataFrame(rows)
    df.to_csv(DATA / "figure_2b_chemical_class.csv", index=False)
    return df


def generate_figure_2() -> None:
    _, retained, pool = load_and_validate()
    trend = compute_publication_trend(pool)
    class_data = compute_chemical_class(retained, pool)
    fig, axes = plt.subplots(
        1, 2, figsize=(7.48, 3.55),
        gridspec_kw={"width_ratios": [1.05, 1.35]},
    )
    ax = axes[0]
    ax.plot(trend["Publication year"], trend["Direct purified-compound"], marker="o", markersize=3.3, linewidth=1.4, color=BLUE, label="Direct purified-compound")
    ax.plot(trend["Publication year"], trend["Contextual"], marker="s", markersize=3.3, linewidth=1.4, linestyle="--", color=ORANGE, label="Contextual")
    ax.set(xlabel="Publication year", ylabel="Unique PMID-linked references", title="Annual retained-reference distribution\n(audited corpus, n = 168)")
    ax.set_xticks([2013, 2015, 2017, 2019, 2021, 2023, 2025])
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, loc="upper left")
    panel_label(ax, "A", x=-0.18, y=1.06)
    sns.despine(ax=ax)

    ax = axes[1]
    plot = class_data.iloc[::-1]
    y = np.arange(len(plot))
    ax.barh(y - 0.18, plot["First-pass row-level annotations"], height=0.34, color=SKY, label="Row-level annotations (n = 141)")
    ax.barh(y + 0.18, plot["Direct unique references"], height=0.34, color=BLUE, label="Direct references (n = 77)")
    ax.set_yticks(y, plot["Chemical class"])
    ax.set(xlabel="Count within stated unit")
    ax.set_title("Chemical-class distribution", loc="left", pad=25)
    ax.set_xlim(left=0)
    ax.legend(
        frameon=False, loc="upper right", bbox_to_anchor=(1.01, 1.16),
        fontsize=5.9, handlelength=1.7, labelspacing=0.20,
    )
    panel_label(ax, "B", x=-0.10, y=1.06)
    sns.despine(ax=ax)
    fig.subplots_adjust(left=0.085, right=0.985, top=0.82, bottom=0.18, wspace=0.40)
    save_figure(fig, "Figure_2")


def compute_reference_distribution(pool: pd.DataFrame) -> pd.DataFrame:
    hierarchy = pool["Primary evidence classification"].value_counts().rename_axis("Reference-level classification").reset_index(name="Unique PMID-linked references")
    hierarchy["Reference-level classification"] = hierarchy["Reference-level classification"].replace({
        "Cell-only isolated monomer evidence": "Cell-only isolated-compound evidence",
        "Animal-validated isolated monomer evidence": "Animal-validated isolated-compound evidence",
    })
    hierarchy["Attribution group"] = np.where(hierarchy["Reference-level classification"].str.contains("Animal-validated|Cell-only", regex=True), "Direct purified-compound evidence", "Contextual evidence")
    hierarchy.to_csv(DATA / "figure_3a_reference_distribution.csv", index=False)
    return hierarchy


def compute_model_features(pool: pd.DataFrame) -> pd.DataFrame:
    direct = pool[pool["Direct isolated monomer experimental evidence"].eq("Yes")].copy()
    text = _combined_text(direct)
    rows = [{"Explicitly recorded model feature": label, "Direct unique references": int(text.map(lambda x, terms=terms: _contains(x, terms)).sum())} for label, terms in MODEL_RULES.items()]
    df = pd.DataFrame(rows)
    df.to_csv(DATA / "figure_3b_model_usage.csv", index=False)
    return df


def generate_figure_3() -> None:
    _, _, pool = load_and_validate()
    hierarchy = compute_reference_distribution(pool)
    model_data = compute_model_features(pool)

    fig, axes = plt.subplots(
        2, 1, figsize=(7.48, 6.40),
        gridspec_kw={"height_ratios": [1.65, 1.0]},
    )

    ax = axes[0]
    plot = hierarchy.sort_values("Unique PMID-linked references")
    colors = plot["Attribution group"].map({
        "Direct purified-compound evidence": BLUE,
        "Contextual evidence": ORANGE,
    })
    y = np.arange(len(plot))
    bars = ax.barh(y, plot["Unique PMID-linked references"], color=colors)
    reference_labels = {
        "Cell-only isolated-compound evidence": "Cell-only isolated-\ncompound evidence",
        "Animal-validated isolated-compound evidence": "Animal-validated isolated-\ncompound evidence",
        "Extract/formula/decoction contextual evidence": "Extract/formula/decoction\ncontextual evidence",
        "Formulation/nanodelivery contextual evidence": "Formulation/nanodelivery\ncontextual evidence",
        "Network pharmacology/in silico evidence": "Network pharmacology /\nin silico evidence",
        "Review/background evidence": "Review/background\nevidence",
        "Derivative/metabolite contextual evidence": "Derivative/metabolite\ncontextual evidence",
        "Multi-compound combination contextual evidence": "Multi-compound combination\ncontextual evidence",
        "Other CRC-relevant contextual evidence": "Other CRC-relevant\ncontextual evidence",
    }
    ax.set_yticks(y, [reference_labels.get(v, v) for v in plot["Reference-level classification"]])
    ax.bar_label(bars, padding=3, fontsize=7)
    ax.set_xlabel("Unique PMID-linked references")
    ax.set_title("Reference-level evidence distribution (n = 168)", loc="left", pad=10)
    ax.set_xlim(0, max(plot["Unique PMID-linked references"]) + 8)
    ax.legend(
        handles=[
            Line2D([0], [0], color=BLUE, lw=6, label="Direct purified-compound"),
            Line2D([0], [0], color=ORANGE, lw=6, label="Contextual"),
        ],
        frameon=False, loc="lower right", fontsize=6.4,
        handlelength=1.7, labelspacing=0.25,
    )
    panel_label(ax, "A", x=-0.43, y=1.02)
    sns.despine(ax=ax)

    ax = axes[1]
    plot = model_data.sort_values("Direct unique references")
    y = np.arange(len(plot))
    bars = ax.barh(y, plot["Direct unique references"], color=GREEN)
    model_labels = {
        "Syngeneic / immune-competent": "Syngeneic /\nimmune-competent",
        "AOM/DSS or inherited tumor model": "AOM/DSS or inherited\ntumor model",
        "Conventional xenograft": "Conventional xenograft",
        "Organ-specific / liver-metastasis": "Organ-specific /\nliver-metastasis",
        "Patient-derived CRC organoid": "Patient-derived\nCRC organoid",
    }
    ax.set_yticks(y, [model_labels.get(v, v) for v in plot["Explicitly recorded model feature"]])
    ax.bar_label(bars, padding=3, fontsize=7)
    ax.set_xlabel("Direct unique references (non-exclusive counts)")
    ax.set_title("Selected model features in the direct evidence set", loc="left", pad=10)
    ax.set_xlim(0, max(1, plot["Direct unique references"].max() + 2))
    panel_label(ax, "B", x=-0.43, y=1.02)
    sns.despine(ax=ax)

    fig.subplots_adjust(left=0.42, right=0.965, top=0.94, bottom=0.09, hspace=0.62)
    save_figure(fig, "Figure_3")

def compute_compound_counts(pool: pd.DataFrame) -> pd.DataFrame:
    direct = pool[pool["Direct isolated monomer experimental evidence"].eq("Yes")].copy()
    compound_col = direct["Monomer(s) mapped in source table"]
    rows = []
    for compound in COMPOUNDS_WITH_COMPARATOR:
        subset = direct[_compound_mask(compound_col, compound)]
        rows.append({
            "Compound": compound,
            "Plot label": f"{compound}\u2020" if compound == COMPARATOR else compound,
            "Animal-validated references": int(subset["Animal-validated monomer evidence"].eq("Yes").sum()),
            "Cell-only references": int(subset["Cell-only monomer evidence"].eq("Yes").sum()),
            "TCM-source status": "Related natural-product comparator" if compound == COMPARATOR else "Core problem-oriented candidate",
        })
    df = pd.DataFrame(rows)
    df.to_csv(DATA / "figure_4a_compound_evidence_levels.csv", index=False)
    return df


def compute_domain_matrix(pool: pd.DataFrame) -> pd.DataFrame:
    direct = pool[pool["Direct isolated monomer experimental evidence"].eq("Yes")].copy()
    compound_col = direct["Monomer(s) mapped in source table"]
    combined = _combined_text(direct)
    rows = []
    for compound in CORE_COMPOUNDS:
        subset = direct[_compound_mask(compound_col, compound)]
        subset_text = combined.loc[subset.index]
        for domain, terms in DOMAIN_RULES.items():
            rows.append({"Compound": compound, "CRC problem domain": domain.replace("\n", " "), "Direct unique references": int(subset_text.map(lambda x, terms=terms: _contains(x, terms)).sum())})
    df = pd.DataFrame(rows)
    df.to_csv(DATA / "figure_4b_compound_domain_heatmap.csv", index=False)
    return df


def generate_figure_4() -> None:
    _, _, pool = load_and_validate()
    compound_data = compute_compound_counts(pool)
    domain_data = compute_domain_matrix(pool)
    fig, axes = plt.subplots(
        2, 1, figsize=(7.48, 6.65),
        gridspec_kw={"height_ratios": [1.05, 1.0]},
    )
    ax = axes[0]
    plot = compound_data.iloc[::-1]
    ax.barh(plot["Plot label"], plot["Cell-only references"], color=SKY, label="Cell-only")
    ax.barh(plot["Plot label"], plot["Animal-validated references"], left=plot["Cell-only references"], color=GREEN, label="Animal-validated")
    ax.set(xlabel="Direct unique references", title="Evidence categories across representative compounds")
    ax.legend(
        frameon=False, loc="lower right", bbox_to_anchor=(0.99, 0.02),
        ncol=1, fontsize=6.2, handlelength=1.7, labelspacing=0.20,
    )
    panel_label(ax, "A", x=-0.08, y=1.05)
    sns.despine(ax=ax)

    ax = axes[1]
    matrix = domain_data.pivot(index="Compound", columns="CRC problem domain", values="Direct unique references").reindex(CORE_COMPOUNDS)
    matrix = matrix[[name.replace("\n", " ") for name in DOMAIN_RULES]]
    sns.heatmap(
        matrix, annot=True, fmt="d", cmap="Blues", linewidths=0.6, linecolor="white",
        cbar_kws={"label": "Direct references", "shrink": 0.88, "pad": 0.02}, ax=ax,
        annot_kws={"fontsize": 7},
    )
    ax.set(xlabel="CRC problem domain", ylabel="", title="Problem-compound distribution within the direct evidence set")
    ax.set_xticklabels([
        "Resistance /\nsensitization",
        "Metastatic /\nimmune niche",
        "Inflammation-associated\ntumorigenesis",
    ], rotation=0, ha="center", fontsize=7.0)
    ax.tick_params(axis="y", labelsize=7, pad=2)
    panel_label(ax, "B", x=-0.08, y=1.05)
    fig.subplots_adjust(left=0.20, right=0.94, top=0.93, bottom=0.10, hspace=0.58)
    save_figure(fig, "Figure_4")


def generate_treatment_context_figure() -> None:
    _, _, pool = load_and_validate()
    direct = pool[pool["Direct isolated monomer experimental evidence"].eq("Yes")].copy()
    compound_col = direct["Monomer(s) mapped in source table"]
    combined = _combined_text(direct)
    rows = []
    for compound in COMPOUNDS_WITH_COMPARATOR:
        subset = direct[_compound_mask(compound_col, compound)]
        subset_text = combined.loc[subset.index]
        for treatment, terms in TREATMENT_RULES.items():
            hit = subset_text.map(lambda x, terms=terms: _contains(x, terms))
            hit_subset = subset.loc[hit.index[hit]]
            total = len(hit_subset)
            animal = int(hit_subset["Animal-validated monomer evidence"].eq("Yes").sum())
            rows.append({
                "Compound": compound,
                "Treatment context": treatment.replace("\n", " "),
                "Direct unique references": total,
                "Animal-validated references": animal,
                "Animal-validation proportion": animal / total if total else np.nan,
            })
    treatment_data = pd.DataFrame(rows)
    treatment_data.to_csv(DATA / "figure_5_treatment_context_bubbles.csv", index=False)

    x_names = list(TREATMENT_RULES)
    y_names = CORE_COMPOUNDS + [f"{COMPARATOR}\u2020"]
    fig, ax = plt.subplots(figsize=(7.48, 4.15))
    for row in treatment_data.itertuples(index=False):
        y_label = f"{row.Compound}\u2020" if row.Compound == COMPARATOR else row.Compound
        x = x_names.index(next(name for name in x_names if name.replace("\n", " ") == row[1]))
        y = y_names.index(y_label)
        total = row[2]
        if total > 0:
            ax.scatter(x, y, s=80 + 90 * total, c=[row[4]], cmap="viridis", vmin=0, vmax=1, edgecolor="#333333", linewidth=0.6)
            ax.text(x, y, str(total), ha="center", va="center", fontsize=7)
    ax.set_xticks(range(len(x_names)), x_names)
    ax.set_yticks(range(len(y_names)), y_names)
    ax.set_xlim(-0.6, len(x_names) - 0.4); ax.set_ylim(len(y_names) - 0.5, -0.5)
    ax.set(xlabel="Treatment-defined context", ylabel="", title="Treatment-context evidence for representative compounds")
    sm = mpl.cm.ScalarMappable(cmap="viridis", norm=mpl.colors.Normalize(vmin=0, vmax=1)); sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02, fraction=0.035)
    cbar.set_label("Animal-validation proportion")
    sns.despine(ax=ax)
    fig.tight_layout()
    save_figure(fig, "Figure_5")


def generate_evidence_gap_figure() -> None:
    _, _, pool = load_and_validate()
    direct = pool[pool["Direct isolated monomer experimental evidence"].eq("Yes")].copy()
    text = _combined_text(direct)
    rows = []
    for feature, terms in POSITIVE_DESIGN_FEATURE_RULES.items():
        recorded = int(text.map(lambda x, terms=terms: _contains(x, terms)).sum())
        if recorded > 0:
            rows.append({"Design feature explicitly captured": feature, "Direct unique references": recorded})
    feature_data = pd.DataFrame(rows).sort_values("Direct unique references", ascending=True)
    feature_data.to_csv(DATA / "figure_6_selected_design_features.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.48, 3.65))
    bars = ax.barh(feature_data["Design feature explicitly captured"], feature_data["Direct unique references"], color=GREEN)
    ax.bar_label(bars, padding=2, fontsize=7)
    ax.set(xlabel="Direct unique references (positive records only)", title="Selected clinically relevant design features explicitly captured\nin the audited direct evidence set")
    ax.set_xlim(0, max(12, feature_data["Direct unique references"].max() + 3))
    sns.despine(ax=ax)
    fig.tight_layout()
    save_figure(fig, "Figure_6")


def generate_all() -> None:
    set_style()
    generate_evidence_curation_flow()
    generate_figure_2()
    generate_figure_3()
    generate_figure_4()
    generate_treatment_context_figure()
    generate_evidence_gap_figure()
    print(f"Figures written to: {FIGURES}")
    print(f"Cleaned figure data written to: {DATA}")


if __name__ == "__main__":
    generate_all()



