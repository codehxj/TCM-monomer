from pathlib import Path

import pandas as pd


PAPER = Path(__file__).resolve().parent
INPUT = PAPER / "supplementary_evidence_table.csv"

# Manual reference-level corrections verified against locally stored article text.
DERIVATIVE_CONTEXT = {
    "32592323": "The tested antitumour agent was a beta-carboline-based quinone derivative, not isolated beta-lapachone.",
    "33153029": "The tested agents were oleoyl hybrids rather than isolated pinocembrin.",
    "37595397": "The tested agents were newly synthesized polyphenols inspired by natural scaffolds, not isolated honokiol.",
    "37982821": "The tested agent was ethylcoprostanol, a beta-sitosterol metabolite, not isolated beta-sitosterol.",
    "38759254": "The tested agents were honokiol thioether derivatives rather than isolated honokiol.",
    "40968164": "The isolated test compound mapped to kaempferol was kaempferol 3-O-rutinoside, not kaempferol aglycone.",
}
PREPARATION_CONTEXT = {
    "36102033": "The tested CRC-relevant material was cocoa pod husk-derived material rather than isolated catechin.",
    "40940406": "The tested material was an Aspergillus niger extract containing pachymic acid, not isolated pachymic acid.",
}
COMBINATION_CONTEXT = {
    "35966396": "The intervention was curcumin in combination with silibinin; it supports combination context only.",
}
ANIMAL_VALIDATED = {
    "35593352",
    "35990576",
    "37326338",
    "37749082",
    "39987121",
    "40184790",
    "40318528",
    "40242441",
}
CELL_VALIDATED = {
    "36355520": "The isolated emodin and chrysophanol cell-treatment results are direct; docking and ADME/Tox components remain hypothesis-generating.",
    "37939611": "Isovitexin was separately administered to CRC cells; the Herba Patriniae intervention remains contextual.",
    "39271360": "Luteolin was directly tested in CRC cells; network-derived target interpretation remains hypothesis-generating.",
    "40255473": "Tanshinone IIA was directly tested in HCT116 and SW480 cells with pathway-agonist reversal; network/docking components remain hypothesis-generating.",
}


def set_context(row: pd.Series, category: str, note: str) -> pd.Series:
    row["Primary evidence classification"] = category
    row["Direct isolated monomer experimental evidence"] = "No"
    row["Animal-validated monomer evidence"] = "No"
    row["Cell-only monomer evidence"] = "No"
    row["Use in review"] = "Context only; " + note
    return row


def correct(row: pd.Series) -> pd.Series:
    pmid = str(row["PMID"])
    if pmid in DERIVATIVE_CONTEXT:
        return set_context(row, "Derivative/metabolite contextual evidence", DERIVATIVE_CONTEXT[pmid])
    if pmid in PREPARATION_CONTEXT:
        row["Extract/formula/decoction contextual evidence"] = "Yes"
        return set_context(row, "Extract/formula/decoction contextual evidence", PREPARATION_CONTEXT[pmid])
    if pmid in COMBINATION_CONTEXT:
        return set_context(row, "Multi-compound combination contextual evidence", COMBINATION_CONTEXT[pmid])
    if pmid in ANIMAL_VALIDATED:
        row["Primary evidence classification"] = "Animal-validated isolated monomer evidence"
        row["Direct isolated monomer experimental evidence"] = "Yes"
        row["Animal-validated monomer evidence"] = "Yes"
        row["Cell-only monomer evidence"] = "No"
        row["Use in review"] = "Direct isolated-monomer experiment with animal evidence verified from the local article record."
        row["Review/background evidence"] = "No"
    if pmid in CELL_VALIDATED:
        row["Primary evidence classification"] = "Cell-only isolated monomer evidence"
        row["Direct isolated monomer experimental evidence"] = "Yes"
        row["Animal-validated monomer evidence"] = "No"
        row["Cell-only monomer evidence"] = "Yes"
        row["Use in review"] = "Direct experimental component only; " + CELL_VALIDATED[pmid]
        if pmid == "36355520":
            row["Monomer(s) mapped in source table"] = "Chrysophanol; Emodin"
        if pmid == "39271360":
            row["Models recorded"] = "CRC cells assessed by CCK-8 and RT-qPCR (verified in local article record)"
    if pmid == "38157250":
        row["Monomer(s) mapped in source table"] = "Apigenin; Luteolin"
        row["Use in review"] = "Direct cell evidence only; luteolin was experimentally assessed as an apigenin metabolite in the local article record."
    return row


def main() -> None:
    table = pd.read_csv(INPUT, dtype={"PMID": str})
    table = table.apply(correct, axis=1)
    summary = table["Primary evidence classification"].value_counts().rename_axis("Classification").reset_index(name="References")
    table.to_csv(PAPER / "supplementary_evidence_table.csv", index=False, encoding="utf-8-sig")
    with pd.ExcelWriter(PAPER / "supplementary_evidence_table.xlsx", engine="openpyxl") as writer:
        table.to_excel(writer, sheet_name="Table S1 evidence map", index=False)
        summary.to_excel(writer, sheet_name="Classification summary", index=False)
    direct = table["Direct isolated monomer experimental evidence"].eq("Yes").sum()
    animal = table["Animal-validated monomer evidence"].eq("Yes").sum()
    cell = table["Cell-only monomer evidence"].eq("Yes").sum()
    print(f"retained={len(table)} direct={direct} animal={animal} cell_only={cell}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
