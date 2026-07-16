"""Shared qPCR table constants and small text helpers."""

from __future__ import annotations

import re

import pandas as pd


SAMPLE_SIZE_COLUMN = "Sample size (n)"
REFERENCE_SOURCE_COLUMN = "Reference source"
MEAN_RQ_COLUMN = "MEAN RQ"
SEM_COLUMN = "SEM"
SAMPLE_ID_COLUMN = "Sample ID"
ANIMAL_ID_COLUMN = "Animal ID"
INDIVIDUAL_RQ_COLUMN = "Individual RQ"
CONTROL_COMPOUNDS = {"NRT", "NTC", "QC"}
PLOTDATA_SHEET_PREFIX = "plotdata-"
REFQC_SHEET_PREFIX = "refqc-"
REFERENCE_GENE_COLUMN = "Reference gene"
MEAN_CT_COLUMN = "Mean CT"


def clean_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value).replace("\n", " ")).strip()


def is_blank(value: object) -> bool:
    return value is None or pd.isna(value) or clean_text(value) == ""


def is_plotdata_sheet_name(sheet_name: object) -> bool:
    return clean_text(sheet_name).lower().startswith(PLOTDATA_SHEET_PREFIX)


def is_refqc_sheet_name(sheet_name: object) -> bool:
    return clean_text(sheet_name).lower().startswith(REFQC_SHEET_PREFIX)


def summary_columns() -> list[str]:
    return [
        "Group",
        "Compound ID",
        REFERENCE_SOURCE_COLUMN,
        SAMPLE_SIZE_COLUMN,
        MEAN_RQ_COLUMN,
        SEM_COLUMN,
    ]


def result_columns() -> list[str]:
    return summary_columns() + [
        SAMPLE_ID_COLUMN,
        ANIMAL_ID_COLUMN,
        INDIVIDUAL_RQ_COLUMN,
    ]
