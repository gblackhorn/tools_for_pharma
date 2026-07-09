"""Small text helpers shared across workflows."""

from __future__ import annotations

import re

import pandas as pd


def clean_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value).replace("\n", " ")).strip()


def sanitize_filename(value: object) -> str:
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", clean_text(value)).strip("_")
    return name or "plot"
