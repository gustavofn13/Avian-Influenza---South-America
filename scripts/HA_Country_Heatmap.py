#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HA subtype distribution across South American countries
South American AIV Dataset - reads from AP3 (Unified sheet)

Standard header shared with the timeline/host heatmap scripts.
Text in SVG is preserved as editable text (svg.fonttype = 'none').
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import numpy as np
from pathlib import Path

# ==========================================================
# STANDARD LOAD BLOCK (shared across all figure scripts)
# ==========================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_FILE = PROJECT_ROOT / "data" / "AP3_subtype_counts_oficial.xlsx"
OUTPUT_DIR = PROJECT_ROOT / "output"
SHEET = "Unified"
TARGET_COUNTRIES = ["Chile", "Argentina", "Brazil", "Peru", "Colombia", "Bolivia"]

OUTPUT_DIR.mkdir(exist_ok=True)
df = pd.read_excel(INPUT_FILE, sheet_name=SHEET)

df["Country"] = (
    df["Location"]
    .astype(str)
    .str.strip()
    .str.split("/").str[0]
    .str.strip()
)
df = df[df["Country"].isin(TARGET_COUNTRIES)].copy()

df["HA"] = (
    df["Subtype"]
    .astype(str)
    .str.upper()
    .str.extract(r"(H\d{1,2})")[0]
)
VALID_HA = [f"H{i}" for i in range(1, 17)]
df = df[df["HA"].isin(VALID_HA)].copy()

# ==========================================================
# EDITABLE-TEXT SVG SETTINGS
# ==========================================================

mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"

# ==========================================================
# BUILD MATRIX
# ==========================================================

country_order = (
    df["Country"]
    .value_counts()
    .reindex(TARGET_COUNTRIES)
    .sort_values(ascending=False)
    .index
    .tolist()
)

heatmap_data = pd.crosstab(df["Country"], df["HA"])
heatmap_data = heatmap_data.reindex(
    index=country_order, columns=VALID_HA, fill_value=0
)

print("Heatmap matrix (matches AP3):\n")
print(heatmap_data)

print("\nCountry totals (must match AP3):")
print(df["Country"].value_counts().reindex(TARGET_COUNTRIES))
print(f"Grand total: {len(df)}")

# ==========================================================
# PLOT
# ==========================================================

sns.set_theme(style="white", context="paper")
plt.rcParams.update({"font.size": 14})

annot_data = heatmap_data.astype(str).mask(heatmap_data == 0, "")
mask = heatmap_data.astype(float).replace(0, np.nan).isnull()
vmax = max(1, heatmap_data.to_numpy().max())

fig, ax = plt.subplots(figsize=(9, 5))
ax.set_facecolor("#D9D9D9")

sns.heatmap(
    heatmap_data,
    mask=mask,
    annot=annot_data,
    fmt="",
    cmap="YlOrRd",
    square=True,
    linewidths=0.4,
    linecolor="white",
    vmin=1, vmax=vmax,
    cbar_kws={"label": "Number of detections"},
    ax=ax,
)

ax.set_title(
    "Distribution of HA subtypes across South American countries",
    fontsize=13, fontweight="bold", pad=12,
)
ax.set_xlabel("HA subtype")
ax.set_ylabel("Country")
plt.xticks(rotation=0)
plt.yticks(rotation=0)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "HA_Country_Heatmap.svg", format="svg", bbox_inches="tight")
plt.savefig(OUTPUT_DIR / "HA_Country_Heatmap.png", dpi=600, bbox_inches="tight")
plt.close()

print("\nSaved: HA_Country_Heatmap.svg / .png")
