#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HA x NA subtype matrix
South American AIV Dataset - reads from AP3 (Unified sheet)

Standard header shared with the timeline / country heatmap / host heatmap scripts.
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

# ==========================================================
# EDITABLE-TEXT SVG SETTINGS
# ==========================================================

mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"

# ==========================================================
# EXTRACT HA AND NA
# --------------------------------------------
# Keep only records where BOTH HA and NA are characterized.
# Records with only HA (e.g., "H3" without NA typing) are
# excluded from this matrix because they cannot be placed
# on a HA x NA cell. This is documented in the figure legend.
# ==========================================================

df["HA"] = df["Subtype"].astype(str).str.upper().str.extract(r"(H\d{1,2})")[0]
df["NA"] = df["Subtype"].astype(str).str.upper().str.extract(r"(N\d{1,2})")[0]

df = df.dropna(subset=["HA", "NA"]).copy()

HA_ORDER = [f"H{i}" for i in range(1, 17)]
NA_ORDER = [f"N{i}" for i in range(1, 10)]

df = df[df["HA"].isin(HA_ORDER) & df["NA"].isin(NA_ORDER)].copy()

# ==========================================================
# BUILD MATRIX
# ==========================================================

heatmap_data = pd.crosstab(df["HA"], df["NA"])
heatmap_data = heatmap_data.reindex(index=HA_ORDER, columns=NA_ORDER, fill_value=0)

# Add row and column totals
heatmap_data["Total"] = heatmap_data.sum(axis=1)
totals_row = heatmap_data.sum(axis=0)
totals_row.name = "Total"
heatmap_data = pd.concat([heatmap_data, totals_row.to_frame().T])

print("HA x NA matrix (matches AP3, HA+NA both characterized):\n")
print(heatmap_data)

print(f"\nGrand total (HA+NA complete pairs): {len(df)}")
print(f"HA subtypes present: {df['HA'].nunique()}")
print(f"NA subtypes present: {df['NA'].nunique()}")

# ==========================================================
# PLOT
# ==========================================================

sns.set_theme(style="white", context="paper")
plt.rcParams.update({"font.size": 14})

# Mask cells with 0 (grey background), keep Total row/column visible
mask = heatmap_data.astype(float).replace(0, np.nan).isnull()

# Annotation: show number, empty string for 0
annot_data = heatmap_data.astype(str).mask(heatmap_data == 0, "")

# Use max of non-Total cells for color scale, so Totals don't wash out the palette
inner_max = heatmap_data.iloc[:-1, :-1].to_numpy().max()
vmax = max(1, inner_max)

fig, ax = plt.subplots(figsize=(10, 10))
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
    "HA × NA subtype matrix in South American avian influenza viruses",
    fontsize=13, fontweight="bold", pad=12,
)
ax.set_xlabel("NA subtype")
ax.set_ylabel("HA subtype")
plt.xticks(rotation=0)
plt.yticks(rotation=0)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "HA_NA_Matrix.svg", format="svg", bbox_inches="tight")
plt.savefig(OUTPUT_DIR / "HA_NA_Matrix.png", dpi=600, bbox_inches="tight")
plt.close()

print("\nSaved: HA_NA_Matrix.svg / .png")
