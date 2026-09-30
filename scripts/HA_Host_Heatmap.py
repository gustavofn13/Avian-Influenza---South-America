#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HA subtype distribution across ecological host groups
South American AIV Dataset - reads from AP3 (Unified sheet)

Standard header shared with the timeline/country heatmap scripts.
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
# STANDARD PLOT STYLE (shared across all heatmap/matrix scripts)
# --------------------------------------------
# Figure size scales with the matrix shape (CELL_SIZE inches per
# row/column), so every heatmap keeps the same cell proportions
# regardless of how many rows/columns it has.
# ==========================================================

CELL_SIZE = 0.55       # inches per heatmap cell
PAD_W = 3.0             # extra width for y-axis labels + colorbar
PAD_H = 2.0             # extra height for title + x-axis tick labels
FONT_SIZE_BASE = 10
FONT_SIZE_TITLE = 13
FONT_SIZE_ANNOT = 10

# ==========================================================
# HOST -> ECOLOGICAL GROUP MAPPING
# --------------------------------------------
# Categorization is done by keyword matching on the Host field
# (case-insensitive). Matched in order: Waterfowl, Gulls/Seabirds,
# Shorebirds, Terrestrial/Domestic. Everything else falls into
# "Environmental / Unspecified" (mainly environmental samples,
# generic "Avian" entries, etc.).
# ==========================================================

def categorize_host(host):
    if pd.isna(host):
        return "Environmental / Unspecified"
    s = str(host).lower()

    waterfowl = [
        "duck", "anas", "anseriformes", "anserini", "anatidae",
        "mallard", "teal", "pintail", "shoveler", "wigeon",
        "dabbling", "anser", "cygnus", "swan", "goose", "waterfowl",
        "dendrocygna", "cairina",
    ]
    gulls_seabirds = [
        "gull", "larus", "tern", "sterna", "skimmer", "rynchops",
        "laridae", "stercorariidae", "skua", "shearwater",
        "albatross", "petrel", "pygoscelis", "penguin", "leucophaeus",
    ]
    shorebirds = [
        "sandpiper", "calidris", "arenaria", "turnstone",
        "plover", "charadrius", "oystercatcher", "haematopus",
        "curlew", "godwit", "whimbrel", "numenius", "stilt",
        "himantopus", "snipe", "shorebird", "wader", "actitis",
        "catoptrophorus", "tringa",
    ]
    terrestrial_domestic = [
        "chicken", "gallus", "turkey", "meleagris", "domestic",
        "poultry", "pheasant", "partridge", "quail", "guineafowl",
        "numida", "fulica", "gallinula", "coot", "moorhen",
    ]

    for group, keywords in [
        ("Waterfowl", waterfowl),
        ("Gulls / Seabirds", gulls_seabirds),
        ("Shorebirds", shorebirds),
        ("Terrestrial / Domestic", terrestrial_domestic),
    ]:
        if any(k in s for k in keywords):
            return group

    return "Environmental / Unspecified"

df["EcoGroup"] = df["Host"].apply(categorize_host)

# ==========================================================
# BUILD MATRIX
# ==========================================================

ECO_ORDER = [
    "Waterfowl",
    "Gulls / Seabirds",
    "Shorebirds",
    "Terrestrial / Domestic",
    "Environmental / Unspecified",
]

heatmap_data = pd.crosstab(df["EcoGroup"], df["HA"])
heatmap_data = heatmap_data.reindex(
    index=ECO_ORDER, columns=VALID_HA, fill_value=0
)

print("Heatmap matrix (matches AP3):\n")
print(heatmap_data)

print("\nEco group totals:")
print(df["EcoGroup"].value_counts().reindex(ECO_ORDER))
print(f"Grand total: {len(df)}")

# ==========================================================
# PLOT
# ==========================================================

sns.set_theme(style="white", context="paper")
plt.rcParams.update({"font.size": FONT_SIZE_BASE})

annot_data = heatmap_data.astype(str).mask(heatmap_data == 0, "")
mask = heatmap_data.astype(float).replace(0, np.nan).isnull()
vmax = max(1, heatmap_data.to_numpy().max())

n_rows, n_cols = heatmap_data.shape
fig, ax = plt.subplots(figsize=(n_cols * CELL_SIZE + PAD_W, n_rows * CELL_SIZE + PAD_H))
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
    annot_kws={"size": FONT_SIZE_ANNOT, "weight": "bold"},
    ax=ax,
)

ax.set_title(
    "Distribution of HA subtypes across ecological host groups",
    fontsize=FONT_SIZE_TITLE, fontweight="bold", pad=12,
)
ax.set_xlabel("HA subtype")
ax.set_ylabel("Ecological group")
plt.xticks(rotation=0)
plt.yticks(rotation=0)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "HA_Host_Heatmap.svg", format="svg", bbox_inches="tight")
plt.savefig(OUTPUT_DIR / "HA_Host_Heatmap.png", dpi=600, bbox_inches="tight")
plt.close()

print("\nSaved: HA_Host_Heatmap.svg / .png")
