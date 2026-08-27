#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Temporal distribution of AIV subtypes in South America (2000-2024)
South American AIV Dataset - reads from AP3 (Unified sheet)

Standard header shared with the country/host heatmap scripts:
  - loads AP3 from a single INPUT_FILE
  - normalizes 'Country' (before first '/', stripped)
  - restricts to the 6 target countries
  - saves SVG with editable text (svg.fonttype = 'none')
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
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

df["Year"] = pd.to_numeric(df["Year"], errors="coerce")
df = df.dropna(subset=["Year"])
df["Year"] = df["Year"].astype(int)

# ==========================================================
# EDITABLE-TEXT SVG SETTINGS
# Keeps text as text in the SVG (not converted to paths),
# so it can be edited in Affinity/Inkscape/Illustrator.
# ==========================================================

mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["pdf.fonttype"] = 42        # TrueType, editable in PDF too
mpl.rcParams["font.family"] = "sans-serif"

# ==========================================================
# AGGREGATE
# ==========================================================

df_agg = (
    df.groupby(["Country", "Year"])
      .size()
      .reset_index(name="Quantidade")
)

print("Totals per country (must match AP3):")
print(df_agg.groupby("Country")["Quantidade"].sum().reindex(TARGET_COUNTRIES))
print()

df_agg["Country"] = pd.Categorical(
    df_agg["Country"], categories=TARGET_COUNTRIES, ordered=True
)

# ==========================================================
# PLOT
# ==========================================================

country_colors = {
    "Chile":     "#E32219",
    "Argentina": "#FF9326",
    "Brazil":    "#3F9492",
    "Peru":      "#C9227A",
    "Colombia":  "#4A148C",
    "Bolivia":   "#5D4037",
}

fig, ax = plt.subplots(figsize=(14, 7), facecolor="white")
ax.set_facecolor("white")

for i, _ in enumerate(TARGET_COUNTRIES):
    ax.axhline(y=i, color="#E0E0E0", linestyle="-", linewidth=1, zorder=1)

ax.scatter(
    df_agg["Year"],
    df_agg["Country"].cat.codes,
    s=df_agg["Quantidade"] * 40,
    c=df_agg["Country"].map(country_colors),
    alpha=0.9,
    zorder=2,
)

for _, row in df_agg.iterrows():
    y_pos = TARGET_COUNTRIES.index(row["Country"])
    ax.text(
        x=row["Year"], y=y_pos - 0.22,
        s=str(row["Quantidade"]),
        ha="center", va="center",
        fontsize=8.5, color="#333333",
    )

ax.set_title(
    "Temporal distribution of AIV subtypes in South America (2000-2024)",
    fontsize=14, fontweight="bold", loc="left", pad=15, color="#111111",
)
ax.set_xlabel("Year", fontsize=11, fontweight="bold", labelpad=10)
ax.set_xlim(2000.5, 2024.5)
ax.set_xticks(range(2000, 2025))
ax.tick_params(axis="x", labelsize=9.5, colors="#333333")

ax.set_yticks(range(len(TARGET_COUNTRIES)))
ax.set_yticklabels(TARGET_COUNTRIES, fontsize=11, fontweight="bold", color="#222222")
ax.invert_yaxis()

for spine in ["top", "left", "right"]:
    ax.spines[spine].set_visible(False)
ax.spines["bottom"].set_color("#CCCCCC")
ax.spines["bottom"].set_linewidth(1)
ax.yaxis.set_ticks_position("none")

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "timeline_bubble.svg", format="svg", bbox_inches="tight")
plt.savefig(OUTPUT_DIR / "timeline_bubble.png", dpi=600, bbox_inches="tight")
plt.close()

print("Saved: timeline_bubble.svg / .png")
