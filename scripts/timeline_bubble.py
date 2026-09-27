#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Temporal distribution of AIV subtype records in South America (2000-2024)
South American AIV Dataset - reads from AP3 (Unified sheet)

Two-panel temporal figure (replaces the previous bubble timeline):
  (A) Inter-annual: subtype records per YEAR, stacked by country.
  (B) Intra-annual: subtype records per MONTH of collection, stacked by
      country, over austral-season context bands (summer/autumn/winter/spring),
      with per-season totals.

Standard header shared with the country/host heatmap scripts:
  - loads AP3 from a single INPUT_FILE
  - normalizes 'Country' (before first '/', stripped)
  - restricts to the 6 target countries
  - saves SVG with editable text (svg.fonttype = 'none')
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# ==========================================================
# STANDARD LOAD BLOCK (shared across all figure scripts)
# ==========================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_FILE = PROJECT_ROOT / "data" / "AP3_subtype_counts_oficial.xlsx"
OUTPUT_DIR = PROJECT_ROOT / "output"
SHEET = "Unified"
TARGET_COUNTRIES = ["Chile", "Argentina", "Brazil", "Peru", "Colombia", "Bolivia"]

df = pd.read_excel(INPUT_FILE, sheet_name=SHEET)

df["Country"] = (
    df["Location"].astype(str).str.strip().str.split("/").str[0].str.strip()
)
df = df[df["Country"].isin(TARGET_COUNTRIES)].copy()

df["Year"]  = pd.to_numeric(df["Year"], errors="coerce")
df["Month"] = pd.to_numeric(df["Month"], errors="coerce")

# ==========================================================
# EDITABLE-TEXT SVG SETTINGS
# ==========================================================

mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"

# ==========================================================
# STYLE  (pastel palette)
#   Chile      -> pale blue
#   Argentina  -> stronger blue
#   Brazil     -> pastel green
#   Peru       -> pastel red
#   Colombia   -> pastel yellow
#   Bolivia    -> pastel taupe
# ==========================================================

COUNTRY_COLORS = {
    "Chile":     "#B3D3E8",
    "Argentina": "#5B79B0",
    "Brazil":    "#9BCDB2",
    "Peru":      "#E39A9A",
    "Colombia":  "#ECD489",
    "Bolivia":   "#C7B39B",
}

YEARS  = list(range(2000, 2025))
MONTHS = list(range(1, 13))
MONTH_LABELS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

# Austral seasons (Southern Hemisphere): month groups + soft pastel tint
SEASON_MONTHS = {
    "Summer": [12, 1, 2],
    "Autumn": [3, 4, 5],
    "Winter": [6, 7, 8],
    "Spring": [9, 10, 11],
}
SEASON_TINT = {
    "Summer": "#FBE7DA",
    "Autumn": "#FBF0D6",
    "Winter": "#DEEBF3",
    "Spring": "#E4F0DE",
}
# contiguous shaded spans on the Jan..Dec axis (label, x0, x1, season)
SEASON_BANDS = [
    ("Summer", 0.5,  2.5),   # Jan-Feb
    ("Autumn", 2.5,  5.5),   # Mar-May
    ("Winter", 5.5,  8.5),   # Jun-Aug
    ("Spring", 8.5,  11.5),  # Sep-Nov
    ("Summer", 11.5, 12.5),  # Dec (summer tail)
]

# ==========================================================
# AGGREGATE
# ==========================================================

yr = df.dropna(subset=["Year"]).copy();  yr["Year"]  = yr["Year"].astype(int)
mo = df.dropna(subset=["Month"]).copy(); mo["Month"] = mo["Month"].astype(int)

A = (pd.crosstab(yr["Year"], yr["Country"])
       .reindex(index=YEARS, columns=TARGET_COUNTRIES, fill_value=0))
B = (pd.crosstab(mo["Month"], mo["Country"])
       .reindex(index=MONTHS, columns=TARGET_COUNTRIES, fill_value=0))

month_tot = B.sum(axis=1)
season_tot = {s: int(month_tot.reindex(m).sum()) for s, m in SEASON_MONTHS.items()}

print("Records per year (per country) — must match AP3:")
print(A.sum(axis=1).to_string())
print(f"\nTotal records with year:  {int(A.values.sum())} / {len(df)}")
print(f"Total records with month: {int(B.values.sum())} / {len(df)} "
      f"({B.values.sum()/len(df)*100:.0f}%)")
print("\nRecords per austral season:")
for s in ["Summer", "Autumn", "Winter", "Spring"]:
    print(f"  {s}: {season_tot[s]}")

# ==========================================================
# PLOT
# ==========================================================

fig, (axA, axB) = plt.subplots(
    2, 1, figsize=(17, 11.5),
    gridspec_kw={"height_ratios": [1.4, 1.0], "hspace": 0.30},
    facecolor="white",
)

# ---------- Panel A: per year, stacked by country ----------
bottom = np.zeros(len(YEARS))
for c in TARGET_COUNTRIES:
    axA.bar(YEARS, A[c].values, bottom=bottom, width=0.78,
            color=COUNTRY_COLORS[c], edgecolor="white", linewidth=0.4, label=c)
    bottom += A[c].values

totalsA = A.sum(axis=1).values
for x, t in zip(YEARS, totalsA):
    if t > 0:
        axA.text(x, t + 0.8, str(int(t)), ha="center", va="bottom",
                 fontsize=12, color="#333333")

gap_top = totalsA.max() * 0.45
axA.annotate("no records\n2020–2021", xy=(2020.5, 1.5), xytext=(2020.5, gap_top),
             ha="center", va="center", fontsize=12, color="#8A8A8A",
             arrowprops=dict(arrowstyle="-", color="#CCCCCC", lw=1.2))

axA.set_title("A. Subtype records per year", fontsize=17, fontweight="bold",
              loc="left", pad=12, color="#111111")
axA.set_ylabel("Records", fontsize=16, fontweight="bold")
axA.set_xlim(1999.4, 2024.6)
axA.set_xticks(YEARS)
axA.set_xticklabels([str(y) for y in YEARS], rotation=90, fontsize=12, color="#333333")
axA.set_ylim(0, totalsA.max() * 1.15)
axA.tick_params(axis="y", labelsize=13)
for s in ["top", "right"]:
    axA.spines[s].set_visible(False)
axA.spines["left"].set_color("#CCCCCC")
axA.spines["bottom"].set_color("#CCCCCC")
axA.legend(title="Country", ncol=3, frameon=False, fontsize=13,
           title_fontsize=14, loc="upper left", bbox_to_anchor=(0.005, 0.99),
           columnspacing=1.4, handletextpad=0.6)

# ---------- Panel B: per month, stacked by country, season bands + totals ----------
band_top = month_tot.max() * 1.32
for name, x0, x1 in SEASON_BANDS:
    axB.axvspan(x0, x1, color=SEASON_TINT[name], alpha=1.0, zorder=0)
    if x1 - x0 > 1:   # label only the main (multi-month) bands
        axB.text((x0 + x1) / 2, band_top * 0.985,
                 f"{name}\nn = {season_tot[name]}",
                 ha="center", va="top", fontsize=13, color="#6E6E6E",
                 style="italic", zorder=1, linespacing=1.25)

bottom = np.zeros(len(MONTHS))
for c in TARGET_COUNTRIES:
    axB.bar(MONTHS, B[c].values, bottom=bottom, width=0.72,
            color=COUNTRY_COLORS[c], edgecolor="white", linewidth=0.4,
            zorder=3, label=c)
    bottom += B[c].values

totalsB = month_tot.values
for x, t in zip(MONTHS, totalsB):
    if t > 0:
        axB.text(x, t + 0.6, str(int(t)), ha="center", va="bottom",
                 fontsize=12, color="#333333", zorder=4)

axB.set_title("B. Subtype records per month of collection  "
              f"(n = {int(totalsB.sum())} with month; austral seasons shaded)",
              fontsize=17, fontweight="bold", loc="left", pad=12, color="#111111")
axB.set_ylabel("Records", fontsize=16, fontweight="bold")
axB.set_xlabel("Month of collection", fontsize=17, fontweight="bold", labelpad=10)
axB.set_xlim(0.5, 12.5)
axB.set_xticks(MONTHS)
axB.set_xticklabels(MONTH_LABELS, fontsize=13, color="#333333")
axB.set_ylim(0, band_top)
axB.tick_params(axis="y", labelsize=13)
for s in ["top", "right"]:
    axB.spines[s].set_visible(False)
axB.spines["left"].set_color("#CCCCCC")
axB.spines["bottom"].set_color("#CCCCCC")

fig.suptitle("Temporal distribution of AIV subtype records in South America (2000–2024)",
             fontsize=20, fontweight="bold", x=0.012, ha="left", y=0.995)

plt.tight_layout()

plt.savefig(
    OUTPUT_DIR / "timeline_bubble.svg",
    format="svg",
    bbox_inches="tight"
)

plt.savefig(
    OUTPUT_DIR / "timeline_bubble.png",
    dpi=600,
    bbox_inches="tight"
)

plt.close()

print("\nSaved: timeline_bubble.svg / .png")