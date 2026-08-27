# Eco-Epidemiology and Evolution of AIV in South America (2000–2024)

Reproducibility package for the figures of the systematic-review manuscript:

**"Eco-Epidemiological and Evolutionary Landscape of Avian Influenza Viruses in South America: A Continental Synthesis (2000–2024)"**

This repository contains the curated dataset (Appendix 3) and the scripts that generate four of the five main figures directly from it. Running the four scripts on the provided AP3 file reproduces the published figures.


---

## Contents

```
aiv_south_america/
├── data/
│   └── AP3_subtype_counts_oficial.xlsx    # curated master dataset (source of truth)
├── scripts/
│   ├── timeline_bubble.py                 # temporal distribution figure
│   ├── HA_Country_Heatmap.py              # HA subtype × country heatmap
│   ├── HA_Host_Heatmap.py                 # HA subtype × ecological group heatmap
│   └── HA_NA_Matrix.py                    # HA × NA subtype matrix
├── output/                                # figures land here (SVG + PNG)
├── requirements.txt                       # exact package versions
├── run_all.sh                             # runs all four scripts in sequence
└── README.md                              # this file
```

The `map_R/` folder is separate — it holds the R pipeline for the sampling effort map, with its own AP6 dataset and README.

---

## What the scripts do

All four scripts share the same load block and read from **`data/AP3_subtype_counts_oficial.xlsx`**, sheet `Unified`. This ensures that any update to AP3 propagates identically to all figures.

Standard processing across scripts:
- `Location` column is stripped of leading/trailing whitespace, split on `/`, and only the first token (country) is retained
- Records are filtered to the six target countries: Chile, Argentina, Brazil, Peru, Colombia, Bolivia
- `Subtype` is parsed with regex `H\d{1,2}` for HA and `N\d{1,2}` for NA
- SVG output uses editable text (`svg.fonttype = "none"`) — text remains as `<text>` elements, not paths, so it can be edited in Affinity Designer, Inkscape or Illustrator

### Expected totals from AP3

| Country   | Records (HA) | Records (HA+NA) |
|-----------|--------------|-----------------|
| Chile     | 208          | 207             |
| Argentina | 64           | 64              |
| Brazil    | 28           | 18              |
| Peru      | 25           | 25              |
| Colombia  | 7            | 7               |
| Bolivia   | 1            | 1               |
| **Total** | **333**      | **322**         |

The HA × NA matrix uses only records where both subtypes are characterized (n=322). Records with only HA (e.g., "H3" without NA typing — 11 records, all from Brazil) are excluded from this figure because they cannot be placed on a HA × NA cell.

---

## Requirements

- Python 3.10 or newer
- The packages listed in `requirements.txt`

Install into a virtual environment (recommended):

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

If starting from scratch on WSL, see "First-time setup" below.

---

## How to run

**All four figures in one command:**

```bash
bash run_all.sh
```

**Or one at a time:**

```bash
python scripts/timeline_bubble.py
python scripts/HA_Country_Heatmap.py
python scripts/HA_Host_Heatmap.py
python scripts/HA_NA_Matrix.py
```

Each script prints a totals table to the terminal for verification (it should match the tables above) and saves the figure as SVG and PNG in the current directory.

---

## Modifying the figures

**To edit visuals (colors, labels, font)** — open the `.svg` in Affinity Designer / Inkscape / Illustrator. All text is editable.

**To update the underlying data** — edit `data/AP3_subtype_counts_oficial.xlsx` and re-run the scripts. The figures will regenerate consistently.

**To add new countries** — modify `TARGET_COUNTRIES` at the top of each script.

---

## First-time setup on WSL (Ubuntu)

If your WSL environment is not yet set up for Python data work, run this once:

```bash
# 1. Install Python and venv support
sudo apt update
sudo apt install -y python3 python3-venv python3-pip

# 2. Clone or copy this repository into your WSL home
cd ~
# (place the aiv_south_america/ folder here)
cd aiv_south_america

# 3. Create and activate the virtual environment
python3 -m venv .venv
source .venv/bin/activate

# 4. Install exact package versions
pip install -r requirements.txt

# 5. Run
bash run_all.sh
```

After the first setup, each new session only needs:

```bash
cd ~/aiv_south_america
source .venv/bin/activate
bash run_all.sh
```

---

## Data source

The `Unified` sheet of AP3 consolidates two independent sources:
- **Sequences** (aba `Sequences`): 308 records retrieved from GISAID EpiFlu (accession IDs preserved)
- **Manuscript** (aba `Manuscript`): 123 records extracted from the 73 studies included in the systematic review (PRISMA-registered, OSF: [protocol DOI])

After deduplication and consolidation, the `Unified` sheet contains 335 rows. Of these:
- **333** enter the country/host/timeline figures (1 record from Ecuador and 1 residual header row are excluded)
- **322** enter the HA × NA matrix (only records with both HA and NA characterized)

---

## Citation

If you reuse this data or code, please cite:

Fenner GO, [coauthors], de Araujo J. *Eco-Epidemiological and Evolutionary Landscape of Avian Influenza Viruses in South America: A Continental Synthesis (2000–2024)*. [Journal, year]. DOI: [to be added upon acceptance]

---

## Contact

Gustavo de Oliveira Fenner — Laboratório de Pesquisas em Vírus Emergentes (LPVE), Instituto de Ciências Biomédicas, Universidade de São Paulo. gustavo.fenner@usp.br
