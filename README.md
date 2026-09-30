# Eco-Epidemiology and Evolution of AIV in South America (2000–2024)

Reproducibility package for the figures of the systematic-review manuscript:

**"Eco-Epidemiological and Evolutionary Landscape of Avian Influenza Viruses in South America (2000–2024)"**

This repository contains the curated dataset (Appendix 3), the scripts that generate the figures directly from it, the BEAST XML files used for the molecular clock / BSSVS analyses, and the phylogroup classification pipeline. Running `run_all.sh` on the provided AP3 file reproduces the published figures and tables.


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
├── xmls/                                  # BEAST XML files, one folder per gene (PB2/PB1/PA/NP/MP)
│   └── <GENE>/
│       ├── <GENE>_MOLCLOCK_SKYRIDE.xml    # molecular clock, Bayesian Skyride
│       ├── <GENE>_MOLCLOCK_COA.xml        # molecular clock, constant coalescent
│       ├── <GENE>_BSSVS_NORTH_HOST.xml    # BSSVS, North clade, host predictor
│       ├── <GENE>_BSSVS_NORTH_LOCAL.xml   # BSSVS, North clade, location predictor
│       ├── <GENE>_BSSVS_SOUTH_HOST.xml    # BSSVS, South clade, host predictor
│       └── <GENE>_BSSVS_SOUTH_LOCAL.xml   # BSSVS, South clade, location (see script docstring)
├── Phylogroup_github/                     # phylogroup classification pipeline
│   ├── phylogroup_counts/                 # per-gene N/S clade assignment lists + table builder
│   │   └── build_phylogroup_tables.py
│   └── phylogroup_tanglegram/             # NP-backbone tanglegram figure
│       └── make_tanglegram_figure.py
├── requirements.txt                       # exact package versions
├── run_all.sh                             # runs the full pipeline in sequence
└── README.md                              # this file
```


---

## What the scripts do

All four scripts share the same load block and read from **`data/AP3_subtype_counts_oficial.xlsx`**, sheet `Unified`. This ensures that any update to AP3 propagates identically to all figures.

Standard processing across scripts:
- `Location` column is stripped of leading/trailing whitespace, split on `/`, and only the first token (country) is retained
- Records are filtered to the six target countries: Chile, Argentina, Brazil, Peru, Colombia, Bolivia
- `Subtype` is parsed with regex `H\d{1,2}` for HA and `N\d{1,2}` for NA
- SVG output uses editable text (`svg.fonttype = "none"`) — text remains as `<text>` elements, not paths, so it can be edited in Affinity Designer, Inkscape or Illustrator

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

**Everything in one command** (four AP3 figures + phylogroup tables + tanglegram figure):

```bash
bash run_all.sh
```

**Or one step at a time:**

```bash
python3 scripts/timeline_bubble.py
python3 scripts/HA_Country_Heatmap.py
python3 scripts/HA_Host_Heatmap.py
python3 scripts/HA_NA_Matrix.py

# phylogroup pipeline
cd Phylogroup_github/phylogroup_counts
python3 build_phylogroup_tables.py --ap3 ../../data/AP3_subtype_counts_oficial.xlsx --out .
cp Table_B_new.csv ../phylogroup_tanglegram/Table_B_new.csv
cd ../phylogroup_tanglegram
python3 make_tanglegram_figure.py
```

(On some systems the interpreter is called `python` instead of `python3` — use whichever one exists in your environment. `run_all.sh` auto-detects this.)

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
- **Manuscript** (aba `Manuscript`): 123 records extracted from the 47 studies included in the systematic review

---

## Citation

If you reuse this data or code, please cite:

Fenner GO, [coauthors], de Araujo J. *Eco-Epidemiological and Evolutionary Landscape of Avian Influenza Viruses in South America [Journal, year]. DOI: [to be added upon acceptance]

---

## Contact

Gustavo de Oliveira Fenner — Laboratório de Pesquisas em Vírus Emergentes (LPVE), Instituto de Ciências Biomédicas, Universidade de São Paulo. gustavo_fenner@icb.usp.br
