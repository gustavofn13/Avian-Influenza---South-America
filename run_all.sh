#!/bin/bash
# Regenerates all four figures from the AP3 dataset.
# Run from the project root: bash run_all.sh
 
set -e   # stop on first error

# Ubuntu/WSL usually only ships "python3"; other setups may only have "python".
if command -v python3 >/dev/null 2>&1; then
    PYTHON=python3
else
    PYTHON=python
fi

echo "==> Generating temporal distribution figure..."
$PYTHON scripts/timeline_bubble.py

echo ""
echo "==> Generating HA x country heatmap..."
$PYTHON scripts/HA_Country_Heatmap.py

echo ""
echo "==> Generating HA x ecological group heatmap..."
$PYTHON scripts/HA_Host_Heatmap.py

echo ""
echo "==> Generating HA x NA subtype matrix..."
$PYTHON scripts/HA_NA_Matrix.py

echo ""
echo "==> Building phylogroup tables (Table A / Table B)..."
(cd Phylogroup_github/phylogroup_counts && $PYTHON build_phylogroup_tables.py \
    --ap3 ../../data/AP3_subtype_counts_oficial.xlsx --out .)

echo ""
echo "==> Generating tanglegram figure..."
cp Phylogroup_github/phylogroup_counts/Table_B_new.csv Phylogroup_github/phylogroup_tanglegram/Table_B_new.csv
(cd Phylogroup_github/phylogroup_tanglegram && $PYTHON make_tanglegram_figure.py)

echo ""
echo "==> All figures generated successfully."
echo "    Check timeline_bubble.svg, HA_Country_Heatmap.svg,"
echo "    HA_Host_Heatmap.svg, HA_NA_Matrix.svg, and"
echo "    Phylogroup_github/phylogroup_tanglegram/tanglegram_figure.pdf"
 