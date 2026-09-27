Run all · SH
#!/bin/bash
# Regenerates all four figures from the AP3 dataset.
# Run from the project root: bash run_all.sh
 
set -e   # stop on first error
 
echo "==> Generating temporal distribution figure..."
python scripts/timeline_bubble.py
 
echo ""
echo "==> Generating HA x country heatmap..."
python scripts/HA_Country_Heatmap.py
 
echo ""
echo "==> Generating HA x ecological group heatmap..."
python scripts/HA_Host_Heatmap.py
 
echo ""
echo "==> Generating HA x NA subtype matrix..."
python scripts/HA_NA_Matrix.py
 
echo ""
echo "==> All figures generated successfully."
echo "    Check timeline_bubble.svg, HA_Country_Heatmap.svg,"
echo "    HA_Host_Heatmap.svg, and HA_NA_Matrix.svg"
 