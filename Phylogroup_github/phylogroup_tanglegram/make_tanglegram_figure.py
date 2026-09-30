#!/usr/bin/env python3
"""
Tanglegram exploration figure (A4, single page)
================================================
Top: 4 small tanglegrams, NP vs each of the other genes (PB2, PB1, PA, MP).
     Branches/links colored by phylogroup (per-gene S/N assignment).
Bottom: two compatibility heatmaps — host x phylogroup and country x
        phylogroup, percentages computed over the full 257-isolate set,
        phylogroup taken from the NP (backbone) column.

Inputs : treefiles/{1_PB2,2_PB1,3_PA,5_NP,7_MP}_INTERSECT.treefile
         Table_B_new.csv
Output : tanglegram_figure.pdf / .png  (A4 portrait)
"""

import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.collections import LineCollection
from Bio import Phylo

# ══════════════════════════════════════════════════════════════════════════════
# CONFIG
# ══════════════════════════════════════════════════════════════════════════════

TABLE_B_FILE = "Table_B_new.csv"

TREE_FILES = {
    "NP":  "treefiles/5_NP_INTERSECT.treefile",
    "PB2": "treefiles/1_PB2_INTERSECT.treefile",
    "PB1": "treefiles/2_PB1_INTERSECT.treefile",
    "PA":  "treefiles/3_PA_INTERSECT.treefile",
    "MP":  "treefiles/7_MP_INTERSECT.treefile",
}
PAIRS = [("NP", "PB2"), ("NP", "PB1"), ("NP", "PA"), ("NP", "MP")]

OUTPUT_PREFIX = "tanglegram_figure"
DPI = 300

# Phylogroup colors — SWAPPED relative to the concordance figure:
# North = green, South = red/maroon
COL_NORTH = "#4DB27C"
COL_SOUTH = "#AC1A40"
COL_MISS  = "#bdbdbd"
PHYLO_COLOR = {"N": COL_NORTH, "S": COL_SOUTH, "": COL_MISS}

FONT = "DejaVu Serif"
FS_TITLE   = 15
FS_SUBTITLE = 10
FS_TABLE_LABEL = 8.5
FS_TABLE_CELL  = 8
FS_LEGEND = 9

LW_TREE = 0.35
LW_LINK = 0.4
ALPHA_LINK = 0.55

HOST_ORDER = [
    "Waterfowl", "Shorebirds", "Gulls/Seabirds",
    "Terrestrial/Domestic", "Environmental/Other",
]
COUNTRY_ORDER = ["Chile", "Argentina", "Peru", "Brazil", "Colombia", "Bolivia", "Unknown"]

# ══════════════════════════════════════════════════════════════════════════════
# helpers
# ══════════════════════════════════════════════════════════════════════════════

def extract_epi(label):
    m = re.search(r'EPI_ISL_(\d+)', label)
    return m.group(1) if m else None


def normalize_name(label):
    return label.replace('_', ' ').replace('\xa0', ' ').strip().lower()


def canonical_id(tip_label):
    epi = extract_epi(tip_label)
    if epi:
        return "EPI_" + epi
    name_field = tip_label.split('|')[1] if '|' in tip_label else tip_label
    return "NAME_" + normalize_name(name_field)


def get_node_coords(clade, x_start=0.0, counter=None):
    if counter is None:
        counter = [0]
    coords = {}
    bl = clade.branch_length or 0.0
    x = x_start + bl
    if clade.is_terminal():
        y = counter[0]
        counter[0] += 1
        coords[clade] = (x, y)
    else:
        child_ys = []
        for child in clade.clades:
            coords.update(get_node_coords(child, x, counter))
            child_ys.append(coords[child][1])
        coords[clade] = (x, np.mean(child_ys))
    return coords


def build_tree_lines(tree, coords, tip_color, default_color):
    h_int, h_term, h_term_colors, v_segs = [], [], [], []

    def walk(clade):
        if clade.is_terminal():
            return
        px, py = coords[clade]
        child_ys = []
        for child in clade.clades:
            cx, cy = coords[child]
            if child.is_terminal():
                h_term.append([(px, cy), (cx, cy)])
                h_term_colors.append(tip_color.get(child, default_color))
            else:
                h_int.append([(px, cy), (cx, cy)])
            child_ys.append(cy)
            walk(child)
        v_segs.append([(px, min(child_ys)), (px, max(child_ys))])

    root = tree.root
    rx, ry = coords[root]
    h_int.append([(0, ry), (rx, ry)])
    walk(root)
    return h_int, h_term, h_term_colors, v_segs


def load_tree(path):
    tree = Phylo.read(path, "newick")
    tree.root_at_midpoint()
    coords = get_node_coords(tree.root)
    tips = tree.get_terminals()
    max_x = max(coords[t][0] for t in tips)
    return tree, coords, tips, max_x


# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Load Table B
# ══════════════════════════════════════════════════════════════════════════════

print("Loading Table B …")
tb = pd.read_csv(TABLE_B_FILE, dtype=str).fillna('')
tb['_id'] = tb.apply(
    lambda r: "EPI_" + r['EPI_ISL'].strip() if r['EPI_ISL'].strip()
    else "NAME_" + normalize_name(r['isolate']),
    axis=1,
)
row_by_id = {r['_id']: r for _, r in tb.iterrows()}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Load all 5 trees
# ══════════════════════════════════════════════════════════════════════════════

trees = {}
for gene, path in TREE_FILES.items():
    print(f"Loading {gene} tree …")
    trees[gene] = load_tree(path)
    print(f"   {len(trees[gene][2])} tips")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Draw figure (A4 portrait)
# ══════════════════════════════════════════════════════════════════════════════

print("Drawing figure …")

fig = plt.figure(figsize=(8.27, 11.69))
fig.patch.set_facecolor('white')

gs_title  = gridspec.GridSpec(1, 1, figure=fig, left=0.05, right=0.95, top=0.99, bottom=0.905)
gs_tangle = gridspec.GridSpec(1, 1, figure=fig, left=0.03, right=0.99, top=0.885, bottom=0.335)
gs_tables = gridspec.GridSpec(1, 1, figure=fig, left=0.17, right=0.97, top=0.295, bottom=0.065)

# ── Title ────────────────────────────────────────────────────────────────────
ax_title = fig.add_subplot(gs_title[0])
ax_title.axis('off')
ax_title.text(0.5, 0.7, "Tanglegram: NP backbone vs. other internal genes",
              ha='center', va='center', fontsize=FS_TITLE, fontweight='bold',
              fontname=FONT, transform=ax_title.transAxes)
ax_title.text(0.5, 0.05,
              "Branches colored by per-gene phylogroup assignment  "
              "(North = green, South = maroon)",
              ha='center', va='center', fontsize=FS_SUBTITLE, fontname=FONT,
              transform=ax_title.transAxes, color='#444444')

# ── Tanglegrams (2x2) ────────────────────────────────────────────────────────
tangle_gs = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs_tangle[0],
                                              wspace=0.25, hspace=0.35)

GAP = 0.25

for k, (left_gene, right_gene) in enumerate(PAIRS):
    ax = fig.add_subplot(tangle_gs[k // 2, k % 2])

    l_tree, l_coords, l_tips, l_max_x = trees[left_gene]
    r_tree, r_coords, r_tips, r_max_x = trees[right_gene]

    # id -> tip for each side
    l_id_tip = {canonical_id(t.name): t for t in l_tips}
    r_id_tip = {canonical_id(t.name): t for t in r_tips}

    def state_color(gene, tip_id):
        row = row_by_id.get(tip_id)
        val = row[gene] if row is not None else ''
        return PHYLO_COLOR.get(val, COL_MISS)

    l_tip_color = {t: state_color(left_gene, canonical_id(t.name)) for t in l_tips}
    r_tip_color = {t: state_color(right_gene, canonical_id(t.name)) for t in r_tips}

    l_h_int, l_h_term, l_h_term_col, l_v = build_tree_lines(l_tree, l_coords, l_tip_color, '#999999')
    r_h_int, r_h_term, r_h_term_col, r_v = build_tree_lines(r_tree, r_coords, r_tip_color, '#999999')

    def t_left(x):
        return (x / l_max_x if l_max_x > 0 else x) - 1.0

    def t_right(x):
        return GAP + (1.0 - (x / r_max_x if r_max_x > 0 else x))

    l_int_n  = [[(t_left(x1), y1), (t_left(x2), y2)] for (x1, y1), (x2, y2) in l_h_int]
    l_term_n = [[(t_left(x1), y1), (t_left(x2), y2)] for (x1, y1), (x2, y2) in l_h_term]
    l_v_n    = [[(t_left(x1), y1), (t_left(x2), y2)] for (x1, y1), (x2, y2) in l_v]

    r_int_n  = [[(t_right(x1), y1), (t_right(x2), y2)] for (x1, y1), (x2, y2) in r_h_int]
    r_term_n = [[(t_right(x1), y1), (t_right(x2), y2)] for (x1, y1), (x2, y2) in r_h_term]
    r_v_n    = [[(t_right(x1), y1), (t_right(x2), y2)] for (x1, y1), (x2, y2) in r_v]

    ax.add_collection(LineCollection(l_int_n, colors='#999999', linewidths=LW_TREE, zorder=2))
    ax.add_collection(LineCollection(l_v_n,   colors='#999999', linewidths=LW_TREE, zorder=2))
    ax.add_collection(LineCollection(l_term_n, colors=l_h_term_col, linewidths=LW_TREE, zorder=3))

    ax.add_collection(LineCollection(r_int_n, colors='#999999', linewidths=LW_TREE, zorder=2))
    ax.add_collection(LineCollection(r_v_n,   colors='#999999', linewidths=LW_TREE, zorder=2))
    ax.add_collection(LineCollection(r_term_n, colors=r_h_term_col, linewidths=LW_TREE, zorder=3))

    # connecting links, colored by the LEFT (NP) gene's phylogroup
    common_ids = set(l_id_tip) & set(r_id_tip)
    link_segs, link_colors = [], []
    for tid in common_ids:
        lt = l_id_tip[tid]
        rt = r_id_tip[tid]
        ly = l_coords[lt][1]
        ry = r_coords[rt][1]
        link_segs.append([(0.0, ly), (GAP, ry)])
        link_colors.append(state_color(left_gene, tid))
    ax.add_collection(LineCollection(link_segs, colors=link_colors,
                                      linewidths=LW_LINK, alpha=ALPHA_LINK, zorder=1))

    n_l = len(l_tips)
    n_r = len(r_tips)
    ax.set_xlim(-1.05, GAP + 1.05)
    ax.set_ylim(-3, max(n_l, n_r) + 2)
    ax.invert_yaxis()
    ax.axis('off')
    ax.set_title(f"{left_gene} vs {right_gene}", fontname=FONT,
                 fontsize=FS_SUBTITLE + 1, fontweight='bold', pad=3)

# ── Compatibility tables ─────────────────────────────────────────────────────
table_gs = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_tables[0], wspace=0.6)

total_n = len(tb)


def compat_matrix(category_col, order, fill_value="Unknown"):
    cats = tb[category_col].replace('', fill_value)
    counts = np.zeros((len(order), 2))  # columns: South, North
    for i, cat in enumerate(order):
        counts[i, 0] = ((cats == cat) & (tb['NP'] == 'S')).sum()
        counts[i, 1] = ((cats == cat) & (tb['NP'] == 'N')).sum()
    pct = counts / total_n * 100
    return counts, pct


def draw_compat_table(ax, title, order, counts, pct):
    ax.imshow(pct, aspect='auto', cmap='YlOrRd', vmin=0, vmax=pct.max() * 1.1,
              interpolation='nearest')
    for i in range(len(order)):
        for j in range(2):
            txt_color = 'white' if pct[i, j] > pct.max() * 0.55 else 'black'
            ax.text(j, i, f"{int(counts[i, j])}\n({pct[i, j]:.1f}%)",
                    ha='center', va='center', fontsize=FS_TABLE_CELL,
                    color=txt_color, fontname=FONT)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['South', 'North'], fontsize=FS_TABLE_LABEL, fontweight='bold', fontname=FONT)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=FS_TABLE_LABEL, fontname=FONT)
    ax.xaxis.set_ticks_position('top')
    ax.set_title(title, fontsize=FS_SUBTITLE + 1, fontweight='bold', fontname=FONT, pad=28)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xticks(np.arange(-0.5, 2, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(order), 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=1.5)
    ax.tick_params(which='minor', length=0)


ax_host = fig.add_subplot(table_gs[0])
host_counts, host_pct = compat_matrix('host_group', HOST_ORDER)
draw_compat_table(ax_host, f"Host group vs. phylogroup (NP)\n% of {total_n} isolates",
                   HOST_ORDER, host_counts, host_pct)

ax_ctry = fig.add_subplot(table_gs[1])
ctry_counts, ctry_pct = compat_matrix('country', COUNTRY_ORDER)
draw_compat_table(ax_ctry, f"Country vs. phylogroup (NP)\n% of {total_n} isolates",
                   COUNTRY_ORDER, ctry_counts, ctry_pct)

# ── Legend ───────────────────────────────────────────────────────────────────
legend_patches = [
    mpatches.Patch(color=COL_SOUTH, label='South phylogroup'),
    mpatches.Patch(color=COL_NORTH, label='North phylogroup'),
    mpatches.Patch(color=COL_MISS,  label='Not classified'),
]
fig.legend(handles=legend_patches, loc='lower center', ncol=3,
           bbox_to_anchor=(0.5, 0.005), frameon=False,
           prop={'family': FONT, 'size': FS_LEGEND})

# ── Save ─────────────────────────────────────────────────────────────────────
for ext in ('pdf', 'png'):
    out = f"{OUTPUT_PREFIX}.{ext}"
    plt.savefig(out, dpi=DPI, facecolor='white')
    print(f"Saved: {out}")

print("Done.")
