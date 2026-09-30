#!/usr/bin/env python3
"""
build_phylogroup_tables.py
==========================
Builds Table A (phylogroup composition by country and by host group) and
Table B (per-isolate cross-gene phylogroup matrix) from FigTree-exported
clade assignment lists and the AP3 metadata spreadsheet.

Inputs
------
  GENE_NORTE.txt / GENE_SUL.txt  – one file per gene per phylogroup,
      each line in the format exported by FigTree:
          EPI ISL 123456|A/species/Country/ID/Year/Date|Year
      or, for GenBank-only isolates:
          A/species/Country/ID/Year
  Appendix_3_-_subtype_counts_oficial.xlsx  – AP3 metadata (host, location)
  Table_B.csv  – previous cross-gene table (used for host_group / country
      when AP3 metadata is unavailable; optional)

Outputs
-------
  Table_A_country.csv   – phylogroup composition per gene × country (S / N)
  Table_A_hostgroup.csv – phylogroup composition per gene × host group
  Table_B_new.csv       – per-isolate cross-gene matrix with class assignment
  phylogroup_report.txt – QC report (conflicts, missing files, counts)

Usage
-----
  python3 build_phylogroup_tables.py                      # defaults below
  python3 build_phylogroup_tables.py --dir ./assignments --ap3 AP3.xlsx
"""

import re
import csv
import sys
import argparse
import collections
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

try:
    import openpyxl
except ImportError:
    sys.exit("Missing dependency: pip install openpyxl")

# ── Configuration ─────────────────────────────────────────────────────────────
GENES = ["PB2", "PB1", "PA", "NP", "MP"]

# Expected file names; script searches the input directory for these.
FILE_PATTERNS = {
    gene: {"N": f"{gene}_NORTE.txt", "S": f"{gene}_SUL.txt"}
    for gene in GENES
}

# Host-group assignment rules (substring match on host string, lower-case).
# Ordered by priority – first match wins.
HOST_GROUP_RULES = [
    ("Waterfowl",           ["duck", "teal", "pintail", "pochard", "shoveler",
                              "wigeon", "gadwall", "mallard", "goose", "swan",
                              "merganser", "scoter", "eider", "bufflehead",
                              "goldeneye", "canvasback", "redhead", "ruddy",
                              "whistling duck", "comb duck", "muscovy",
                              "brazilian teal", "ringed teal", "silver teal",
                              "white-cheeked pintail", "yellow-billed pintail",
                              "yellow-billed teal", "cinnamon teal",
                              "rosy-billed pochard"]),
    ("Shorebirds",          ["sandpiper", "plover", "turnstone", "dunlin",
                              "knot", "dowitcher", "snipe", "curlew", "godwit",
                              "phalarope", "willet", "oystercatcher",
                              "black-necked stilt", "white-backed stilt",
                              "stilt", "avocet", "sanderling", "calidris",
                              "charadri", "whimbrel", "lapwing"]),
    ("Gulls/Seabirds",      ["gull", "tern", "skimmer", "petrel", "albatross",
                              "booby", "frigatebird", "pelican", "cormorant",
                              "gannet", "shearwater", "skua", "egret",
                              "heron", "ibis", "penguin", "moorhen",
                              "coot", "seagull", "kelp gull", "dolphin gull"]),
    ("Terrestrial/Domestic",["chicken", "turkey", "guineafowl", "meleagris",
                              "gallopavo", "gallus", "pheasant", "quail",
                              "partridge", "pigeon", "dove", "tinamou",
                              "parrot", "rapa nui"]),
    ("Environmental/Other", ["environment", "water", "sediment",
                              "fecal", "cloacal", "nasal", "wild bird"]),
]

# ── Parsing helpers ───────────────────────────────────────────────────────────
def parse_epi(line: str):
    """Return canonical EPI_ISL_XXXXXX or None."""
    m = re.search(r'EPI[\s_]*ISL[\s_]*(\d+)', line, re.I)
    return f"EPI_ISL_{m.group(1)}" if m else None


def parse_genbank_id(line: str):
    """Return GB_<strain_code> for GenBank-only entries."""
    h = re.sub(r'\|.*', '', line.strip().lstrip('('))
    h = re.sub(r'\(.*?\)', '', h)
    h = re.sub(r'GenBank.*', '', h, flags=re.I)
    h = re.sub(r'[\s\xa0_]+', '_', h).strip('_').lower()
    parts = [p for p in h.split('/') if p and len(p) > 2]
    if len(parts) >= 2:
        return f"GB_{parts[-2]}"
    return f"GB_{h}" if h else None


def get_id(line: str):
    """Extract canonical sequence ID from one text-file line."""
    line = line.strip()
    if not line or line.startswith('#'):
        return None
    epi = parse_epi(line)
    if epi:
        return epi
    if line.startswith('A/') or line.startswith('(A/'):
        return parse_genbank_id(line)
    return None


def load_assignment_file(path: Path):
    """Return set of canonical IDs in the file."""
    ids = set()
    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            i = get_id(line)
            if i:
                ids.add(i)
    return ids


# ── Metadata loading ──────────────────────────────────────────────────────────
def load_ap3(path: Path, sheet: str = "Sequences"):
    """Return dict epi → {host, location, subtype, year}."""
    wb  = openpyxl.load_workbook(path, data_only=True)
    ws  = wb[sheet]
    meta = {}
    for row in list(ws.iter_rows(values_only=True))[1:]:
        src = str(row[0]).strip() if row[0] else ''
        if not src:
            continue
        info = {
            'host':     str(row[1]).strip() if row[1] else '',
            'subtype':  str(row[2]).strip() if row[2] else '',
            'location': str(row[4]).strip() if row[4] else '',
            'year':     str(row[6]).strip() if row[6] else '',
            'isolate':  str(row[5]).strip() if row[5] else '',
        }
        m = re.search(r'EPI_ISL_(\d+)', src)
        if m:
            meta[f"EPI_ISL_{m.group(1)}"] = info
        else:
            gb = re.search(r'GenBank[:\s]+([\w\.]+)', src, re.I)
            if gb:
                nn = parse_genbank_id(info['isolate'])
                if nn:
                    meta[nn] = info
    return meta


def load_table_b_csv(path: Path):
    """Load previous Table B for host_group / country fallback."""
    out = {}
    with open(path, encoding='utf-8', errors='replace') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            epi = row.get('EPI_ISL', '').strip()
            if epi:
                key = f"EPI_ISL_{epi}" if not epi.startswith('EPI') else epi
                out[key] = {
                    'host_group': row.get('host_group', ''),
                    'country':    row.get('country', ''),
                    'host':       row.get('host', ''),
                    'isolate':    row.get('isolate', ''),
                }
    return out


def assign_host_group(host_str: str) -> str:
    h = host_str.lower()
    for group, keywords in HOST_GROUP_RULES:
        if any(k in h for k in keywords):
            return group
    return "Environmental/Other"


def parse_country(location: str) -> str:
    """Extract clean country name from AP3 location field."""
    loc = location.strip().lstrip()
    # Remove leading spaces / ISO codes
    loc = re.sub(r'^\s+', '', loc)
    # Take first token if comma-separated
    loc = loc.split(',')[0].strip()
    return loc if loc else 'Unknown'


# ── Core logic ────────────────────────────────────────────────────────────────
def build_assignments(input_dir: Path, report_lines: list):
    """
    Returns dict: epi → {gene: 'S'|'N'|'?'} and QC report lines.
    '?' means the ID appeared in BOTH files for that gene (conflict).
    """
    assignments = collections.defaultdict(dict)  # epi → {gene: pg}
    report_lines.append("\n" + "="*60)
    report_lines.append("PHYLOGROUP FILE LOADING")
    report_lines.append("="*60)

    for gene in GENES:
        n_file = input_dir / FILE_PATTERNS[gene]["N"]
        s_file = input_dir / FILE_PATTERNS[gene]["S"]

        n_ids = load_assignment_file(n_file) if n_file.exists() else set()
        s_ids = load_assignment_file(s_file) if s_file.exists() else set()

        if not n_file.exists():
            report_lines.append(f"\n  {gene}: ⚠️  NORTE file missing ({n_file.name})")
        if not s_file.exists():
            report_lines.append(f"\n  {gene}: ⚠️  SUL file missing ({s_file.name})")

        # Detect contamination: if most of N is inside S, S file is suspect
        overlap = n_ids & s_ids
        if overlap and len(overlap) / max(len(n_ids), 1) > 0.5:
            report_lines.append(
                f"\n  {gene}: ⚠️  SUL file appears contaminated — "
                f"{len(overlap)}/{len(n_ids)} NORTE IDs also present in SUL. "
                f"Treating SUL as (SUL_file − NORTE_file) = "
                f"{len(s_ids - n_ids)} unique S assignments."
            )
            s_ids = s_ids - n_ids  # strip contamination

        conflicts = n_ids & s_ids
        if conflicts:
            report_lines.append(
                f"\n  {gene}: ⚠️  {len(conflicts)} IDs in BOTH N and S "
                f"(marked as '?' — topologically impossible, check source files)"
            )
            for c in sorted(conflicts)[:5]:
                report_lines.append(f"      {c}")
            if len(conflicts) > 5:
                report_lines.append(f"      ... (+{len(conflicts)-5} more)")

        for i in n_ids - s_ids:
            assignments[i][gene] = 'N'
        for i in s_ids - n_ids:
            assignments[i][gene] = 'S'
        for i in conflicts:
            assignments[i][gene] = '?'

        report_lines.append(
            f"\n  {gene}: {len(n_ids)} N  |  {len(s_ids)} S  "
            f"|  {len(conflicts)} conflicts  |  "
            f"{len(n_ids)+len(s_ids)-len(conflicts)} total assigned"
        )

    return dict(assignments)


def classify(gene_calls: dict) -> str:
    """
    Classify an isolate based on its gene-level phylogroup calls.
    gene_calls: {gene: 'S'|'N'|'?'}
    """
    vals = [v for v in gene_calls.values() if v in ('S', 'N')]
    if not vals:
        return 'Unassigned'
    unique = set(vals)
    if len(unique) == 1:
        return 'Consistent-South' if 'S' in unique else 'Consistent-North'
    return 'Reassortant'


# ── Table building ────────────────────────────────────────────────────────────
def build_table_b(assignments, meta_ap3, meta_csv, report_lines):
    """
    Returns list of dicts — one per isolate.
    """
    all_ids = set(assignments.keys())
    # Add isolates from meta that have ≥1 assignment but may not be in assignments
    rows = []
    report_lines.append("\n" + "="*60)
    report_lines.append("TABLE B CONSTRUCTION")
    report_lines.append("="*60)

    no_meta = 0
    for epi in sorted(all_ids):
        # Prefer AP3 metadata, fall back to CSV
        m = meta_ap3.get(epi) or {}
        c = meta_csv.get(epi) or {}

        host     = m.get('host')     or c.get('host', '')
        location = m.get('location') or c.get('country', '')
        isolate  = m.get('isolate')  or c.get('isolate', '')
        subtype  = m.get('subtype', '')
        year     = m.get('year', '')

        if not host and not location:
            no_meta += 1

        host_group = c.get('host_group') or assign_host_group(host)
        country    = parse_country(location) if location else c.get('country', '')

        calls = assignments[epi]
        row = {
            'EPI_ISL':    epi.replace('EPI_ISL_', '') if 'GB_' not in epi else '',
            'isolate':    isolate,
            'year':       year,
            'subtype':    subtype,
            'country':    country,
            'host':       host,
            'host_group': host_group,
        }
        for gene in GENES:
            row[gene] = calls.get(gene, '')
        n_calls = sum(1 for g in GENES if calls.get(g) in ('S', 'N'))
        row['n_genes'] = n_calls
        row['class']   = classify(calls)
        rows.append(row)

    report_lines.append(f"\n  Total isolates: {len(rows)}")
    report_lines.append(f"  Without metadata: {no_meta}")
    cls_ctr = collections.Counter(r['class'] for r in rows)
    for cls, n in sorted(cls_ctr.items()):
        report_lines.append(f"  {cls}: {n}")

    return rows


def build_table_a_by_field(table_b_rows, field, field_label, report_lines):
    """
    Returns list of dicts for Table A grouped by `field` (e.g. 'country' or 'host_group').
    """
    groups = sorted(set(r[field] for r in table_b_rows if r[field]))
    out = []
    for grp in groups:
        grp_rows = [r for r in table_b_rows if r[field] == grp]
        row = {field_label: grp}
        total_s = total_n = 0
        for gene in GENES:
            s = sum(1 for r in grp_rows if r[gene] == 'S')
            n = sum(1 for r in grp_rows if r[gene] == 'N')
            row[f"{gene}_S"] = s
            row[f"{gene}_N"] = n
            total_s += s
            total_n += n
        row['Total_S'] = total_s
        row['Total_N'] = total_n
        out.append(row)

    # TOTAL row
    tot = {field_label: 'TOTAL'}
    for gene in GENES:
        tot[f"{gene}_S"] = sum(r[f"{gene}_S"] for r in out)
        tot[f"{gene}_N"] = sum(r[f"{gene}_N"] for r in out)
    tot['Total_S'] = sum(r['Total_S'] for r in out)
    tot['Total_N'] = sum(r['Total_N'] for r in out)
    out.append(tot)

    return out


# ── Output helpers ────────────────────────────────────────────────────────────
def write_csv(path: Path, rows: list, fieldnames: list):
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    print(f"  ✓ {path}")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dir",     default=".",
                   help="Directory containing GENE_NORTE/SUL.txt files")
    p.add_argument("--ap3",     default="Appendix_3_-_subtype_counts_oficial.xlsx",
                   help="AP3 metadata spreadsheet")
    p.add_argument("--ap3_sheet", default="Sequences")
    p.add_argument("--table_b_csv", default="Table_B.csv",
                   help="Previous Table B CSV (used as metadata fallback)")
    p.add_argument("--out",     default=".",
                   help="Output directory")
    args = p.parse_args()

    input_dir  = Path(args.dir)
    output_dir = Path(args.out)
    output_dir.mkdir(parents=True, exist_ok=True)

    report_lines = ["PHYLOGROUP TABLE BUILD REPORT", "="*60]

    # ── Load metadata ─────────────────────────────────────────────────────────
    report_lines.append("\nMETADATA SOURCES")
    meta_ap3, meta_csv = {}, {}

    ap3_path = Path(args.ap3)
    if ap3_path.exists():
        meta_ap3 = load_ap3(ap3_path, args.ap3_sheet)
        report_lines.append(f"  AP3 xlsx: {len(meta_ap3)} records  ({ap3_path.name})")
    else:
        report_lines.append(f"  ⚠️  AP3 xlsx not found: {ap3_path}")

    csv_path = Path(args.table_b_csv)
    if csv_path.exists():
        meta_csv = load_table_b_csv(csv_path)
        report_lines.append(f"  Table_B CSV: {len(meta_csv)} records  ({csv_path.name})")
    else:
        report_lines.append(f"  ℹ️  Table_B CSV not found: {csv_path} (optional)")

    # ── Build assignments ─────────────────────────────────────────────────────
    assignments = build_assignments(input_dir, report_lines)

    # ── Build Table B ─────────────────────────────────────────────────────────
    table_b = build_table_b(assignments, meta_ap3, meta_csv, report_lines)

    # ── Build Table A ─────────────────────────────────────────────────────────
    table_a_country   = build_table_a_by_field(table_b, 'country',    'Country',    report_lines)
    table_a_hostgroup = build_table_a_by_field(table_b, 'host_group', 'Host_group', report_lines)

    # ── Write outputs ─────────────────────────────────────────────────────────
    report_lines.append("\n" + "="*60)
    report_lines.append("OUTPUT FILES")
    report_lines.append("="*60)

    # Table B
    b_fields = (['EPI_ISL', 'isolate', 'year', 'subtype', 'country',
                  'host', 'host_group'] + GENES + ['n_genes', 'class'])
    write_csv(output_dir / "Table_B_new.csv", table_b, b_fields)

    # Table A — country
    a_country_fields = ['Country']
    for g in GENES:
        a_country_fields += [f"{g}_S", f"{g}_N"]
    a_country_fields += ['Total_S', 'Total_N']
    write_csv(output_dir / "Table_A_country.csv", table_a_country, a_country_fields)

    # Table A — host group
    a_host_fields = ['Host_group']
    for g in GENES:
        a_host_fields += [f"{g}_S", f"{g}_N"]
    a_host_fields += ['Total_S', 'Total_N']
    write_csv(output_dir / "Table_A_hostgroup.csv", table_a_hostgroup, a_host_fields)

    # Report
    report_path = output_dir / "phylogroup_report.txt"
    report_lines.append(f"\n  ✓ {output_dir / 'Table_B_new.csv'}")
    report_lines.append(f"  ✓ {output_dir / 'Table_A_country.csv'}")
    report_lines.append(f"  ✓ {output_dir / 'Table_A_hostgroup.csv'}")
    report_lines.append(f"  ✓ {report_path}")

    report_text = "\n".join(report_lines)
    print(report_text)
    report_path.write_text(report_text, encoding='utf-8')


if __name__ == "__main__":
    main()
