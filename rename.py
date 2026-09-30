#!/usr/bin/env python3
"""
Padroniza os nomes dos arquivos .xml dentro de xmls/<GENE>/.

Cada gene tem dois tipos de analise:

  - Molecular clock (relogio molecular): marcador SKYRIDE ou COA
        -> {GENE}_MOLCLOCK_{SKYRIDE|COA}.xml

  - BSSVS (coalescente), dividido por regiao (north/south) e por
    preditor (host/local)
        -> {GENE}_BSSVS_{NORTH|SOUTH}_{HOST|LOCAL}.xml

Uso:
    python rename.py            # apenas mostra o mapeamento (dry-run)
    python rename.py --apply    # renomeia os arquivos de fato
"""
import argparse
import re
import sys
from pathlib import Path

XMLS_DIR = Path(__file__).resolve().parent / "xmls"


def classify(stem: str):
    """Retorna o novo 'stem' (sem extensao) para um nome de arquivo antigo,
    ou None se o arquivo nao se encaixar nos padroes conhecidos."""
    tokens = re.split(r"[^A-Za-z0-9]+", stem.upper())
    tokens = [t for t in tokens if t]

    north = "NORTH" in tokens
    south = "SOUTH" in tokens

    host = any(t.startswith("HOST") for t in tokens)
    local = any(t.startswith("LOC") for t in tokens)  # LOC, LOCAL, LOCALS

    skyride = "SKYRIDE" in tokens
    coa = "COA" in tokens

    if north or south:
        region = "NORTH" if north else "SOUTH"
        if host:
            predictor = "HOST"
        elif local:
            predictor = "LOCAL"
        else:
            return None
        return f"BSSVS_{region}_{predictor}"

    if skyride or coa:
        marker = "SKYRIDE" if skyride else "COA"
        return f"MOLCLOCK_{marker}"

    return None


def build_rename_plan(xmls_dir: Path):
    plan = []
    unmatched = []

    for gene_dir in sorted(p for p in xmls_dir.iterdir() if p.is_dir()):
        gene = gene_dir.name.upper()
        for xml_file in sorted(gene_dir.iterdir()):
            if xml_file.suffix.lower() != ".xml":
                continue

            new_stem = classify(xml_file.stem)
            if new_stem is None:
                unmatched.append(xml_file)
                continue

            new_name = f"{gene}_{new_stem}.xml"
            new_path = xml_file.with_name(new_name)
            plan.append((xml_file, new_path))

    return plan, unmatched


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Executa os renames. Sem essa flag, apenas mostra o que seria feito.",
    )
    args = parser.parse_args()

    if not XMLS_DIR.is_dir():
        sys.exit(f"Diretorio nao encontrado: {XMLS_DIR}")

    plan, unmatched = build_rename_plan(XMLS_DIR)

    # detecta colisoes de nome de destino dentro da mesma pasta
    seen = {}
    conflicts = []
    for old_path, new_path in plan:
        key = new_path
        if key in seen and seen[key] != old_path:
            conflicts.append((seen[key], old_path, new_path))
        seen[key] = old_path

    for old_path, new_path in plan:
        marker = "==" if old_path.resolve() == new_path.resolve() else "->"
        print(f"{old_path.relative_to(XMLS_DIR)} {marker} {new_path.name}")

    if unmatched:
        print("\nArquivos nao reconhecidos (nao foram renomeados):")
        for f in unmatched:
            print(f"  {f.relative_to(XMLS_DIR)}")

    if conflicts:
        print("\nCONFLITOS: mais de um arquivo mapeia para o mesmo nome final:")
        for old_a, old_b, new_path in conflicts:
            print(f"  {old_a.name} e {old_b.name} -> {new_path.name}")
        if args.apply:
            sys.exit("\nAbortando: resolva os conflitos antes de aplicar o rename.")

    if not args.apply:
        print("\n(dry-run) nenhuma alteracao foi feita. Use --apply para renomear.")
        return

    for old_path, new_path in plan:
        if old_path.resolve() == new_path.resolve():
            continue
        if new_path.exists():
            sys.exit(f"Destino ja existe, abortando: {new_path}")
        old_path.rename(new_path)

    print("\nRename concluido.")


if __name__ == "__main__":
    main()
