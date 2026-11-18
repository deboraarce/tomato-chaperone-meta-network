"""Systematic search of tomato abiotic stress RNA-Seq datasets in NCBI GEO.

Arce DP & Pratta GR (2026)

Queries GEO DataSets (Entrez) for Solanum lycopersicum high-throughput
sequencing studies on abiotic stress, and writes the records found with a
preliminary classification of available processed files. Records are then
screened manually (tissue, stress type, number of biological replicates,
inclusion and exclusion reasons) following PRISMA.

Requirements: biopython
Usage (from the repository root):
    NCBI_EMAIL=you@institution.org python scripts/00_geo_systematic_search.py
"""

from __future__ import annotations

import csv
import json
import os
import time
from datetime import datetime
from pathlib import Path

from Bio import Entrez

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

OUTPUT_DIR = Path("results/systematic_review")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
RAW_RESULTS_PATH = OUTPUT_DIR / "geo_raw_results.json"
TABLE_PATH = OUTPUT_DIR / "geo_results_table.csv"
ELIGIBILITY_TABLE_PATH = OUTPUT_DIR / "geo_processed_data_table.csv"
DESEQ_INPUT_TABLE_PATH = OUTPUT_DIR / "geo_deseq_input_table.csv"
LOG_PATH = OUTPUT_DIR / "search_log.txt"

NCBI_EMAIL = os.getenv("NCBI_EMAIL")
if not NCBI_EMAIL:
    raise SystemExit("Set the NCBI_EMAIL environment variable (required by NCBI Entrez).")
Entrez.email = NCBI_EMAIL
Entrez.tool = "prisma_tomato_abiotic_rnaseq_review"

DB = "gds"
RETMAX = 500
STRESS_TERMS = [
    "drought",
    "water deficit",
    "salt",
    "salinity",
    "heat",
    "high temperature",
    "cold",
    "chilling",
    "freezing",
    "flooding",
    "waterlogging",
    "hypoxia",
    "osmotic",
    "cadmium",
    "heavy metal",
    "UV",
    "ozone",
]

QUERY = (
    '("Solanum lycopersicum"[Organism]) '
    f"AND ({' OR '.join(STRESS_TERMS)}) "
    'AND ("expression profiling by high throughput sequencing"[DataSet Type])'
)

PROCESSED_DATA_TERMS = (
    "count", "counts", "fpkm", "tpm", "rpkm", "normalized",
    "normalised", "matrix", "quantification", "dge",
)
RAW_DATA_TERMS = ("sra", "fastq", "raw", "reads", "bam", "sam")


# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------

def run_search(query: str, db: str = DB, retmax: int = RETMAX) -> dict:
    handle = Entrez.esearch(db=db, term=query, retmax=retmax)
    try:
        return Entrez.read(handle)
    finally:
        handle.close()


def fetch_summaries(id_list: list[str], db: str = DB, batch_size: int = 100) -> list[dict]:
    summaries: list[dict] = []
    for start in range(0, len(id_list), batch_size):
        batch = id_list[start : start + batch_size]
        handle = Entrez.esummary(db=db, id=",".join(batch))
        try:
            result = Entrez.read(handle)
        finally:
            handle.close()
        summaries.extend(result)
        time.sleep(0.4)
    return summaries


def parse_summary(item: dict) -> dict:
    """Record fields plus empty columns for manual screening."""
    return {
        "GEO_accession": item.get("Accession", ""),
        "title": item.get("title", ""),
        "GPL_platform": item.get("GPL", ""),
        "n_samples": item.get("n_samples", ""),
        "PDAT_fecha_publicacion": item.get("PDAT", ""),
        "summary": (item.get("summary", "") or "")[:300],
        "tejido": "",
        "tipo_estres": "",
        "n_replicas_biologicas": "",
        "archivo_procesado_disponible": "",
        "incluido_en_revision": "",
        "motivo_exclusion": "",
    }


def detect_processed_files(item: dict) -> dict:
    """Keyword-based classification of available processed or raw data."""
    ftp_link = item.get("FTPLink", "") or ""
    supp_file = item.get("suppFile", "") or ""
    title = item.get("title", "") or ""
    summary = item.get("summary", "") or ""

    metadata_text = f"{title} {summary} {supp_file} {ftp_link}".lower()

    has_processed_signal = any(term in metadata_text for term in PROCESSED_DATA_TERMS)
    has_raw_signal = any(term in metadata_text for term in RAW_DATA_TERMS)

    if has_processed_signal and not has_raw_signal:
        status = "procesados"
    elif has_raw_signal and not has_processed_signal:
        status = "crudos"
    else:
        status = "revision_manual"

    is_count_like = any(term in metadata_text for term in ("count", "counts", "dge"))
    is_normalized_like = any(term in metadata_text for term in ("tpm", "fpkm", "rpkm", "normalized", "normalised"))
    has_matrix_hint = "matrix" in metadata_text

    if is_count_like or (has_matrix_hint and not is_normalized_like):
        de_seq_input = "raw count matrix"
        is_deseq_valid = True
    elif is_normalized_like:
        de_seq_input = "normalized matrix"
        is_deseq_valid = False
    elif has_processed_signal:
        de_seq_input = "processed file, type to be checked manually"
        is_deseq_valid = False
    else:
        de_seq_input = "no count or normalized table detected"
        is_deseq_valid = False

    return {
        "GEO_accession": item.get("Accession", ""),
        "title": title,
        "suppFile": supp_file,
        "FTPLink": ftp_link,
        "processed_files_found": "matrix" if has_matrix_hint else "",
        "raw_files_found": "raw" if has_raw_signal else "",
        "status": status,
        "deSeq_input_type": de_seq_input,
        "is_deseq_valid": is_deseq_valid,
        "notes": "SRA/FASTQ/BAM/SAM files not used as input",
    }


def save_results(raw_record: dict, rows: list[dict]) -> None:
    with RAW_RESULTS_PATH.open("w", encoding="utf-8") as file_obj:
        json.dump(raw_record, file_obj, indent=2, default=str)

    fieldnames = list(rows[0].keys()) if rows else list(parse_summary({}).keys())
    with TABLE_PATH.open("w", newline="", encoding="utf-8") as file_obj:
        writer = csv.DictWriter(file_obj, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def save_eligibility_table(rows: list[dict]) -> None:
    fieldnames = [
        "GEO_accession", "title", "suppFile", "FTPLink",
        "processed_files_found", "raw_files_found", "status",
        "deSeq_input_type", "is_deseq_valid", "notes",
    ]
    with ELIGIBILITY_TABLE_PATH.open("w", newline="", encoding="utf-8") as file_obj:
        writer = csv.DictWriter(file_obj, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def save_deseq_input_table(rows: list[dict]) -> None:
    fieldnames = [
        "GEO_accession", "title", "suppFile", "FTPLink",
        "processed_files_found", "deSeq_input_type", "notes",
    ]
    with DESEQ_INPUT_TABLE_PATH.open("w", newline="", encoding="utf-8") as file_obj:
        writer = csv.DictWriter(file_obj, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_log(total_hits: int, retrieved: int) -> None:
    with LOG_PATH.open("a", encoding="utf-8") as file_obj:
        file_obj.write(
            f"[{datetime.now().isoformat()}] DB={DB} | "
            f"Query={QUERY} | Count={total_hits} | Retrieved={retrieved}\n"
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print(f"Query:\n{QUERY}\n")

    raw_record = run_search(QUERY)
    id_list = raw_record.get("IdList", [])
    total_hits = int(raw_record.get("Count", 0))

    print(f"Records found: {total_hits}")
    print(f"Records retrieved (retmax={RETMAX}): {len(id_list)}")

    if not id_list:
        print("No records found.")
        return

    summaries = fetch_summaries(id_list)
    rows = [parse_summary(item) for item in summaries]
    save_results(raw_record, rows)

    eligibility_rows = [detect_processed_files(item) for item in summaries]
    processed_rows = [row for row in eligibility_rows if row["status"] == "procesados"]
    de_seq_rows = [row for row in eligibility_rows if row["is_deseq_valid"]]
    save_eligibility_table(eligibility_rows)
    save_deseq_input_table(de_seq_rows)
    write_log(total_hits=total_hits, retrieved=len(id_list))

    print(f"{len(rows)} records written to {TABLE_PATH.name}")
    print(f"{len(processed_rows)} records with processed data in {ELIGIBILITY_TABLE_PATH.name}")
    print(f"{len(de_seq_rows)} records with count tables in {DESEQ_INPUT_TABLE_PATH.name}")


if __name__ == "__main__":
    main()
