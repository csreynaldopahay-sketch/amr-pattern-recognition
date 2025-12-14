import copy
import csv
import os
import re
from collections import Counter
from typing import Dict, List, Tuple

MAJORITY_THRESHOLD = 0.5
MISSING_THRESHOLD = 0.5

STATUS_MAP = {"S": 0, "I": 1, "R": 2}

METADATA_COLUMNS = [
    "isolate_code",
    "species",
    "region",
    "site",
    "environment",
    "sampling_source",
    "esbl",
]

SPECIES_CORRECTIONS = {
    "pseudomoans": "Pseudomonas",
    "vibrio cholarae": "Vibrio Cholerae",
}


def standardize_species(name: str) -> str:
    cleaned = re.sub(r"\s+", " ", name or "").strip()
    titled = cleaned.title() if cleaned else ""
    normalized = titled.lower()
    return SPECIES_CORRECTIONS.get(normalized, titled)


def standardize_antibiotic_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "", name or "").upper()


def normalize_esbl(value: str) -> str:
    token = (value or "").strip().upper()
    if not token:
        return ""
    if token in {"POSITIVE", "POS", "+"}:
        return "POS"
    if token in {"NEGATIVE", "NEG", "-"}:
        return "NEG"
    return token


def clean_result(value: str) -> str:
    if not value:
        return ""
    statuses = {ch for ch in value.upper() if ch in STATUS_MAP}
    if len(statuses) == 1:
        return statuses.pop()
    return ""


def parse_metadata_from_filename(path: str) -> Tuple[str, str]:
    base = os.path.basename(path)
    patterns = [
        r"1NET_P2-AMR_(.+?) - Copy - (.+)\.csv",
        r"1NET_P2-AMR_(.+?) - (.+)\.csv",
    ]
    for pattern in patterns:
        match = re.match(pattern, base)
        if match:
            region = match.group(1).strip()
            site = os.path.splitext(match.group(2).strip())[0]
            return region, site
    return "Unknown", os.path.splitext(base)[0]


def read_antibiotic_mapping(reader: csv.reader) -> List[Tuple[str, int]]:
    header = None
    for row in reader:
        if len(row) > 2 and row[2].strip().upper() == "CODE":
            header = row
            break
    if header is None:
        return []
    antibiotic_row = next(reader, [])
    mic_row = next(reader, [])
    mapping: List[Tuple[str, int]] = []
    current_ab = ""
    for idx, name in enumerate(antibiotic_row):
        if name.strip():
            current_ab = standardize_antibiotic_name(name)
        mic_token = mic_row[idx].strip().upper() if idx < len(mic_row) else ""
        if mic_token == "INT." and current_ab:
            mapping.append((current_ab, idx))
    return mapping


def load_records(path: str) -> List[Dict[str, str]]:
    region, site = parse_metadata_from_filename(path)
    records: List[Dict[str, str]] = []
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        mapping = read_antibiotic_mapping(reader)
        for row in reader:
            if not any(cell.strip() for cell in row):
                continue
            isolate_code = (row[2] if len(row) > 2 else "").strip()
            species = standardize_species(row[3] if len(row) > 3 else "")
            if not isolate_code and not species:
                continue
            record: Dict[str, str] = {
                "isolate_code": isolate_code,
                "species": species,
                "region": region,
                "site": site,
                "environment": "Unknown",
                "sampling_source": "Unknown",
                "esbl": normalize_esbl(row[4] if len(row) > 4 else ""),
            }
            for antibiotic, idx in mapping:
                value = row[idx] if idx < len(row) else ""
                record[antibiotic] = clean_result(value)
            records.append(record)
    return records


def deduplicate(records: List[Dict[str, str]]) -> List[Dict[str, str]]:
    seen = set()
    result = []
    for record in records:
        key = record.get("isolate_code")
        if key and key in seen:
            continue
        if key:
            seen.add(key)
        result.append(record)
    return result


def compute_coverage(records: List[Dict[str, str]], antibiotics: List[str]) -> Dict[str, float]:
    total = len(records) or 1
    coverage = {}
    for ab in antibiotics:
        observed = sum(1 for r in records if r.get(ab))
        coverage[ab] = observed / total
    return coverage


def filter_by_missing(records: List[Dict[str, str]], antibiotics: List[str]) -> Tuple[List[Dict[str, str]], int]:
    filtered = []
    dropped = 0
    if not antibiotics:
        return records, dropped
    for record in records:
        missing = sum(1 for ab in antibiotics if not record.get(ab))
        if missing / len(antibiotics) > MISSING_THRESHOLD:
            dropped += 1
            continue
        filtered.append(record)
    return filtered, dropped


def encode_record(record: Dict[str, str], antibiotics: List[str]) -> Dict[str, object]:
    encoded: Dict[str, object] = {key: record.get(key, "") for key in METADATA_COLUMNS}
    num_tested = 0
    num_resistant = 0
    for ab in antibiotics:
        status = record.get(ab)
        if not status:
            encoded[ab] = ""
            encoded[f"{ab}_binary"] = ""
            continue
        encoded[ab] = STATUS_MAP[status]
        encoded[f"{ab}_binary"] = 1 if status == "R" else 0
        num_tested += 1
        if status == "R":
            num_resistant += 1
    encoded["num_antibiotics_tested"] = num_tested
    encoded["num_resistant"] = num_resistant
    encoded["MAR_index"] = round(num_resistant / num_tested, 4) if num_tested else ""
    encoded["MDR_flag"] = 1 if num_resistant >= 3 else 0
    return encoded


def write_csv(path: str, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_summary(
    path: str,
    total_records: int,
    kept_records: int,
    dropped_missing: int,
    duplicates_removed: int,
    coverage: Dict[str, float],
    kept_antibiotics: List[str],
) -> None:
    lines = [
        f"Total isolates ingested: {total_records}",
        f"Duplicates removed: {duplicates_removed}",
        f"Isolates removed for missing data: {dropped_missing}",
        f"Antibiotics retained (>= {MAJORITY_THRESHOLD:.0%} coverage): {', '.join(kept_antibiotics)}",
        "Antibiotic coverage:",
    ]
    for ab in sorted(coverage):
        lines.append(f"  - {ab}: {coverage[ab]:.2%}")
    with open(path, "w") as handle:
        handle.write("\n".join(lines))


def main() -> None:
    csv_files = [f for f in os.listdir(".") if f.lower().endswith(".csv")]
    ingested: List[Dict[str, str]] = []
    for file in csv_files:
        ingested.extend(load_records(file))
    duplicates_removed = len(ingested)
    all_records = deduplicate(ingested)
    duplicates_removed -= len(all_records)
    raw_records = copy.deepcopy(all_records)

    antibiotics = sorted(
        {key for record in all_records for key in record.keys() if key not in METADATA_COLUMNS}
    )
    coverage = compute_coverage(all_records, antibiotics)
    kept_antibiotics = [ab for ab, cov in coverage.items() if cov >= MAJORITY_THRESHOLD]

    filtered_records, dropped_missing = filter_by_missing(all_records, kept_antibiotics)

    processed_records = [encode_record(record, kept_antibiotics) for record in filtered_records]

    os.makedirs("data/processed", exist_ok=True)

    raw_fieldnames = METADATA_COLUMNS + antibiotics
    write_csv("data/processed/unified_raw.csv", raw_records, raw_fieldnames)

    processed_fieldnames = METADATA_COLUMNS + kept_antibiotics + [f"{ab}_binary" for ab in kept_antibiotics]
    processed_fieldnames += ["num_antibiotics_tested", "num_resistant", "MAR_index", "MDR_flag"]
    write_csv("data/processed/analysis_ready.csv", processed_records, processed_fieldnames)

    feature_fieldnames = ["isolate_code"] + kept_antibiotics
    feature_rows = [
        {"isolate_code": row["isolate_code"], **{ab: row[ab] for ab in kept_antibiotics}}
        for row in processed_records
    ]
    write_csv("data/processed/feature_matrix.csv", feature_rows, feature_fieldnames)

    metadata_fieldnames = METADATA_COLUMNS + ["num_antibiotics_tested", "num_resistant", "MAR_index", "MDR_flag"]
    metadata_rows = [
        {key: row.get(key, "") for key in metadata_fieldnames} for row in processed_records
    ]
    write_csv("data/processed/metadata.csv", metadata_rows, metadata_fieldnames)

    write_summary(
        "data/processed/preprocessing_summary.txt",
        total_records=len(raw_records) + duplicates_removed,
        kept_records=len(processed_records),
        dropped_missing=dropped_missing,
        duplicates_removed=duplicates_removed,
        coverage=coverage,
        kept_antibiotics=kept_antibiotics,
    )


if __name__ == "__main__":
    main()
