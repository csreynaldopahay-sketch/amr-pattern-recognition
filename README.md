## AMR data preprocessing (Phase 1)

This repository bundles the raw AMR CSV exports and a minimal preprocessing script that consolidates and cleans them into analysis-ready datasets.

### How to run

```bash
python preprocess.py
```

The script requires only Python 3 and reads all `*.csv` files in the repository root.

### What the script does

1. **Ingestion** – loads every site/region CSV and adds explicit metadata (Region, Site, Environment, Sampling source).
2. **Cleaning** – standardizes species and antibiotic codes, normalizes ESBL values, removes duplicate isolate codes, and discards impossible AST symbols.
3. **Missing data handling** – keeps antibiotics tested in ≥50% of isolates and drops isolates with >50% missing AST values among the retained drugs. All exclusions are logged in `data/processed/preprocessing_summary.txt`.
4. **Encoding** – converts AST outcomes to numeric form (S=0, I=1, R=2) and adds binary R vs non-R flags.
5. **Feature engineering** – recomputes MAR index and an MDR flag (≥3 resistant antibiotics).

### Outputs (`data/processed/`)

- `unified_raw.csv` – merged deduplicated dataset with letter AST results.
- `analysis_ready.csv` – metadata, encoded AST results, binary resistance flags, MAR, and MDR columns.
- `feature_matrix.csv` – isolate code + numeric resistance fingerprint (retained antibiotics only).
- `metadata.csv` – metadata with MAR/MDR metrics.
- `preprocessing_summary.txt` – counts of exclusions and antibiotic coverage.
