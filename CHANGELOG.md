# Changelog

All notable changes to this project will be documented in this file.

## [1.0.3] - 2026-07-08

### Release Notes

### Pipeline

- [Fixed] `publishDir` targets in `main.nf` now use closures instead of static strings, so per-gene outputs are published to the correct directory.
- [Changed] `--exonic` is now optional; the pipeline runs without it instead of requiring the file.

### Plot

- [Fixed] `aa_to_genomic()` compared strand against `"+"`/`"-"`, but Ensembl REST returns `1`/`-1`. Domain placement could be skipped or misplaced on the affected strand check.
- [Added] Variant track y-axis/legend adapts when a gene has no non-coding or no exonic data, instead of assuming both are always present.
- [Added] `set.seed(42)` before label placement (`ggrepel`) for reproducible plots across runs.
- [Added] Fallback density panel ("No UK Biobank variants") for genes with too few UK Biobank points to draw a density curve.
- [Changed] Legend/axis text sizes and domain legend label truncation retuned for readability with long auto-fetched domain names.
- [Changed] Genomic/non-coding/SV/exonic inputs are parsed as semicolon-delimited, comma-decimal CSVs (`read_csv2`) instead of comma-delimited, matching the actual input file format.
- [Changed] Coding filter now accepts CADD **or** REVEL passing threshold (previously required both), with `af_max` relaxed 0.0001 → 0.001 and `cadd_min` 25 → 20.
- [Changed] UK Biobank/gnomAD density panels normalized to a shared 0 to 1 scale instead of separate per-panel maxima.
- [Changed] Exported `*_genomic_table.csv` / `*_exonic_table.csv` now contain the full assembled variant tables used for plotting, not just the raw filtered columns.


## [1.0.2] - 2026-02-26

### Release Notes

- [Fixed] error: getBM() timeout
- Removed the dependency on Ensembl `biomart`.
- Using an offline `ensembldb` insted of biomart.

## [1.0.1] - 2026-02-19

### Release Notes

- [Fixed] Enhanced logging message if no variants are plotted.
- [Added] Genomic and exonic variants table to the output directory.

## [1.0.0] - 2026-01-31

### Release Notes

- Initial release for testing.
