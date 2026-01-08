# Utils

## QC Table

`qc_table.sh` runs `qc_table.R` to generate the QC summary table from the final Isec barcodes and CellBender-filtered counts.

### Output columns

- `min_umi` / `max_umi`: min/max UMI counts per cell (CellBender-filtered matrix)
- `min_detected_genes` / `max_detected_genes`: min/max detected genes per cell
- `max_percent_mito`: max mitochondrial percent per cell
- `cells_from_cellranger`: Cell Ranger estimated cells
- `cells_filtered_out`: `cells_from_cellranger - cells_passed_qc`
- `cells_passed_qc`: final Isec cell count
- `variable_genes`: set by `VARIABLE_GENES` (default `3000`)

### QC thresholds (current)

- `percent.mt <= 15%` (hard cap)
- Lower bounds: `UMI >= 500`, `genes >= 200`
- Upper bounds per sample:
  - `UMI <= min(MAD_high_log10, Q99.5, 80000)`
  - `genes <= min(MAD_high_log10, Q99.5, 11000)`

These thresholds are applied in `src/5.scdblfinder.R`, and the QC-passed barcodes are carried into Isec.

## QC Plots

`qc_plots.sh` runs `qc_plots.R` to generate per-sample QC plots from CellBender-filtered counts and QC-passed barcodes.

Outputs under `analysis/qc_plots/<SRR>/`:
- `<SRR>_qc_violin_metrics.png` (nCount, nFeature, percent.mt by QC pass/fail)
- `<SRR>_qc_scatter_ncount_nfeature.png` (with QC thresholds)

## Doublet UMAP

`doublet_umap.sh` runs `doublet_umap.R` to create UMAPs colored by scDblFinder class.

Outputs under `analysis/doublet/<SRR>/`:
- `<SRR>_umap_doublet.png`

## SPEC2 Table1 QC

`qc_table_spec2.sh` runs `qc_table_spec2.R` to build a SPEC2-compliant QC summary table.

Output:
- `analysis/Table1_QC_SPEC2.tsv`

## Cell type composition

`celltype_composition.sh` runs `celltype_composition.R` to generate cell type composition plots using the manual annotation map.

Outputs under `analysis/celltype_composition/`:
- `celltype_composition_by_sample.png`
- `celltype_composition_by_group.png`
- `celltype_composition_by_sample.tsv`
- `celltype_composition_by_group.tsv`

## Plot style

All plotting scripts read the same style env vars:
- `PLOT_FONT` (default `Helvetica`)
- `PLOT_BASE_SIZE` (default `8`)
- `PLOT_LABEL_SIZE` (optional; defaults to `PLOT_BASE_SIZE * 0.35`)
- `PLOT_DPI` (default `300`)
