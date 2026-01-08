# EEC scRNA-seq pipeline (PRJNA786266, 6 samples)

This repository contains an end-to-end scRNA-seq processing and annotation pipeline for six 10x Genomics samples from endometrial tissues (EEC vs normal controls). The workflow starts from SRA download and FASTQ conversion, and proceeds through Cell Ranger quantification, CellBender background removal, scDblFinder doublet detection, barcode intersection (Isec), Seurat SCTransform, clustering/UMAP, SingleR annotation, and manual-annotation evidence (FeaturePlot, DotPlot, heatmaps). All plots are generated as PNG only.

The project also includes course-required deliverables: per-sample QC plots, doublet UMAPs, clustering/annotation UMAPs, manual-annotation evidence panels, a QC summary table, and group-level cell-type composition figures. Final reports are provided in both Chinese and English.

## Dataset and groups

Dataset: PRJNA786266 / SRP349751 (6 samples).

| sample_id | group | label in meta |
| --- | --- | --- |
| SRR17165223 | Disease | EEC_3 |
| SRR17165224 | Disease | EEC_2 |
| SRR17165229 | Disease | EEC_1 |
| SRR17165227 | Control | normal_3 |
| SRR17165230 | Control | normal_2 |
| SRR17165231 | Control | normal_1 |

Group labels are derived from `meta/SraRunTable.csv`.

## Pipeline overview

1. SRA download and FASTQ conversion.
2. Cell Ranger quantification (GRCh38 reference).
3. CellBender remove-background.
4. scDblFinder doublet detection with standardized QC.
5. Isec barcode intersection (Cell Ranger filtered + CellBender + scDblFinder singlet + QC-passed).
6. Seurat SCTransform, UMAP, clustering.
7. SingleR annotation.
8. Manual annotation evidence: FeaturePlot panels, DotPlot panels, heatmaps.
9. Group-level composition plots and figure panels.

## QC rules and intersection definition

QC filtering is applied prior to scDblFinder and Isec:

- percent.mt <= 15%
- UMI >= 500
- genes >= 200
- Upper bounds computed per sample:
  - UMI <= min(MAD_high_log10, Q99.5, 80000)
  - genes <= min(MAD_high_log10, Q99.5, 11000)

Isec defines the final retained barcodes as the intersection of:

- Cell Ranger filtered barcodes
- CellBender cell_barcodes
- scDblFinder singlets
- QC-passed barcodes

## Directory layout

- `data/`: raw data downloads and FASTQ files.
- `qc/`: Cell Ranger, CellBender, scDblFinder, and Isec outputs.
- `analysis/`: Seurat objects, UMAPs, SingleR, manual annotation, and figures.
- `src/`: main pipeline scripts (01-10).
- `utils/`: plotting utilities, QC tables, and panel assembly.
- `final_report_en.md`: final report.

## Scripted workflow

Main scripts in `src/` (run in order):

- `src/1.download_data.sh`
- `src/2.convert2fq.sh`
- `src/3.cellranger.sh`
- `src/4.cellbender.sh`
- `src/5.scdblfinder.sh` (calls `src/5.scdblfinder.R`)
- `src/6.isec.sh` (calls `src/6.isec.R`)
- `src/7.seurat_sct.sh` (calls `src/7.seurat_sct.R`)
- `src/8.umap_clustering.sh` (calls `src/8.umap_clustering.R`)
- `src/9.singler_annot.sh` (calls `src/9.singler_annot.R`)
- `src/10.manual_annotation.sh` (calls `src/10.manual_annotation.R`)

Manual-annotation plotting can be run without recomputing markers by setting:

```
MANUAL_PLOT_ONLY=1 bash src/10.manual_annotation.sh
```

Utilities in `utils/` generate QC plots, doublet UMAPs, QC tables, composition plots, and composite panels:

- `utils/qc_plots.sh`
- `utils/doublet_umap.sh`
- `utils/qc_table_spec2.sh`
- `utils/celltype_composition.sh`
- `utils/make_panels.sh`

## Key outputs

Per-sample outputs:

- Cell Ranger outputs: `qc/cellranger/<SRR>/outs/`
- CellBender outputs: `qc/cellbender/<SRR>/`
- scDblFinder outputs: `qc/scdblfinder/<SRR>/`
- Isec outputs: `qc/isec/<SRR>/`
- Seurat objects: `analysis/seurat/<SRR>/`
- UMAP/clustering: `analysis/umap_clustering/<SRR>/`
- SingleR annotations: `analysis/singler/<SRR>/`
- Manual annotation evidence: `analysis/manual_annotation/<SRR>/`
- Heatmaps: `analysis/heatmap_evidence/<SRR>/`

Project-level outputs:

- QC summary table: `analysis/Table1_QC_SPEC2.tsv`
- QC panels: `analysis/figures/panels/panel_qc_violin.png`, `analysis/figures/panels/panel_qc_scatter.png`
- Grouped panels: `analysis/figures/panels/` (Disease vs Control)
- Composition plots: `analysis/celltype_composition/`
- Reports: `final_report_en.md`

All generated figures are PNG only.

## Manual annotation and evidence requirements

Manual annotation follows a marker-based evidence chain and is required to satisfy course specifications:

- Cluster marker statistics (FindAllMarkers).
- FeaturePlot marker panels covering major lineages.
- DotPlot marker panels for cluster-level comparison.
- Heatmaps for top markers per cluster and marker panels by manual cell type.
- A cluster-to-cell-type mapping table.

These outputs support and refine SingleR labels and are included in the final reports.

## Software and references

- sra-tools (prefetch, fasterq-dump)
- Cell Ranger (GRCh38 reference: `refdata-gex-GRCh38-2024-A`)
- CellBender remove-background (model=full, epochs=150, fpr=0.01)
- R: Seurat, SingleR, scDblFinder, celldex, ggplot2, dplyr, patchwork

## Notes

- SingleR uses celldex references and may require internet access on first run to download caches.
- Figures are designed for course reporting and are assembled into Disease/Control panels in `analysis/figures/panels/`.
- The QC summary table is presented in compact form in the main report, with full metrics in the appendix.
