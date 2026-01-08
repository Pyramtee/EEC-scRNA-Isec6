#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd -P)"

out_dir="${PROJECT_ROOT}/analysis/figures/panels"
mkdir -p "${out_dir}"

if ! command -v ffmpeg >/dev/null 2>&1; then
  echo "[ERROR] ffmpeg not found" >&2
  exit 1
fi

tile_filter() {
  local w="$1"
  local h="$2"
  local gap="$3"
  echo "scale=${w}:${h}:force_original_aspect_ratio=decrease,pad=${w}:${h}:(ow-iw)/2:(oh-ih)/2:color=white,pad=${w}+${gap}:${h}+${gap}:${gap}/2:${gap}/2:color=white"
}

get_dim() {
  ffprobe -v error -select_streams v:0 -show_entries stream=width,height -of csv=p=0 "$1"
}

panel_6_grid() {
  local out="$1"
  local w="$2"
  local h="$3"
  local gap="$4"
  shift 4
  if [[ $# -ne 6 ]]; then
    echo "[ERROR] panel_6_grid needs 6 inputs" >&2
    exit 1
  fi
  local i0="$1" i1="$2" i2="$3" i3="$4" i4="$5" i5="$6"
  if [[ "${w}" == "auto" || "${h}" == "auto" ]]; then
    IFS=',' read -r w h <<<"$(get_dim "${i0}")"
  fi
  local tf
  tf="$(tile_filter "${w}" "${h}" "${gap}")"

  ffmpeg -y \
    -i "$i0" -i "$i1" -i "$i2" -i "$i3" -i "$i4" -i "$i5" \
    -filter_complex "\
[0:v]${tf}[v0];\
[1:v]${tf}[v1];\
[2:v]${tf}[v2];\
[3:v]${tf}[v3];\
[4:v]${tf}[v4];\
[5:v]${tf}[v5];\
[v0][v1][v2]hstack=inputs=3[top];\
[v3][v4][v5]hstack=inputs=3[bottom];\
[top][bottom]vstack=inputs=2[out]" \
    -map "[out]" "$out"
}

panel_3_row() {
  local out="$1"
  local w="$2"
  local h="$3"
  local gap="$4"
  shift 4
  if [[ $# -ne 3 ]]; then
    echo "[ERROR] panel_3_row needs 3 inputs" >&2
    exit 1
  fi
  local i0="$1" i1="$2" i2="$3"
  if [[ "${w}" == "auto" || "${h}" == "auto" ]]; then
    IFS=',' read -r w h <<<"$(get_dim "${i0}")"
  fi
  local tf
  tf="$(tile_filter "${w}" "${h}" "${gap}")"
  ffmpeg -y \
    -i "$i0" -i "$i1" -i "$i2" \
    -filter_complex "\
[0:v]${tf}[v0];\
[1:v]${tf}[v1];\
[2:v]${tf}[v2];\
[v0][v1][v2]hstack=inputs=3[out]" \
    -map "[out]" "$out"
}

panel_3_col() {
  local out="$1"
  local w="$2"
  local h="$3"
  local gap="$4"
  shift 4
  if [[ $# -ne 3 ]]; then
    echo "[ERROR] panel_3_col needs 3 inputs" >&2
    exit 1
  fi
  local i0="$1" i1="$2" i2="$3"
  if [[ "${w}" == "auto" || "${h}" == "auto" ]]; then
    IFS=',' read -r w h <<<"$(get_dim "${i0}")"
  fi
  local tf
  tf="$(tile_filter "${w}" "${h}" "${gap}")"
  ffmpeg -y \
    -i "$i0" -i "$i1" -i "$i2" \
    -filter_complex "\
[0:v]${tf}[v0];\
[1:v]${tf}[v1];\
[2:v]${tf}[v2];\
[v0][v1][v2]vstack=inputs=3[out]" \
    -map "[out]" "$out"
}

samples=(SRR17165223 SRR17165224 SRR17165227 SRR17165229 SRR17165230 SRR17165231)
disease=(SRR17165223 SRR17165224 SRR17165229)
control=(SRR17165227 SRR17165230 SRR17165231)

gap=40
gap_zero=0
gap_wide=160

panel_6_grid "${out_dir}/panel_qc_violin.png" auto auto "${gap}" \
  "${PROJECT_ROOT}/analysis/qc_plots/${samples[0]}/${samples[0]}_qc_violin_metrics.png" \
  "${PROJECT_ROOT}/analysis/qc_plots/${samples[1]}/${samples[1]}_qc_violin_metrics.png" \
  "${PROJECT_ROOT}/analysis/qc_plots/${samples[2]}/${samples[2]}_qc_violin_metrics.png" \
  "${PROJECT_ROOT}/analysis/qc_plots/${samples[3]}/${samples[3]}_qc_violin_metrics.png" \
  "${PROJECT_ROOT}/analysis/qc_plots/${samples[4]}/${samples[4]}_qc_violin_metrics.png" \
  "${PROJECT_ROOT}/analysis/qc_plots/${samples[5]}/${samples[5]}_qc_violin_metrics.png"

panel_6_grid "${out_dir}/panel_qc_scatter.png" auto auto "${gap}" \
  "${PROJECT_ROOT}/analysis/qc_plots/${samples[0]}/${samples[0]}_qc_scatter_ncount_nfeature.png" \
  "${PROJECT_ROOT}/analysis/qc_plots/${samples[1]}/${samples[1]}_qc_scatter_ncount_nfeature.png" \
  "${PROJECT_ROOT}/analysis/qc_plots/${samples[2]}/${samples[2]}_qc_scatter_ncount_nfeature.png" \
  "${PROJECT_ROOT}/analysis/qc_plots/${samples[3]}/${samples[3]}_qc_scatter_ncount_nfeature.png" \
  "${PROJECT_ROOT}/analysis/qc_plots/${samples[4]}/${samples[4]}_qc_scatter_ncount_nfeature.png" \
  "${PROJECT_ROOT}/analysis/qc_plots/${samples[5]}/${samples[5]}_qc_scatter_ncount_nfeature.png"

panel_3_row "${out_dir}/panel_doublet_umap_disease.png" auto auto "${gap}" \
  "${PROJECT_ROOT}/analysis/doublet/${disease[0]}/${disease[0]}_umap_doublet.png" \
  "${PROJECT_ROOT}/analysis/doublet/${disease[1]}/${disease[1]}_umap_doublet.png" \
  "${PROJECT_ROOT}/analysis/doublet/${disease[2]}/${disease[2]}_umap_doublet.png"
panel_3_row "${out_dir}/panel_doublet_umap_control.png" auto auto "${gap}" \
  "${PROJECT_ROOT}/analysis/doublet/${control[0]}/${control[0]}_umap_doublet.png" \
  "${PROJECT_ROOT}/analysis/doublet/${control[1]}/${control[1]}_umap_doublet.png" \
  "${PROJECT_ROOT}/analysis/doublet/${control[2]}/${control[2]}_umap_doublet.png"

panel_3_row "${out_dir}/panel_cluster_umap_disease.png" auto auto "${gap_wide}" \
  "${PROJECT_ROOT}/analysis/umap_clustering/${disease[0]}/${disease[0]}_umap_clusters.png" \
  "${PROJECT_ROOT}/analysis/umap_clustering/${disease[1]}/${disease[1]}_umap_clusters.png" \
  "${PROJECT_ROOT}/analysis/umap_clustering/${disease[2]}/${disease[2]}_umap_clusters.png"
panel_3_row "${out_dir}/panel_cluster_umap_control.png" auto auto "${gap_wide}" \
  "${PROJECT_ROOT}/analysis/umap_clustering/${control[0]}/${control[0]}_umap_clusters.png" \
  "${PROJECT_ROOT}/analysis/umap_clustering/${control[1]}/${control[1]}_umap_clusters.png" \
  "${PROJECT_ROOT}/analysis/umap_clustering/${control[2]}/${control[2]}_umap_clusters.png"

panel_3_row "${out_dir}/panel_singler_umap_disease.png" auto auto "${gap}" \
  "${PROJECT_ROOT}/analysis/singler/${disease[0]}/${disease[0]}_umap_singler_refined.png" \
  "${PROJECT_ROOT}/analysis/singler/${disease[1]}/${disease[1]}_umap_singler_refined.png" \
  "${PROJECT_ROOT}/analysis/singler/${disease[2]}/${disease[2]}_umap_singler_refined.png"
panel_3_row "${out_dir}/panel_singler_umap_control.png" auto auto "${gap}" \
  "${PROJECT_ROOT}/analysis/singler/${control[0]}/${control[0]}_umap_singler_refined.png" \
  "${PROJECT_ROOT}/analysis/singler/${control[1]}/${control[1]}_umap_singler_refined.png" \
  "${PROJECT_ROOT}/analysis/singler/${control[2]}/${control[2]}_umap_singler_refined.png"

panel_3_row "${out_dir}/panel_manual_umap_disease.png" auto auto "${gap}" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[0]}/ManualAnnotation/umap_manual_celltype.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[1]}/ManualAnnotation/umap_manual_celltype.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[2]}/ManualAnnotation/umap_manual_celltype.png"
panel_3_row "${out_dir}/panel_manual_umap_control.png" auto auto "${gap}" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[0]}/ManualAnnotation/umap_manual_celltype.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[1]}/ManualAnnotation/umap_manual_celltype.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[2]}/ManualAnnotation/umap_manual_celltype.png"

panel_3_col "${out_dir}/panel_featureplot_lineage_disease.png" auto auto "${gap_zero}" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[0]}/ManualAnnotation/featureplot_lineage_panels.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[1]}/ManualAnnotation/featureplot_lineage_panels.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[2]}/ManualAnnotation/featureplot_lineage_panels.png"
panel_3_col "${out_dir}/panel_featureplot_lineage_control.png" auto auto "${gap_zero}" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[0]}/ManualAnnotation/featureplot_lineage_panels.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[1]}/ManualAnnotation/featureplot_lineage_panels.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[2]}/ManualAnnotation/featureplot_lineage_panels.png"

panel_3_col "${out_dir}/panel_featureplot_immune_disease.png" auto auto "${gap_zero}" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[0]}/ManualAnnotation/featureplot_immune_panels.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[1]}/ManualAnnotation/featureplot_immune_panels.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[2]}/ManualAnnotation/featureplot_immune_panels.png"
panel_3_col "${out_dir}/panel_featureplot_immune_control.png" auto auto "${gap_zero}" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[0]}/ManualAnnotation/featureplot_immune_panels.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[1]}/ManualAnnotation/featureplot_immune_panels.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[2]}/ManualAnnotation/featureplot_immune_panels.png"

panel_3_col "${out_dir}/panel_featureplot_cycle_disease.png" auto auto "${gap_zero}" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[0]}/ManualAnnotation/featureplot_cycle_panels.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[1]}/ManualAnnotation/featureplot_cycle_panels.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[2]}/ManualAnnotation/featureplot_cycle_panels.png"
panel_3_col "${out_dir}/panel_featureplot_cycle_control.png" auto auto "${gap_zero}" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[0]}/ManualAnnotation/featureplot_cycle_panels.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[1]}/ManualAnnotation/featureplot_cycle_panels.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[2]}/ManualAnnotation/featureplot_cycle_panels.png"

panel_3_col "${out_dir}/panel_dotplot_markers_disease.png" auto auto "${gap_zero}" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[0]}/ManualAnnotation/dotplot_marker_panels.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[1]}/ManualAnnotation/dotplot_marker_panels.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[2]}/ManualAnnotation/dotplot_marker_panels.png"
panel_3_col "${out_dir}/panel_dotplot_markers_control.png" auto auto "${gap_zero}" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[0]}/ManualAnnotation/dotplot_marker_panels.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[1]}/ManualAnnotation/dotplot_marker_panels.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[2]}/ManualAnnotation/dotplot_marker_panels.png"

panel_3_col "${out_dir}/panel_heatmap_top_markers_disease.png" auto auto "${gap}" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[0]}/HeatmapEvidence/heatmap_top_markers_per_cluster.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[1]}/HeatmapEvidence/heatmap_top_markers_per_cluster.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[2]}/HeatmapEvidence/heatmap_top_markers_per_cluster.png"
panel_3_col "${out_dir}/panel_heatmap_top_markers_control.png" auto auto "${gap}" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[0]}/HeatmapEvidence/heatmap_top_markers_per_cluster.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[1]}/HeatmapEvidence/heatmap_top_markers_per_cluster.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[2]}/HeatmapEvidence/heatmap_top_markers_per_cluster.png"

panel_3_col "${out_dir}/panel_heatmap_marker_panel_disease.png" auto auto "${gap}" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[0]}/HeatmapEvidence/heatmap_celltype_marker_panel.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[1]}/HeatmapEvidence/heatmap_celltype_marker_panel.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${disease[2]}/HeatmapEvidence/heatmap_celltype_marker_panel.png"
panel_3_col "${out_dir}/panel_heatmap_marker_panel_control.png" auto auto "${gap}" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[0]}/HeatmapEvidence/heatmap_celltype_marker_panel.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[1]}/HeatmapEvidence/heatmap_celltype_marker_panel.png" \
  "${PROJECT_ROOT}/analysis/manual_annotation/${control[2]}/HeatmapEvidence/heatmap_celltype_marker_panel.png"

echo "[INFO] Panels saved to ${out_dir}"
