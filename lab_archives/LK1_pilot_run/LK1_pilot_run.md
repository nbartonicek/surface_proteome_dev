# LK1 pilot run

**Run:** `260423_VH01624_453_222HWMYNX`
**Dates of work:** 2026-04-28 to 2026-05-20 (a few follow-ups up to 2026-06-15)
**Status:** superseded - LK1 was re-sequenced as `260522_VH01624_461_222JLJVNX`. Everything
below is the first pass on the pilot data and should not be quoted as a current result.
**Written:** 2026-07-31, retroactively, from the scripts and the output files still on disk.

This is the first run of the project and most of the pipeline was worked out on it. A lot of
the scripts that produced these outputs have since been renumbered or replaced, so where the
version that actually ran is now in `backup/`, that is the copy frozen alongside this report.

---

## 1. Run, libraries and samples

Three libraries off one GEM well:

| Library | FASTQ prefix | R1 reads |
|---|---|---|
| GEX | `LK1-GEX_S5_R1_001.fastq.gz` | 475,071,355 |
| ADT (TotalSeq-A) | `LK1-TotalSeqA-ADT_S6_R1_001.fastq.gz` | 348,013,031 |
| HTO (TotalSeq-A) | `LK1-TotalSeqA-HTO_S7_R1_001.fastq.gz` | 15,864,696 |

Read counts are from `0.count_reads.sh` (line count / 4 on R1), written to
`raw/260423_VH01624_453_222HWMYNX/fastq_read_counts.tsv`. The HTO library is ~30x shallower
than ADT, which is the expected ratio but worth remembering when reading the demux figures.

Four samples hashed into the well, from the `hto_names` map in `6a.initial_QC.R`:

| HTO barcode | Sample |
|---|---|
| `HTO1-GTCAACTCTTTAGCG` | MOLM13 (AML cell line) |
| `HTO2-TGATGGCCTATTGGG` | HBDN206-MNpCT |
| `HTO3-TTCCGCCTCTCTTTG` | HBDN392-AML-MDS |
| `HTO4-AGTAAGTTCAGCGTA` | HBDN501-AML-KMT2A |

Note the naming drift: the HTO map writes them with underscores (`HBDN206_MNpCT`), the
figures and later metadata use hyphens (`HBDN206-MNpCT`). Same samples.

## 2. Where the files are

**Raw (as written in the scripts, on the cluster):**
`/scratch/users/nbartonicek/projects/amgen/raw/260423_VH01624_453_222HWMYNX/` with
`Sample_LK1-GEX/`, `Sample_LK1-TotalSeqA-ADT/`, `Sample_LK1-TotalSeqA-HTO/`.

On `bioinf_scratch` that directory now holds only `fastq_read_counts.tsv` - the pilot FASTQs
are no longer there. They are either still on cluster `/scratch` or have been cleaned up;
this was not checked when writing this report.

**Reference:** `annotation/refdata-gex-GRCh38-2024-A`

**Outputs**, all under `results/<dir>/260423_VH01624_453_222HWMYNX/`:

| Directory | Dates | What is in it |
|---|---|---|
| `cite_seq_count/` | 04-28 to 05-05 | ADT and HTO count matrices, three versions (raw, emptydrops-whitelisted, `_raw`) |
| `seurat_demux/` | 04-28 to 05-03 | first Seurat object, HTO demux, ADT QC, BoneMarrowMap projection |
| `emptydrops/` | 05-01 | called barcodes and the HTO/scDblFinder-filtered object |
| `scDblFinder/` | 05-03 | doublet scores and the filtered/nonfiltered objects |
| `genotype_demux/` | 05-03 to 05-04 | cellSNP pileup and vireo donor calls |
| `demux_comparison/` | 05-04 to 05-05 | HTO vs vireo agreement, the parameter sweep |
| `seurat_annotated/` | 05-05 to 06-15 | annotated object, projection figures, copykat, numbat |
| `data_integration/` | 05-06 to 05-20 | DSB vs CLR comparison, composition by sample |

## 3. What was run, in order

### 3.1 Read counting - `0.count_reads.sh`
SLURM array (3 tasks, one per library), `rhel_short`, 2 h, 4 GB. Counts R1 lines and divides
by 4. Trivial but it is the denominator for everything downstream.

### 3.2 Cell Ranger
**The pilot Cell Ranger invocation is not preserved.** `1.cellranger.sh` at the top level now
points at `260528`/`LK2-GEX`, and `1a.cellranger_resequence.sh` is the re-sequencing version
(which needs Cell Ranger 10.0.0). `6a.initial_QC.R` reads from
`../results/cellranger/260423_VH01624_453_222HWMYNX/LK1-GEX/outs/filtered_feature_bc_matrix`,
and **that directory does not exist on the mount** - only `results/cellranger_withbam/` does,
and it has no `260423` subdirectory. So the pilot count matrices have either been removed or
live on cluster `/scratch`. The commented-out module line in `1.cellranger.sh` is
`cellranger/7.2.0-gcc-13.2.0`, which is the best available guess at the version used, but it
is a guess.

### 3.3 CITE-seq-Count (ADT and HTO)
Run three times against different whitelists, which is why `cite_seq_count/260423.../` has
three sets of counts:

- `adt_counts_LK1` / `hto_counts_LK1` (04-28) - first pass
- `adt_counts_LK1_emptydrops` / `hto_counts_LK1_emptydrops` (05-01) - re-run against the
  emptyDrops barcode list once that existed
- `adt_counts_LK1_raw` (05-05) - raw, unwhitelisted, needed for DSB background estimation

Each has `read_count/`, `umi_count/`, `uncorrected_cells/dense_umis.tsv`, `unmapped.csv` and
`run_report.yaml`. The `unmapped.csv` files are worth keeping - the ADT unmapped-read problem
was only properly chased down much later (`7g.CITE_unmapped_QC.R`, July).

### 3.4 Initial QC and HTO demultiplexing - `6a.initial_QC.R`
Also `backup_4.initial_QC.R`, `backup_6.demux.R` and `backup_8.create_seurat_object.R`, which
are the earlier splits of the same work.

Parameters actually used:

- `CreateSeuratObject(min.cells = 3, min.features = 50)`
- `percent.mt` from `^MT-`, `percent.ribo` from `^RP[SL]`
- HTO normalised CLR, `margin = 2`
- `HTODemux(positive.quantile = 0.80)`
- RNA: `ScaleData(vars.to.regress = "percent.mt")`, PCA, `dims = 1:30`,
  `FindClusters(resolution = 0.5)`, UMAP on 1:30
- ADT normalised CLR, `margin = 2`

The 0.80 positive quantile is the value that later got swept - `demux_comparison/` carries
`_all_positive_quantiles` and `_positive_quantile_0.99` variants of every figure.

Barcode handling worth noting: the HTO and ADT barcode lists come out of CITE-seq-Count
without the GEM-well suffix, so the script does `paste0(hto_barcodes, "-1")` before matching
to the Cell Ranger barcodes. The same class of problem reappeared much later with long-read
BAMs in the mitochondrial work.

### 3.5 emptyDrops and scDblFinder
`emptydrops/260423.../LK1-GEX/` (05-01) and `scDblFinder/260423.../LK1-GEX/` (05-03). The
emptyDrops barcode list (`emptydrops_barcodes.txt`) is what the second CITE-seq-Count run was
whitelisted against. Both a filtered and a nonfiltered Seurat object were kept from
scDblFinder, which is the right call and was not always done later.

### 3.6 Genotype demultiplexing and the HTO comparison
cellSNP pileup (05-03) then vireo (05-04), both under `genotype_demux/260423.../LK1-GEX/`.
`demux_comparison/` (05-05) is the first version of the HTO-vs-vireo agreement analysis that
later became its own benchmarking stream. Every figure exists in a default version and an
`_all_positive_quantiles` version.

### 3.7 CITE-seq QC and DSB - `backup_7a.CITE_qc.R`, `7b.CITE_vs_DSB_comparison.R`, `7c.evaluate_DSB.R`
DSB needs background droplets, which is what the `_raw` CITE-seq-Count run was for. `7c` is
the DSB evaluation on this run; `7b` is the head-to-head against CLR on HBDN206-MNpCT.

### 3.8 Annotation - `backup_10.annotate.R`, `backup_10a.annotate.R`, `backup_10b.annotate_DSB.R`
BoneMarrowMap Symphony projection. Two objects came out of this:
`LK1_GEX_HTO_ADT_BoneMarrowMap_projected.rds` (05-03, ADT/CLR) and
`LK1_GEX_CITE_DSB_BoneMarrowMap_projected.rds` (05-06, DSB), plus
`demux_singlets_annotated_seurat.rds` which is the one later scripts read.

### 3.9 DSB vs CLR and composition - `12.data_integration.R`, `backup_10c.compare_DSB_CLR.R`
This produced `data_integration/260423.../` and the marker-specificity CSVs.
**Caveat:** the frozen `12.data_integration.R` was last edited 2026-06-09, but these outputs
are dated 05-06 to 05-20. The figures were made by an earlier version of the script, which is
not recoverable (`scripts/` was not under git until 2026-07-31).

### 3.10 CNV calling - `9.copykat.R`, `9a.copykat_perSample.R`
`9.copykat.R` ran CopyKAT pooled (05-08, one `..._CNA_raw_results_gene_by_cell.txt` at the
top of `copykat_annotated/`), then `9a.copykat_perSample.R` re-ran it per donor (05-12), which
is what the three donor subdirectories are. Numbat calls were integrated into the object on
05-08 (`seurat_annotated/260423.../numbat/`), though the Numbat run scripts themselves are
dated June and belong to the later benchmarking report.

### 3.11 Cancer cell identification and surfaceome - `12.identify_cancer_cells.R`, `14.surfaceome_trial.R`
First attempt at calling malignant cells from the CNV output, and the first surfaceome pass.
`14.surfaceome_trial.R` is dated 05-25, right at the point the re-sequenced run arrived, so it
is the boundary of the pilot work.

### 3.12 One-off - `add_cite_for_hiep.R`
Added the CITE data to an object for a colleague (05-08). Included for completeness, no
analytical result.

---

## 4. Key figures

Paths are relative to `results/`. Dates are file modification times.

### Demultiplexing and first-pass QC - `seurat_demux/260423_VH01624_453_222HWMYNX/`

| Figure | Date | What it shows |
|---|---|---|
| `01_initial_qc_violin.pdf` | 04-28 | nFeature_RNA, nCount_RNA and percent.mt split by `hash.ID`. The standard first look - use it to check no hashed sample is obviously degraded relative to the others. |
| `02_hto_ridgeplots.pdf` | 04-28 | CLR-normalised HTO signal per hashtag. This is the figure that tells you whether hashing worked at all; clean bimodality per HTO is what you want. |
| `03_umap_demux.pdf` | 04-28 | Three UMAPs side by side: HTO global classification (singlet/doublet/negative), assigned `hash.ID`, and top HTO by raw count. The third panel is the sanity check on the second. |
| `04_umap_qc_features.pdf` | 04-28 | UMI counts, detected genes and mitochondrial % on the UMAP. Localised high-MT regions are the dying-cell clusters. |
| `05_hto_featureplots.pdf` | 04-28 | Each HTO on the UMAP. Confirms the hashed samples occupy distinct territory rather than smearing. |
| `06_ADT_QC.pdf` | 04-29 | ADT totals and features per cell. |
| `07_top_ADT_markers.pdf` | 04-29 | Feature plots of the top 5 ADT markers by mean CLR plus biotin. Biotin was added deliberately here and is the start of the whole biotin thread that runs through to July. |
| `08_top_ADT_cluster_markers_dotplot.pdf` | 04-29 | Top ADT markers per RNA cluster. |
| `10_top10_discriminative_ADT_plus_biotin_UMAP.pdf` | 04-29 | Most discriminative ADTs with biotin overlaid. |
| `11_biotin_expression_by_cluster_violin.pdf` | 04-29 | Biotin by cluster - the first look at the biotin heterogeneity later split high vs low. |
| `09_mapping_error_QC.pdf` | 05-03 | BoneMarrowMap mapping error. The QC that decides which cells get an annotation at all. |
| `10_projected_predicted_celltypes_pass_only.pdf` | 05-03 | Predicted cell types on the reference UMAP, cells passing mapping QC only. |
| `11_projected_pseudotime_pass_only.pdf` | 05-03 | Projected pseudotime. |
| `12_projection_by_sampleID.pdf` | 05-03 | Projection split by sample - where each donor sits on the reference. |
| `13.broad_celltype_composition.pdf`, `14.specific_celltype_composition.pdf` | 05-03 | Composition at both annotation levels. |
| `15.heatmap_cell_type_broad.pdf`, `15a.heatmap_cell_type_all.pdf` | 05-03 | Cell type heatmaps. |

### HTO vs vireo - `demux_comparison/260423_VH01624_453_222HWMYNX/LK1-GEX/`
All dated 05-05. Each exists in a default version and an `_all_positive_quantiles` version;
the sweep is the point of the figure set.

| Figure | What it shows |
|---|---|
| `01_demux_efficiency_barplot.pdf` | Singlet/doublet/negative rates for HTO and vireo side by side. The headline comparison. |
| `02_hto_vireo_doublet_comparison.pdf` | Do the two methods call the same cells doublets. |
| `03_hto_vireo_donor_heatmap.pdf` | Cross-tabulation of HTO sample against vireo donor. Off-diagonal mass is the disagreement. |
| `04_upset_HTO_vireo_overlap.pdf` | UpSet of the call overlap. |
| `05_HTO_QC.pdf`, `05_HTO_QC_threshold_independent.pdf` | HTO QC with and without a threshold applied. |
| `06_HTO_ridge.pdf`, `06b_HTO_ratio_thresholds_threshold_independent.pdf` | Ridge plots and the ratio thresholds. |
| `07_UMAP_Vireo_and_HTO_all_positive_quantiles.pdf` | Both calls on the UMAP across the quantile sweep. |

### Annotation and projection - `seurat_annotated/260423_VH01624_453_222HWMYNX/projectionFigures/`

| Figure | Date | What it shows |
|---|---|---|
| `01_mapping_error_QC.pdf` | 05-06 | Mapping error distribution, re-done on the DSB object. |
| `02_projected_predicted_celltypes_pass_only.pdf` | 05-06 | Predicted cell types, passing cells. |
| `03_projected_samples_pass_only.pdf` | 05-06 | Samples on the projected UMAP. |
| `04_projected_pseudotime_pass_only.pdf` | 05-06 | Pseudotime. |
| `05_projection_by_sampleID_BMM_helper_manual_split.pdf` | 05-06 | Per-sample projection using the BoneMarrowMap helper. |
| `06_ADT_CLR_top_marker_projection.pdf` / `06_CITE_DSB_top_marker_projection.pdf` | 05-06 | The same top markers under CLR and under DSB. **This pair is the CLR-vs-DSB comparison in its most direct form.** |
| `07_CITE_DSB_selected_markers_manual_projectedUMAP.pdf` | 05-06 | Hand-picked markers on the projected UMAP. |
| `density_<sample>_projectedUMAP.pdf` (4 files) | 05-05 | Per-sample cell density on the reference: MOLM13, HBDN206-MNpCT, HBDN392-AML-MDS, HBDN501-AML-KMT2A. MOLM13 should be tight and restricted - it is the cell line, so it is the positive control for the projection behaving sensibly. |
| `DSB_vs_CLR_HBDN206-MNpCT/01_marker_mean_shift_DSB_minus_CLR.pdf` | 05-20 | Per-marker mean shift, DSB minus CLR, for one donor. Quantifies what the eyeball comparison above shows. |

### DSB quality - `data_integration/260423_VH01624_453_222HWMYNX/`

| Figure | Date | What it shows |
|---|---|---|
| `01_background_suppression.pdf` | 05-20 | How much background DSB removes. The main argument for using DSB. |
| `02_expected_marker_specificity.pdf`, `02_expected_marker_percent_specificity.pdf` | 05-20 | Whether expected markers land on the expected cell types. |
| `03_marker_leakage_diffusion.pdf` | 05-20 | Marker signal leaking into populations that should be negative. |
| `04_expected_vs_other_violin.pdf` | 05-20 | Expected vs other cell types per marker. |
| `05_featureplots_CLR_vs_DSB_expected_markers.pdf` | 05-20 | Side-by-side feature plots, the visual summary of the whole comparison. |
| `07_broad_celltype_composition_by_sample.pdf`, `08_detailed_celltype_composition_by_sample.pdf` | 05-06 | Composition by sample at both levels. |
| `09_DSB_UMAP.pdf`, `09_DSB_UMAP_by_sample.pdf`, `10_DSB_UMAP_by_sample_variable_proteins.pdf` | 05-06 | UMAP computed on the DSB protein space. |
| `11_copykat_by_sample_and_broad_annotation.pdf` | 05-08 | CopyKAT aneuploid/diploid call against sample and broad cell type. |

Accompanying CSVs (`01_background_suppression.csv`, `02_expected_marker_specificity.csv`,
`03_marker_leakage_diffusion.csv`, `06_summary_metrics.csv`) carry the numbers behind these.

### CNV - `seurat_annotated/260423_VH01624_453_222HWMYNX/`

Per donor under `copykat_annotated/<donor>/` (05-12), for HBDN206-MNpCT, HBDN392-AML-MDS and
HBDN501-AML-KMT2A:

| Figure | What it shows |
|---|---|
| `01_copykat_call_composition_by_sampleID.pdf` | Aneuploid/diploid proportions. |
| `02_copykat_calls_on_projected_RNA_UMAP.pdf`, `02a_...` | Calls on the projected UMAP. |
| `03_copykat_by_BoneMarrowMap_broad_celltype.pdf`, `03a_...counts...` | Calls against annotation, proportions and counts. The check that "aneuploid" is not just landing on one cell type for technical reasons. |
| `LK1_copykat_copykat_heatmap.jpeg` | The CopyKAT CNV heatmap itself. |

Numbat integration under `numbat/` (05-08), 15 figures. The informative ones are
`09_numbat_compartment_on_umap.pdf` (tumour/normal compartment),
`10b_numbat_clone_on_umap_split_by_sample.pdf` (clone assignment per sample) and
`11_numbat_compartment_by_sample_and_broad_annotation.pdf` (compartment against annotation).

`numbat/CITE_HSC_Numbat_analysis/` (06-15) is later work (`7e.1.find_cancer_CITE_markers_run1.R`)
that used the pilot object after the re-sequenced run had already arrived. It is listed here
because it lives under this run's directory, but it belongs with the CITE benchmarking report:
HSC marker presence, the cancer-vs-normal-by-Numbat volcano, validation in normal-sample HSCs,
and the BAFF-R RNA/protein correlation.

---

## 5. Tools

Cluster side, from the scripts:

| Tool | Version | Source |
|---|---|---|
| Cell Ranger | 7.2.0 (uncertain, see 3.2) | commented module line in `1.cellranger.sh` |
| Reference | `refdata-gex-GRCh38-2024-A` | `annotation/` |
| CITE-seq-Count | not recorded | `3.hto_count.sbatch`, `4.cite_seq_count.sbatch` |
| cellsnp-lite, vireo | not recorded | `3.cellSNP.sh` |

R side. **No `sessionInfo()` was captured at the time**, and `scripts/` was not under version
control until 2026-07-31, so the exact versions used in April/May are not recoverable. The
analysis ran on the laptop against the `bioinf_scratch` mount, and the current laptop session
is:

```
R version 4.5.1 (2025-06-13), aarch64-apple-darwin20, macOS Tahoe 26.5.2
Seurat 5.4.0, SeuratObject 5.3.0, Matrix 1.7-4
harmony 1.2.4, scIntegrationMetrics 1.2.0
edgeR 4.6.3, limma 3.64.3, DESeq2 1.48.2
AUCell 1.30.1, UCell 2.12.0, fgsea 1.34.2, msigdbr 25.1.1, decoupleR 2.14.0
rhdf5 2.52.1, data.table 1.18.2.1, tidyverse 2.0.0, patchwork 1.3.2
```

Treat that as "the environment as it stands now", not as the environment that produced these
figures. R 4.5.1 was released 2025-06-13 so it is plausible it was in use in April 2026, but
the package versions will have moved.

CopyKAT and dsb are not in the list above because they were not attached in that session;
both were used here.

---

## 6. Caveats

1. **The run is superseded.** LK1 was re-sequenced as `260522_VH01624_461_222JLJVNX`. Nothing
   here should be presented as a current result.
2. **The pilot Cell Ranger output is not on the mount.** `6a.initial_QC.R` reads
   `../results/cellranger/260423.../LK1-GEX/outs/filtered_feature_bc_matrix` and that path does
   not exist. Check cluster `/scratch` before assuming it was deleted.
3. **The pilot FASTQs are not on the mount either** - `raw/260423.../` holds only the read
   count table.
4. **Dates are file modification times.** `scripts/` had no git history before 2026-07-31, so
   there is no commit trail for this period. The output dates are reliable; the script dates
   only tell you when a file was last edited.
5. **`12.data_integration.R` post-dates its own output** by about three weeks. The
   `data_integration/260423.../` figures were made by a version of that script that no longer
   exists.
6. **`13.CITE_residuals.R` writes to `results/residual_cite_analysis`, which does not exist.**
   The residual CITE analysis is dated 05-06 and belongs to this period, but no output
   survives on the mount. It is not included in the frozen scripts for that reason - see the
   CITE benchmarking report instead.
7. **`numbat/CITE_HSC_Numbat_analysis/` is not pilot-era** (06-15) even though it sits under
   this run.

## 7. Frozen scripts

`scripts/` next to this report holds the versions as they stand on 2026-07-31. Files prefixed
`backup_` came from `scripts/backup/` and are the ones that actually ran during the pilot;
the un-prefixed ones are the current top-level versions, which in some cases have been edited
since (see caveat 5).

| File | Section |
|---|---|
| `0.count_reads.sh` | 3.1 |
| `6a.initial_QC.R`, `backup_4.initial_QC.R`, `backup_6.demux.R`, `backup_8.create_seurat_object.R` | 3.4 |
| `backup_7a.CITE_qc.R`, `7b.CITE_vs_DSB_comparison.R`, `7c.evaluate_DSB.R` | 3.7 |
| `backup_10.annotate.R`, `backup_10a.annotate.R`, `backup_10b.annotate_DSB.R` | 3.8 |
| `12.data_integration.R`, `backup_10c.compare_DSB_CLR.R`, `backup_11.combine_cite.R` | 3.9 |
| `9.copykat.R`, `9a.copykat_perSample.R`, `backup_11a.example_cite_copykat.R` | 3.10 |
| `12.identify_cancer_cells.R`, `14.surfaceome_trial.R` | 3.11 |
| `add_cite_for_hiep.R` | 3.12 |
