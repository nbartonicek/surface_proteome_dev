---
title: "LK1 pilot run"
subtitle: "260423_VH01624_453_222HWMYNX"
date: "2026-07-31"
output:
  github_document:
    toc: true
    html_preview: false
params:
  run_id: "260423_VH01624_453_222HWMYNX"
  proj: "/Volumes/bioinf_scratch/users/nbartonicek/projects/amgen"
  sample_short: "LK1"
  sample_name: "LK1-GEX"
---

<!--
# ------------------------------------------------------------------
# LK1 pilot run - step 15 of 15
#
# Renders the lab archive report for the pilot run from the outputs of
# steps 01-14. Everything that can be read off disk is read off disk, so
# the report cannot drift from the results: tables come from the summary
# CSVs each step wrote, the figure index is built by listing the PDFs that
# are actually there, and the timeline is built from file modification
# times. Prose that cannot be derived (why a thing was done, what a figure
# means) is written out below and marked as such.
#
# Missing inputs are skipped with a note rather than erroring - the pilot
# predates the nextflow pipeline and some of its outputs have since been
# cleaned up, so a partial render is the expected case, not a failure.
#
# Render with (needs pandoc, e.g. from RStudio):
#   Rscript -e 'rmarkdown::render("15.LK1_pilot_run_report.Rmd",
#                                 output_file = "LK1_pilot_run.md")'
#
# Or without pandoc - the output is markdown anyway, so knit is enough:
#   Rscript -e 'knitr::knit("15.LK1_pilot_run_report.Rmd",
#                           output = "LK1_pilot_run.md")'
# knit() does not process the YAML header, so set params by hand first if you
# need anything other than the defaults above.
#
# Written 2026-07-31 for the retroactive lab archive.
# ------------------------------------------------------------------
-->



# LK1 pilot run

**Run:** `` 260423_VH01624_453_222HWMYNX ``
**Status:** superseded - LK1 was re-sequenced as `260522_VH01624_461_222JLJVNX`. Everything
below is the first pass on the pilot data and should not be quoted as a current result.
**Report generated:** 2026-07-31 by `15.LK1_pilot_run_report.Rmd` from the outputs of
steps 01-14 in this folder.

Found **1266 output files** for this run across **8 results directories**, dated 2026-04-28 to 2026-06-15.

This is the first run of the project and most of the pipeline was worked out on it. Several
of the scripts that produced these outputs have since been renumbered or replaced; the copies
in this folder are the versions that ran, renumbered into execution order.

# 1. Run, libraries and samples


Table: R1 read counts per library (from 01.count_reads.sh)

|fastq                                |reads       |library          |
|:------------------------------------|:-----------|:----------------|
|LK1-GEX_S5_R1_001.fastq.gz           |475,071,355 |GEX              |
|LK1-TotalSeqA-ADT_S6_R1_001.fastq.gz |348,013,031 |ADT (TotalSeq-A) |
|LK1-TotalSeqA-HTO_S7_R1_001.fastq.gz |15,864,696  |HTO (TotalSeq-A) |

Read counts are line count / 4 on R1. The HTO library is far shallower than ADT, which is the
expected ratio but worth remembering when reading the demux figures.

Four samples were hashed into the well. The map is hardcoded in
`02.initial_QC_and_HTO_demux.R`:

| HTO barcode | Sample |
|---|---|
| `HTO1-GTCAACTCTTTAGCG` | MOLM13 (AML cell line) |
| `HTO2-TGATGGCCTATTGGG` | HBDN206-MNpCT |
| `HTO3-TTCCGCCTCTCTTTG` | HBDN392-AML-MDS |
| `HTO4-AGTAAGTTCAGCGTA` | HBDN501-AML-KMT2A |

Note the naming drift: the HTO map writes them with underscores (`HBDN206_MNpCT`), the
figures and later metadata use hyphens (`HBDN206-MNpCT`). Same samples.


Table: HTO demultiplexing outcome (demux_summary.csv)

|HTO_classification.global |hash.ID           | n_cells|
|:-------------------------|:-----------------|-------:|
|Doublet                   |Doublet           |    5813|
|Negative                  |Negative          |    3455|
|Singlet                   |HBDN392-AML-MDS   |    3281|
|Singlet                   |HBDN206-MNpCT     |    6813|
|Singlet                   |HBDN501-AML-KMT2A |    1731|
|Singlet                   |MOLM13            |    1996|




Table: Per-sample QC medians (sample_qc_summary.csv)

|hash.ID           |HTO_classification.global | n_cells| median_nCount_RNA| median_nFeature_RNA| median_percent_mt| median_percent_ribo| median_nCount_HTO| median_nFeature_HTO|
|:-----------------|:-------------------------|-------:|-----------------:|-------------------:|-----------------:|-------------------:|-----------------:|-------------------:|
|Doublet           |Doublet                   |    5813|              4723|                2225|          2.994497|            21.36240|                72|                   3|
|HBDN392-AML-MDS   |Singlet                   |    3281|              3622|                1816|          1.885681|            17.17778|                50|                   2|
|HBDN206-MNpCT     |Singlet                   |    6813|              3928|                2013|          2.598923|            20.77890|                42|                   2|
|HBDN501-AML-KMT2A |Singlet                   |    1731|              2468|                1352|          3.079803|            19.38263|                49|                   3|
|MOLM13            |Singlet                   |    1996|             20793|                5092|          3.374237|            25.66427|                60|                   3|
|Negative          |Negative                  |    3455|              4179|                2017|          2.639442|            20.50888|                 8|                   2|

# 2. Where the files are

**Raw**, as written in the scripts (cluster paths):
`/scratch/users/nbartonicek/projects/amgen/raw/`` 260423_VH01624_453_222HWMYNX ``/` with `Sample_LK1-GEX/`,
`Sample_LK1-TotalSeqA-ADT/`, `Sample_LK1-TotalSeqA-HTO/`.

On `bioinf_scratch`, `raw/260423_VH01624_453_222HWMYNX/` currently holds: `fastq_read_counts.tsv`.

**Reference:** `annotation/refdata-gex-GRCh38-2024-A`


Table: Output directories under results/<dir>/260423_VH01624_453_222HWMYNX/

|dir              |first      |last       | files|
|:----------------|:----------|:----------|-----:|
|cite_seq_count   |2026-04-28 |2026-05-05 |    48|
|seurat_demux     |2026-04-28 |2026-05-11 |   453|
|emptydrops       |2026-05-01 |2026-05-01 |     4|
|genotype_demux   |2026-05-03 |2026-05-04 |    12|
|scDblFinder      |2026-05-03 |2026-05-03 |     5|
|demux_comparison |2026-05-04 |2026-05-05 |    32|
|seurat_annotated |2026-05-05 |2026-06-15 |   195|
|data_integration |2026-05-06 |2026-06-01 |   517|

# 3. Timeline

Built from file modification times, so it is when outputs landed, not when jobs started.


Table: What was written, by day

|date       |directories                                                      | files|
|:----------|:----------------------------------------------------------------|-----:|
|2026-04-28 |cite_seq_count, seurat_demux                                     |    27|
|2026-04-29 |seurat_demux                                                     |     7|
|2026-05-01 |cite_seq_count, emptydrops                                       |    22|
|2026-05-03 |cite_seq_count, genotype_demux, scDblFinder, seurat_demux        |    47|
|2026-05-04 |cite_seq_count, demux_comparison, genotype_demux, seurat_demux   |    47|
|2026-05-05 |cite_seq_count, demux_comparison, seurat_annotated, seurat_demux |   396|
|2026-05-06 |data_integration, seurat_annotated, seurat_demux                 |   131|
|2026-05-07 |data_integration                                                 |    12|
|2026-05-08 |data_integration, seurat_annotated                               |    33|
|2026-05-11 |seurat_demux                                                     |    17|
|2026-05-12 |seurat_annotated                                                 |    37|
|2026-05-19 |seurat_annotated                                                 |     1|
|2026-05-20 |data_integration, seurat_annotated                               |    26|
|2026-05-21 |data_integration, seurat_annotated                               |    76|
|2026-05-25 |data_integration                                                 |    38|
|2026-05-27 |data_integration                                                 |   207|
|2026-06-01 |data_integration                                                 |   128|
|2026-06-09 |seurat_annotated                                                 |     2|
|2026-06-15 |seurat_annotated                                                 |    12|

# 4. What was run, in order

The scripts in this folder are numbered in execution order. Three steps that ran on the
pilot have **no surviving script**: Cell Ranger, emptyDrops, and cellSNP/vireo. The current
top-level versions of those all point at later runs, so they are not included here.

| Step | Script | What it does |
|---|---|---|
| - | *(not preserved)* | Cell Ranger count. `02` reads `../results/cellranger/<run>/LK1-GEX/outs/filtered_feature_bc_matrix`. |
| 01 | `01.count_reads.sh` | R1 read counts per library. |
| 02 | `02.initial_QC_and_HTO_demux.R` | Seurat object, HTO + ADT attached, CLR, `HTODemux(positive.quantile = 0.80)`. |
| - | *(not preserved)* | emptyDrops. Its barcode list is what CITE-seq-Count was re-run against. |
| - | *(not preserved)* | cellSNP pileup and vireo donor calling. |
| 03 | `03.demux_vireo_and_doublets.R` | Vireo calls onto the object, HTO cross-tab, doublet detection. |
| 04 | `04.CITE_qc.R` | CITE QC on emptyDrops cells: raw ADT, CLR, marker distributions, overlays. |
| 05 | `05.evaluate_DSB.R` | DSB using raw background droplets, with IgG-control subsampling. |
| 06 | `06.combine_cite.R` | Folds CITE assays into one object, writes the `analysis_bundle`. |
| 07 | `07.annotate_reference_setup.R` | BoneMarrowMap Symphony reference, uwot path, first projection. |
| 08 | `08.annotate_ADT_CLR.R` | Projection on the ADT/CLR object. |
| 09 | `09.annotate_CITE_DSB.R` | Same projection on the CITE_DSB object. |
| 10 | `10.compare_DSB_vs_CLR_one_sample.R` | DSB vs CLR head-to-head, HBDN206-MNpCT. |
| 11 | `11.data_integration_DSB_vs_CLR.R` | Background suppression, marker specificity, leakage, composition. |
| 12 | `12.copykat_pooled.R` | CopyKAT, all cells pooled. |
| 13 | `13.copykat_per_sample.R` | CopyKAT per donor. |
| 14 | `14.surfaceome_trial.R` | First surfaceome pass. Exploratory, saves nothing. |
| 15 | `15.LK1_pilot_run_report.Rmd` | This report. |

Key parameters actually used in step 02: `CreateSeuratObject(min.cells = 3, min.features = 50)`,
`percent.mt` from `^MT-`, HTO and ADT normalised CLR with `margin = 2`,
`HTODemux(positive.quantile = 0.80)`, RNA scaled regressing `percent.mt`, PCA and UMAP on
`dims = 1:30`, `FindClusters(resolution = 0.5)`.

The 0.80 positive quantile is the value that later got swept - `demux_comparison/` carries
`_all_positive_quantiles` and `_positive_quantile_0.99` variants of every figure.

# 5. Results

## 5.1 Demultiplexing: HTO vs vireo


Table: Demux efficiency, HTO vs vireo

|method        | total_cells| singlets| doublets| negatives_or_unassigned| singlet_fraction| doublet_fraction| negative_or_unassigned_fraction|
|:-------------|-----------:|--------:|--------:|-----------------------:|----------------:|----------------:|-------------------------------:|
|HTO           |       23520|    11457|      778|                   11285|        0.4871173|        0.0330782|                       0.4798044|
|cellSNP_vireo |       23520|    22405|     1115|                       0|        0.9525935|        0.0474065|                       0.0000000|




Table: HTO / vireo agreement

| total_cells| both_singlet| HTO_doublet_only| Vireo_doublet_only| both_doublet| HTO_singlet_Vireo_unassigned| HTO_negative_Vireo_singlet| both_singlet_fraction| both_doublet_fraction|
|-----------:|------------:|----------------:|------------------:|------------:|----------------------------:|--------------------------:|---------------------:|---------------------:|
|       23520|        10982|              538|                875|          240|                            0|                      10885|             0.4669218|             0.0102041|




Table: HTO sample against vireo donor

|hto_donor         |vireo_donor | n_cells|
|:-----------------|:-----------|-------:|
|HBDN392-AML-MDS   |donor3      |    2974|
|HBDN392-AML-MDS   |unassigned  |     108|
|HBDN392-AML-MDS   |donor1      |      85|
|HBDN392-AML-MDS   |donor0      |      18|
|HBDN392-AML-MDS   |donor2      |       9|
|HBDN501-AML-KMT2A |donor2      |    1450|
|HBDN501-AML-KMT2A |donor1      |     106|
|HBDN501-AML-KMT2A |unassigned  |      86|
|HBDN501-AML-KMT2A |donor3      |      67|
|HBDN501-AML-KMT2A |donor0      |      16|
|HBDN206-MNpCT     |donor0      |    3664|
|HBDN206-MNpCT     |donor1      |     205|
|HBDN206-MNpCT     |unassigned  |     101|
|HBDN206-MNpCT     |donor3      |      70|
|HBDN206-MNpCT     |donor2      |      36|
|MOLM13            |donor1      |    1912|
|MOLM13            |unassigned  |      47|
|MOLM13            |donor3      |      14|
|MOLM13            |donor0      |       7|
|MOLM13            |donor2      |       7|

## 5.2 Cell type composition


Table: Cell type composition by sample

|sampleID      |predicted_CellType          | n_cells| total_cells|  fraction|    percent|
|:-------------|:---------------------------|-------:|-----------:|---------:|----------:|
|HBDN206-MNpCT |HSC                         |    1789|        5697| 0.3140249| 31.4024925|
|HBDN206-MNpCT |MPP-MyLy                    |    1375|        5697| 0.2413551| 24.1355099|
|HBDN206-MNpCT |MPP-MkEry                   |    1172|        5697| 0.2057223| 20.5722310|
|HBDN206-MNpCT |CD14 Mono                   |     282|        5697| 0.0494997|  4.9499737|
|HBDN206-MNpCT |CD8 Effector Memory 2       |     184|        5697| 0.0322977|  3.2297701|
|HBDN206-MNpCT |LMPP                        |     164|        5697| 0.0287871|  2.8787081|
|HBDN206-MNpCT |CD4 Naive                   |      93|        5697| 0.0163244|  1.6324381|
|HBDN206-MNpCT |CD4 Central Memory          |      86|        5697| 0.0150957|  1.5095664|
|HBDN206-MNpCT |NK                          |      57|        5697| 0.0100053|  1.0005266|
|HBDN206-MNpCT |Orthochromatic Erythroblast |      56|        5697| 0.0098297|  0.9829735|
|HBDN206-MNpCT |CD8 Effector Memory 1       |      50|        5697| 0.0087765|  0.8776549|
|HBDN206-MNpCT |Mature B                    |      48|        5697| 0.0084255|  0.8425487|
|HBDN206-MNpCT |Early ProMono               |      41|        5697| 0.0071968|  0.7196770|
|HBDN206-MNpCT |CD4 Effector Memory         |      38|        5697| 0.0066702|  0.6670177|
|HBDN206-MNpCT |MEP                         |      36|        5697| 0.0063191|  0.6319115|
|HBDN206-MNpCT |Late ProMono                |      30|        5697| 0.0052659|  0.5265929|
|HBDN206-MNpCT |Megakaryocyte Precursor     |      29|        5697| 0.0050904|  0.5090398|
|HBDN206-MNpCT |CD8 Central Memory          |      25|        5697| 0.0043883|  0.4388275|
|HBDN206-MNpCT |Cycling Progenitor          |      21|        5697| 0.0036862|  0.3686151|
|HBDN206-MNpCT |GMP-Mono                    |      17|        5697| 0.0029840|  0.2984027|
|HBDN206-MNpCT |Early GMP                   |      13|        5697| 0.0022819|  0.2281903|
|HBDN206-MNpCT |NK CD56high                 |      13|        5697| 0.0022819|  0.2281903|
|HBDN206-MNpCT |CD4 Regulatory              |      11|        5697| 0.0019308|  0.1930841|
|HBDN206-MNpCT |CD8 Tissue Resident Memory  |      11|        5697| 0.0019308|  0.1930841|
|HBDN206-MNpCT |CFU-E                       |       8|        5697| 0.0014042|  0.1404248|
|HBDN206-MNpCT |Immature B                  |       8|        5697| 0.0014042|  0.1404248|
|HBDN206-MNpCT |BFU-E                       |       7|        5697| 0.0012287|  0.1228717|
|HBDN206-MNpCT |CD16 Mono                   |       7|        5697| 0.0012287|  0.1228717|
|HBDN206-MNpCT |Polychromatic Erythroblast  |       7|        5697| 0.0012287|  0.1228717|
|HBDN206-MNpCT |GMP-Neut                    |       4|        5697| 0.0007021|  0.0702124|
|HBDN206-MNpCT |CD8 Naive                   |       3|        5697| 0.0005266|  0.0526593|
|HBDN206-MNpCT |CLP                         |       2|        5697| 0.0003511|  0.0351062|
|HBDN206-MNpCT |Pro-B VDJ                   |       2|        5697| 0.0003511|  0.0351062|
|HBDN206-MNpCT |Basophilic Erythroblast     |       1|        5697| 0.0001755|  0.0175531|
|HBDN206-MNpCT |EoBasoMast Precursor        |       1|        5697| 0.0001755|  0.0175531|
|HBDN206-MNpCT |GMP-Cycle                   |       1|        5697| 0.0001755|  0.0175531|
|HBDN206-MNpCT |Plasma Cell                 |       1|        5697| 0.0001755|  0.0175531|
|HBDN206-MNpCT |Pre-pDC                     |       1|        5697| 0.0001755|  0.0175531|
|HBDN206-MNpCT |Pro-Erythroblast            |       1|        5697| 0.0001755|  0.0175531|
|HBDN206-MNpCT |T Proliferating             |       1|        5697| 0.0001755|  0.0175531|




Table: Lineage composition by sample

|sampleID          |predicted_Lineage | n_cells|    percent|
|:-----------------|:-----------------|-------:|----------:|
|HBDN206-MNpCT     |Lymphoid          |      58|  1.0180797|
|HBDN206-MNpCT     |Myeloid / DC      |       1|  0.0175531|
|HBDN206-MNpCT     |Other             |    5374| 94.3303493|
|HBDN206-MNpCT     |Stem / progenitor |     264|  4.6340179|
|HBDN392-AML-MDS   |Lymphoid          |      90|  1.5728766|
|HBDN392-AML-MDS   |Myeloid / DC      |     617| 10.7829430|
|HBDN392-AML-MDS   |Other             |    4708| 82.2789235|
|HBDN392-AML-MDS   |Stem / progenitor |     307|  5.3652569|
|HBDN501-AML-KMT2A |Lymphoid          |      22|  0.7848733|
|HBDN501-AML-KMT2A |Other             |    1791| 63.8958259|
|HBDN501-AML-KMT2A |Stem / progenitor |     990| 35.3193007|
|MOLM13            |Myeloid / DC      |       8|  0.1427042|
|MOLM13            |Other             |    4836| 86.2647164|
|MOLM13            |Stem / progenitor |     762| 13.5925794|

## 5.3 DSB vs CLR

*Not available on disk.*



Table: Background suppression

|marker               |  mean_CLR|   mean_DSB| delta_DSB_minus_CLR|
|:--------------------|---------:|----------:|-------------------:|
|IgM                  | 2.8168731| -0.1249694|          -2.9418425|
|IgG-Fc               | 2.6873402| -0.0377656|          -2.7251058|
|Mouse-IgG2a          | 0.5459849|  0.0542048|          -0.4917802|
|Rat-IgG1-l           | 0.3367272|  0.1759701|          -0.1607571|
|Rat-IgG2c-k          | 0.2259942|  0.0989940|          -0.1270002|
|FCER1-alpha          | 0.1691688|  0.0823210|          -0.0868477|
|Mouse-IgG2b          | 0.1457389|  0.1309483|          -0.0147906|
|Armenian-Hamster-IgG | 0.1067859|  0.1315473|           0.0247614|
|Rat-IgG2a-k          | 0.1058658|  0.1784481|           0.0725823|
|Mouse-IgG1           | 0.0900038|  0.1976486|           0.1076448|
|Rat-IgG2b            | 0.0376943|  0.1582002|           0.1205059|
|Rat-IgG1-k           | 0.0618827|  0.1836524|           0.1217696|




Table: Expected marker specificity

|marker |expected_broad    | expected_mean_CLR| expected_mean_DSB| other_median_CLR| other_median_DSB| specificity_CLR| specificity_DSB| expected_pct_positive_CLR| expected_pct_positive_DSB| other_pct_positive_CLR| other_pct_positive_DSB| pct_specificity_CLR| pct_specificity_DSB| delta_specificity_DSB_minus_CLR| delta_pct_specificity_DSB_minus_CLR|
|:------|:-----------------|-----------------:|-----------------:|----------------:|----------------:|---------------:|---------------:|-------------------------:|-------------------------:|----------------------:|----------------------:|-------------------:|-------------------:|-------------------------------:|-----------------------------------:|
|CD8    |Lymphoid          |         1.0951425|        15.6111854|        0.0000000|       -0.1631998|       1.0951425|      15.7743852|                  61.45340|                  61.45340|               8.937921|               8.937921|          52.5154758|          52.5154758|                      14.6792427|                            0.000000|
|CD11b  |Myeloid / DC      |         1.8591633|         9.3511384|        0.0000000|        0.0931207|       1.8591633|       9.2580177|                  98.89503|                  98.61878|              27.015483|              62.270867|          71.8795444|          36.3479178|                       7.3988544|                          -35.531627|
|CD4    |Lymphoid          |         0.7939794|         4.6337958|        0.0000000|        0.1094111|       0.7939794|       4.5243847|                  57.81991|                  66.19273|              44.913987|              60.284218|          12.9059187|           5.9085146|                       3.7304053|                           -6.997404|
|CD19   |Lymphoid          |         0.2405145|         2.8677163|        0.0000000|       -0.0869771|       0.2405145|       2.9546933|                  25.43444|                  36.80885|              13.219895|              30.572177|          12.2145439|           6.2366702|                       2.7141789|                           -5.977874|
|CD133  |Stem / progenitor |         2.1438014|         3.7335483|        0.8987586|        0.5035712|       1.2450427|       3.2299771|                 100.00000|                  99.50313|              97.928994|              71.745562|           2.0710059|          27.7575703|                       1.9849344|                           25.686564|
|CD34   |Stem / progenitor |         4.2281987|         4.5069696|        0.9300531|        0.1510431|       3.2981456|       4.3559264|                  99.97840|                  99.65435|              96.745562|              60.798817|           3.2328349|          38.8555364|                       1.0577808|                           35.622701|
|CD20   |Lymphoid          |         0.0448804|         0.9003657|        0.0000000|        0.0872981|       0.0448804|       0.8130676|                  17.06161|                  63.98104|               5.422588|              64.285714|          11.6390235|          -0.3046716|                       0.7681873|                          -11.943695|
|CD14   |Myeloid / DC      |         0.7182740|         0.8920683|        0.4792751|        0.1233388|       0.2389989|       0.7687294|                  95.02762|                  76.51934|              89.250756|              56.807261|           5.7768679|          19.7120759|                       0.5297306|                           13.935208|
|CD163  |Myeloid / DC      |         1.6092160|         0.5717646|        1.4400683|        0.0495966|       0.1691477|       0.5221680|                 100.00000|                  81.21547|              99.839829|              53.425877|           0.1601708|          27.7895931|                       0.3530203|                           27.629422|
|CD33   |Myeloid / DC      |         1.3813027|         1.1280148|        0.6163499|        0.0929426|       0.7649528|       1.0350722|                 100.00000|                  96.68508|              89.980424|              59.156433|          10.0195764|          37.5286493|                       0.2701194|                           27.509073|

## 5.4 CNV calls


Table: CopyKAT calls - HBDN206-MNpCT

|sampleID      |copykat_call | n_cells|   percent|
|:-------------|:------------|-------:|---------:|
|HBDN206-MNpCT |aneuploid    |    4619| 81.077760|
|HBDN206-MNpCT |diploid      |     774| 13.586098|
|HBDN206-MNpCT |not.defined  |     304|  5.336142|




Table: CopyKAT calls - HBDN501-AML-KMT2A

|sampleID          |copykat_call | n_cells|  percent|
|:-----------------|:------------|-------:|--------:|
|HBDN501-AML-KMT2A |aneuploid    |    1965| 70.10346|
|HBDN501-AML-KMT2A |diploid      |     498| 17.76668|
|HBDN501-AML-KMT2A |not.defined  |     340| 12.12986|




Table: Numbat compartment by sample and broad annotation

|sample        |predicted_CellType_Broad |numbat_compartment |    n|     percent|
|:-------------|:------------------------|:------------------|----:|-----------:|
|HBDN206-MNpCT |B                        |tumor              |    1|   1.7857143|
|HBDN206-MNpCT |B                        |normal             |   18|  32.1428571|
|HBDN206-MNpCT |B                        |not_called         |   37|  66.0714286|
|HBDN206-MNpCT |CD4 Memory T             |tumor              |    1|   0.7407407|
|HBDN206-MNpCT |CD4 Memory T             |normal             |   27|  20.0000000|
|HBDN206-MNpCT |CD4 Memory T             |not_called         |  107|  79.2592593|
|HBDN206-MNpCT |CD8 Memory T             |normal             |   54|  19.9261993|
|HBDN206-MNpCT |CD8 Memory T             |not_called         |  217|  80.0738007|
|HBDN206-MNpCT |Cycling Progenitor       |tumor              |   20|  95.2380952|
|HBDN206-MNpCT |Cycling Progenitor       |not_called         |    1|   4.7619048|
|HBDN206-MNpCT |Early Erythroid          |tumor              |    8|  80.0000000|
|HBDN206-MNpCT |Early Erythroid          |not_called         |    2|  20.0000000|
|HBDN206-MNpCT |Early GMP                |tumor              |    5|  38.4615385|
|HBDN206-MNpCT |Early GMP                |not_called         |    8|  61.5384615|
|HBDN206-MNpCT |Early Lymphoid           |not_called         |    2| 100.0000000|
|HBDN206-MNpCT |EoBasoMast Precursor     |not_called         |    1| 100.0000000|
|HBDN206-MNpCT |HSC MPP                  |tumor              | 3920|  90.4059041|
|HBDN206-MNpCT |HSC MPP                  |normal             |    7|   0.1614391|
|HBDN206-MNpCT |HSC MPP                  |not_called         |  409|   9.4326568|
|HBDN206-MNpCT |LMPP                     |tumor              |  148|  90.2439024|
|HBDN206-MNpCT |LMPP                     |not_called         |   16|   9.7560976|
|HBDN206-MNpCT |Late Erythroid           |tumor              |    5|   7.9365079|
|HBDN206-MNpCT |Late Erythroid           |not_called         |   58|  92.0634921|
|HBDN206-MNpCT |Late GMP                 |tumor              |   16|  72.7272727|
|HBDN206-MNpCT |Late GMP                 |not_called         |    6|  27.2727273|
|HBDN206-MNpCT |MEP                      |tumor              |   15|  34.8837209|
|HBDN206-MNpCT |MEP                      |not_called         |   28|  65.1162791|
|HBDN206-MNpCT |Megakaryocyte Precursor  |tumor              |    7|  24.1379310|
|HBDN206-MNpCT |Megakaryocyte Precursor  |not_called         |   22|  75.8620690|
|HBDN206-MNpCT |Monocyte                 |tumor              |  122|  42.2145329|
|HBDN206-MNpCT |Monocyte                 |normal             |    1|   0.3460208|
|HBDN206-MNpCT |Monocyte                 |not_called         |  166|  57.4394464|
|HBDN206-MNpCT |NK                       |tumor              |    1|   1.4285714|
|HBDN206-MNpCT |NK                       |normal             |    9|  12.8571429|
|HBDN206-MNpCT |NK                       |not_called         |   60|  85.7142857|
|HBDN206-MNpCT |Naive T                  |normal             |    6|   6.2500000|
|HBDN206-MNpCT |Naive T                  |not_called         |   90|  93.7500000|
|HBDN206-MNpCT |Plasma Cell              |normal             |    1| 100.0000000|
|HBDN206-MNpCT |Pro-B                    |not_called         |    2| 100.0000000|
|HBDN206-MNpCT |Pro-Monocyte             |tumor              |   28|  39.4366197|

# 6. Figures

Every PDF, PNG and JPEG under this run's output directories, with the date it was written.
Descriptions are hand-written for the figures that carry the argument; the rest are listed so
the index is complete.




### `data_integration/260423_VH01624_453_222HWMYNX`

| Figure | Date | What it shows |
|---|---|---|
| `01_background_suppression.pdf` | 2026-05-20 | How much background DSB removes. The main argument for using DSB. |
| `02_expected_marker_percent_specificity.pdf` | 2026-05-20 |  |
| `02_expected_marker_specificity.pdf` | 2026-05-20 | Whether expected markers land on the expected cell types. |
| `03_marker_leakage_diffusion.pdf` | 2026-05-20 | Marker signal leaking into populations that should be negative. |
| `04_expected_vs_other_violin.pdf` | 2026-05-20 |  |
| `05_featureplots_CLR_vs_DSB_expected_markers.pdf` | 2026-05-20 | Side-by-side feature plots. Visual summary of the whole comparison. |
| `07_broad_celltype_composition_by_sample.pdf` | 2026-05-06 |  |
| `08_detailed_celltype_composition_by_sample.pdf` | 2026-05-06 |  |
| `09_DSB_UMAP.pdf` | 2026-05-06 |  |
| `09_DSB_UMAP_by_sample.pdf` | 2026-05-06 |  |
| `10_DSB_UMAP_by_sample_variable_proteins.pdf` | 2026-05-06 |  |
| `11_copykat_by_sample_and_broad_annotation.pdf` | 2026-05-07 |  |


### `data_integration/260423_VH01624_453_222HWMYNX/09_WNN_RNA_DSB`

| Figure | Date | What it shows |
|---|---|---|
| `01_WNN_RNA_DSB_clusters.pdf` | 2026-05-06 |  |
| `02_WNN_RNA_DSB_clusters_by_sample.pdf` | 2026-05-06 |  |
| `03_WNN_RNA_DSB_broad_celltypes_by_sample.pdf` | 2026-05-06 |  |


### `data_integration/260423_VH01624_453_222HWMYNX/09_WNN_RNA_DSB/CITE_marker_featureplots_by_sample`

| Figure | Date | What it shows |
|---|---|---|
| `CITE_markers_HBDN206.MNpCT_CD8.Memory.T.pdf` | 2026-05-06 |  |
| `CITE_markers_HBDN206.MNpCT_HSC.MPP.pdf` | 2026-05-06 |  |
| `CITE_markers_HBDN392.AML.MDS_Monocyte.pdf` | 2026-05-06 |  |
| `CITE_markers_HBDN392.AML.MDS_cDC.pdf` | 2026-05-06 |  |
| `CITE_markers_HBDN392.AML.MDS_pDC.pdf` | 2026-05-06 |  |
| `CITE_markers_HBDN501.AML.KMT2A_Monocyte.pdf` | 2026-05-06 |  |
| `CITE_markers_MOLM13_Cycling.Progenitor.pdf` | 2026-05-06 |  |
| `CITE_markers_MOLM13_Late.GMP.pdf` | 2026-05-06 |  |
| `CITE_markers_MOLM13_Pro.Monocyte.pdf` | 2026-05-06 |  |
| `DotPlot_CITE_markers_HBDN206.MNpCT_CD8.Memory.T.pdf` | 2026-05-06 |  |
| `DotPlot_CITE_markers_HBDN206.MNpCT_HSC.MPP.pdf` | 2026-05-06 |  |
| `DotPlot_CITE_markers_HBDN392.AML.MDS_Monocyte.pdf` | 2026-05-06 |  |
| `DotPlot_CITE_markers_HBDN392.AML.MDS_cDC.pdf` | 2026-05-06 |  |
| `DotPlot_CITE_markers_HBDN392.AML.MDS_pDC.pdf` | 2026-05-06 |  |
| `DotPlot_CITE_markers_HBDN501.AML.KMT2A_Monocyte.pdf` | 2026-05-06 |  |
| `DotPlot_CITE_markers_MOLM13_Cycling.Progenitor.pdf` | 2026-05-06 |  |
| `DotPlot_CITE_markers_MOLM13_Late.GMP.pdf` | 2026-05-06 |  |
| `DotPlot_CITE_markers_MOLM13_Pro.Monocyte.pdf` | 2026-05-06 |  |


### `data_integration/260423_VH01624_453_222HWMYNX/cellranger_downsampled_DE_and_clustering/DE_results`

| Figure | Date | What it shows |
|---|---|---|
| `cellranger_downsampled_DE_total.pdf` | 2026-06-01 |  |
| `cellranger_downsampled_DE_up_down.pdf` | 2026-06-01 |  |


### `data_integration/260423_VH01624_453_222HWMYNX/cellranger_downsampled_DE_and_clustering/HBDN501_downsampled_umaps`

| Figure | Date | What it shows |
|---|---|---|
| `HBDN501_cluster_separation_metrics.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_010_HBDN501-AML-KMT2A_UMAP_broad_celltypes.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_010_HBDN501-AML-KMT2A_UMAP_clusters.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_010_HBDN501-AML-KMT2A_UMAP_combined.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_020_HBDN501-AML-KMT2A_UMAP_broad_celltypes.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_020_HBDN501-AML-KMT2A_UMAP_clusters.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_020_HBDN501-AML-KMT2A_UMAP_combined.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_030_HBDN501-AML-KMT2A_UMAP_broad_celltypes.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_030_HBDN501-AML-KMT2A_UMAP_clusters.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_030_HBDN501-AML-KMT2A_UMAP_combined.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_040_HBDN501-AML-KMT2A_UMAP_broad_celltypes.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_040_HBDN501-AML-KMT2A_UMAP_clusters.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_040_HBDN501-AML-KMT2A_UMAP_combined.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_050_HBDN501-AML-KMT2A_UMAP_broad_celltypes.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_050_HBDN501-AML-KMT2A_UMAP_clusters.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_050_HBDN501-AML-KMT2A_UMAP_combined.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_060_HBDN501-AML-KMT2A_UMAP_broad_celltypes.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_060_HBDN501-AML-KMT2A_UMAP_clusters.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_060_HBDN501-AML-KMT2A_UMAP_combined.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_070_HBDN501-AML-KMT2A_UMAP_broad_celltypes.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_070_HBDN501-AML-KMT2A_UMAP_clusters.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_070_HBDN501-AML-KMT2A_UMAP_combined.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_080_HBDN501-AML-KMT2A_UMAP_broad_celltypes.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_080_HBDN501-AML-KMT2A_UMAP_clusters.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_080_HBDN501-AML-KMT2A_UMAP_combined.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_090_HBDN501-AML-KMT2A_UMAP_broad_celltypes.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_090_HBDN501-AML-KMT2A_UMAP_clusters.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_090_HBDN501-AML-KMT2A_UMAP_combined.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_100_HBDN501-AML-KMT2A_UMAP_broad_celltypes.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_100_HBDN501-AML-KMT2A_UMAP_clusters.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_100_HBDN501-AML-KMT2A_UMAP_combined.pdf` | 2026-06-01 |  |


### `data_integration/260423_VH01624_453_222HWMYNX/cellranger_downsampled_DE_and_clustering/annotation_robustness`

| Figure | Date | What it shows |
|---|---|---|
| `LK1-GEX_frac_010_mapping_error_QC.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_010_projected_predicted_celltypes.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_020_mapping_error_QC.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_020_projected_predicted_celltypes.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_030_mapping_error_QC.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_030_projected_predicted_celltypes.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_040_mapping_error_QC.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_040_projected_predicted_celltypes.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_050_mapping_error_QC.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_050_projected_predicted_celltypes.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_060_mapping_error_QC.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_060_projected_predicted_celltypes.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_070_mapping_error_QC.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_070_projected_predicted_celltypes.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_080_mapping_error_QC.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_080_projected_predicted_celltypes.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_090_mapping_error_QC.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_090_projected_predicted_celltypes.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_100_mapping_error_QC.pdf` | 2026-06-01 |  |
| `LK1-GEX_frac_100_projected_predicted_celltypes.pdf` | 2026-06-01 |  |
| `annotation_accuracy_across_downsample_depths.pdf` | 2026-06-01 |  |
| `annotation_confusion_broad_across_depths.pdf` | 2026-06-01 |  |
| `annotation_confusion_fine_across_depths.pdf` | 2026-06-01 |  |


### `data_integration/260423_VH01624_453_222HWMYNX/projection_umaps`

| Figure | Date | What it shows |
|---|---|---|
| `01_projected_predicted_celltypes_broad.pdf` | 2026-05-08 |  |
| `02_projected_predicted_celltypes_granulated.pdf` | 2026-05-08 |  |
| `03_projected_samples.pdf` | 2026-05-08 |  |
| `04_projected_pseudotime.pdf` | 2026-05-08 |  |
| `05_broad_annotations_split_by_sample.pdf` | 2026-05-08 |  |
| `06_detailed_annotations_split_by_sample.pdf` | 2026-05-08 |  |
| `06a_detailed_annotations_split_by_sample_empty.pdf` | 2026-05-07 |  |
| `07_broad_celltype_composition_by_sample.pdf` | 2026-05-06 |  |
| `09_copykat_on_joint_umap.pdf` | 2026-05-07 |  |
| `10_copykat_on_umap_split_by_sample.pdf` | 2026-05-07 |  |
| `11.top3_DSB_CITE_markers_per_broad_annotation_cutoff_clustered_3.pdf` | 2026-05-08 |  |
| `12.top3_DSB_CITE_markers_per_detailed_annotation_clustered_cutoff_3.pdf` | 2026-05-08 |  |
| `ADT_raw_marker_read_distributions_Tier_1_highest.pdf` | 2026-05-07 |  |
| `ADT_raw_marker_read_distributions_Tier_2_middle.pdf` | 2026-05-07 |  |
| `ADT_raw_marker_read_distributions_Tier_3_lowest.pdf` | 2026-05-07 |  |
| `ADT_raw_marker_read_distributions_density.pdf` | 2026-05-07 |  |
| `top3_DSB_CITE_markers_per_broad_annotation_cutoff_2.pdf` | 2026-05-07 |  |
| `top3_DSB_CITE_markers_per_broad_annotation_cutoff_3.pdf` | 2026-05-08 |  |
| `top3_DSB_CITE_markers_per_broad_annotation_cutoff_3.v2.pdf` | 2026-05-07 |  |
| `top3_DSB_CITE_markers_per_broad_annotation_manual_dotplot.pdf` | 2026-05-07 |  |


### `data_integration/260423_VH01624_453_222HWMYNX/projection_umaps/HSC_MPP_copykat_DSB_HBDN206-MNpCT`

| Figure | Date | What it shows |
|---|---|---|
| `01_HSC_MPP_DSB_volcano.pdf` | 2026-05-20 |  |
| `02_HSC_MPP_top_DSB_featureplots.pdf` | 2026-05-20 |  |


### `data_integration/260423_VH01624_453_222HWMYNX/projection_umaps/HSC_MPP_numbat_DSB_HBDN206-MNpCT`

| Figure | Date | What it shows |
|---|---|---|
| `01_HSC_MPP_DSB_numbat_volcano.pdf` | 2026-05-20 |  |
| `02_all_cells_numbat_malignancy_umap.pdf` | 2026-05-20 |  |
| `03_HSC_MPP_numbat_malignancy_umap.pdf` | 2026-05-20 |  |
| `04_HSC_MPP_top_DSB_featureplots_numbat.pdf` | 2026-05-20 |  |
| `05_HSC_MPP_top_DSB_violins_numbat.pdf` | 2026-05-20 |  |


### `data_integration/260423_VH01624_453_222HWMYNX/projection_umaps/numbat_copykat_comparison_HBDN206-MNpCT`

| Figure | Date | What it shows |
|---|---|---|
| `04_binary_overlap_heatmap.pdf` | 2026-05-21 |  |
| `05_binary_overlap_classes.pdf` | 2026-05-21 |  |
| `06_celltype_by_binary_overlap.pdf` | 2026-05-21 |  |
| `07_binary_overlap_umap.pdf` | 2026-05-21 |  |
| `09_overlap_heatmap.pdf` | 2026-05-21 |  |
| `10_overlap_class_barplot.pdf` | 2026-05-21 |  |
| `11_celltype_composition_by_overlap.pdf` | 2026-05-21 |  |
| `12_celltype_discordance_enrichment.pdf` | 2026-05-21 |  |
| `12b_discordant_direction_by_celltype.pdf` | 2026-05-21 |  |
| `13_umap_numbat_copykat_overlap_celltype.pdf` | 2026-05-21 |  |
| `16_numbat_clone_umap.pdf` | 2026-05-21 |  |
| `17_numbat_clone_split_by_overlap.pdf` | 2026-05-21 |  |
| `18_numbat_only_cells_umap.pdf` | 2026-05-21 |  |
| `20_clone_overlap_barplot.pdf` | 2026-05-21 |  |


### `data_integration/260423_VH01624_453_222HWMYNX/projection_umaps/numbat_copykat_comparison_HBDN206-MNpCT/clone_1_inspection`

| Figure | Date | What it shows |
|---|---|---|
| `04_clone_1_inspection_umaps.pdf` | 2026-05-21 |  |
| `05_clone_1_qc_violin.pdf` | 2026-05-21 |  |
| `08_possible_confidence_featureplots.pdf` | 2026-05-21 |  |


### `data_integration/260423_VH01624_453_222HWMYNX/projection_umaps/numbat_copykat_comparison_HBDN501-AML-KMT2A`

| Figure | Date | What it shows |
|---|---|---|
| `04_binary_overlap_heatmap.pdf` | 2026-05-21 |  |
| `05_binary_overlap_classes.pdf` | 2026-05-21 |  |
| `06_celltype_by_binary_overlap.pdf` | 2026-05-21 |  |
| `07_binary_overlap_umap.pdf` | 2026-05-21 |  |
| `16_numbat_clone_umap.pdf` | 2026-05-21 |  |
| `17_numbat_clone_split_by_overlap.pdf` | 2026-05-21 |  |
| `18_numbat_only_cells_umap.pdf` | 2026-05-21 |  |
| `20_clone_overlap_barplot.pdf` | 2026-05-21 |  |


### `data_integration/260423_VH01624_453_222HWMYNX/projection_umaps/numbat_copykat_comparison_noDoublet_HBDN206-MNpCT`

| Figure | Date | What it shows |
|---|---|---|
| `04_binary_overlap_heatmap.pdf` | 2026-05-25 |  |
| `05_binary_overlap_classes.pdf` | 2026-05-25 |  |
| `06_celltype_by_binary_overlap.pdf` | 2026-05-25 |  |
| `07_binary_overlap_umap.pdf` | 2026-05-25 |  |
| `16_numbat_clone_umap.pdf` | 2026-05-25 |  |
| `17_numbat_clone_split_by_overlap.pdf` | 2026-05-25 |  |
| `18_numbat_only_cells_umap.pdf` | 2026-05-25 |  |
| `20_clone_overlap_barplot.pdf` | 2026-05-25 |  |


### `data_integration/260423_VH01624_453_222HWMYNX/projection_umaps/numbat_copykat_comparison_noDoublet_HBDN206-MNpCT/clone_1_inspection`

| Figure | Date | What it shows |
|---|---|---|
| `04_clone_1_inspection_umaps.pdf` | 2026-05-25 |  |
| `05_clone_1_qc_violin.pdf` | 2026-05-25 |  |
| `08_possible_confidence_featureplots.pdf` | 2026-05-25 |  |


### `data_integration/260423_VH01624_453_222HWMYNX/read_downsampled_pseudobulk_DE`

| Figure | Date | What it shows |
|---|---|---|
| `read_downsampled_pseudobulk_DE_total.pdf` | 2026-05-27 |  |
| `read_downsampled_pseudobulk_DE_total.png` | 2026-05-27 |  |
| `read_downsampled_pseudobulk_DE_up_down.pdf` | 2026-05-27 |  |
| `read_downsampled_pseudobulk_DE_up_down.png` | 2026-05-27 |  |


### `demux_comparison/260423_VH01624_453_222HWMYNX/LK1-GEX`

| Figure | Date | What it shows |
|---|---|---|
| `01_demux_efficiency_barplot.pdf` | 2026-05-04 | Singlet/doublet/negative rates for HTO and vireo side by side. The headline comparison. |
| `01_demux_efficiency_barplot_all_positive_quantiles.pdf` | 2026-05-05 |  |
| `02_hto_vireo_doublet_comparison.pdf` | 2026-05-04 | Whether the two methods call the same cells doublets. |
| `02_hto_vireo_doublet_comparison_all_positive_quantiles.pdf` | 2026-05-05 |  |
| `03_hto_vireo_donor_heatmap.pdf` | 2026-05-04 | HTO sample against vireo donor. Off-diagonal mass is the disagreement. |
| `03_hto_vireo_donor_heatmap_all_positive_quantiles.pdf` | 2026-05-05 |  |
| `04_upset_HTO_vireo_overlap.pdf` | 2026-05-04 | UpSet of the call overlap. |
| `04_upset_HTO_vireo_overlap_all_positive_quantiles.pdf` | 2026-05-05 |  |
| `05_HTO_QC.pdf` | 2026-05-04 |  |
| `05_HTO_QC_threshold_independent.pdf` | 2026-05-05 |  |
| `06_HTO_ridge.pdf` | 2026-05-04 |  |
| `06_HTO_ridge_threshold_independent.pdf` | 2026-05-05 |  |
| `06b_HTO_ratio_thresholds_threshold_independent.pdf` | 2026-05-05 |  |
| `07_UMAP_Vireo_and_HTO_all_positive_quantiles.pdf` | 2026-05-05 |  |


### `emptydrops/260423_VH01624_453_222HWMYNX/LK1-GEX`

| Figure | Date | What it shows |
|---|---|---|
| `HTO_QC.pdf` | 2026-05-01 |  |


### `genotype_demux/260423_VH01624_453_222HWMYNX/LK1-GEX/vireo`

| Figure | Date | What it shows |
|---|---|---|
| `fig_GT_distance_estimated.pdf` | 2026-05-04 |  |


### `scDblFinder/260423_VH01624_453_222HWMYNX/LK1-GEX`

| Figure | Date | What it shows |
|---|---|---|
| `RNA_QC_before_filtering.pdf` | 2026-05-03 |  |
| `scDblFinder_QC.pdf` | 2026-05-03 |  |


### `seurat_annotated/260423_VH01624_453_222HWMYNX`

| Figure | Date | What it shows |
|---|---|---|
| `00_BoneMarrowMap_reference_celltype.pdf` | 2026-05-05 |  |
| `00a_BoneMarrowMap_reference_broad_celltype.pdf` | 2026-05-05 |  |
| `01_background_suppression.pdf` | 2026-05-20 | How much background DSB removes. The main argument for using DSB. |
| `02_expected_marker_specificity.pdf` | 2026-05-20 | Whether expected markers land on the expected cell types. |
| `08_celltype_composition_by_sampleID_stacked_bar.pdf` | 2026-05-05 |  |
| `09_celltype_counts_by_sampleID_stacked_bar.pdf` | 2026-05-05 |  |
| `11_lineage_composition_by_sampleID.pdf` | 2026-05-05 |  |
| `13_specific_celltype_composition_by_sampleID.pdf` | 2026-05-05 |  |


### `seurat_annotated/260423_VH01624_453_222HWMYNX/copykat_annotated/HBDN206-MNpCT`

| Figure | Date | What it shows |
|---|---|---|
| `01_copykat_call_composition_by_sampleID.pdf` | 2026-05-12 | CopyKAT aneuploid/diploid proportions. |
| `02_copykat_calls_on_projected_RNA_UMAP.pdf` | 2026-05-12 | CopyKAT calls on the projected UMAP. |
| `02a_copykat_calls_on_projected_RNA_UMAP.pdf` | 2026-05-12 |  |
| `03_copykat_by_BoneMarrowMap_broad_celltype.pdf` | 2026-05-12 | Calls against annotation. Checks 'aneuploid' is not landing on one cell type for technical reasons. |
| `03a_copykat_counts_by_BoneMarrowMap_broad_celltype.pdf` | 2026-05-12 |  |
| `LK1_copykat_copykat_heatmap.jpeg` | 2026-05-12 | The CopyKAT CNV heatmap itself. |


### `seurat_annotated/260423_VH01624_453_222HWMYNX/copykat_annotated/HBDN392-AML-MDS`

| Figure | Date | What it shows |
|---|---|---|
| `01_copykat_call_composition_by_sampleID.pdf` | 2026-05-12 | CopyKAT aneuploid/diploid proportions. |
| `02_copykat_calls_on_projected_RNA_UMAP.pdf` | 2026-05-12 | CopyKAT calls on the projected UMAP. |
| `02a_copykat_calls_on_projected_RNA_UMAP.pdf` | 2026-05-12 |  |
| `03_copykat_by_BoneMarrowMap_broad_celltype.pdf` | 2026-05-12 | Calls against annotation. Checks 'aneuploid' is not landing on one cell type for technical reasons. |
| `03a_copykat_counts_by_BoneMarrowMap_broad_celltype.pdf` | 2026-05-12 |  |
| `LK1_copykat_copykat_heatmap.jpeg` | 2026-05-12 | The CopyKAT CNV heatmap itself. |


### `seurat_annotated/260423_VH01624_453_222HWMYNX/copykat_annotated/HBDN501-AML-KMT2A`

| Figure | Date | What it shows |
|---|---|---|
| `01_copykat_call_composition_by_sampleID.pdf` | 2026-05-12 | CopyKAT aneuploid/diploid proportions. |
| `02_copykat_calls_on_projected_RNA_UMAP.pdf` | 2026-05-12 | CopyKAT calls on the projected UMAP. |
| `02a_copykat_calls_on_projected_RNA_UMAP.pdf` | 2026-05-12 |  |
| `03_copykat_by_BoneMarrowMap_broad_celltype.pdf` | 2026-05-12 | Calls against annotation. Checks 'aneuploid' is not landing on one cell type for technical reasons. |
| `03a_copykat_counts_by_BoneMarrowMap_broad_celltype.pdf` | 2026-05-12 |  |
| `LK1_copykat_copykat_heatmap.jpeg` | 2026-05-12 | The CopyKAT CNV heatmap itself. |


### `seurat_annotated/260423_VH01624_453_222HWMYNX/numbat`

| Figure | Date | What it shows |
|---|---|---|
| `09_numbat_compartment_on_umap.pdf` | 2026-05-08 | Numbat tumour/normal compartment on the UMAP. |
| `09_numbat_compartment_on_umap_projected.pdf` | 2026-05-08 |  |
| `09_numbat_on_joint_umap.pdf` | 2026-05-08 |  |
| `09a_numbat_on_joint_umap_projected.pdf` | 2026-05-08 |  |
| `10_numbat_compartment_on_harmony_projected_split_by_sample.pdf` | 2026-05-08 |  |
| `10_numbat_compartment_on_umap_projected_split_by_sample.pdf` | 2026-05-08 |  |
| `10_numbat_compartment_on_umap_split_by_sample.pdf` | 2026-05-08 |  |
| `10_numbat_on_umap_split_by_sample.pdf` | 2026-05-08 |  |
| `10a_numbat_on_umap_split_by_sample_projected.pdf` | 2026-05-08 |  |
| `10b_numbat_clone_on_harmony_projected_split_by_sample.pdf` | 2026-05-08 |  |
| `10b_numbat_clone_on_umap_projected_split_by_sample.pdf` | 2026-05-08 |  |
| `10b_numbat_clone_on_umap_split_by_sample.pdf` | 2026-05-08 | Numbat clone assignment per sample. |
| `11_numbat_by_sample_and_broad_annotation.pdf` | 2026-05-08 |  |
| `11_numbat_compartment_by_sample_and_broad_annotation.pdf` | 2026-05-08 | Compartment against annotation - checks the call is not just tracking one cell type. |
| `12_numbat_clone_by_sample_and_broad_annotation.pdf` | 2026-05-08 |  |


### `seurat_annotated/260423_VH01624_453_222HWMYNX/numbat/CITE_HSC_Numbat_analysis`

| Figure | Date | What it shows |
|---|---|---|
| `01_HSC_marker_presence_vs_other_cells.pdf` | 2026-06-15 |  |
| `02_volcano_HSC_CITE_DE_Cancer_vs_Normal_by_Numbat.pdf` | 2026-06-15 |  |
| `03_validation_candidate_markers_in_normal_sample_HSC.pdf` | 2026-06-15 |  |
| `04_heatmap_candidate_CITE_markers_HSC.pdf` | 2026-06-15 |  |
| `05_featureplots_candidate_HSC_CITE_markers.pdf` | 2026-06-15 |  |
| `06_BAFFR_mean_RNA_vs_mean_protein_by_celltype.pdf` | 2026-06-15 |  |
| `07_selected_CITE_markers_by_celltype_manual_dotplot.pdf` | 2026-06-15 |  |


### `seurat_annotated/260423_VH01624_453_222HWMYNX/numbat/LK1_HBDN206-MNpCT/numbat_final`

| Figure | Date | What it shows |
|---|---|---|
| `bulk_clones_1.png` | 2026-05-06 |  |
| `bulk_clones_final.png` | 2026-05-06 |  |
| `bulk_subtrees_1.png` | 2026-05-06 |  |
| `panel_1.png` | 2026-05-06 |  |


### `seurat_annotated/260423_VH01624_453_222HWMYNX/numbat/LK1_HBDN501-AML-KMT2A/numbat_final`

| Figure | Date | What it shows |
|---|---|---|
| `bulk_clones_1.png` | 2026-05-06 |  |
| `bulk_clones_final.png` | 2026-05-06 |  |
| `bulk_subtrees_1.png` | 2026-05-06 |  |
| `panel_1.png` | 2026-05-06 |  |


### `seurat_annotated/260423_VH01624_453_222HWMYNX/projectionFigures`

| Figure | Date | What it shows |
|---|---|---|
| `01_mapping_error_QC.pdf` | 2026-05-05 | Mapping error distribution, re-done on the DSB object. |
| `02_projected_predicted_celltypes_pass_only.pdf` | 2026-05-05 |  |
| `03_projected_samples_pass_only.pdf` | 2026-05-05 |  |
| `04_projected_pseudotime_pass_only.pdf` | 2026-05-05 |  |
| `05_projection_by_sampleID_BMM_helper_manual_split.pdf` | 2026-05-05 |  |
| `06_ADT_CLR_top_marker_projection.pdf` | 2026-05-05 | Top markers under CLR. Pair with the DSB version below - this is the CLR-vs-DSB comparison in its most direct form. |
| `06_CITE_DSB_top_marker_projection.pdf` | 2026-05-05 | The same top markers under DSB. |
| `07_ADT_CLR_selected_markers_manual_projectedUMAP.pdf` | 2026-05-05 |  |
| `07_CITE_DSB_selected_markers_manual_projectedUMAP.pdf` | 2026-05-05 |  |
| `07_CITE_DSB_selected_markers_ordered_projectedUMAP.pdf` | 2026-05-05 |  |
| `HBDN206-MNpCTdensity_HBDN206-MNpCT_projectedUMAP.pdf` | 2026-05-05 |  |
| `HBDN392-AML-MDSdensity_HBDN392-AML-MDS_projectedUMAP.pdf` | 2026-05-05 |  |
| `HBDN501-AML-KMT2Adensity_HBDN501-AML-KMT2A_projectedUMAP.pdf` | 2026-05-05 |  |
| `MOLM13density_MOLM13_projectedUMAP.pdf` | 2026-05-05 |  |
| `density_HBDN206-MNpCT_projectedUMAP.pdf` | 2026-05-05 |  |
| `density_HBDN392-AML-MDS_projectedUMAP.pdf` | 2026-05-05 |  |
| `density_HBDN501-AML-KMT2A_projectedUMAP.pdf` | 2026-05-05 |  |
| `density_MOLM13_projectedUMAP.pdf` | 2026-05-05 |  |


### `seurat_annotated/260423_VH01624_453_222HWMYNX/projectionFigures/DSB_vs_CLR_HBDN206-MNpCT`

| Figure | Date | What it shows |
|---|---|---|
| `01_marker_mean_shift_DSB_minus_CLR.pdf` | 2026-05-20 | Per-marker mean shift, DSB minus CLR. Quantifies what the eyeball comparison shows. |


### `seurat_demux/260423_VH01624_453_222HWMYNX`

| Figure | Date | What it shows |
|---|---|---|
| `01_initial_qc_violin.pdf` | 2026-04-28 | nFeature_RNA, nCount_RNA and percent.mt by hash.ID. First look - checks no hashed sample is degraded relative to the others. |
| `02_hto_ridgeplots.pdf` | 2026-04-28 | CLR-normalised HTO signal per hashtag. Tells you whether hashing worked at all; clean bimodality per HTO is what you want. |
| `03_umap_demux.pdf` | 2026-04-28 | HTO global classification, assigned hash.ID, and top HTO by raw count. The third panel is the sanity check on the second. |
| `04_umap_qc_features.pdf` | 2026-04-28 | UMI counts, detected genes and mitochondrial % on the UMAP. Localised high-MT regions are the dying-cell clusters. |
| `05_hto_featureplots.pdf` | 2026-04-28 | Each HTO on the UMAP. Confirms hashed samples occupy distinct territory rather than smearing. |
| `06_ADT_QC.pdf` | 2026-04-29 | ADT totals and features per cell. |
| `07_top_ADT_markers.pdf` | 2026-04-29 | Top 5 ADT markers by mean CLR plus biotin. Biotin was added deliberately - this is the start of the biotin thread that runs to July. |
| `08_top_ADT_cluster_markers_dotplot.pdf` | 2026-04-29 | Top ADT markers per RNA cluster. |
| `09_mapping_error_QC.pdf` | 2026-05-03 | BoneMarrowMap mapping error. Decides which cells get an annotation at all. |
| `10_projected_predicted_celltypes_pass_only.pdf` | 2026-05-03 |  |
| `10_top10_discriminative_ADT_plus_biotin_UMAP.pdf` | 2026-04-29 | Most discriminative ADTs with biotin overlaid. |
| `10a.database_umap.pdf` | 2026-05-03 |  |
| `11_biotin_expression_by_cluster_violin.pdf` | 2026-04-29 | Biotin by cluster - first look at the heterogeneity later split high vs low. |
| `11_projected_pseudotime_pass_only.pdf` | 2026-05-03 |  |
| `12_projection_by_sampleID.pdf` | 2026-05-03 |  |
| `13.broad_celltype_composition.pdf` | 2026-05-03 |  |
| `14.specific_celltype_composition.pdf` | 2026-05-03 |  |
| `15.heatmap_cell_type_broad.pdf` | 2026-05-03 |  |
| `15a.heatmap_cell_type_all.pdf` | 2026-05-03 |  |
| `celltype_composition_by_sampleID_stacked_bar.pdf` | 2026-05-03 |  |
| `celltype_counts_by_sampleID_stacked_bar.pdf` | 2026-05-03 |  |


### `seurat_demux/260423_VH01624_453_222HWMYNX/CITE_Late_GMP_HBDN206_MNpCT`

| Figure | Date | What it shows |
|---|---|---|
| `01_LateGMP_ADT_subclusters.pdf` | 2026-05-04 |  |
| `02_LateGMP_CopyKAT_on_ADT_UMAP.pdf` | 2026-05-04 |  |
| `03_top10_abundant_CITE_markers_LateGMP.pdf` | 2026-05-04 |  |
| `03_top12_abundant_CITE_markers_LateGMP.pdf` | 2026-05-04 |  |
| `04_top10_variable_CITE_markers_LateGMP.pdf` | 2026-05-04 |  |
| `04_top12_variable_CITE_markers_LateGMP.pdf` | 2026-05-04 |  |
| `05_top_separating_CITE_markers_on_ADT_UMAP.pdf` | 2026-05-04 |  |
| `06_top_separating_CITE_markers_heatmap.pdf` | 2026-05-04 |  |
| `07_top_separating_CITE_markers_dotplot.pdf` | 2026-05-04 |  |
| `08_LateGMP_ADT_subcluster_CopyKAT_composition.pdf` | 2026-05-04 |  |


### `seurat_demux/260423_VH01624_453_222HWMYNX/CITE_refinement`

| Figure | Date | What it shows |
|---|---|---|
| `01_ADT_UMAP_BoneMarrowMap_broad.pdf` | 2026-05-03 |  |
| `02_ADT_UMAP_ADT_clusters.pdf` | 2026-05-03 |  |
| `03_ADT_UMAP_sampleID.pdf` | 2026-05-03 |  |
| `04_WNN_UMAP_BoneMarrowMap_broad.pdf` | 2026-05-03 |  |
| `05_WNN_UMAP_WNN_clusters.pdf` | 2026-05-03 |  |
| `06_WNN_UMAP_sampleID.pdf` | 2026-05-03 |  |
| `07_WNN_canonical_ADT_markers.pdf` | 2026-05-03 |  |
| `08_WNN_cluster_composition_by_sampleID.pdf` | 2026-05-03 |  |
| `09_BoneMarrowMapBroad_vs_WNN_heatmap.pdf` | 2026-05-03 |  |
| `10_top_ADT_markers_by_WNN_cluster_heatmap.pdf` | 2026-05-03 |  |
| `11_within_lineage_WNN_refinement.pdf` | 2026-05-03 |  |


### `seurat_demux/260423_VH01624_453_222HWMYNX/copykat`

| Figure | Date | What it shows |
|---|---|---|
| `01_copykat_call_composition_by_sampleID.pdf` | 2026-05-11 | CopyKAT aneuploid/diploid proportions. |
| `02_copykat_calls_on_RNA_UMAP.pdf` | 2026-05-04 |  |
| `02_copykat_calls_on_projected_RNA_UMAP.pdf` | 2026-05-11 | CopyKAT calls on the projected UMAP. |
| `02a_copykat_calls_on_projected_RNA_UMAP.pdf` | 2026-05-11 |  |
| `03_copykat_by_BoneMarrowMap_broad_celltype.pdf` | 2026-05-11 | Calls against annotation. Checks 'aneuploid' is not landing on one cell type for technical reasons. |
| `03a_copykat_counts_by_BoneMarrowMap_broad_celltype.pdf` | 2026-05-11 |  |
| `LK1_copykat_copykat_heatmap.jpeg` | 2026-05-11 | The CopyKAT CNV heatmap itself. |
| `LK1_copykat_copykat_with_genes_heatmap.pdf` | 2026-05-11 |  |


### `seurat_demux/260423_VH01624_453_222HWMYNX/numbat/LK1_HBDN206-MNpCT/numbat_final`

| Figure | Date | What it shows |
|---|---|---|
| `exp_roll_clust.png` | 2026-05-05 |  |


### `seurat_demux/260423_VH01624_453_222HWMYNX/numbat/LK1_HBDN392-AML-MDS/numbat_final`

| Figure | Date | What it shows |
|---|---|---|
| `exp_roll_clust.png` | 2026-05-05 |  |


### `seurat_demux/260423_VH01624_453_222HWMYNX/numbat/LK1_HBDN501-AML-KMT2A/numbat_final`

| Figure | Date | What it shows |
|---|---|---|
| `exp_roll_clust.png` | 2026-05-05 |  |


### `seurat_demux/260423_VH01624_453_222HWMYNX/numbat/LK1_MOLM13/numbat_final`

| Figure | Date | What it shows |
|---|---|---|
| `exp_roll_clust.png` | 2026-05-05 |  |


Total: **303 figures**, of which 44 are described above.

# 7. Tools

Cluster side, from the scripts:

| Tool | Version | Source |
|---|---|---|
| Cell Ranger | 7.2.0 (uncertain) | commented module line in `scripts/1.cellranger.sh` |
| Reference | `refdata-gex-GRCh38-2024-A` | `annotation/` |
| CITE-seq-Count | not recorded | `scripts/3.hto_count.sbatch`, `4.cite_seq_count.sbatch` |
| cellsnp-lite, vireo | not recorded | `scripts/3.cellSNP.sh` |

R side. **No `sessionInfo()` was captured at the time**, and `scripts/` was not under version
control until 2026-07-31, so the exact versions used in April and May 2026 are not
recoverable. The analysis ran on the laptop against the `bioinf_scratch` mount. The session
that renders this report is:

```
R version 4.5.1 (2025-06-13) 
Platform: aarch64-apple-darwin20 
Running under: macOS Tahoe 26.5.2 

tidyverse 2.0.0 
```

Treat that as the environment as it stands now, not as the environment that produced these
figures. CopyKAT and dsb were both used here but are not attached when this report renders.

# 8. Caveats

1. **The run is superseded.** LK1 was re-sequenced as `260522_VH01624_461_222JLJVNX`.
2. **The pilot Cell Ranger output is not on the mount.** Step 02 reads
   `../results/cellranger/<run>/LK1-GEX/outs/filtered_feature_bc_matrix` and that path does
   not exist. Check cluster `/scratch` before assuming it was deleted.
3. **The pilot FASTQs are not on the mount either** - only the read count table survives.
4. **Dates are file modification times.** `scripts/` had no git history before 2026-07-31, so
   there is no commit trail for this period. Output dates are reliable; script dates only say
   when a file was last edited.
5. **Step 11 post-dates its own output** by about three weeks (`12.data_integration.R` was
   last edited 2026-06-09, the `data_integration/` figures are 05-06 to 05-20). Those figures
   were made by a version of that script that no longer exists.
6. **`13.CITE_residuals.R` is not included.** It belongs to this period but writes to
   `results/residual_cite_analysis`, which does not exist on the mount, so nothing can be
   shown for it. See the CITE benchmarking report.
7. **`numbat/CITE_HSC_Numbat_analysis/` is not pilot-era** (2026-06-15) even though it sits
   under this run's directory. It belongs with the CITE benchmarking report.


**Scripts dropped when collating this folder**, and why:

| Dropped | Reason |
|---|---|
| `backup/8.create_seurat_object.R` | Byte-identical to `scripts/6a.initial_QC.R` (step 02). |
| `backup/10b.annotate_DSB.R` | Byte-identical to `scripts/12.identify_cancer_cells.R` (step 09). |
| `backup/4.initial_QC.R` | Earlier, shorter version of step 02. |
| `backup/10c.compare_DSB_CLR.R` | Superseded by step 11, which covers the same comparison in more detail. |
| `backup/11a.example_cite_copykat.R` | Exploratory scratch, not part of the pipeline. |
| `add_cite_for_hiep.R` | One-off to hand CITE data to a colleague, no analytical result. |

Note that `scripts/12.identify_cancer_cells.R` does not identify cancer cells - it is the
BoneMarrowMap projection on the DSB object, byte-identical to `backup/10b.annotate_DSB.R`.
The name is misleading and it is called `09.annotate_CITE_DSB.R` here.
