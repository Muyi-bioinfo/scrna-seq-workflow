# scrna-seq-workflow

English | [简体中文](README_CN.md)

> A complete single-cell RNA-seq analysis pipeline — from FASTQ quantification to downstream analysis.
> Built as a self-learning project from a public tutorial + 10x Genomics official documentation (no real project data involved), adapted to Seurat v5.
> The architecture follows mainstream scRNA-seq projects on GitHub: one script per analysis stage + centralized config + a sequential runner (cf. [hossainlab/sc-workflow](https://github.com/hossainlab/sc-workflow)). A comparison with workflow engines such as Snakemake/Nextflow can be found in docs/01.

## Project story (how this repo came about)

1. **Learned from a tutorial**: started from a public Seurat v4 tutorial covering single-sample basics through multi-sample integration — 10 teaching scripts in total (the original learning track stays on my machine, not part of this repo)
2. **Adapted to a new major version**: the tutorial targets Seurat v4; my actual environment is Seurat v5 + Linux. Every script was rewritten line by line — the v4→v5 lessons are distilled into the script comments and the [key pitfalls](#key-pitfalls) table below
3. **Filled in the upstream gap**: the tutorial starts from the expression matrix; I added the cellranger quantification step myself (verified end-to-end on 10x tiny test data)
4. **Extended beyond the main line**: self-studied additional integration methods (RPCA / FastMNN / scVI / Liger) and gene-set scoring (AddModuleScore / UCell / AUCell)
5. **Re-engineered into a pipeline**: after surveying mainstream GitHub project structures, reorganized everything into staged scripts + centralized config + a resumable runner (script granularity follows [hossainlab/sc-workflow](https://github.com/hossainlab/sc-workflow): packaged by analysis stage, not one script per operation)

## Repository layout

```
scrna-seq-workflow/
├── README.md                    # Project overview (this file, English)
├── README_CN.md                 # Chinese version of the overview
├── environment.yaml             # mamba environment spec (rebuild the scrna env in one command)
├── run_pipeline.sh              # Orchestrator: sequential execution + resume + per-step logs
├── config/
│   ├── config.yaml              # Single-sample pipeline config (default)
│   ├── config.multi.yaml        # Multi-sample config (auto-loaded by --mode multi)
│   ├── sample_sheet.csv         # Sample sheet (sample, group, fastqs, matrix)
│   ├── cellranger.yaml          # cellranger tool-level config (decoupled from R; reusable across projects)
│   ├── trajectory_presets/      # Lineage presets for 07 trajectory (root + cell subset + rationale)
│   └── annotation_presets/      # Annotation marker presets for 03 (species folders + aliases + sources)
├── bash/                        # Upstream shell scripts (cellranger is a shell tool — language-agnostic)
│   └── 00_run_cellranger.sh     # cellranger count in batch (reads config + sample sheet)
├── R/                           # R pipeline scripts (staged; each stage can rerun or be submitted alone)
│   ├── 01_load_qc.R             # Load + QC (single/multi-sample branch on config automatically)
│   ├── 02_preprocess_cluster.R  # Preprocessing + batch correction + UMAP + clustering
│   ├── 03_annotate.R            # Cell annotation + optional subclustering (config switch)
│   ├── 04_multi_group_plot.R    # Multi-group visualization
│   ├── 05_diff_gsea.R           # Exploratory DE + GSEA (cell-level statistics)
│   ├── 06_pseudobulk_de.R       # Confirmatory DE: pseudobulk + DESeq2 (sample-aware)
│   ├── 07_gene_set_score.R      # Advanced: gene-set scoring
│   ├── 08_trajectory_analysis.R # Advanced: trajectory (monocle3)
│   └── 09_cell_communication.R  # Advanced: cell-cell communication (CellChat)
├── Python/                      # Python track (planned: parallel scanpy version, see Python/README.md)
├── docs/                        # Pipeline docs / cellranger guide / R package compilation notes
├── utils/                       # Shared utilities: utils.R (config/archiving/QC/preprocess/clustering)
│                                #   + save_fig.R (figure saving, cairo_pdf with CJK support)
├── data/                        # Data layering (symlinks/large files, not in git):
│   ├── raw/                     #   raw FASTQs (read-only, cellranger input)
│   ├── matrices/                #   public datasets / prequantified matrices (PBMC3k, GSE96583)
│   └── reference/               #   references:
│       ├── refdata/             #     cellranger reference genomes (tiny demo reference in place)
│       └── gene_sets/           #     MSigDB gene-set gmt files
└── output/<run>/                # Each run is fully self-contained (not in git):
                                 #   per-step result rds + figures/ + logs/ (per-step logs)
                                 #   + config_used.yaml (parameter snapshot of this run, for traceability)
```

## Pipeline overview

![Pipeline overview](docs/assets/pipeline_overview.png)

<!-- vector version + generator script: docs/assets/pipeline_overview.svg · docs/assets/make_pipeline_figure.py -->

## Results showcase

Real outputs of the run03 validation batch (GSE96583, 24,644 cells after QC).

**Multi-sample integration — before vs after harmony**

| Before integration | After harmony |
|---|---|
| ![before](docs/assets/showcase/integration_before.png) | ![after](docs/assets/showcase/integration_after_harmony.png) |

Before correction, STIM and CTRL cells separate inside shared clusters; after harmony, corresponding cell types from the two groups align, enabling cross-group comparison.

**Cell-type annotation (13 types)**

![Annotated UMAP](docs/assets/showcase/dimplot_annotated.png)

The cluster→cell-type mapping is machine-checked against marker evidence (`check_mapping_evidence()` in step 03), so a silent cluster-ID renumbering cannot mislabel cells.

**GSEA — IFN-β stimulation response (CD14 monocytes, STIM vs CTRL)**

![GSEA](docs/assets/showcase/gsea_hallmark_cd14mono.png)

Hallmark GSEA ranks interferon-α response as the top activated pathway (GeneRatio ≈ 0.8), with interferon-γ response also strongly enriched — consistent with the IFN-β stimulation design of GSE96583.

**Trajectory — T-cell lineage (monocle3)**

![Pseudotime](docs/assets/showcase/trajectory_pseudotime.png)

Pseudotime on the T-lineage subset (12,828 cells, root = CD4 Naive T, lineage preset from `config/trajectory_presets/`); 0 Inf values — subsetting to a connected lineage avoids the disconnected-graph Inf problem.

## Quick start

```bash
# 1. Activate the environment (create it first: mamba env create -f environment.yaml)
mamba activate scrna        # R 4.5.3 + Seurat 5.5.1 + presto acceleration

# 2. Upstream quantification (starting from FASTQ; requires the reference genome
#    and the fastqs column of the sample sheet)
bash bash/00_run_cellranger.sh
#    Tool-level parameters: config/cellranger.yaml (reference genome / hardware resources)
#    Reference download & FASTQ layout: docs/02_cellranger_guide.md (reference ~11GB)
#    ⚠️ No FASTQs? Skip this step — the example datasets ship as ready-made matrices
#    (see Example data below) and the pipeline starts at step 3

# 3. Single-sample pipeline (PBMC3k example, default config: config/config.yaml)
bash run_pipeline.sh

# 4. Multi-sample pipeline (STIM/CTRL; auto-loads config/config.multi.yaml)
bash run_pipeline.sh --mode multi

# 5. Append the advanced modules
bash run_pipeline.sh --mode multi --with-advanced

# 6. Flexible execution (same-run resume: steps with existing outputs are skipped automatically)
bash run_pipeline.sh --only 02            # rerun preprocessing+clustering only (common after changing resolution)
bash run_pipeline.sh --from 02 --to 03    # run steps 02–03
bash run_pipeline.sh --force              # force-rerun every step of the current run
Rscript R/02_preprocess_cluster.R   # run a single stage directly (e.g. submit to an HPC node)

# 7. Run management (the `batch` field in the config): every run archives its results
#    and logs under its own directory
#    Same run name = resume from breakpoint; a new run name (e.g. batch: run02) = fresh
#    run, previous results are kept for comparison
#    Each run auto-saves a config_used.yaml parameter snapshot, for full traceability
#    output/run01/ (results + logs/ + config_used.yaml) ← one complete run
```

### Example data

| Dataset | Source | Used by |
|---|---|---|
| PBMC 3k (hg19) | 10x Genomics website | single-sample pipeline |
| STIM/CTRL matrices | GSE96583 (IFN-β stimulated PBMC) | multi-sample integration + DE + GSEA + pseudobulk DE |
| pbmc_1k_v3 FASTQ + GRCh38-2024-A reference | 10x Genomics website | real cellranger practice |

Download locations and directory layout: data-sources section at the end of [docs/01_pipeline_overview.md](docs/01_pipeline_overview.md).
Note: the GSE96583 GEO submission has mitochondrial genes filtered out, so `percent.mt` degrades to 0 for that dataset (a real cellranger run is unaffected).

### Sample sheet (multi-sample mode)

`config/sample_sheet.csv` is the entry point for plugging your own data into the multi-sample pipeline — every downstream comparison starts from it:

| Column | Consumed by | Role |
|---|---|---|
| `sample` | 01 (orig.ident), 06 (default `sample_id`) | sample name = the biological-replicate unit of pseudobulk DE |
| `group` | 02 integration, 04 group plots, 05 DE+GSEA, 06 (`condition_col`) | the grouping unit of **all downstream comparisons** — values must match `pseudobulk$test_group` / `reference_group` in the config exactly |
| `fastqs` | 00 cellranger | upstream quantification input (either this or `matrix`) |
| `matrix` | 01 | 10x matrix directory/h5 (either this or `fastqs`) |

Two things to keep in mind:

- A wrong `group` value does **not** raise an error — it silently produces wrong comparisons. After switching datasets, check the `group` column against `test_group`/`reference_group` in the config first.
- With pooled matrices (each row is a mix of many donors, e.g. GSE96583), the `sample` column no longer represents biological replicates — configure `multi$cell_metadata` so step 01 joins donors by (group, barcode) into `sample_id`.

## Outputs

Every run archives under `output/<batch>/<step>/` — result `.rds` + `figures/` + a `.done` resume marker; each run also keeps per-step `logs/` and a `config_used.yaml` parameter snapshot. What each step produces (real files from the run03 validation batch):

| Step | Key results | Representative figures |
|---|---|---|
| 00 cellranger | `<sample>/outs/` filtered matrix + `web_summary.html` QC report | — |
| 01 load + QC | `seurat_qc.rds` | `qc_violin_before` |
| 02 preprocess + cluster | `seurat_clustered.rds` | `dimplot_clusters`, `umap_before_integration`, `dimplot_by_group_harmony`, `elbow` |
| 03 annotate (+ subcluster) | `seurat_annotated.rds`, `all_markers_table.rds`, `top10_markers.csv`, `annotation_evidence.csv` | `dimplot_annotated`, `dotplot_markers_by_cluster`, `featureplot_celltypes/` |
| 04 multi-group plot | figures only | `barplot_cell_proportion`, `dotplot_grouped`, `heatmap_top5_markers` |
| 05 DE + GSEA | `split_markers.rds`, `diff_stim_vs_ctrl.rds` | `compareCluster_dotplot`, `gsea_kegg_dotplot`, `gsea_hallmark_dotplot` |
| 06 pseudobulk DE | `summary_degs.csv`, `by_celltype/<type>/` (`deseq2_results.tsv` + `significant_degs.tsv`) | `figures/<type>/` (`volcano`, `MAplot`, `PCA`, `heatmap`), `summary_05_vs_06` |
| 07 gene-set scoring | `seurat_scored.rds` | `featureplot_<pathway>` (one per signature) |
| 08 trajectory | `trajectory_cds.rds`, `trajectory_resolved.yaml` (parameter snapshot) | `trajectory_pseudotime`, `trajectory_celltype`, `trajectory_by_group` |
| 09 cell communication | `cellchat.rds` (single) / `cellchat_list.rds` (multi) | `compare_interaction_counts`, `diff_network_all`, `ranknet_pathway_comparison` |

## Key pitfalls

| Pitfall | Root cause | Fix |
|----|------|------|
| `FindAllMarkers` fails after `merge()` | Seurat v5 merge creates split layers | `JoinLayers()` |
| `JoinLayers` errors after loading a v4-format rds | the Assay class of a v4 object is "Assay" | `scobj[["RNA"]] <- as(scobj[["RNA"]], "Assay5")` |
| Cell Ranger 9 complains about a missing `--create-bam` | the flag is mandatory since Cell Ranger 9 | pass `--create-bam=false` explicitly |
| scale.data grows explosively | scaling all genes = a dense genes × cells matrix | scale HVGs only (default); clear the slot before saving |
| Top markers lack the canonical markers | the v5 logFC formula change inflates logFC for lowly expressed genes | prefilter with `annotate$marker_pct` before ranking (default 0.5) |
| GSEA results not reproducible | GSEA estimates p-values by random permutation | `set.seed(123)` |
| Annotation map mixed up between single/multi-sample runs | different datasets have different cluster counts and cell types | separate config profiles: config.yaml / config.multi.yaml |
| cluster_map maps to the wrong clusters after a rerun (e.g. CD14 Mono swapped with naive T) | Louvain cluster IDs are arbitrary labels from community-discovery order; any upstream change re-shuffles all IDs while the clusters themselves stay the same (run02 incident, 2026-08-21) | ① machine check via `check_mapping_evidence()` in 03: a wrong mapping raises a warning + an evidence table `annotation_evidence.csv` ② the clustering chain is fully seeded (`set.seed` inside `cluster_cells()` + `random.seed=42`): same input ⇒ same IDs ③ verification reruns use a new run name, keeping the old run for comparison |
| irlba segfaults (monocle3 `preprocess_cds`, and step 02 `RunPCA` after the 2026-09-01 environment refresh) | irlba + multithreaded OpenBLAS race (r-irlba 2.3.7 + OpenBLAS 0.3.33) | run_pipeline.sh pins `OPENBLAS_NUM_THREADS=1` for every R step (Sys.setenv inside R does not work — it must be set at shell level) |
| `order_cells` errors with "root_pr_nodes or root_cells must be provided" | non-interactive Rscript cannot open the interactive root picker (only RStudio does) | lineage presets specify the root explicitly (config/trajectory_presets/) |
| ~30% Inf pseudotime when running all cell types together | no continuous trajectory exists across heterogeneous cell types — the graph is disconnected | subset to a lineage (e.g. T cells only); verified Inf = 0 on run02 |
| FeaturePlot inconsistent with the DimPlot layout | `umap_naive` hijacks the default dimrec (DefaultDimReduc matches by name), so FeaturePlot without an explicit reduction draws on the pre-correction embedding | always pass an explicit `reduction = "umap"` |
| Group-level DE looks inflated (thousands of "significant" STIM-vs-CTRL genes) | `FindAllMarkers`/`FindMarkers` treat each cell as an independent observation, but cells from the same donor are correlated — pseudoreplication (GSE96583 is 8 paired donors) | confirm with `06_pseudobulk_de.R`: aggregate raw counts per donor × cell type, DESeq2 with the donor as blocking factor; the same-threshold comparison lands in `summary_degs.csv` |

## Repository scope & data availability

- Only **code and docs** are committed (`R/`, `bash/`, `Python/`, `utils/`, `config/`, `docs/`, `run_pipeline.sh`, `environment.yaml`, `LICENSE`, READMEs)
- Data and results (`data/`, `output/`, symlinks) are git-ignored
- Data sources and how to obtain them: end of [docs/01_pipeline_overview.md](docs/01_pipeline_overview.md)

## Environment

- Full dependency list in [environment.yaml](environment.yaml); rebuild with one command:
  `mamba env create -f environment.yaml`
- Core versions: R 4.5.3, Seurat 5.5.1, SeuratObject 5.4.0, monocle3 1.4.27, harmony, scDblFinder, yaml
- ⚠️ 4 packages must be installed from GitHub (not available via conda, see notes at the top of the yaml): SeuratWrappers, CellChat (required) + presto, scCustomize (optional)
- cellranger 9.0.1 (standalone software, not an R package; see docs/02_cellranger_guide.md)

### Advanced-module dependencies (check before running 08/09)

```bash
# 08 trajectory: monocle3 installed and fully validated (run02 batch, T-cell lineage preset)
#   install: mamba install -c conda-forge -c bioconda r-monocle3
#   note: requires single-threaded BLAS (handled automatically by run_pipeline.sh; see Key pitfalls)

# 09 cell-cell communication: CellChat 2.2.0.9001 installed (compiled from GitHub source; removed from CRAN)
#   install (for reference): devtools::install_github("jinworks/CellChat")
#   compilation notes: docs/03_r_package_compile.md (conda env + R 4.5 headers)
```

## References & acknowledgments

- [hossainlab/sc-workflow](https://github.com/hossainlab/sc-workflow) — the architecture reference: script granularity (one script per analysis stage) + a sequential runner
- [Seurat PBMC 3k tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) — the learning starting point (v4); [v5 vignettes](https://satijalab.org/seurat/) for the migration
- Core tools: [Seurat](https://github.com/satijalab/seurat) · [harmony](https://github.com/immunogenomics/harmony) · [monocle3](https://github.com/cole-trapnell-lab/monocle3) · [CellChat](https://github.com/jinworks/CellChat)

## License

MIT — see [LICENSE](LICENSE). This project is for educational and research purposes only.
