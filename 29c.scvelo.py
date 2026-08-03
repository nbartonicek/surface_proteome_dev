#!/usr/bin/env python3
# 29c.scvelo.py
# RNA velocity analysis for LK2 (AML, human GRCh38).
# Prerequisites: run 29a (velocyto) and 29b (export) first.
# conda activate scvelo

import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scvelo as scv
from scipy import io
from scipy.sparse import csr_matrix

RESULTS = "/scratch/users/nbartonicek/projects/amgen/results_nf/260528_VH01624_464_222K7VKNX"
SAMPLE  = "LK2-GEX"
VEL_DIR = os.path.join(RESULTS, "29_velocity", SAMPLE)
LOOM    = os.path.join(RESULTS, "01_cellranger", SAMPLE, SAMPLE, "outs", "velocyto", f"{SAMPLE}.loom")
FIG_DIR = os.path.join(VEL_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

scv.settings.figdir = FIG_DIR
scv.settings.set_figure_params(dpi=120, frameon=False)

# ── 1. Build AnnData from exported counts ───────────────────────────────────
X = io.mmread(os.path.join(VEL_DIR, "counts.mtx")).T.tocsr()
meta = pd.read_csv(os.path.join(VEL_DIR, "metadata.csv"), index_col="barcode")
genes = pd.read_csv(os.path.join(VEL_DIR, "gene_names.csv"), header=None)[0].tolist()

adata = ad.AnnData(X=X, obs=meta, var=pd.DataFrame(index=genes))

# Projected UMAP from BoneMarrowMap
adata.obsm["X_umap"] = adata.obs[["UMAP1_projected", "UMAP2_projected"]].to_numpy()
adata.obs.drop(columns=["UMAP1_projected", "UMAP2_projected"], inplace=True)

# ── 2. Merge with velocyto loom (spliced / unspliced counts) ────────────────
ldata = scv.read(LOOM, cache=True)

# Velocyto barcodes are formatted as <sample>:<barcode>x — strip prefix/suffix
ldata.obs.index = [bc.split(":")[1].rstrip("x") + "-1" for bc in ldata.obs.index]
ldata.var_names_make_unique()

adata = scv.utils.merge(adata, ldata)

# ── 3. Pre-processing ────────────────────────────────────────────────────────
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=3000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# ── 4. Velocity (dynamical model — more accurate than stochastic for AML) ───
scv.tl.recover_dynamics(adata, n_jobs=4)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)

# ── 5. Velocity confidence & latent time ────────────────────────────────────
scv.tl.velocity_confidence(adata)
scv.tl.latent_time(adata)

# ── 6. Plots ─────────────────────────────────────────────────────────────────
color_by = ["predicted_CellType_Broad", "sample_name", "numbat_compartment",
            "copykat_call", "latent_time", "velocity_confidence"]

for c in color_by:
    if c not in adata.obs.columns:
        continue
    scv.pl.umap(adata, color=c, legend_loc="on data", title=c,
                save=f"_{c}.pdf")

scv.pl.velocity_embedding_stream(
    adata, basis="umap", color="predicted_CellType_Broad",
    legend_loc="right margin", title="RNA velocity (dynamical)",
    save="_velocity_stream.pdf"
)

scv.pl.velocity_embedding_grid(
    adata, basis="umap", color="predicted_CellType_Broad",
    scale=0.3, title="RNA velocity grid",
    save="_velocity_grid.pdf"
)

# Latent time coloured by compartment (tumor vs normal)
scv.pl.scatter(
    adata, color="latent_time", color_map="gnuplot",
    size=80, title="Latent time",
    save="_latent_time.pdf"
)

# ── 7. Top velocity genes ────────────────────────────────────────────────────
scv.tl.rank_velocity_genes(adata, groupby="predicted_CellType_Broad", min_corr=0.3)
top_genes = scv.DataFrame(adata.uns["rank_velocity_genes"]["names"]).head(5)
print("Top velocity genes per cell type:")
print(top_genes)
top_genes.to_csv(os.path.join(VEL_DIR, "top_velocity_genes.csv"))

# ── 8. Save ──────────────────────────────────────────────────────────────────
adata.write(os.path.join(VEL_DIR, "LK2_velocity.h5ad"))
print(f"Done. Results in {VEL_DIR}")
