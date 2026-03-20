# Dockerfile for celltypist GenePattern Module
# Automated cell type annotation for single-cell RNA-seq data using CellTypist.
#
# Resource requirements: 4 CPU cores, 16 GB memory
# Supported input formats: .h5ad, .csv, .tsv, .txt, .tab, .mtx
# Wrapper script: wrapper.py (Python 3)
#
# Parsed imports from wrapper.py:
#   Standard library: argparse, sys, os, logging, pathlib
#   Third-party: anndata, pandas, scanpy, celltypist, leidenalg

FROM python:3.11-slim

# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------
LABEL maintainer="GenePattern"
LABEL module.name="celltypist"
LABEL module.version="1.0.0"
LABEL module.language="python"
LABEL description="Automated cell type annotation for scRNA-seq using CellTypist"

# ---------------------------------------------------------------------------
# Working directory
# ---------------------------------------------------------------------------
WORKDIR /module

# ---------------------------------------------------------------------------
# System dependencies
#
# libhdf5-dev  — required by h5py (AnnData .h5ad I/O)
# libopenblas-dev — fast linear algebra for numpy/scipy
# gcc / g++ / python3-dev — compile native extensions (e.g. igraph, leidenalg)
# pkg-config   — used by several C-extension builds
# ca-certificates, wget, curl — network access for model downloads at run-time
# git          — occasionally needed by pip source installs
# ---------------------------------------------------------------------------
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        wget \
        git \
        gcc \
        g++ \
        python3-dev \
        libhdf5-dev \
        pkg-config \
        libopenblas-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# Python dependencies
#
# celltypist — core annotation library (pulls in numpy, scipy, scikit-learn)
# anndata    — AnnData object I/O (.h5ad support)
# scanpy     — single-cell analysis (PCA, neighbors, UMAP, Leiden pre-processing)
# pandas     — tabular data (CSV / TSV loading)
# leidenalg  — Leiden clustering (required by scanpy majority-voting workflow)
# python-igraph — graph library dependency of leidenalg
# matplotlib — figure generation for optional UMAP / dot-plot output
# ---------------------------------------------------------------------------
RUN pip install --no-cache-dir \
        "celltypist>=1.6.0" \
        "anndata>=0.10.0" \
        "scanpy>=1.9.0" \
        "pandas>=1.5.0" \
        "leidenalg>=0.10.0" \
        "python-igraph>=0.11.0" \
        "matplotlib>=3.7.0" \
        "scipy>=1.10.0" \
        "numpy>=1.24.0" \
        "h5py>=3.8.0"

# ---------------------------------------------------------------------------
# Copy module files
# ---------------------------------------------------------------------------
COPY wrapper.py /module/

# Ensure the wrapper is executable
RUN chmod +x /module/wrapper.py

# ---------------------------------------------------------------------------
# Environment variables
# ---------------------------------------------------------------------------
ENV MODULE_CPU_CORES=4
ENV MODULE_MEMORY=16GB
ENV MODULE_NAME=celltypist
# Disable interactive matplotlib backend (no display available in container)
ENV MPLBACKEND=Agg
# Prevent Python from writing .pyc files and buffering stdout
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# ---------------------------------------------------------------------------
# Default command — drop to bash so the GenePattern framework can invoke
# the wrapper with its own argument string
# ---------------------------------------------------------------------------
CMD ["/bin/bash"]
