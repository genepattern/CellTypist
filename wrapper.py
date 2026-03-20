#!/usr/bin/env python3
"""
Wrapper script for the celltypist GenePattern module.
Performs automated cell type annotation using pre-trained or custom
CellTypist models on single-cell RNA-seq expression data.

Supported input formats: AnnData (.h5ad), CSV, TSV, 10x MTX
Output: predicted labels, probability matrix, decision matrix,
        annotated AnnData, and optional UMAP/dot-plot figures.
"""

import argparse
import sys
import os
import logging
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
VALID_INPUT_EXTENSIONS = {".h5ad", ".csv", ".tsv", ".txt", ".tab", ".mtx"}
VALID_PREDICTION_MODES = {"best match", "prob match"}
DEFAULT_BUILTIN_MODEL = "Immune_All_Low.pkl"
DEFAULT_OUTPUT_PREFIX = "celltypist_output"
DEFAULT_PREDICTION_MODE = "best match"
DEFAULT_PROB_THRESHOLD = 0.5
DEFAULT_MIN_PROPORTION = 0.0


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
def setup_logging() -> None:
    """Configure root logger to write INFO+ messages to stdout."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
        stream=sys.stdout,
    )


# ---------------------------------------------------------------------------
# Argument Parsing
# ---------------------------------------------------------------------------
def parse_arguments() -> argparse.Namespace:
    """Parse GenePattern-style command-line arguments.

    Flag names use dots (e.g. --input.file) to match the GenePattern manifest
    commandLine / prefix_when_specified values exactly. The ``dest`` attribute
    maps each flag to a Python-safe underscore identifier.
    """
    parser = argparse.ArgumentParser(
        description=(
            "CellTypist GenePattern wrapper — automated cell type annotation "
            "for single-cell RNA-seq data."
        )
    )

    # ---- Required -----------------------------------------------------------
    parser.add_argument(
        "--input.file",
        dest="input_file",
        required=True,
        metavar="INPUT_FILE",
        help=(
            "Input expression data file. Supported formats: AnnData (.h5ad), "
            "comma-separated (.csv), tab-separated (.tsv/.txt/.tab), or 10x "
            "Genomics sparse matrix (.mtx). For .h5ad files the expression "
            "matrix must be log1p-normalised to 10,000 counts per cell. For "
            "CSV/TSV/MTX files raw counts are accepted and normalised "
            "internally. Matrix must be cells (rows) x genes (columns) unless "
            "--transpose.input is set to true."
        ),
    )

    # ---- MTX companions -----------------------------------------------------
    parser.add_argument(
        "--gene.file",
        dest="gene_file",
        default=None,
        metavar="GENE_FILE",
        help=(
            "Path to gene names file. Required when --input.file is in 10x "
            "MTX format. One gene name per line corresponding to matrix rows."
        ),
    )
    parser.add_argument(
        "--cell.file",
        dest="cell_file",
        default=None,
        metavar="CELL_FILE",
        help=(
            "Path to cell barcodes file. Required when --input.file is in 10x "
            "MTX format. One barcode per line corresponding to matrix columns."
        ),
    )

    # ---- Model selection ----------------------------------------------------
    parser.add_argument(
        "--builtin.model",
        dest="builtin_model",
        default=DEFAULT_BUILTIN_MODEL,
        metavar="MODEL_NAME",
        help=(
            "Built-in pre-trained CellTypist model. Ignored when "
            "--custom.model.file is supplied. Examples: "
            "Immune_All_Low.pkl (~31 broad immune types, default), "
            "Immune_All_High.pkl (~98 fine-grained immune subtypes), "
            "or tissue-specific models (Gut_High, Lung_High, Brain_High, etc.)."
        ),
    )
    parser.add_argument(
        "--custom.model.file",
        dest="custom_model_file",
        default=None,
        metavar="CUSTOM_MODEL_FILE",
        help=(
            "Optional user-supplied CellTypist model (.pkl). When provided "
            "this takes precedence over --builtin.model. Custom models can be "
            "trained with celltypist.train() on any labelled reference dataset."
        ),
    )

    # ---- Matrix orientation -------------------------------------------------
    parser.add_argument(
        "--transpose.input",
        dest="transpose_input",
        default="false",
        choices=["true", "false"],
        metavar="TRANSPOSE_INPUT",
        help=(
            "Set to true if the input matrix is genes (rows) x cells (columns) "
            "rather than the default cells x genes orientation."
        ),
    )

    # ---- Prediction settings ------------------------------------------------
    parser.add_argument(
        "--prediction.mode",
        dest="prediction_mode",
        default=DEFAULT_PREDICTION_MODE,
        metavar="PREDICTION_MODE",
        help=(
            "Prediction mode for cell type assignment. "
            "'best match' assigns each cell exactly one label (highest "
            "probability). 'prob match' assigns one or more pipe-separated "
            "labels for every type whose probability exceeds "
            "--probability.threshold; cells below threshold are 'Unassigned'."
        ),
    )
    parser.add_argument(
        "--probability.threshold",
        dest="probability_threshold",
        type=float,
        default=DEFAULT_PROB_THRESHOLD,
        metavar="PROBABILITY_THRESHOLD",
        help=(
            "Probability cutoff used in 'prob match' mode (0.0–1.0). "
            "Cell types below this threshold are excluded from multi-label "
            "assignment. Has no effect in 'best match' mode. Default: 0.5."
        ),
    )

    # ---- Majority voting ----------------------------------------------------
    parser.add_argument(
        "--majority.voting",
        dest="majority_voting",
        default="false",
        choices=["true", "false"],
        metavar="MAJORITY_VOTING",
        help=(
            "Enable majority-vote post-processing. Leiden over-clustering "
            "(resolution=5) is applied and each sub-cluster receives the "
            "plurality cell type. Sub-clusters where the plurality comprises "
            "less than --min.proportion of cells are labelled 'Heterogeneous'."
        ),
    )
    parser.add_argument(
        "--min.proportion",
        dest="min_proportion",
        type=float,
        default=DEFAULT_MIN_PROPORTION,
        metavar="MIN_PROPORTION",
        help=(
            "Minimum fraction of the dominant cell type within a sub-cluster "
            "for majority voting (0.0–1.0). Sub-clusters below this fraction "
            "are labelled 'Heterogeneous'. Only used when --majority.voting "
            "is true. Default: 0.0 (no minimum)."
        ),
    )

    # ---- Output settings ----------------------------------------------------
    parser.add_argument(
        "--output.prefix",
        dest="output_prefix",
        default=DEFAULT_OUTPUT_PREFIX,
        metavar="OUTPUT_PREFIX",
        help=(
            "Prefix string for all output file names. Outputs produced: "
            "{prefix}_predicted_labels.csv, {prefix}_probability_matrix.csv, "
            "{prefix}_decision_matrix.csv, {prefix}_annotated.h5ad. "
            "Default: 'celltypist_output'."
        ),
    )
    parser.add_argument(
        "--plot.results",
        dest="plot_results",
        default="false",
        choices=["true", "false"],
        metavar="PLOT_RESULTS",
        help=(
            "Generate UMAP and dot-plot visualisations. Plots are saved as "
            "PDF and PNG. If UMAP coordinates are absent they are computed "
            "automatically using scanpy."
        ),
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Input Validation
# ---------------------------------------------------------------------------
def validate_inputs(args: argparse.Namespace) -> None:
    """Validate all user-supplied parameters; exit with code 1 on failure."""

    # --- input.file ----------------------------------------------------------
    input_path = Path(args.input_file)
    if not input_path.is_file():
        logging.error("Input file not found: %s", args.input_file)
        sys.exit(1)

    ext = input_path.suffix.lower()
    if ext not in VALID_INPUT_EXTENSIONS:
        logging.error(
            "Unsupported input file extension '%s'. Supported: %s",
            ext,
            ", ".join(sorted(VALID_INPUT_EXTENSIONS)),
        )
        sys.exit(1)

    # --- MTX companion files -------------------------------------------------
    if ext == ".mtx":
        if not args.gene_file:
            logging.error(
                "--gene.file is required when --input.file is in MTX format."
            )
            sys.exit(1)
        if not args.cell_file:
            logging.error(
                "--cell.file is required when --input.file is in MTX format."
            )
            sys.exit(1)
        if not Path(args.gene_file).is_file():
            logging.error("Gene names file not found: %s", args.gene_file)
            sys.exit(1)
        if not Path(args.cell_file).is_file():
            logging.error("Cell barcodes file not found: %s", args.cell_file)
            sys.exit(1)

    # --- custom model --------------------------------------------------------
    if args.custom_model_file and not Path(args.custom_model_file).is_file():
        logging.error(
            "Custom model file not found: %s", args.custom_model_file
        )
        sys.exit(1)

    # --- prediction mode -----------------------------------------------------
    if args.prediction_mode not in VALID_PREDICTION_MODES:
        logging.error(
            "Invalid --prediction.mode '%s'. Must be one of: %s",
            args.prediction_mode,
            ", ".join(sorted(VALID_PREDICTION_MODES)),
        )
        sys.exit(1)

    # --- probability threshold -----------------------------------------------
    if not (0.0 <= args.probability_threshold <= 1.0):
        logging.error(
            "--probability.threshold must be between 0.0 and 1.0 (got %s).",
            args.probability_threshold,
        )
        sys.exit(1)

    # --- min proportion ------------------------------------------------------
    if not (0.0 <= args.min_proportion <= 1.0):
        logging.error(
            "--min.proportion must be between 0.0 and 1.0 (got %s).",
            args.min_proportion,
        )
        sys.exit(1)

    logging.info("Input validation passed.")


# ---------------------------------------------------------------------------
# Data Loading
# ---------------------------------------------------------------------------
def load_input_data(args: argparse.Namespace):
    """Load the input expression matrix into an AnnData object."""
    import anndata
    import pandas as pd
    import scanpy as sc

    input_path = Path(args.input_file)
    ext = input_path.suffix.lower()
    transpose = args.transpose_input.lower() == "true"

    logging.info("Loading input data: %s", input_path)

    if ext == ".h5ad":
        adata = anndata.read_h5ad(str(input_path))

    elif ext in {".csv", ".tsv", ".txt", ".tab"}:
        sep = "," if ext == ".csv" else "\t"
        df = pd.read_csv(str(input_path), sep=sep, index_col=0)
        if transpose:
            logging.info("Transposing input matrix (genes x cells → cells x genes).")
            df = df.T
        adata = anndata.AnnData(
            X=df.values,
            obs=pd.DataFrame(index=df.index.astype(str)),
            var=pd.DataFrame(index=df.columns.astype(str)),
        )

    else:  # .mtx
        mat = sc.read_mtx(str(input_path))
        genes = pd.read_csv(args.gene_file, header=None, sep="\t")
        cells = pd.read_csv(args.cell_file, header=None, sep="\t")
        # MTX files are genes × cells; transpose to cells × genes
        adata = anndata.AnnData(
            X=mat.X.T,
            obs=pd.DataFrame(index=cells.iloc[:, 0].astype(str)),
            var=pd.DataFrame(index=genes.iloc[:, 0].astype(str)),
        )
        if transpose:
            logging.info(
                "Transposing MTX-derived matrix as requested by --transpose.input."
            )
            adata = adata.T

    logging.info(
        "Data loaded successfully: %d cells × %d genes.",
        adata.n_obs,
        adata.n_vars,
    )
    return adata


# ---------------------------------------------------------------------------
# Model Loading
# ---------------------------------------------------------------------------
def load_model(args: argparse.Namespace):
    """Return a CellTypist Model object from custom file or built-in name."""
    from celltypist import models

    if args.custom_model_file:
        logging.info("Loading custom model: %s", args.custom_model_file)
        return models.Model.load(args.custom_model_file)

    model_name = args.builtin_model or DEFAULT_BUILTIN_MODEL
    logging.info("Downloading / loading built-in model: %s", model_name)
    models.download_models(force_update=False, model=model_name)
    return models.Model.load(model=model_name)


# ---------------------------------------------------------------------------
# Majority-Voting Preprocessing
# ---------------------------------------------------------------------------
def prepare_over_clustering(adata):
    """Compute PCA → neighbors → Leiden (res=5) for majority-vote smoothing."""
    import scanpy as sc

    logging.info(
        "Preparing over-clustering for majority voting (Leiden resolution=5)."
    )

    if "X_pca" not in adata.obsm:
        logging.info("PCA coordinates not found — computing PCA.")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(
            adata, min_mean=0.0125, max_mean=3, min_disp=0.5
        )
        sc.pp.pca(adata, svd_solver="arpack")

    if "neighbors" not in adata.uns:
        logging.info("Neighbor graph not found — computing neighbors.")
        sc.pp.neighbors(adata)

    sc.tl.leiden(adata, resolution=5, key_added="over_clustering")
    logging.info("Over-clustering complete.")
    return "over_clustering"


# ---------------------------------------------------------------------------
# Output Saving
# ---------------------------------------------------------------------------
def save_outputs(predictions, out_prefix: str) -> None:
    """Persist all tabular and AnnData outputs."""
    labels_path = f"{out_prefix}_predicted_labels.csv"
    prob_path = f"{out_prefix}_probability_matrix.csv"
    decision_path = f"{out_prefix}_decision_matrix.csv"
    h5ad_path = f"{out_prefix}_annotated.h5ad"

    predictions.predicted_labels.to_csv(labels_path)
    logging.info("Predicted labels saved: %s", labels_path)

    predictions.probability_matrix.to_csv(prob_path)
    logging.info("Probability matrix saved: %s", prob_path)

    predictions.decision_matrix.to_csv(decision_path)
    logging.info("Decision matrix saved: %s", decision_path)

    adata_out = predictions.to_adata()
    adata_out.write_h5ad(h5ad_path)
    logging.info("Annotated AnnData saved: %s", h5ad_path)

    return adata_out


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def generate_plots(adata_out, predictions, out_prefix: str) -> None:
    """Produce UMAP and dot-plot figures (non-fatal on failure)."""
    import scanpy as sc
    import celltypist

    logging.info("Generating visualisation plots.")

    try:
        if "X_umap" not in adata_out.obsm:
            logging.info("UMAP coordinates not found — computing UMAP.")
            if "X_pca" not in adata_out.obsm:
                sc.pp.normalize_total(adata_out, target_sum=1e4)
                sc.pp.log1p(adata_out)
                sc.pp.highly_variable_genes(
                    adata_out, min_mean=0.0125, max_mean=3, min_disp=0.5
                )
                sc.pp.pca(adata_out, svd_solver="arpack")
            if "neighbors" not in adata_out.uns:
                sc.pp.neighbors(adata_out)
            sc.tl.umap(adata_out)

        celltypist.plot.predicted_cell_types(
            adata_out,
            predictions,
            basis="umap",
            show=False,
            return_fig=False,
            prefix=out_prefix,
            format=["pdf", "png"],
        )
        logging.info("Plots saved with prefix: %s", out_prefix)

    except Exception as exc:  # non-fatal — annotation outputs already written
        logging.warning(
            "Plot generation encountered an error (non-fatal): %s", exc
        )


# ---------------------------------------------------------------------------
# Core Execution
# ---------------------------------------------------------------------------
def run_celltypist(args: argparse.Namespace) -> None:
    """Orchestrate data loading, annotation, and output generation."""
    import celltypist

    out_prefix = args.output_prefix or DEFAULT_OUTPUT_PREFIX
    use_majority_voting = args.majority_voting.lower() == "true"

    # 1. Load data
    adata = load_input_data(args)

    # 2. Load model
    model = load_model(args)

    # 3. Prepare over-clustering if majority voting is requested
    over_clustering = None
    if use_majority_voting:
        over_clustering = prepare_over_clustering(adata)

    # 4. Annotate
    logging.info(
        "Running CellTypist annotation (mode: '%s', majority_voting: %s).",
        args.prediction_mode,
        use_majority_voting,
    )
    predictions = celltypist.annotate(
        adata,
        model=model,
        majority_voting=use_majority_voting,
        over_clustering=over_clustering,
        min_prop=args.min_proportion,
        p_thres=args.probability_threshold,
        mode=args.prediction_mode,
    )

    # 5. Save tabular outputs and annotated AnnData
    adata_out = save_outputs(predictions, out_prefix)

    # 6. Optional plots
    if args.plot_results.lower() == "true":
        generate_plots(adata_out, predictions, out_prefix)

    logging.info("CellTypist annotation completed successfully.")


# ---------------------------------------------------------------------------
# Entry Point
# ---------------------------------------------------------------------------
def main() -> None:
    setup_logging()
    args = parse_arguments()
    validate_inputs(args)

    try:
        run_celltypist(args)
    except Exception as exc:
        logging.error("CellTypist execution failed: %s", exc)
        sys.exit(2)

    sys.exit(0)


if __name__ == "__main__":
    main()
