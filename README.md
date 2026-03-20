# CellTypist (v1)

**Description**: Automated cell type annotation for single-cell RNA-seq data using pre-trained machine learning models from the CellTypist framework.
**Authors**: Domínguez Conde et al.; Wellcome Sanger Institute / Broad Institute of MIT and Harvard
**Contact**: https://groups.google.com/forum/#!forum/genepattern-help
**Algorithm Version**: CellTypist v1.6.x

## Summary

CellTypist is an automated cell type annotation tool designed for single-cell RNA sequencing (scRNA-seq) data. Given a query expression matrix and a pre-trained reference model, CellTypist assigns a cell type label — drawn from the model's vocabulary — to each individual cell. This enables rapid, reproducible cell type identification across diverse tissues, organisms, and experimental contexts without requiring extensive manual curation.

### How Does CellTypist Work?

At its core, CellTypist uses a logistic regression classifier trained on large, curated reference atlases (e.g., the Human Cell Atlas). For each cell in the query dataset, the model computes a probability score for every cell type in its vocabulary. Depending on the prediction mode selected, a cell is either assigned:
- The single **best-matching** cell type (highest probability), or
- **Multiple labels** (all cell types exceeding a user-defined probability threshold), which is useful for capturing transitional or ambiguous cell states.

An optional **majority voting** post-processing step further refines predictions by leveraging transcriptional similarity among neighboring cells. In this mode, CellTypist performs over-clustering of the query data and assigns the plurality cell type within each sub-cluster as the consensus label — effectively smoothing out noisy single-cell predictions.

### Where Does CellTypist Fit?

CellTypist is typically applied after initial preprocessing, dimensionality reduction, and clustering steps in a standard scRNA-seq workflow (e.g., using Scanpy or Seurat). It is particularly well-suited for:

- **Immune profiling**: Annotating complex immune landscapes in blood, tissue biopsies, and tumors using high-resolution immune models.
- **Cross-dataset integration**: Rapidly annotating newly generated data using atlas-scale reference models.
- **Large-scale studies**: Efficiently processing hundreds of thousands of cells in a single run.
- **Multi-tissue annotation**: Leveraging tissue-specific models (gut, lung, brain, skin, thymus, fetal) for targeted annotation.

For users new to scRNA-seq cell type annotation, CellTypist provides an accessible entry point: simply supply your expression matrix and choose a built-in model — no programming experience is required when using the GenePattern interface.

## References

1. Domínguez Conde C, Xu C, Jarvis LB, et al. **Cross-tissue immune cell analysis reveals tissue-specific features in humans.** *Science*. 2022;376(6594):eabl5197. doi:[10.1126/science.abl5197](https://doi.org/10.1126/science.abl5197)

2. Xu C, Lopez R, Mereu E, et al. **Probabilistic harmonization and annotation of single‐cell transcriptomics data with deep generative models.** *Molecular Systems Biology*. 2021;17(1):e9620.

3. Wolf FA, Angerer P, Theis FJ. **SCANPY: large-scale single-cell gene expression data analysis.** *Genome Biology*. 2018;19:15. doi:[10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)

4. CellTypist documentation: [https://www.celltypist.org](https://www.celltypist.org)

## Source Links

* [CellTypist GitHub Repository](https://github.com/Teichlab/celltypist)
* [CellTypist Docker Image](https://hub.docker.com/r/genepattern/celltypist)
* [CellTypist Model Zoo](https://www.celltypist.org/models)

## Parameters

| Name | Description | Default Value |
| :--- | :--- | :--- |
| input_file * | The input single-cell RNA-seq expression data file. Accepted formats: AnnData (`.h5ad`), CSV (`.csv`), TSV (`.tsv`, `.txt`, `.tab`), or 10x Genomics sparse matrix (`.mtx`). | — |
| gene_file | Gene names file; required only when input is in 10x MTX format. One gene symbol per line. | None |
| cell_file | Cell barcodes file; required only when input is in 10x MTX format. One barcode per line. | None |
| model | Built-in pre-trained CellTypist model to use. Ignored if `custom_model_file` is provided. | `Immune_All_Low.pkl` |
| custom_model_file | Optional user-supplied custom CellTypist model (`.pkl`). Overrides `model` when provided. | None |
| transpose_input | Set to `true` if input matrix is genes (rows) × cells (columns) instead of cells × genes. | `false` |
| prediction_mode | Assignment mode: `best match` (one label per cell) or `prob match` (multiple labels above threshold). | `best match` |
| probability_threshold | Probability cutoff for `prob match` mode. Cell types below this score are excluded. Range: 0.0–1.0. | `0.5` |
| majority_voting | Enable majority-vote post-processing to smooth predictions using neighborhood sub-clusters. | `false` |
| min_proportion | Minimum fraction of dominant cell type in a sub-cluster for majority voting. Range: 0.0–1.0. | `0` |
| output_prefix | Prefix string for all output file names. | `celltypist_output` |
| plot_results | Generate UMAP and dot plot visualizations of annotation results. | `false` |

\* required

## Input Files

1. **input_file**
   The primary input data file containing single-cell gene expression measurements. CellTypist accepts several formats:

   - **AnnData (`.h5ad`)**: The recommended format. The expression matrix (`adata.X`) should be **log1p-normalized to 10,000 counts per cell** (i.e., log1p(CPM/10)). This is the standard normalization produced by `scanpy.pp.normalize_total(target_sum=1e4)` followed by `scanpy.pp.log1p()`. If a pre-computed neighborhood graph and UMAP coordinates are stored in the AnnData object, they will be used for majority voting and plot generation, respectively.
   - **CSV / TSV (`.csv`, `.tsv`, `.txt`, `.tab`)**: A delimited matrix of raw counts or normalized counts. The matrix should be oriented as **cells (rows) × genes (columns)** unless `transpose_input` is set to `true`. The first column should contain cell identifiers, and the first row should contain gene symbols.
   - **10x Genomics MTX (`.mtx`)**: Sparse matrix format produced by Cell Ranger. When using MTX input, the `gene_file` and `cell_file` parameters must also be provided.

   **Size guidance**: CellTypist can process datasets ranging from a few hundred to over one million cells, though very large datasets may require significant memory (≥32 GB RAM recommended for >500,000 cells).

2. **gene_file** *(MTX input only)*
   A plain-text file containing gene names (one per line) corresponding to the rows of the MTX matrix. Typically named `features.tsv` or `genes.tsv` in 10x Genomics output directories.

3. **cell_file** *(MTX input only)*
   A plain-text file containing cell barcodes (one per line) corresponding to the columns of the MTX matrix. Typically named `barcodes.tsv` in 10x Genomics output directories.

4. **custom_model_file** *(optional)*
   A serialized CellTypist model file in `.pkl` format. Custom models can be trained on any user-labeled reference dataset using the `celltypist.train()` function and saved with `model.write()`. When provided, this file takes precedence over the `model` (built-in model) selection.

## Output Files

1. **`{output_prefix}_predicted_labels.csv`**
   A comma-separated table with one row per cell, containing:
   - `cell_id`: Cell barcode or index.
   - `predicted_labels`: Assigned cell type label(s). In `best match` mode, a single label; in `prob match` mode, pipe-separated labels for all types above the probability threshold.
   - `majority_voting` *(if enabled)*: Consensus label derived from neighborhood sub-cluster voting.
   - `conf_score`: Confidence score (maximum predicted probability) for the top prediction.

2. **`{output_prefix}_probability_matrix.csv`**
   A matrix of predicted probabilities with cells as rows and all cell types in the model vocabulary as columns. Each value represents the probability that a given cell belongs to a given cell type. This file is useful for downstream analyses such as plotting custom heatmaps or applying alternative thresholds.

3. **`{output_prefix}_decision_matrix.csv`**
   A matrix of raw decision function scores (prior to softmax normalization) with the same dimensions as the probability matrix. These scores reflect the raw logistic regression output and can be used for advanced scoring comparisons.

4. **`{output_prefix}_annotated.h5ad`**
   An updated AnnData object (`.h5ad`) with all annotation results embedded in `adata.obs`. This file is ready for immediate use in downstream Scanpy workflows, including differential expression analysis, trajectory inference, and further visualization.

5. **`{output_prefix}_umap.pdf` / `{output_prefix}_umap.png`** *(if `plot_results = true`)*
   UMAP scatter plots colored by predicted cell type labels and confidence scores. If UMAP coordinates are absent from the input AnnData, they are computed automatically.

6. **`{output_prefix}_dotplot.pdf` / `{output_prefix}_dotplot.png`** *(if `plot_results = true`)*
   Dot plots showing the distribution of predicted probabilities or label frequencies across cell types, providing a compact summary of annotation results.

## Example Data

Input:
[Example AnnData file (PBMC 3k)](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) — a standard 10x Genomics PBMC dataset commonly used for benchmarking single-cell workflows.

Output:
Example output files (predicted labels, probability matrix, and annotated `.h5ad`) are available in the [CellTypist GitHub repository tutorials](https://github.com/Teichlab/celltypist/tree/main/docs/notebooks).

## Requirements

- **Platform**: GenePattern server (cloud or local installation)
- **Docker Image**: `genepattern/celltypist:latest`
  - Python 3.8+
  - CellTypist ≥ 1.6
  - Scanpy ≥ 1.9
  - NumPy, Pandas, scikit-learn, anndata
- **Memory**: Minimum 8 GB RAM; 32 GB+ recommended for datasets with >500,000 cells.
- **Operating System**: Linux (via Docker container; no local installation required when running on GenePattern).
- **Input Normalization**: For `.h5ad` input, expression values must be log1p-normalized to 10,000 counts per cell prior to submission. CSV/TSV inputs with raw counts are normalized automatically.

## License

CellTypist is released under the [MIT License](https://github.com/Teichlab/celltypist/blob/main/LICENSE). The GenePattern module wrapper is provided for research use under the same terms.

## Version Comments

| Version | Release Date | Description |
| :--- | :--- | :--- |
| 1 | 2024-06-01 | Initial GenePattern module release. Supports AnnData, CSV, TSV, and 10x MTX input formats; built-in and custom model selection; best match and prob match prediction modes; optional majority voting and visualization outputs. |
