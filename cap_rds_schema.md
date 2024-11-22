# RDS file schema

## AnnData fields mapping to Seurat5

In the table below `seurat_obj` represents a Seurat5 object, the result of the conversion. It consists of only one assay, which represents an `Assay5` object and is defined as `assay <- seurat_obs[["RNA"]]`

| AnnData Component  | Seurat Component                       |
|--------------------|----------------------------------------|
| `adata.X`          | `assay@layers$data`                    |
| `adata.var`        | `assay@meta.data`                      |
| `adata.raw.X`*     | `assay@layers$counts`                  |
| `adata.raw.var`*   | -                                      |
| `adata.obs`        | `seurat_obj@meta.data`                 |
| `adata.uns`        | `seurat_obj@misc`                      |
| `adata.obsm`       | `seurat_obj@reductions`                |
`*` - if it exists

## Naming

The Seurat5 RDS schema (rules for field naming) mirrors the input AnnData file schema. So, if one uses an AnnData file from the [Cell Annotation Platform](https://celltype.info/) the same naming conventions will be applied to the Seurat5 object (see [Cap-AnnData schema](https://github.com/cellannotation/cell-annotation-schema/blob/main/docs/cap_anndata_schema.md)). The only exception is the names of keys in mapping fields like `.uns` and `.obsm`. For those fields, all underscores `_` will be replaced with dots `.`, and if the key starts with `X_`, this prefix will be removed. 

Examples:
 
| Field name in AnnData      | Field name in Seurat object            |
|----------------------------|----------------------------------------|
| `adata.uns["schema_v"]`    | `seurat_obj@misc$schema.v`             |
| `adata.obsm["X_tsne"]`     | `assay@reductions$tsne`                |
| `adata.obsm["X_tsne_2"]`   | `assay@reductions$tsne.2`              |
