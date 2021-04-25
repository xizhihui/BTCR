# BTCR

Simple Stats and Plots for BCR and TCR.


## usage

* load data

```R
library(Seurat)
library(BTCR)

paths <- c(
    "S1" = "cellranger_vdj/S1/outs",
    "S2" = "cellranger_vdj/S2/outs",
    "A1" = "cellranger_vdj/A1/outs",
    "A2" = "cellranger_vdj/A2/outs"
)
scrna <- readRDS("scRNA_seurat_obj_integrated_these_four_samples.rds")

# for single sample: loadSingle
btcr <- loadMultiple(paths, groups = c("S", "S", "A", "A"), CR = "TCR")
```

* barcode correction

If vdj samples are less than samples in scRNA, barcodes in `btcr` should be corrected.

```R
all_samples <- unique(scrna$orig.ident)
# S1, S2, S3, A1, A2, A3

for (sample in names(paths)) {
    sample_idx <- which(sample == all_samples)
    
    bc_idx <- btcr$contigs$sample == sample
    btcr$contigs$barcode[bc_idx] <- gsub("\\d$", sample_idx, btcr$contigs$barcode[bc_idx])
    
    bc_idx <- btcr$merge$sample == sample
    btcr$merge$barcode[bc_idx] <- gsub("\\d$", sample_idx, btcr$merge$barcode[bc_idx])
}
```

* preprocess

```R
# not all cells in btcr are in scrna, maybe some have been filtered.
btcr <- filterBtcr(btcr, scrna = scrna)

btcr <- refreshId(btcr, by = "cdr3s_nt")
btcr <- assignStatus(btcr, idcol = "clonotype_newid", scrna = scrna, group.by = "orig.ident")

# use VDJdb to annotate the TCRs. BCR is not supported.
btcr <- annotateTCR(btcr, latest = TRUE)
```

* visulizations

```R
# barplot for CDR3 sequence usage
plotCdr3(btcr, group = "sample", topn = 10, dtype = "percentage")


# barplot for V J C gene usage
plotGeneBar(btcr, group = "sample", gene = "j", dtype = "percentage")


# heatmap for V-J gene pair usage
plotVJpair(btcr, group = "sample", topn = 30, ncol = 2, dtype = "percentage")


# barplot for TCR/BCR-detected cells in scrna
if ("celltype" %in% colnames(scrna@meta.data)) {
    # count percentage in all cells
    plotCellCR(btcr, group = "celltype", scrna = scrna)
    
    # count percentage in only T cells
    tcells <- subset(scrna, subset = celltype %in% c("CD8 T", "CD4 T"))
    plotCellCR(btcr, group = "celltype", scrna = tcells)
}


# tsne/umap for TCR/BCR-detected cells in scrna
# only show TCR/BCR or not
plotCloneDim(btcr, scrna = tcells, colorby = "CR", reduction = "tsne")
# show top frequency cdr3s_aa
plotCloneDim(btcr, scrna = tcells, colorby = "cdr3s_aa", reduction = "tsne")


# barplot for clonal status: clonal, in-cluster clonal, or cross-cluster clonal
plotCloneState(btcr, scrna = tcells, group.by = "celltype", split.by = "orig.ident", status = "status2")


# check the TRA/TRB antigen species
plotAnnoDim(btcr, scrna, dtype = "TRA")
```


## reference

* Plots and Stats are learnt from: [Single-cell landscape of immunological responses in patients with COVID-19](https://www.nature.com/articles/s41590-020-0762-x)
