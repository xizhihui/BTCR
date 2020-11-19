#' load data
#'
#' load Cellranger VDJ data, including BCR and TCR.
#'
#' @importFrom magrittr %>%
#' @export
#' @param path Path to the Cellranger VDJ outs directory.
#' @param sample Name for this sample.
#' @param CR TCR or BCR.
#' @return A BTCR list contains "clonotypes", "contigs" and "merge" data.
#' @details The "merge" element in BTCR list object is a combination of "clonotypes" and "contigs". In "merge", each gene isoforms are pasted togother so that every cell has just one line.
loadSingle <- function(path, sample, CR = c("BCR", "TCR")) {
    CR <- match.arg(CR)
    origin_file <- c("clonotypes.csv", "filtered_contig_annotations.csv")

    clonotypes <- read.csv(file.path(path, origin_file[1]))
    contigs <- read.csv(file.path(path, origin_file[2]))
    contigs$clonotype_id <- contigs$raw_clonotype_id

    warning("Non-productive contigs are removed.")
    contigs <- contigs[contigs$productive == "True", ]

    contigs_cols <- c("barcode", "chain", "v_gene", "d_gene", "j_gene", "c_gene", "clonotype_id")
    clonotypes_cols <- c("clonotype_id", "cdr3s_aa", "cdr3s_nt")

    data <- merge(contigs[, contigs_cols], clonotypes[, clonotypes_cols], by = "clonotype_id")
    data_list <- lapply(split(data, data$barcode), function(dt) {
        dt2 <- dt[1, ]
        chain <- unique(dt$chain)
        #browser()
        if (CR == "BCR") {
            dt2$category <- dplyr::case_when(
                "Multi" %in% chain ~ "Multi-chain",
                ("IGH" %in% chain) && any(c("IGL", "IGK") %in% chain) ~ "Productive",
                any(c("IGL", "IGK") %in% chain) ~ "Light-chain only",
                "IGH" %in% chain ~ "Heavy-chain only"
            )
        } else {
            dt2$category <- dplyr::case_when(
                "Multi" %in% chain ~ "Multi-chain",
                ("TRA" %in% chain) && ("TRB" %in% chain) ~ "Productive",
                "TRA" %in% chain ~ "α-chain only",
                "TRB" %in% chain ~ "β-chain only"
            )
        }
        dt2$chain <- paste0(chain, collapse = ";")
        for (gene in c("v_gene", "d_gene", "j_gene", "c_gene")) {
            dt2[, gene] <- paste0(unique(dt[, gene]), collapse = ";")
        }
        dt2
    })
    data <- do.call(rbind, data_list)

    final <- list(
        clonotypes = clonotypes,
        contigs = contigs,
        merge = data
    )
    final <- lapply(final, function(dt) {
        dt$sample <- sample
        dt
    })
    attributes(final)$type <- CR
    final
}

#' @rdname loadSingle
#' @export
#' @inheritParams loadSingle
#' @param paths Named paths to Cellranger VDJ outs directory.
#' @param groups Groups corresponded to each sample.
loadMultiple <- function(paths, groups = NULL, CR = c("BCR", "TCR")) {
    CR <- match.arg(CR)
    if (is.null(names(paths))) {
        samples <- setNames(nm = basename(paths))
        message("Paths are not named. Use it's basename instead.")
        stopifnot(length(unique(samples)) != length(paths))
        names(paths) <- samples
    } else {
        samples <- setNames(nm = names(paths))
    }

    if (!is.null(groups) && length(groups) != length(samples)) {
        stop("Groups must be the same length to samples(paths).")
    } else if (is.null(groups)) {
        groups <- samples
        message("No groups provided. Use samples instead.")
    } else {
        names(groups) <- samples
    }

    final <- list() # clonotypes, contigs, merges
    datalist <- lapply(samples, function(sample) {
        data <- loadSingle(paths[sample], sample = sample, CR = CR)
        data <- lapply(data, function(dt) {
            idx <- which(sample == samples)
            #dt$sample <- sample
            dt$group <- groups[sample]

            if ("barcode" %in% colnames(dt)) {
                dt$barcode <- paste0(dt$barcode, "_", idx)
                message("Add sample idx suffix to barcode.")
            }
            dt$clonotype_newid <- paste0(dt$clonotype_id, "_", dt$sample)
            dt
        })
        data
    })

    for (nm in names(datalist[[1]])) {
        final[[nm]] <- do.call(rbind, lapply(datalist, function(x) x[[nm]]))
    }

    attributes(final)$type <- CR

    final
}

#' filter data
#'
#' Remove the Non-productive chains. If scrna is provided, cells not in scrna would be filtered either.
#'
#' @export
#' @param btcr The BTCR object.
#' @param scrna The Single Cell RNAseq Seurat object corresponded to the BCR or TCR. Used for filtering the barcodes.
#' @return A BTCR filtered by non-productive chains and invalid barcodes.
filterBtcr <- function(btcr, scrna = NULL) {
    btcr$merge <- btcr$merge %>% dplyr::filter(category == "Productive")
    btcr$contigs <- btcr$contigs %>% dplyr::filter(barcode %in% btcr$merge$barcode)
    if (!is.null(scrna)) {
        cells <- colnames(scrna)
        if (sum(btcr$merge$barcode %in% cells) == 0) {
            warning("There are no barcodes overlapped in scrna. SKip this filtering.")
        } else {
            btcr$merge <- dplyr::filter(btcr$merge, barcode %in% cells)
            btcr$contigs <- dplyr::filter(btcr$contigs, barcode %in% cells)
        }
    } else {
        warning("No scrna provided. The downstream plots would not be calculated based on the identified cells. \nFor example, T/B cells that failed to pass the quality control would be included in calculation.")
    }
    btcr
}

#' create the final clonotype id
#'
#' @importFrom magrittr %>%
#' @param by use the "cdr3s_nt" or "cdr3s_aa" to create id. Default is nt.
#' @param btcr The BTCR object.
#' @return A new BTCR object with "final_clonotype_id" column.
#' @export
refreshId <- function(btcr, by = c("cdr3s_nt", "cdr3s_aa")) {
    by <- match.arg(by)

    clonotypes <- btcr$merge[, c("barcode", "cdr3s_nt", "cdr3s_aa")]

    # reorder by frequency
    clonotypes_freq <- as.data.frame(table(clonotypes[, by]), stringsAsFactors = FALSE) %>%
        dplyr::arrange(dplyr::desc(Freq))
    # include all clonotypes
    all_clonotypes <- unique(c(clonotypes_freq$Var1, btcr$clonotypes[, by]))

    # make new id
    btcr$merge$final_clonotype_id <- paste0("clonotype", match(btcr$merge[, by], all_clonotypes))
    btcr$clonotypes$final_clonotype_id <- paste0("clonotype", match(btcr$clonotypes[, by], all_clonotypes))
    btcr$contigs$final_clonotype_id <- btcr$merge$final_clonotype_id[match(btcr$contigs$barcode, btcr$merge$barcode)]

    # use for-loop instead of lapply to retain the attributess.
    for (nm in names(btcr)) {
        btcr[[nm]] <- btcr[[nm]][order(gsub("clonotype", "", btcr[[nm]]$final_clonotype_id)), ]
    }
    btcr
}

#' assign status: clonal or not
#'
#' @importFrom magrittr %>%
#' @param idcol Column name used as the clonotype id.
#' @param scrna The Single Cell RNAseq Seurat object corresponded to the BCR or TCR.
#' @param group.by The meta.data column in scrna used to group the cells.
#' @return A new BTCR object with "status" and "status2"(if scrna and group.by provided).
#' @export
#' @details "status" is used to indicate whether a cell is clonal or not, "status2" is used to indicate their groups relationship.
#'   For each clonotypes, if there are at least two cells with it, these cells would be clonals.
#'   For each clonals, if they are from different cell group, they are cross-cluster; otherwise, they are within-cluster.
assignStatus <- function(btcr, idcol = "clonotype_id", scrna = NULL, group.by = NULL) {
    merged_clonotypes <- btcr$merge
    stopifnot(idcol %in% colnames(btcr$merge))
    # for each clonotypes,
    # 1. if there are more than two cells in it, it would be clonals;
    # 2. if these cells are from different cell group, they are cross-cluster; or they are within-cluster
    clonotypes <- as.data.frame(table(merged_clonotypes[, idcol]))
    clonals <- clonotypes$Var1[clonotypes$Freq > 1]
    merged_clonotypes$status <- ifelse(merged_clonotypes[, idcol] %in% clonals, "Clonal", "No clonal")

    if (is.null(scrna) || is.null(group.by)) {
        warning("scrna or group.by is NULL. Skip the Cross/Within Clonal assign.")
    } else {
        merged_clonotypes$status2 <- ifelse(merged_clonotypes[, idcol] %in% clonals, "Clonal", "No clonal")
        cellinfo <- scrna@meta.data[merged_clonotypes$barcode, group.by, drop = FALSE]
        clonals_within <- c()
        clonals_cross <- c()
        for (clonal in clonals) {
            cell_in_clonal <- merged_clonotypes$barcode[merged_clonotypes[, idcol] == clonal]
            cell_group <- unique(cellinfo[cell_in_clonal, group.by])
            if (length(cell_group) > 1) clonals_cross <- c(clonals_cross, cell_in_clonal)
            else clonals_within <- c(clonals_within, cell_in_clonal)
        }
        merged_clonotypes$status2[merged_clonotypes$barcode %in% clonals_cross] <- "Cross Clonal"
        merged_clonotypes$status2[merged_clonotypes$barcode %in% clonals_within] <- "Within Clonal"
    }
    btcr$merge <- merged_clonotypes
    btcr
}

#' calculate the clonotype overlaps in each group
#'
#' @param btcr A BTCR object.
#' @param scrna scrna The Single Cell RNAseq Seurat object corresponded to the BCR or TCR.
#' @param groupby Overlaps groups. The meta.data column used to group cells/barcodes, always be clusters.
#' @param splitby Sample groups. Overlaps won't be calculated cross groups. The meta.data column used to group cells/barcodes, always be sample groups.
#' @param id Which clonotype id used to calculate overlaps.
#' @return A BTCR object with "overlaps" and "id" attributes.
#' @export
calcOverlap <- function(btcr, scrna, groupby, splitby = NULL, id = c("final_clonotype_id", "clonotype_id", "clonotype_newid")) {
    id <- match.arg(id)
    # values is used for the groups to calculate the overlaps; if NULL, use all groups.

    stopifnot(groupby %in% colnames(scrna@meta.data))
    groups <- unique(as.character(scrna@meta.data[, groupby]))

    if (!is.null(splitby)) {
        stopifnot(splitby %in% colnames(scrna@meta.data))
        barcodes <- split(colnames(scrna), scrna@meta.data[, splitby])
    } else {
        barcodes <- list(all = colnames(scrna))
    }
    group_split_barcodes <- lapply(barcodes, function(barcode) {
        group_bc <- split(barcode, scrna@meta.data[barcode, groupby])
        group_bc <- group_bc[groups]
        group_clonotypes <- lapply(group_bc, function(grbc) {
            grbc <- intersect(grbc, btcr$merge$barcode)
            unique(btcr$merge[match(grbc, btcr$merge$barcode), id])
        })
        overlap <- gplots::venn(group_clonotypes, show.plot = FALSE)
        overlap <- attributes(overlap)$intersections
        overlap
    })

    attributes(btcr)$overlaps <- group_split_barcodes
    attributes(btcr)$id <- id
    attributes(btcr)$group <- groupby
    attributes(btcr)$split <- splitby
    btcr
}

#' save the BTCR data to a csv file
#'
#' @param btcr A BTCR object.
#' @param outfile The path to file to write in. It should be with a "csv" suffix.
#' @export
saveBtcr <- function(btcr, outfile) {
    dt <- btcr$merge
    if (!endsWith(outfile, "csv")) {
        warning("Results will saved as comma seperated values file. But outfile has no 'csv' suffix")
    }
    overlaps <- attributes(btcr)$overlaps

    if (!is.null(overlaps)) {
        dt$overlaps <- "unique"
        groupby <- attributes(btcr)$group
        splitby <- attributes(btcr)$split
        id <- attributes(btcr)$id
        message(groupby, splitby, id)

        groups <- names(overlaps)
        for (grp in groups) {
            curoverlaps <- overlaps[[grp]]
            innerov <- names(curoverlaps)
            for (ov in innerov) {
                if (length(groups) == 1 && groups == "all") {
                    dt[dt[, id] %in% curoverlaps[[ov]], "overlaps"] <- ov
                } else {
                    dt[dt[, splitby] == grp & dt[, id] %in% curoverlaps[[ov]], "overlaps"] <- ov
                }
            }
        }
    }
    write.csv(dt, file = outfile, quote = TRUE, row.names = FALSE)
}