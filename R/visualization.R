#' barplot of cells with TCR or BCR percentages in grouped T cells or B cells
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise
#' @importFrom scales percent
#' @import ggplot2
#' @param btcr A BTCR object.
#' @param scrna The Single Cell RNAseq Seurat object corresponded to the BCR or TCR.
#' @param group The meta.data column in scrna to group cells.
#' @return A ggplot object.
#' @export
plotCellCR <- function(btcr, scrna, group) {
    stopifnot(group %in% colnames(scrna@meta.data))
    data <- data.frame(
        row.names = colnames(scrna),
        group = scrna@meta.data[, group],
        CR = FALSE
    )
    cr_cells <- btcr$merge$barcode %>% unique() %>% intersect(rownames(data))
    data[cr_cells, "CR"] <- TRUE

    ylab <- attributes(btcr)$type

    data <- group_by(data, group) %>% summarise(percent = sum(CR) / length(CR))
    gg <- ggplot(data, aes(x = group, y = percent)) +
        geom_bar(stat = "identity", fill = "#8BC0C9") +
        scale_y_continuous(labels = percent, limits = c(0, round(max(data$percent + 0.1), digits = 2))) +
        theme_classic() +
        labs(x = "", y = ylab)
    gg
}

#' barplot of distributions of clone status in grouped T cells or B cells.
#'
#' @importFrom magrittr %>%
#' @importFrom reshape2 melt dcast
#' @importFrom scales percent_format hue_pal
#' @import ggplot2
#' @param btcr A BTCR object.
#' @param scrna The Single Cell RNAseq Seurat object corresponded to the BCR or TCR.
#' @param group.by The meta.data column in scrna to group cells, like clusters.
#' @param split.by The meta.data column in scrna to group cells, like sample groups.
#' @param status The clonal status used.
#' @return ggplot.
#' @export
plotCloneState <- function(btcr, scrna, group.by, split.by = NULL, status = c("status", "status2")) {
    status <- match.arg(status)
    stopifnot(group.by %in% colnames(scrna@meta.data))

    scrna$status <- "Non-VDJ"
    scrna$status[btcr$merge$barcode] <- btcr$merge[, status]

    status_uniq <- unique(btcr$merge[, status])
    colormap <- setNames(hue_pal()(length(status_uniq)), nm = status_uniq)
    colormap <- c("Non-VDJ" = "gray", colormap)

    data <- data.frame(
        x = scrna@meta.data[, group.by],
        fill = scrna@meta.data$status,
        row.names = colnames(scrna)
    )

    reformat <- function(dt) {
        dt <- dcast(dt, x ~ fill)
        dt[, -1] <- dt[, -1] / sum(dt[, -1])
        melt(dt, variable.name = "fill", value.name = "percent")
    }

    if (!is.null(split.by)) {
        stopifnot(split.by %in% colnames(scrna@meta.data))
        data$grp <- scrna@meta.data[, split.by]
        data2 <- split(data, data$grp)
        data2 <- do.call(rbind, lapply(data2, reformat))
        data2$grp <- gsub("\\..*$", "", rownames(data2))
    } else {
        data2 <- reformat(data)
    }

    gg <- ggplot(data2, aes(x = x, y = percent, fill = fill)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = colormap, name = "Clone state") +
        scale_y_continuous(labels = percent_format()) +
        theme_classic() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_blank(),
            strip.text = element_text(size = 12)
        ) + labs(y = "Distribution of clone status", x = "")
    if (!is.null(split.by)) gg <- gg + facet_wrap(~ grp, nrow = 1)
    gg
}


#' barplot with points and error bars based on mulitple samples data.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise mutate n arrange top_n desc
#' @import ggplot2
#' @param data The data frame to plot.
#' @param x The column used as X axis.
#' @param y THe column used as Y axis.
#' @param fill THe column used as groups.
#' @param point Whether to plot points.
#' @param errbar Wheter to plot error bars.
#' @export
#' @return ggplot.
barPointError <- function(data, x, y, fill, point = TRUE, errbar = TRUE) {
    se <- function(x) { sd(x) / sqrt(length(x)) }
    newdata <- data[, c(x, y, fill)]
    colnames(newdata) <- c("x", "y", "fill")

    data_mean <- group_by(newdata, fill, x) %>%
        summarise(mean = mean(y), se = se(y))
    #data_mean$se[is.na(data_mean$se)] <- 0

    dodge <- position_dodge(width = 0.8)
    gg <- ggplot(data_mean, aes(x = x, y = mean, fill = fill)) +
        geom_bar(position = dodge, stat = "identity", width = 0.7)

    if (point) {
        gg <- gg + geom_point(data = newdata, aes(y = y), position = dodge, show.legend = FALSE)
    }
    if (errbar) {
        gg <- gg + geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = dodge, width = 0.25)
    }

    gg <- gg + theme_classic(base_size = 12) + labs(x = "", y = y, fill = "")
    gg
}

#' barplot of each VDJ gene usages
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise mutate n case_when
#' @importFrom patchwork wrap_plots
#' @importFrom scales percent
#' @import ggplot2
#' @param btcr A BTCR object.
#' @param group Column in btcr$contigs to group each VDJ genes, always be sample or group.
#' @param sample The sample column in btcr$contigs.
#' @param gene Plot v_gene, j_gene or c_gene.
#' @param dtype Plot barcode usage counts or percentages.
#' @param isoform If set TRUE, plot the gene isoforms like TRAV-xxx. Or just plot the main gene like TRA, TRB (isoforms are summed.).
#' @param combine Wrap plots togother. Or return the plot list.
#' @param rotate Rotate the x labels.
#' @param ncol When plot isoform, plot ncol in the layout.
#' @return A ggplot or ggplot list or patchwork wrapped plot.
#' @export
plotGeneBar <- function(btcr, group, sample = "sample", gene = c("v", "j", "c"),
                        dtype = c("percentage", "count"), isoform = TRUE,
                        combine = TRUE, rotate = TRUE, ncol = 1) {
    CR <- attributes(btcr)$type
    dtype <- match.arg(dtype)
    gene <- match.arg(gene)
    genecol <- paste0(match.arg(gene), "_gene")
    stopifnot(group %in% colnames(btcr$contigs))
    stopifnot(sample %in% colnames(btcr$contigs))

    # BCR:
    #   v_gene: IG[HKL]V, use 3 chars
    #   j_gene: IG[HKL]J, use 3 chars
    #   c_gene: IG[KL]C, IGH[ADGEM], use 4 chars
    # TCR:
    #   v_gene: TR[AB]V, use 3 chars
    #   j_gene: TR[AB]J, use 3 chars
    #   c_gene: TR[AB]C[12], use 4 chars
    pre_nchar <- case_when(
        gene %in% c("v", "j") ~ 3,
        gene == "c" ~ 4
    )

    # we must include "sample" inside to calculate gene usage in each sample.
    # group means: sample-group
    contigs <- btcr$contigs[, c("barcode", genecol, sample, group)]
    colnames(contigs) <- c("barcode", "gene", "sample", "group")

    # need points? if a group has multiple samples, need.
    show_point <- length(unique(contigs$sample)) != length(unique(contigs$group))

    # calculate gene isoform usage percentage in single sample.
    if (isoform) {
        count_list <- split(contigs, substr(contigs$gene, start = 1, stop = pre_nchar))
    } else {
        contigs$gene <- substr(contigs$gene, start = 1, stop = pre_nchar)
        count_list <- list(gene = contigs)
    }

    count_list <- lapply(count_list, function(innercount) {
        innercount <- innercount %>%
            group_by(group, sample, gene) %>%
            summarise(count = n()) %>%
            mutate(percentage = count / sum(count))
    })

    plist <- list()
    for (gene_name in names(count_list)) {
        plist[[gene_name]] <- barPointError(
            data = count_list[[gene_name]],
            x = "gene",
            y = dtype,
            fill = "group",
            point = show_point,
            errbar = show_point
        ) + ggtitle(gene_name)
        if (dtype == "percentage") plist[[gene_name]] <- plist[[gene_name]] + scale_y_continuous(labels = percent)
        if (rotate) plist[[gene_name]] <- plist[[gene_name]] + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }

    if (combine && isoform) {
        wrap_plots(plist, ncol = ncol, guides = "collect")
    } else if (!isoform) {
        # not isoform
        plist[[1]] + ggtitle("")
    } else {
        plist
    }
}

#' dimension plots of clonotypes
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange desc pull
#' @importFrom scales hue_pal
#' @importFrom Seurat DimPlot
#' @import ggplot2
#' @param btcr A BTCR object.
#' @param scrna The Single Cell RNAseq Seurat object corresponded to the BCR or TCR.
#' @param colorby can be "CR", "cdr3s_aa", "cdr3s_nt". If "CR", only plot whether a cell has a TCR or BCR.
#'   If "cdr3s", use topn or highlight to limit the clonotypes.
#' @param topn Top N barcode frequency cdr3s to plot. Will be ignored if highlight is provided.
#' @param highlight Which cdr3s in colorby to plot.
#' @param ... Arguments for \code{\link[Seurat::DimPlot]{Seurat::DimPlot}}.
#' @return Seurat::DimPlot returns.
#' @export
plotCloneDim <- function(btcr, scrna, colorby = c("CR", "cdr3s_aa", "cdr3s_nt"),
                         topn = 5, highlight = NULL, ...) {
    # topn and highlight are conflicted. Highlight has a higher priority.

    colorby <- match.arg(colorby)
    CR <- attributes(btcr)$type

    if (colorby == "CR") {
        scrna$clonedim <- paste0("No ", CR)
        scrna$clonedim[btcr$merge$barcode] <- CR
        colors <- setNames(
            c("gray", hue_pal()(1)),
            nm = paste0(c("No ", ""), CR)
        )
        DimPlot(scrna, group.by = "clonedim", cols = colors, ...)
    } else {
        scrna$clonedim <- "Others"

        if (is.null(highlight)) {
            # get the topn colorby values
            colorbydata <- as.data.frame(table(btcr$merge[, colorby])) %>%
                arrange(desc(Freq)) %>%
                pull(Var1) %>% .[1:topn]
            # if (length(colorbydata) > topn) {
            #     warning("Too many ", colorby, "s have the same Frequence. Randomly select ", topn, ".")
            #     colorbydata <- colorbydata[1:topn]
            # }
        } else {
            colorbydata <- intersect(unique(btcr$merge[, colorby]), highlight)
            if (length(colorbydata) == 0) {
                stop("Items to highlight are not in ", colorby)
            }
        }

        celldata <- btcr$merge[btcr$merge[, colorby] %in% colorbydata, c("barcode", colorby)]
        celldata <- split(celldata$barcode, celldata[, colorby])

        colors <- setNames(hue_pal()(length(colorbydata)), nm = colorbydata)

        DimPlot(scrna, cells.highlight = celldata, cols.highlight = colors, ...)
    }
}

#' barplot of top n cdr3 sequences usage in each group.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise mutate n arrange top_n desc
#' @importFrom ggplotify as.ggplot
#' @importFrom scales percent
#' @importFrom patchwork wrap_plots
#' @import ggplot2
#' @param btcr A BTCR object.
#' @param group Column in btcr$contigs to group cdr3.
#' @param dtype Plot usage counts or percentages.
#' @param topn Top N cdr3 to plot. Default: 10.
#' @param ncol Ncols in plot layout.
#' @param combine Combine each plots or return plot list.
#' @return patchwork wrapped plots or ggplot list.
#' @export
plotCdr3 <- function(btcr, group, dtype = c("percentage", "count"),
                     topn = 10, ncol = 1, combine = TRUE) {
    # here we don't use facet_grid( chain ~ group)
    stopifnot(group %in% colnames(btcr$contigs))
    dtype <- match.arg(dtype)
    data <- data.frame(
        chain = btcr$contigs$chain,
        cdr3 = btcr$contigs$cdr3,
        group = btcr$contigs[, group]
    )

    data_list <- split(data, data$chain)

    data_list <- lapply(data_list, function(dt) {
        group_by(dt, group, cdr3) %>%
            summarise(count = n()) %>%
            mutate(percentage = count / sum(count))
    })

    data_list <- lapply(data_list, function(dt) {
        # cdr3_grp is used for ordering the bars in facet.
        dt <- arrange(dt, desc(count)) %>%
            mutate(cdr3_grp = paste0(cdr3, "_", group)) %>%
            group_by(group) %>%
            top_n(n = topn, wt = count)
        dt$cdr3_grp <- factor(dt$cdr3_grp, levels = unique(dt$cdr3_grp))
        dt
    })

    plist <- lapply(names(data_list), function(nm) {
        dt <- data_list[[nm]]
        shared_cdr3 <- unique(dt$cdr3[duplicated(dt$cdr3)])
        gg <- ggplot(dt, aes_string(x = "cdr3_grp", y = dtype, fill = "group")) +
            geom_col() +
            scale_x_discrete(labels = function(x) gsub("_.*$", "", x)) +
            scale_fill_discrete(guide = FALSE) +
            facet_wrap(~ group, nrow = 1, scales = "free_x") +
            labs(x = "", y = nm) +
            theme_classic() +
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                #axis.text.x = element_blank(),
                strip.background = element_blank(),
                strip.text = element_text(size = 12)
            )
        if (dtype == "percentage") gg <- gg + scale_y_continuous(labels = percent)
        ggrob <- ggplotGrob(gg)
        idxes <- which(startsWith(ggrob$layout$name, "axis-b"))
        for (idx in idxes) {
            texts <- ggrob$grobs[[idx]]$children$axis$grobs[[2]]$children[[1]]$label
            cols <- rep("gray30", length(texts))
            cols[texts %in% shared_cdr3] <- "red"
            ggrob$grobs[[idx]]$children$axis$grobs[[2]]$children[[1]]$gp$col <- cols
        }
        as.ggplot(ggrob)
    })
    if (combine) {
        wrap_plots(plist, ncol = ncol)
    } else {
        plist
    }
}


#' heatmap of v-j pairs usage
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise mutate n arrange top_n desc
#' @import ggplot2
#' @param btcr A BTCR object.
#' @param group COlumn in btcr$contigs to group genes.
#' @param dtype Plot usage counts or percentages.
#' @param topn Just plot the v-j pairs with top N usage. Default: 10.
#' @param ncol Ncols in plot layout.
#' @param position The position of group labels in plot.
#' @return ggplot
#' @export
plotVJpair <- function(btcr, group, dtype = c("count", "percentage"), topn = 10,
                       ncol = 1, position = c("right", "left", "bottom", "top")) {
    dtype <- match.arg(dtype)
    position <- match.arg(position)
    stopifnot(group %in% colnames(btcr$contigs))

    data <- data.frame(
        v_gene = btcr$contigs$v_gene,
        j_gene = btcr$contigs$j_gene,
        group = btcr$contigs[, group]
    )

    data <- group_by(data, group, v_gene, j_gene) %>%
        summarise(count = n()) %>%
        group_by(group) %>%
        mutate(percentage = count / sum(count)) %>%
        top_n(n = topn, wt = count)

    ggplot(data, aes_string(x = "v_gene", y = "j_gene", fill = dtype)) +
        geom_tile() +
        scale_fill_gradient(low = "#FFF5FA", high = "blue") +
        facet_wrap(~ group, ncol = ncol, strip.position = position) +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            strip.background = element_blank(),
            strip.text.y = element_text(size = 12, angle = 0),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
            axis.ticks = element_blank()
        ) +
        labs(x = "", y = "")
}


#' dimension plots of shared or unique clonotypes in each groups.
#'
#' @importFrom magrittr %>%
#' @importFrom Seurat DimPlot NoLegend
#' @importFrom scales hue_pal
#' @importFrom patchwork wrap_plots
#' @import ggplot2
#' @param btcr A overlaps-calculated BTCR object. see \code{\link{calcOverlap}}.
#' @param scrna The Single Cell RNAseq Seurat object corresponded to the BCR or TCR.
#' @param values The overlapped groups to plot, at least two values in \code{attributes(btcr)$group}.
#' @param dtype Plot the shared clonotypes or the unique clonotypes.
#' @param ncol Ncols in plot layout.
#' @param ... Arguments for \code{\link[Seurat::DimPlot]{Seurat::DimPlot}}.
#' @return Seurat::DimPlot returns.
#' @export
plotSharedCloneDim <- function(btcr, scrna, values, dtype = c("share", "unique"), ncol = 1, ...) {
    dtype <- match.arg(dtype)
    overlaps <- attributes(btcr)$overlaps
    id <- attributes(btcr)$id

    all_values <-  unique(unlist(sapply(overlaps, names)))
    all_values2 <- unique(unlist(sapply(all_values, function(x) unlist(strsplit(x, ":")))))

    values <- intersect(values, all_values2)
    if (length(values) <= 1) {
        # there must be at least 2 values to get overlaps.
        stop("Values are not all in the overlap items. Does values are the subset of groupby items in CalcOverlaps?")
    }

    if (dtype == "share") {
        item <- Filter(function(x) {
            all(sapply(values, grepl, x = x))
        }, all_values)
    } else {
        item <- values
    }

    plist <- lapply(overlaps, function(ov) {
        clono_ids <- ov[item]
        clono_ids2 <- unique(unlist(clono_ids))
        #clono_ids <- unique(unlist(ov[item]))

        if (length(clono_ids) == 0 || is.null(clono_ids2)) {
            cells <- colnames(scrna)
            subtitle <- "0 clonotypes, 0 cells"
            plt <- DimPlot(scrna, cells.highlight = cells, cols.highlight = "gray", ...) + NoLegend()
        } else {
            cells <- sapply(clono_ids, function(cid) {
                btcr$merge$barcode[btcr$merge[, id] %in% cid]
            })
            subtitle <- paste0(length(clono_ids2), " clonotypes, ", length(unique(unlist(cells))), " cells")

            if (dtype == "share") colors <- "#DE2D26"
            else colors <- setNames(hue_pal()(length(cells)), nm = names(cells))
            plt <- DimPlot(scrna, cells.highlight = cells, cols.highlight = colors, ...)
        }

        if (dtype == "share") plt <- plt + NoLegend()
        plt <- plt + labs(subtitle = subtitle)
    })

    plist <- lapply(names(plist), function(nm) {
        plist[[nm]] + labs(title = nm) + theme(
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)
        )
    })
    if (length(plist) == 1) plist[[1]] <- plist[[1]] + labs(title = "")
    wrap_plots(plist, ncol = ncol, guides = "collect")
}


#' plot the antigen species
#'
#' @importFrom Seurat DimPlot
#' @importFrom Seurat DimPlot
#' @param btcr The BTCR object annotated by annotateTCR
#' @param scrna The Single Cell RNAseq Seurat object corresponded to the BCR or TCR.
#' @param dtype "TRA" or "TRB". This function supports TCR only.
#' @param ... Extra params(except "group.by") passed to \code{Seurat::DimPlot}.
#' @return \code{Seurat::DimPlot} returns.
plotAnnoDim <- function(btcr, scrna, dtype = c("TRA", "TRB"), ...) {
    if (missing(scrna)) {
        stop("scrna is required.")
    }
    dtype <- match.arg(dtype)
    CR <- attr(btcr, "type")
    if (CR != "TCR") {
        stop("plotAnnoDim supports TCR only, now it is ", CR)
    }

    annodata <- attr(btcr, "annotation")
    if (is.null(annodata)) {
        stop("There is no annotation in btcr. Please run annotateTCR first.")
    }

    scrna@meta.data$vdjanno <- NA
    scrna@meta.data[annodata[[dtype]]$barcode, "vdjanno"] <- annodata[[dtype]]$antigen.species

    DimPlot(scrna, group.by = "vdjanno", ...)
}