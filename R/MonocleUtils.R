# Utility functions for working with Monocle3

#' Make palette for categorical values
#'
#' Automatically choose palette given the required number of colors
#'
#' @param n Number of required colors (integer).
#'
#' @return A vector of colors
#'
#' @import RColorBrewer
.make_color_palette <- function(n) {
  if (n <= 8) {
    pal <- RColorBrewer::brewer.pal(n, name = "Dark2")
  } else if (n <= 12) {
    pal <- RColorBrewer::brewer.pal(n, name = "Paired")
  } else {
    if (requireNamespace("Polychrome", quietly = TRUE)) {
      if (n <= 31) {
        pal <- Polychrome::glasbey.colors(n + 1)[-1]
      } else if (n <= 36) {
        pal <- Polychrome::palette36.colors(n)
      } else {
        pal <- Polychrome::createPalette(n, c("#ff0000", "#00ff00", "#0000ff"))
      }
    } else {
      set.seed(1)
      pal <- sample(rainbow(n, end = 0.8), n, replace = FALSE)
    }
  }
  unname(pal)
}


#' Add row labels to a pheatmap
#'
#' Selectively add row labels to a pheatmap object with pretty layout
#' Code taken from
#' https://stackoverflow.com/questions/52599180/partial-row-labels-heatmap-r
#' (all credit goes to the original author) with some modifications to styling
#'
#' @param pheatmap A pheatmap object
#' @param kept_labels A pheatmap object
#' @param repel_degree A pheatmap object
#'
#'
#' @return A vector of colors
.add_flag <- function(pheatmap,
                     kept_labels,
                     repel_degree,
                     add_pointer = FALSE) {

  # repel_degree = number within [0, 1], which controls how much
  #                space to allocate for repelling labels.
  ## repel_degree = 0: spread out labels over existing range of kept labels
  ## repel_degree = 1: spread out labels over the full y-axis

  heatmap <- pheatmap$gtable

  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]]

  # keep only labels in kept_labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept_labels, new.label$label, "")

  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel_degree) {
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant

    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }

      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }

    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))

    return(unit(
      seq(
        from = max(selected.range) + k*(max(full.range) - max(selected.range)),
        to = min(selected.range) - k*(min(selected.range) - min(full.range)),
        length.out = sum(d.select)
      ), "npc"
    ))
  }
  new.y.positions <- repelled.y(new.label$y, d.select = new.label$label != "")

  if (add_pointer) {
    new.flag <- grid::segmentsGrob(
        x0 = new.label$x,
        x1 = new.label$x + unit(0.07, "npc"),
        y0 = new.label$y[new.label$label != ""],
        y1 = new.y.positions,
        gp = gpar(lwd = 0.5)
    )
    heatmap <- gtable::gtable_add_grob(x = heatmap, grobs = new.flag, t = 4, l = 4)

    # shift position for selected labels
    new.label$x <- new.label$x + unit(0.1, "npc")
  }

  new.label$y[new.label$label != ""] <- new.y.positions

  # add flag to heatmap

  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label

  # return a copy of the heatmap invisibly
  invisible(heatmap)
}


#' Plot trajectory
#'
#' This function enhances monocle3::plot_cells() including using better color
#' palette and some other features
#'
#' @param cds CDS object.
#' @param color_cells_by Column in `colData(cds)` that dictates colouring
#'   (character).
#' @param legend_loc Location of color legend, can be "on data", "right" or
#'   "none", default "on data" (character).
#' @param color_palette Color palette to use, can be "viridis", name of brewer
#'   pal or a vector of colors, default NULL (character).
#' @param label_graph_nodes Show root, leaf and branching nodes with principal
#'   point names, default "all" (character).
#' @param cell_size Dot size for cells, if NULL set to `50/sqrt(ncol(cds))`,
#'   default NULL (numeric).
#' @param ... Additional options passed to `monocle3::plot_cells()`.
#'
#' @return A ggplot object
plot_trajectory <- function(cds, color_cells_by, legend_loc = "on data", color_palette = NULL, label_graph_nodes = FALSE, cell_size = NULL, ...) {
  if (requireNamespace("monocle3", quietly = TRUE)) {
    color_by <- colData(cds)[[color_cells_by]]
    if (is.numeric(color_by) || color_cells_by %in% c("pseudotime")) {
      if (is.null(color_palette)) {
        scale_color <- viridis::scale_color_viridis(name = color_cells_by)
      } else if (is.character(color_palette) && length(color_palette) == 1) {
        scale_color <- scale_color_distiller(palette = colorRampPalette(color_palette))
      } else {
        scale_color <- scale_color_distiller(palette = color_palette)
      }
    } else {
      n_color <- length(unique(color_by))
      if (is.null(color_palette)) {
        palette <- .make_color_palette(n_color)
        scale_color <- scale_color_manual(values = palette)
      } else if (is.character(color_palette) && length(color_palette) == 1) {
        scale_color <- scale_color_brewer(palette = color_palette)
      } else {
        scale_color <- scale_color_manual(values = color_palette)
      }
    }
    if (is.null(cell_size)) cell_size <- 50 / sqrt(ncol(cds))
    p <- monocle3::plot_cells(
      cds,
      color_cells_by = color_cells_by, label_cell_groups = (legend_loc == "on data"), label_groups_by_cluster = FALSE,
      label_leaves = any(label_graph_nodes %in% c("leaf", "all")),
      label_branch_points = any(label_graph_nodes %in% c("branch", "branching", "all")),
      label_roots = any(label_graph_nodes %in% c("root", "all")),
      label_principal_points = label_graph_nodes == "all", cell_size = cell_size, ...
    ) + scale_color
    if (legend_loc == "none") {
      p <- p + guides(color = "none")
    }
  } else {
    stop("package:monocle3 not available")
  }
  p
}


#' Get principal graph node by group
#'
#' Get the node in principal graph that is the closest to the greatest number of
#' cells in a group. See
#' https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/#order-cells
#'
#' @param cds CDS object.
#' @param group_by Column in `colData(cds)` that define grouping (character).
#' @param group The group of interest (character).
#' @param node_type Restrict to this type of node if specified: "any", "leaf" or
#'   "branch", default "any" (character).
#' @param choose_one Choose one to return if multiple nodes are found, default
#'   TRUE (logical).
#' @param choose_by If `choose_one` is TRUE, order nodes by this criteria. Can
#'   be a column in `colData(cds)` or "n" (the number of cells the node is
#'   closest to). Default "n" (character).
#' @param ascending Order nodes by `choose_by` in ascending order, default FALSE
#'   (logical).
#' @param reduction_method Reduction method, default "UMAP" (character).
#'
#' @return A node in the principal graph
get_principal_node <- function(cds,
                               group_by=NULL,
                               group=NULL,
                               node_type = "any",
                               choose_one = TRUE,
                               choose_by = "n",
                               ascending = FALSE,
                               reduction_method = "UMAP") {
  node_type <- match.arg(node_type, c("any", "leaf", "branch"))
  if (requireNamespace("monocle3", quietly = TRUE)) {
    # `node_names` is a vector of node names
    node_names <- igraph::V(monocle3::principal_graph(cds)[[reduction_method]])$name

    # `node_pool` is a vector of node names restricted by `node_type`
    if (node_type == "leaf") {
      node_pool <- names(c(monocle3:::leaf_nodes(cds), monocle3:::root_nodes(cds)))
    } else if (node_type == "branch") {
      node_pool <- names(monocle3:::branch_nodes(cds))
    } else {
      node_pool <- node_names
    }

    # `closest_vertex` is a named vector, with name being cell_id and value being node index
    closest_vertex <- cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_closest_vertex[colnames(cds), ]

    # `closest_node` is a named vector, with name being cell_id and value being node name
    closest_node <- node_names[closest_vertex]
    names(closest_node) <- names(closest_vertex)

    if (choose_one) {
      if (choose_by %in% colnames(SummarizedExperiment::colData(cds))) {
        metric <- SummarizedExperiment::colData(cds)[[choose_by]]
        names(metric) <- colnames(cds)
      } else {
        metric <- NULL
      }
    }

    if (!is.null(group_by) && !is.null(group)) {
      cell_idx <- which(SummarizedExperiment::colData(cds)[, group_by] == group)

      closest_node <- closest_node[cell_idx]
      if (choose_one) metric <- metric[cell_idx]
    }

    if (choose_one) {
      ncell_by_node <- tapply(rep(1, length(closest_node)), closest_node, FUN = sum, simplify = TRUE)[node_pool]
      top_nodes <- names(ncell_by_node)[order(ncell_by_node, decreasing = TRUE)][1:5]
      top_nodes <- top_nodes[!is.na(top_nodes)]
      if (!is.null(metric)) {
        summarised_metric <- tapply(metric, closest_node, FUN = mean, simplify = TRUE)
        if (ascending) {
          pr_nodes <- names(which.min(summarised_metric[top_nodes]))
        } else {
          pr_nodes <- names(which.max(summarised_metric[top_nodes]))
        }
      } else {
        pr_nodes <- top_nodes[1]
      }
    } else {
      pr_nodes <- intersect(unique(closest_node), node_pool)
    }

    pr_nodes
  } else {
    stop("package:monocle3 not available")
  }
}


#' Plot heatmap alone trajectory
#'
#' Plot a gene (row) by cell (column) heatmap with cells ordered by pseudotime
#'
#' @param cds CDS object.
#' @param genes A vector of genes to plot (verctor of character).
#' @param min_gene_sd Minimum gene standard deviation across cells to remove
#'   non-variable genes, default 0.2 (numeric).
#' @param min_gene_fraction Minimum fraction of cells the included genes are
#'   expressed, default 0.05 (numeric).
#' @param n_cells Number of cells to plot, down-sample to this number if fewer
#'   then `ncol(cds)`, default 500 (integer).
#' @param order_genes_by Method to re-order genes: "ARSA", "GW", "OLO", "TSP",
#'   "HC", "HC_single", "HC_complete", "HC_average", "HC_ward", "Spectral",
#'   "Spectral_norm", "corr", "rank", FALSE to disalbe ordering, default
#'   "Spectral" if `package:seriation` is available otherwise "HC"
#'   (character|bool).
#' @param rank_threshold Take only cells with expression above this threshold
#'   when ranking genes, default 0.95 (numeric),
#' @param scale_by_gene Scale each gene to have zero mean and unit std, default
#'   TRUE (logical),
#' @param z_max Maximum absolute value after scaling to cap extreme values,
#'   default 3 (numeric).
#' @param gene_annot Annotation for genes (rows) passed to
#'   `pheatmap::pheatmap()`, default NULL.
#' @param cell_annot Annotation for cells (columns) passed to
#'   `pheatmap::pheatmap()`, default NULL.
#' @param show_rownames Show row (gene) names in heatmap, default FALSE, can
#'   selectively show row names if a vector of names is given (logical or vector
#'   of character).
#' @param show_arrows If show_rownames == TRUE, show lines linking gene names
#'   and actual row positions, default FALSE (logical).
#' @param n_genes_per_level Max of number of genes per level of a categorical
#'   variable of cell annotation, e.g. list(celltype=20), can only specify one
#'   variable, default NULL (list)
#' @param smooth_heatmap Smooth heatmap using a kernal function, the larger the
#'   value the more smoothing, 0 disable smoothing default 3 (integer).
#' @param palette Palette function for heatmap, default
#'   `colorRampPalette(rev(brewer.pal(n = 100, name = "RdBu")))` (function).
#' @param annotation_colors Passed to pheatmap::pheatmap(), default NA.
#' @param ... Additional options passed to `pheatmap::pheatmap()`.
#'
#' @return A (re-)ordered matrix for plotting
plot_heatmap <- function(cds, genes, min_gene_sd = 0.2, min_gene_fraction=0.05,
                         n_cells = 500, order_genes_by = "Spectral",
                         rank_threshold = NULL, scale_by_gene = TRUE, z_max = 3,
                         gene_annot = NULL, cell_annot = NULL,
                         show_rownames = FALSE, show_arrows = FALSE,
                         n_genes_per_level = NULL,
                         smooth_heatmap = 3,
                         palette=ifelse(
                           scale_by_gene,
                           colorRampPalette(rev(RColorBrewer::brewer.pal(
                             n = 11, name = "RdBu")
                           )),
                           colorRampPalette(
                             RColorBrewer::brewer.pal(n = 9, name = "Reds"))
                         ),
                         annotation_colors = NA,
                         ...) {
  if (requireNamespace("monocle3", quietly = TRUE)) {
    if (class(cell_annot) == "character") {
      cell_mdat <- as.data.frame(SummarizedExperiment::colData(cds)[cell_annot])
    } else if (any(class(cell_annot) == "data.frame")) {
      cell_mdat <- cell_annot
    } else {
      stop("invalid `cell_annot` provided")
    }

    cell_mdat$pseudotime <- monocle3::pseudotime(cds)
    cell_order <- order(cell_mdat$pseudotime)

    X <- monocle3:::normalize_expr_data(cds, norm_method = "log")[genes, ]

    set.seed(1)
    if (ncol(X) > n_cells) {
      cell_idx <- sort(sample(ncol(X), size = n_cells, replace = FALSE))
    } else {
      cell_idx <- seq_len(ncol(X))
    }

    cell_mdat <- cell_mdat[cell_order[cell_idx], ]
    X <- X[, cell_order[cell_idx]]

    if (min_gene_sd > 0) {
      X_sd <- apply(X, 1, sd)
      X <- X[!is.na(X_sd) & (X_sd >= min_gene_sd), ]
    }

    if (min_gene_fraction > 0) {
      X_frac1 <- Matrix::rowMeans(X > 0)
      X <- X[(X_frac1 >= min_gene_fraction) & (X_frac1 <= 0.99), ]
    }

    if (scale_by_gene) {
      X <- t(scale(t(X)))
      X[is.na(X)] <- 0
      X[X > z_max] <- z_max
      X[X < -z_max] <- -z_max
    } else {
      X <- as.matrix(X)
    }
    if (smooth_heatmap[1] > 0) {
      X1 <- cbind(
        matrix(rep(X[, 1], smooth_heatmap[1]), ncol=smooth_heatmap[1]),
        X,
        matrix(rep(X[, ncol(X)], smooth_heatmap[1]), ncol=smooth_heatmap[1])
      )
      kn <- stats::kernel("daniell", smooth_heatmap[1])
      X <- t(stats::kernapply(t(X1), kn))
    }

    rank_threshold <- ifelse(is.null(rank_threshold), 0.95, rank_threshold)
    gene_ranks <- lapply(1:nrow(X), function(i) {
      which(X[i, ] > quantile(X[i, ], rank_threshold))
    })
    names(gene_ranks) <- rownames(X)
    gene_rank <- as.integer(sapply(gene_ranks, median))
    names(gene_rank) <- rownames(X)

    assign_group <- function(x, groups, method="majority") {
      method <- match.arg(method, c("majority", "median"))
      if (method == "majority") {
        grp_name <- names(sort(-table(groups[x]))[1])
      } else if (method == "median") {
        grp_name <- groups[as.integer(median(x))]
      }
      grp_name
    }

    if (class(n_genes_per_level) == "list" && !is.null(cell_annot)) {
      grp_var <- names(n_genes_per_level)[1]
      groups <- cell_mdat[, grp_var]
      max_n <- n_genes_per_level[[1]]
      tmp_df <- data.frame(
        gene = rownames(X),
        rank = seq_len(nrow(X)),
        annot = sapply(gene_ranks, assign_group, groups)
      )
      selected_genes <- dplyr::slice_min(dplyr::group_by(tmp_df, annot), order_by = rank, n = max_n)$gene
      selected_k <- rownames(X) %in% selected_genes
      X <- X[selected_k, ]
      gene_rank <- gene_rank[selected_k]
    }

    if (order_genes_by != FALSE) {
      nr <- nrow(X)
      nc <- ncol(X)
      if (order_genes_by == "corr") {
        corrs <- sapply(1:nr, function(i) cor(X[i, ], 1:nc))
        o <- order(corrs)
      } else if (order_genes_by == "rank") {
        o <- order(gene_rank)
      } else {
        if (requireNamespace("parallelDist", quietly = TRUE)) {
          d <- parallelDist::parDist(X)
        } else {
          d <- dist(X)
        }
        o <- seriation::get_order(
          seriation::seriate(d, method = order_genes_by)
        )
        if (cor(X[o[1], ], seq_len(ncol(X))) > 0) o <- rev(o)
      }
      X <- X[o, ]
      gene_rank <- gene_rank[o]
    }

    if (length(smooth_heatmap) > 1 && smooth_heatmap[2] > 0) {
      X1 <- rbind(
        matrix(
          rep(X[1, ], smooth_heatmap[2]), nrow = smooth_heatmap[2], byrow = TRUE
        ),
        X,
        matrix(
          rep(X[nrow(X), ], smooth_heatmap[2]),
          nrow = smooth_heatmap[2], byrow = TRUE
        )
      )
      kn <- stats::kernel("daniell", smooth_heatmap[2])
      X <- stats::kernapply(X1, kn)
      z_max <- min(3, max(X))
    }

    if (requireNamespace("pheatmap", quietly = TRUE)) {
      if (scale_by_gene) {
        breaks <- seq(-z_max, z_max, length.out = 101)
      } else {
        breaks <- NA
      }
      rownames_to_show <- NULL
      if (class(show_rownames) == "character") {
        rownames_to_show <- show_rownames
        show_rownames <- TRUE
      }
      p <- pheatmap::pheatmap(
        X, cluster_rows = FALSE, cluster_cols = FALSE,
        annotation_col = cell_mdat, show_rownames = show_rownames,
        show_colnames = FALSE, treeheight_row = 0, treeheight_col = 0,
        border_color = NA, color = palette(100), breaks = breaks,
        annotation_colors = annotation_colors,
        silent = TRUE, ...)
      if (show_rownames && !is.na(annotation_colors)[1]) {
          grp_var <- names(annotation_colors)[1]
          groups <- as.character(cell_mdat[, grp_var])
          gene_label_colors <- annotation_colors[[grp_var]][
            sapply(gene_ranks[names(gene_rank)], assign_group, groups)
          ]
          p$gtable$grobs[[which(p$gtable$layout$name == "row_names")]]$gp$col <- gene_label_colors
      }
      if (!is.null(rownames_to_show)) {
        p <- .add_flag(p, rownames_to_show, 1, add_pointer = show_arrows)
      }
      # plot result
      grid::grid.newpage()
      grid::grid.draw(p)

    } else {
      p <- image(
        X, col = palette(100), xaxt = "n", yaxt = "n", useRaster = TRUE)
    }

    invisible(list(X = X, cell_mdat = cell_mdat, plot = p))
  } else {
    stop("package:monocle3 not available")
  }
}
