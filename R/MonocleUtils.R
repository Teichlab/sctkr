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
      node_pool <- names(monocle3:::leaf_nodes(cds))
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
#'   non-variable genes. Default 0.2 (numeric).
#' @param n_cells Number of cells to plot, down-sample to this number if fewer
#'   then `ncol(cds)`, default 500 (integer).
#' @param order_genes_by Method to re-order genes: "ARSA", "GW", "OLO", "TSP",
#'   "HC", "HC_single", "HC_complete", "HC_average", "HC_ward", "Spectral",
#'   "Spectral_norm", "corr" FALSE to disalbe ordering, default "Spectral" if
#'   `package:seriation` is available otherwise "HC" (character|bool).
#' @param rank_threshold Take only cells with expression above this threshold
#'   when ranking genes, default 1.5 if `scale_by_gene = TRUE` else 4 (numeric),
#' @param scale_by_gene Scale each gene to have zero mean and unit std, default
#'   TRUE (logical),
#' @param z_max Maximum absolute value after scaling to cap extreme values,
#'   default 3 (numeric).
#' @param gene_annot Annotation for genes (rows) passed to
#'   `pheatmap::pheatmap()`, default NULL.
#' @param cell_annot Annotation for cells (columns) passed to
#'   `pheatmap::pheatmap()`, default NULL.
#' @param show_rownames Show row (gene) names in heatmap, default FALSE
#'   (logical).
#' @param n_genes_per_level Max of number of genes per level of a categorical
#'   variable of cell annotation, e.g. list(celltype=20), can only specify one
#'   variable, default NULL (list)
#' @param palette Palette function for heatmap, default
#'   `colorRampPalette(rev(brewer.pal(n = 100, name = "RdBu")))` (function).
#' @param ... Additional options passed to `pheatmap::pheatmap()`.
#'
#' @return A (re-)ordered matrix for plotting
plot_heatmap <- function(cds, genes, min_gene_sd = 0.2, min_gene_fraction=0.1,
                         n_cells = 500, order_genes_by = "Spectral",
                         rank_threshold = NULL, scale_by_gene = TRUE, z_max = 3,
                         gene_annot = NULL, cell_annot = NULL, show_rownames = FALSE,
                         n_genes_per_level = NULL,
                         palette=ifelse(
                           scale_by_gene,
                           colorRampPalette(rev(RColorBrewer::brewer.pal(
                             n = 11, name = "RdBu")
                           )),
                           colorRampPalette(
                             RColorBrewer::brewer.pal(n = 9, name = "Reds"))
                         ), ...) {
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

    X <- SummarizedExperiment::assay(cds[genes, ])

    if (ncol(X) > n_cells) {
      cell_idx <- sort(sample(ncol(X), size = n_cells, replace = FALSE))
    } else {
      cell_idx <- seq_len(ncol(X))
    }

    cell_mdat <- cell_mdat[cell_order[cell_idx], ]
    X <- X[, cell_order[cell_idx]]

    X <- log1p(X / (Matrix::rowSums(X) / 1e4))

    if (min_gene_sd > 0) {
      X_sd <- apply(X, 1, sd)
      X <- X[!is.na(X_sd) & (X_sd >= min_gene_sd), ]
    }

    if (min_gene_fraction > 0) {
      X_frac1 <- Matrix::rowMeans(X > 0)
      X <- X[(X_frac1 >= min_gene_fraction) & (X_frac1 <= 0.95), ]
    }

    if (scale_by_gene) {
      X <- t(scale(t(X)))
      X[is.na(X)] <- 0
      X[X > z_max] <- z_max
      X[X < -z_max] <- -z_max
    } else {
      X <- as.matrix(X)
    }

    rank_threshold <- ifelse(is.null(rank_threshold), 0.8, rank_threshold)
    gene_rank <- sapply(1:nrow(X), function(i) {
      median(which(X[i, ] > quantile(X[i, ], rank_threshold)))
    })

    if (class(n_genes_per_level) == "list" && !is.null(cell_annot)) {
      grp_var <- names(n_genes_per_level)[1]
      groups <- cell_mdat[, grp_var]
      max_n <- n_genes_per_level[[1]]
      tmp_df <- data.frame(gene=rownames(X), rank=seq_len(nrow(X)), annot=groups[gene_rank])
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
    }

    if (requireNamespace("pheatmap", quietly = TRUE)) {
      if (scale_by_gene) {
        breaks <- seq(-z_max, z_max, length.out = 101)
      } else {
        breaks <- NA
      }
      p <- pheatmap::pheatmap(
        X, cluster_rows = FALSE, cluster_cols = FALSE,
        annotation_col = cell_mdat, show_rownames = show_rownames,
        show_colnames = FALSE, treeheight_row = 0, treeheight_col = 0,
        border_color = NA, color = palette(100), breaks = breaks, ...)
    } else {
      p <- image(
        X, col = palette(100), xaxt = "n", yaxt = "n", useRaster = TRUE)
    }

    invisible(list(X = X, cell_mdat = cell_mdat, plot = p))
  } else {
    stop("package:monocle3 not available")
  }
}
