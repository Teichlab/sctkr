# Method and initial implementation by Natsuhiko Kumasaka (natsuhiko@github)
# refactored and maintained by Ni Huang (nh3@github)

#' Make cell type by sample count matrix
#'
#' This function makes cell type by sample count matrix from a cell level annotation table
#'
#' @param obs_tbl Cell level annotation table, such as AnnData.obs (data.frame-like object)
#' @param colSample Column in "obs_tbl" specifying sample ID (str)
#' @param colCelltype Column in "obs_tbl" specifying cell type annotation (str)
#'
#' @return matrix containing cell type counts of dimension nSample x nCelltype
#'
#' @import magrittr
.make_count_matrix <- function(obs_tbl, colSample, colCelltype) {
  cnt_mat <- table(obs_tbl[[colSample]], obs_tbl[[colCelltype]]) %>% as.matrix()
  cnt_mat
}

#' Make sample level metadata table
#'
#' This function makes sample level metadata table from a cell level annotation table
#'
#' @param obs_tbl Cell level annotation table, such as AnnData.obs (data.frame-like object)
#' @param colSample Column in "obs_tbl" specifying sample ID (str)
#' @param colVarCats Columns in "obs_tbl" specifying categorical variables to control for (vector of str)
#' @param colVarNums Columns in "obs_tbl" specifying numerical variables to control for (vector of str)
#' @param extra_term Extra term(s) to be added to the model (str or vector of str)
#'
#' @return tibble containing metadata for each sample, one per row
#'
#' @import dplyr
.make_sample_metadata <- function(obs_tbl, colSample, colVarCats = c(), colVarNums = c(), extra_term = NULL) {
  if (!is.null(extra_term)) {
    extra_vars <- unique(unlist(strsplit(extra_term, ":", fixed = T)))
    extra_vars <- extra_vars[extra_vars %in% colnames(obs_tbl)]
  } else {
    extra_vars <- c()
  }
  metadata_tbl <- obs_tbl %>%
    select(one_of(unique(c(colSample, colVarCats, colVarNums, extra_vars)))) %>%
    unique()
  metadata_tbl
}

#' Make input data for model fitting
#'
#' This function makes input data for model fitting by assembling sample level
#' metadata and cell type by sample count matrix
#'
#' @param Y Cell type by sample count matrix (matrix of dimension nSample x nCelltype)
#' @param metadata_tbl Sample level metadata table (data.frame-like object)
#'
#' @return tibble containing counts and metadata for each cell type in each sample
#'
#' @import tibble
#' @import dplyr
.make_input_for_glmer <- function(Y, metadata_tbl, colSample) {
  samples <- rownames(Y)
  celltypes <- colnames(Y)
  nSample <- length(samples)
  nCelltype <- length(celltypes)

  metadata_tbl <- metadata_tbl[match(samples, metadata_tbl[[colSample]]), ]

  input_tbl <- bind_cols(
    metadata_tbl[rep(1:nSample, nCelltype), ],
    enframe(factor(rep(celltypes, rep(nSample, nCelltype))), name = NULL, value = "Celltype")
  )

  input_tbl$Y <- as.vector(Y)
  cat("input prepared\n")

  input_tbl
}

#' Make model formula
#'
#' This function makes a formula for a GLMM with random effects for each
#' specified categorical variables and interactions between each variable and
#' cell type
#'
#' @param colSample Column in "obs_tbl" specifying sample ID (str)
#' @param colVarCats Columns in "obs_tbl" specifying categorical variables to control for (vector of str)
#' @param colVarNums Columns in "obs_tbl" specifying numerical variables to control for (vector of str)
#' @param extra_term Extra term(s) to be added to the model (str or vector of str)
#'
#' @return formula specifying a GLMM model
.make_formula <- function(colSample, colCelltype, colVarCats, colVarNums, extra_term = NULL) {
  if (!is.null(extra_term)) {
    extra_vars <- unique(unlist(strsplit(extra_term, ":", fixed = T)))
    extra_vars <- extra_vars[extra_vars %in% colnames(obs_tbl)]
    extra_vars[extra_vars == colCelltype] <- "Celltype"
    extra_term <- sprintf("(1|%s)", paste(extra_vars, collapse = ":"))
  }
  terms <- c(
    colVarNums,
    sprintf("(1|%s)", c(colSample, colVarCats)),
    sprintf("(%s-1|Celltype)", colVarNums),
    sprintf("(1|%s:Celltype)", c(colSample, colVarCats))
  )
  formula_str <- paste(c("I(Y) ~ (1|Celltype)", terms, extra_term), collapse = "+")
  cat(paste0(formula_str, "\n"))
  as.formula(formula_str)
}

#' Get mean estimate of random effects and LTSR values for all specified variables
#'
#' This function extracts random effect estimates (mean and sd) of specified
#' variables, calculates LTSR using a z-score test and prepare categories in
#' specified orders for plotting
#'
#' @param ranef_tbl Table of random effect estimates (mean and sd) extracted from fitted model
#'     (data.frame-like object)
#' @param vars A list specifying required categorical sample metadata variables and their order,
#'     e.g. list(size=c('small', 'medium', 'big'), price=c('cheap', 'expensive'))
#' @param celltypes If specified, only keep those cell types (vector of str)
#' @param references If specified, fold changes are calculated relative to the specified level
#'     instead of global mean. Must be a list, e.g. list(size='small', price='cheap')
#'
#' @return tibble containing mean random effects and respective LTSR
#'
#' @import dplyr
#' @import tidyr
#' @import stringr
.getCondValLtsr <- function(ranef_tbl, vars, celltypes = NULL, references = NULL) {
  if (!is.list(vars)) stop('"vars" must be a list of which names are co-variables to plot', call. = F)
  if (!is.null(references) && !is.list(references)) stop('"references" must be a list of which names match that of "vars"', call. = F)

  cat_ranef_tbl <- ranef_tbl %>%
    filter(
      grepl(":Celltype$", grpvar)
    ) %>%
    mutate(
      grpvar = factor(sub(":Celltype$", "", grpvar)),
      grp = factor(ifelse(
        grepl(":.*:", grp),
        sapply(str_split(grp, ":"), function(x) paste(paste(x[-length(x)], collapse = ","), x[length(x)], sep = ":")),
        as.character(grp)
      ))
    ) %>%
    separate(
      grp,
      into = c("grpval", "Celltype"), sep = ":"
    )

  num_ranef_tbl <- ranef_tbl %>%
    filter(
      grpvar == "Celltype" & term != "(Intercept)"
    ) %>%
    mutate(
      grpvar = factor(term),
      grpval= "Slope",
      Celltype = grp
    ) %>%
    select (
      -grp
    )

  if (!is.null(references)) {
    cat_ranef_tbl <- bind_rows(
      lapply(names(vars), function(vname) {
        if (vname %in% names(references)) {
          ref <- references[[vname]]
          full_join(
            cat_ranef_tbl %>% filter(grpvar == vname, grpval != ref),
            cat_ranef_tbl %>% filter(grpvar == vname, grpval == ref) %>% select(grpvar, Celltype, condval, condsd),
            by = c("grpvar", "Celltype")
          ) %>%
            mutate(
              condval = condval.x - condval.y,
              condsd = sqrt(condsd.x^2 + condsd.y^2)
            ) %>%
            select(-condval.x, -condval.y, -condsd.x, -condsd.y)
        } else {
          return(cat_ranef_tbl %>% filter(grpvar == vname))
        }
      })
    )
  }

  ranef_tbl <- bind_rows(cat_ranef_tbl, num_ranef_tbl)

  ranef_tbl <- ranef_tbl %>%
    mutate(
      lfsr = pnorm(condval, 0, condsd)
    ) %>%
    mutate(
      lfsr = ifelse(lfsr > 0.5, 1 - lfsr, lfsr)
    ) %>%
    mutate(
      ltsr = 1 - lfsr
    ) %>%
    select(
      grpvar, grpval, Celltype, condval, ltsr
    )

  if (!is.null(celltypes)) ranef_tbl <- ranef_tbl %>% filter(Celltype %in% celltypes)

  ranef_tbl <- ranef_tbl %>%
    filter(
      grpvar %in% names(vars)
    ) %>%
    mutate(
      grpvar = factor(grpvar, levels = names(vars)),
      grpval = factor(grpval, levels = unlist(vars, use.names = F))
    )

  ranef_tbl
}

#' Cell type composition analysis
#'
#' This function performs cell type composition analysis across samples grouped by certain metadata
#' variable while controling for other variables, by fitting a Poisson Generalised Linear Mixed
#' Model
#'
#' @param obs_tbl Cell level annotation table, such as AnnData.obs (data.frame-like object)
#' @param colSample Column in "obs_tbl" specifying sample ID (str)
#' @param colCelltype Column in "obs_tbl" specifying cell type annotation (str)
#' @param colVarCats Columns in "obs_tbl" specifying categorical variables to control for (vector of str)
#' @param colVarNums Columns in "obs_tbl" specifying numerical variables to control for (vector of str)
#' @param extra_term Extra term(s) to be added to the model (str or vector of str)
#' @param save If specified, used as file name prefix to save returned values (str)
#'
#' @return list containing two tables, 1) estimated random effects, 2) estimated explained variance
#'     (in the form of standard deviation)
#'
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @import lme4
#' @import numDeriv
#'
#' @export
CellTypeCompositionAnalysis <- function(obs_tbl, colSample, colCelltype, colVarCats, colVarNums = NULL, extra_term = NULL, save = NULL) {
  metadata_tbl <- .make_sample_metadata(obs_tbl, colSample, colVarCats, colVarNums, extra_term = extra_term)
  Y <- .make_count_matrix(obs_tbl, colSample, colCelltype)

  input_tbl <- .make_input_for_glmer(Y, metadata_tbl, colSample)

  f <- .make_formula(colSample, colCelltype, colVarCats, colVarNums, extra_term = extra_term)
  cat("model constructed\n")
  res.prop <- glmer(
    f,
    data = input_tbl, family = poisson,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )
  cat("model fitted\n")

  # standard errors of standard deviations (squre root of the variance parameters)
  devfun <- update(res.prop, devFunOnly = T)
  pars <- getME(res.prop, c("theta", "fixef"))
  hess <- hessian(devfun, unlist(pars))
  sdse.prop <- data.frame(sd = unlist(pars), se = sqrt(diag(solve(hess))))
  sdse_tbl <- sdse.prop %>%
    rownames_to_column() %>%
    as_tibble()
  # posterior means and their standard deviations
  res.prop.ranef <- ranef(res.prop)
  ranef_tbl <- data.frame(res.prop.ranef) %>% as_tibble()

  if (is.character(save)) {
    output_dir <- ifelse(grepl("[/]$", save), save, dirname(save))
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = T)
    output_sep <- ifelse(grepl("[/.]$", save), "", ".")
    write_tsv(ranef_tbl, paste(save, "ranef.tsv", sep = output_sep))
    write_tsv(sdse_tbl, paste(save, "sdse.tsv", sep = output_sep))
  }

  list(ranef = ranef_tbl, sdse = sdse_tbl)
}

#' Plot variable explained variance
#'
#' This function plots (square root of) variance explained by variables as a forest plot
#'
#' @param sdse_tbl Table containing square-root of estimated explained variance and standard errors
#'     (data.frame-like object)
#' @param colSample Column in "obs_tbl" specifying sample ID (str)
#' @param ci Confidence interval level (numeric between 0 and 1, default 0.95)
#'
#' @return ggplot object of a forest plot
#'
#' @import dplyr
#' @import ggplot2
#'
#' @export
plot_sdse <- function(sdse_tbl, colSample, ci = 0.95, xlim = c(-0.5, 1.5)) {
  n_se <- qnorm(1 - (1 - ci) / 2)
  p_tbl <- sdse_tbl %>%
    filter(
      grepl(":Celltype.\\(Intercept\\)", rowname)
    ) %>%
    mutate(
      rowname = sub(paste0("^", colSample, "$"), "Residual", sub("theta.(.*):Celltype.*", "\\1", rowname))
    ) %>%
    arrange(sd) %>%
    mutate(
      rowname = factor(rowname, levels = c("Residual", rowname[rowname != "Residual"]))
    ) %>%
    mutate(ci_l = sd - se * n_se, ci_h = sd + se * n_se)
  p <- (
    ggplot(p_tbl, aes(x = sd, y = rowname)) +
      geom_point(shape = 18, size = 3) +
      geom_errorbarh(aes(xmin = ci_l, xmax = ci_h), height = 0) +
      geom_vline(xintercept = 0, lty = 2) +
      xlab("Square root of explained variance") +
      coord_cartesian(xlim = xlim) +
      theme_bw() +
      theme(axis.title.y = element_blank())
  )
}

#' Plot estimated random effects of variables
#'
#' This function plots random effects of variables on cell type composition as a dot plot
#'
#' @param ranef_tbl Table of random effect estimates (mean and sd) extracted from fitted model
#'     (data.frame-like object)
#' @param vars A list specifying required categorical sample metadata variables and their order,
#'     e.g. list(size=c('small', 'medium', 'big'), price=c('cheap', 'expensive'))
#' @param celltypes If specified, only keep those cell types (vector of str)
#' @param celltype_order Determine how cell types are ordered. Either "hclust" which orders by
#'     hierarchical clustering, or a vector of cell types. (str or vector of str)
#' @param references If specified, fold changes are calculated relative to the specified level
#'     instead of global mean. Must be a list, e.g. list(size='small', price='cheap')
#' @param maxFC Cap fold change by this value (numeric, default 3)
#' @param LTSR2p Whether to convert LTSR to p value (logical, default FALSE)
#' @param filterLtsr Only keeps cell types of which at least one LTSR is greater or equal to
#'     specified value. (numeric, default 0)
#' @param swap_axes Whether to swap axis when plotting (logical, default FALSE)
#'
#' @return ggplot object of a dot plot
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom scales squish
#'
#' @export
plot_ranef <- function(ranef_tbl, vars, celltypes = NULL, celltype_order = "hclust", references = NULL,
                       maxFC = 3, LTSR2p = F, highlightLtsr = 0.0, filterLtsr = 0.0, swap_axes = F) {
  ranef_tbl <- .getCondValLtsr(ranef_tbl, vars, celltypes = celltypes, references = references)

  condval_mat <- ranef_tbl %>%
    select(
      Celltype, grpval, condval
    ) %>%
    spread(
      "grpval", "condval"
    ) %>%
    column_to_rownames(
      var = "Celltype"
    ) %>%
    as.matrix()
  if (length(celltype_order) == 1 && celltype_order == "hclust") {
    dendy <- hclust(dist(condval_mat))
    ordered_celltype <- rownames(condval_mat)[dendy$ord]
  } else if (!is.null(celltype_order) && length(celltype_order) == dim(condval_mat)[1]) {
    ordered_celltype <- celltype_order
  }

  ranef_tbl <- ranef_tbl %>% mutate(
    Celltype = factor(Celltype, levels = ordered_celltype),
    condval = condval %>% pmin(log(maxFC)) %>% pmax(log(1 / maxFC)),
    ltsr = ltsr %>% pmin(0.9999) %>% pmax(0.5)
  )

  if (swap_axes) {
    ranef_tbl$Celltype <- factor(ranef_tbl$Celltype, levels = rev(levels(ranef_tbl$Celltype)))
    ranef_tbl$grpval <- factor(ranef_tbl$grpval, levels = rev(levels(ranef_tbl$grpval)))
  }

  if (filterLtsr > 0) {
    filtered_celltypes <- ranef_tbl %>%
      group_by(Celltype) %>%
      summarise(maxLtsr = max(ltsr)) %>%
      dplyr::filter(maxLtsr >= filterLtsr) %>%
      select(Celltype) %>%
      unlist(use.names = F)
    ranef_tbl <- ranef_tbl %>% dplyr::filter(Celltype %in% filtered_celltypes)
  }

  geom_dots <- geom_point(
    aes(
      fill = log2(exp(condval)),
      size = -log10(1 - ltsr)
    ),
    color = "white",
    shape = 21
  )

  if (swap_axes) {
    p <- (
      ggplot(ranef_tbl, aes(y = grpval, x = Celltype)) +
        facet_grid(grpvar ~ ., scales = "free_y", space = "free_y", switch = "x") +
        geom_dots
    )
  } else {
    p <- (
      ggplot(ranef_tbl, aes(x = grpval, y = Celltype)) +
        facet_grid(. ~ grpvar, scales = "free_x", space = "free_x", switch = "x") +
        geom_dots
    )
  }

  p <- (
    p + scale_fill_distiller(
      palette = "RdBu",
      limits = log2(c(1 / maxFC, maxFC)),
      breaks = log2(c(1 / maxFC, maxFC)),
      labels = c(paste0("1/", maxFC), maxFC),
      oob = squish,
      guide = guide_colorbar(
        title = "Fold change", title.position = "top", direction = "horizontal",
        barwidth = 5, barheight = 0.75, raster = F, order = 1)
    )
    + scale_size(
      limits = -log10(1 - c(0.5, 0.9999)),
      breaks = -log10(1 - c(0.5, 0.9, 0.99, 0.999, 0.9999)),
      range = c(0.5, 9),
      labels = ifelse(
        rep(LTSR2p, 5),
        c("0.5", "0.1", "0.01", "0.001", "<0.0001"),
        c("0.5", "0.9", "0.99", "0.999", ">0.9999")
      ),
      guide = guide_legend(
        title = ifelse(LTSR2p, "p", "LTSR"), reverse = T, order = 2,
        override.aes = list(fill = "black", color = "white")
      )
    )
    + theme_bw()
    + theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      strip.placement = "outside",
      strip.background = element_blank(),
      legend.spacing.y = unit(0.5, "line")
    )
  )

  if (highlightLtsr > 0) {
    p <- (
      p + geom_point(
        aes(
          color = ltsr > highlightLtsr,
          alpha = ltsr > highlightLtsr
        ),
        shape = 21, size = 9
      )
      + scale_color_manual(
        label = c(
          "", ifelse(LTSR2p, paste0("p < ", 1 - highlightLtsr), paste0("LTSR > ", highlightLtsr))
        ),
        values = c("white", "red"),
        guide = guide_legend(title = NULL, override.aes = list(size = 9), reverse = T, order = 3)
      )
      + scale_alpha_manual(values = c(0, 1), guide = F)
    )
  }

  p
}

#' Theme that mimics dotplots generated by Natsuhiko's code
#'
#' This function sets a number of options in ggplot2::theme() to mimic the style
#' of dotplots generated by Natsuhiko's code
#'
#' @param ... options passed to ggplot2::theme()
#'
#' @importFrom ggplot2 theme
#'
#' @export
theme_ctca <- function(...) {
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_blank(),
    strip.text = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    ...
  )
}
