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
#'
#' @return tibble containing metadata for each sample, one per row
#' 
#' @import dplyr
.make_sample_metadata <- function(obs_tbl, colSample, colVarCats=c(), colVarNums=c()) {
    metadata_tbl <- obs_tbl %>% select(one_of(c(colSample, colVarCats, colVarNums))) %>% unique()
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
.make_input_for_glmer <- function(Y, metadata_tbl) {
    samples <- rownames(Y)
    celltypes <- colnames(Y)
    nSample <- length(samples)
    nCelltype <- length(celltypes)

    metadata_tbl <- metadata_tbl[match(metadata_tbl[[colSample]], samples), ]

    input_tbl <- bind_cols(
        metadata_tbl[rep(1:nSample, nCelltype), ],
        enframe(factor(rep(celltypes, rep(nSample, nCelltype))), name=NULL, value='Celltype')
    )

    input_tbl$Y <- as.vector(Y)
    cat('input prepared\n')

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
.make_formula <- function(colSample, colVarCats, colVarNums, extra_term=NULL) {
    terms <- c(
        colVarNums,
        sprintf('(1|%s)', c(colSample, colVarCats)),
        sprintf('(%s-1|Celltype)', colVarNums),
        sprintf('(1|%s:Celltype)', c(colSample, colVarCats))
    )
    formula_str <- paste(c('I(c(Y)) ~ (1|Celltype)', terms, extra_term), collapse='+')
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
#'
#' @return tibble containing mean random effects and respective LTSR
#' 
#' @import dplyr
#' @import tidyr
.getCondValLtsr <- function(ranef_tbl, vars=NULL, celltypes=NULL) {
    ranef_tbl <- ranef_tbl %>% filter(
        grepl(':Celltype$', grpvar)
    ) %>% mutate(
        grpvar=factor(sub(':Celltype$', '', grpvar)),
        grp=factor(ifelse(grepl(':.*:', grp), sub(':', ',', grp), as.character(grp)))
    ) %>% separate(
        grp, into=c('grpval', 'Celltype'), sep=':'
    ) %>% mutate(
        lfsr=pnorm(condval, 0, condsd)
    ) %>% mutate(
        lfsr=ifelse(lfsr>0.5, 1-lfsr, lfsr)
    ) %>% mutate(
        ltsr=1-lfsr
    ) %>% select(
        grpvar, grpval, Celltype, condval, ltsr
    )

    if (is.list(vars)) {
        ranef_tbl <- ranef_tbl %>% filter(
            grpvar %in% names(vars)
        ) %>% mutate(
            grpvar=factor(grpvar, levels=names(vars)),
            grpval=factor(grpval, levels=unlist(vars, use.names=F))
        )
    } else if (!is.null(vars)) {
        ranef_tbl <- ranef_tbl %>% filter(grpvar %in% vars)
    }

    if (!is.null(celltypes)) ranef_tbl <- ranef_tbl %>% filter(Celltype %in% celltypes)

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
CellTypeCompositionAnalysis <- function(
        obs_tbl, colSample, colCelltype, colVarCats, colVarNums, extra_term=NULL, save=NULL) {
    metadata_tbl <- .make_sample_metadata(obs_tbl, colSample, colVarCats, colVarNums)
    Y <- .make_count_matrix(obs_tbl, colSample, colCelltype)

    input_tbl <- .make_input_for_glmer(Y, metadata_tbl)

    f <- .make_formula(colSample, colVarCats, colVarNums, extra_term=extra_term)
    cat('model constructed\n')
    res.prop <- glmer(
        f, data=input_tbl, family=poisson,
        control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
    )
    cat('model fitted\n')

    # standard errors of standard deviations (squre root of the variance parameters)
    devfun <- update(res.prop, devFunOnly=T)
    pars <- getME(res.prop, c('theta', 'fixef'))
    hess <- hessian(devfun, unlist(pars))
    sdse.prop <- data.frame(sd=unlist(pars), se=sqrt(diag(solve(hess))))
    sdse_tbl <- sdse.prop %>% rownames_to_column() %>% as_tibble()
    # posterior means and their standard deviations
    res.prop.ranef <- ranef(res.prop)
    ranef_tbl <- data.frame(res.prop.ranef) %>% as_tibble()

    if (is.character(save)) {
        output_dir <- ifelse(grepl('[/]$', save), save, dirname(save))
        if (!dir.exists(output_dir)) dir.create(output_dir, recursive=T)
        output_sep <- ifelse(grepl('[/.]$', save), '', '.')
        write_tsv(ranef_tbl, paste(save, 'ranef.tsv', sep=output_sep))
        write_tsv(sdse_tbl, paste(save, 'sdse.tsv', sep=output_sep))
    }

    list(ranef=ranef_tbl, sdse=sdse_tbl)
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
plot_sdse <- function(sdse_tbl, colSample, ci=0.95) {
    n_se <- qnorm(1-(1-ci)/2)
    p_tbl <- sdse_tbl %>% filter(
        grepl(':Celltype.\\(Intercept\\)', rowname)
    ) %>% mutate(
        rowname=sub(paste0('^', colSample, '$'), 'Residual', sub('theta.(.*):Celltype.*', '\\1', rowname))
    ) %>% arrange(sd) %>% mutate(
        rowname=factor(rowname, levels=c('Residual', rowname[rowname != 'Residual']))
    ) %>% mutate(ci_l=sd-se*n_se, ci_h=sd+se*n_se)
    p <- (
        ggplot(p_tbl, aes(x=sd, y=rowname)) +
        geom_point(shape=18, size=3) +
        geom_errorbarh(aes(xmin=ci_l, xmax=ci_h), height=0) +
        geom_vline(xintercept=0, lty=2) +
        xlab('Square root of explained variance') +
        coord_cartesian(xlim=c(-0.5, 1.5)) +
        theme_bw() +
        theme(axis.title.y=element_blank())
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
#' @param maxFC Cap fold change by this value (numeric, default 3)
#'
#' @return ggplot object of a dot plot
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom scales squish
#' 
#' @export
plot_ranef <- function(ranef_tbl, vars=NULL, celltypes=NULL, maxFC=3) {
    ranef_tbl <- .getCondValLtsr(ranef_tbl, vars=vars, celltypes=celltypes)

    condval_mat <- ranef_tbl %>% select(
        Celltype, grpval, condval
    ) %>% spread(
        'grpval', 'condval'
    ) %>% column_to_rownames(
        var='Celltype'
    ) %>% as.matrix()
    dendy <- hclust(dist(condval_mat))
    ordered_celltype = rownames(condval_mat)[dendy$ord]

    ranef_tbl <- ranef_tbl %>% mutate(
        Celltype=factor(Celltype, levels=ordered_celltype),
        condval=condval %>% pmin(log(maxFC)) %>% pmax(log(1/maxFC)),
        ltsr=ltsr %>% pmin(0.9999) %>% pmax(0.5)
    )

    p <- (
        ggplot(ranef_tbl, aes(x=grpval, y=Celltype)) +
        facet_grid(.~grpvar, scales='free_x', space='free_x', switch='x') +
        geom_point(aes(color=log2(exp(condval)), size=-log10(1-ltsr))) +
        scale_color_distiller(
            palette='RdBu', limits=log2(c(1/maxFC, maxFC)), oob=squish,
            guide=guide_colorbar(title='Log2FC', barwidth=1)) +
        scale_size(
            limits=-log10(1-c(0.5, 0.9999)),
            breaks=-log10(1-c(0.5, 0.9, 0.99, 0.999, 0.9999)),
            labels=c('0.5', '0.9', '0.99', '0.999', '>0.9999'),
            guide=guide_legend(title='LTSR')) +
        theme_bw() +
        theme(
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
            axis.title.x=element_blank(),
            strip.placement='outside',
            strip.background=element_blank()
        )
    )
}
