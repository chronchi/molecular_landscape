# Author: Carlos Ronchi
# Date: 2021.07.22
# Email: carlos.ronchi@epfl.ch
# Description: A set of helper functions that are reused across the different
#   documents in order to decluter them.


#' Violin plot of pathway score given cohort and subgroups
#'
#' @param sum_exp Summarized experiment object. Used to fetch the scores
#'   that are saved on the col data slot
#' @param which_score A string. Which pathway score to select and plot 
#'   the violin plots.
#' @param clinical_variable A string. With which clinical variable
#'   to subgroup
#' @param color_by A string. A clinical variable to color the dots.
#' @param base_size An integer. Argument for ggplot2::theme_bw
#' @param point_size An integer. Argument for ggplot2::geom_point
#' @param title_plot. An integer Argument for ggplot2::labs
#' @return A plot. 
#' @examples
#' \dontrun{
#' plot_scores_vs_clinics(
#'     sum_exp  = metabric,
#'     which_score = "set_erpr", 
#'     clinical_variable = "er_status",
#'     color_by = "er_status"
#' )
#' }
plot_scores_vs_clinics <- function(
    sum_exp, which_score, clinical_variable, color_by, 
    base_size = 10,
    point_size = 1,
    title_plot = ""
){
    
    SummarizedExperiment::colData(sum_exp) %>% data.frame %>% 
        ggplot2::ggplot(
            aes_string(x = clinical_variable, y = which_score, color = color_by)
        ) + 
        ggplot2::geom_violin(color = "black") + 
        ggplot2::geom_jitter(size = point_size) +
        ggplot2::labs(title = title_plot, y = "Score") +
        ggplot2::theme_bw(base_size = base_size)
    
}

#' Umap projection plot from the specified cohort
#'
#' @inheritParams plot_scores_vs_clinics
#' @param cohort A string. Name of the cohort that is being passed. It is
#'   used for the title.
#' @param umap_projection A matrix. Output from uwot::umap.
#' @param alpha_val A double. Specified transparency for ggplot2::geom_point.
#' @return A plot. 
#' @examples
#' \dontrun{
#' plot_umap(
#'     sum_exp  = metabric,
#'     cohort = "METABRIC",
#'     color_by = "er_status",
#'     umap_projection = umap_metabric
#' )
#' }
plot_umap <- function(
    sum_exp, 
    cohort,
    color_by,
    umap_projection,
    base_size = 10,
    point_size = 1,
    alpha_val = 0.8
){
    
    umap_embedding <- umap_projection %>% data.frame %>%
        `colnames<-`(c("umap1", "umap2")) %>%
        tibble::rownames_to_column(var = "sample_name") %>%
        dplyr::inner_join(
            .,
            colData(sum_exp) %>% data.frame,
            by = "sample_name"
        )
    
    ggplot2::ggplot(
        umap_embedding,
        aes_string(x = "umap1", y = "umap2", color = color_by)
    ) + 
        ggplot2::geom_point(size = point_size, alpha = alpha_val) + 
        ggplot2::labs(title = toupper(cohort)) +
        ggplot2::theme_bw(base_size = base_size)
        
}

#' Forest plot of the fitted models
#'
#' @param fit A model fit from survival package. 
#' @param cohort A string. Name of cohort that will be used for title
#'    purposes.
#' @param names_coefficients A named list. Used for converting the name
#'   of the coefficients to something more readable in the forest plot.
#' @param patients A string. Specifies which patients were selected to
#'   calculate the model fit.
#' @param type_survival A string. Either OS (overall survival) or 
#'   RFS (recurrence free survival).
#' @param name_signature A string. Name of the score pathway used to 
#'   test in the survival fit.
#' @param path_to_save A string. Path to save the forest plot.
#' @param width An integer. Width of the plot
#' @param height An integer. Height of the plot
#' @param ... Other options for the function forestplot::forestplot
#' @examples
#' \dontrun{
#' forest_plot_fits(
#'     fit,
#'     cohort = "metabric",
#'     names_coefficients = c(metabric = c(
#'         "age" = "Age", "set_erpr" = "SET ER/PR"
#'     ),
#'     patients = "Endo Only, ER+ BC",
#'     type_survival = "OS",
#'     name_signature = "set_erpr",
#'     path_to_save = "results/plots",
#'     clip = c(0.2, 2)
#' )
#' }
forest_plot_fits <- function(
    fit, cohort, names_coefficients, 
    patients, type_survival,
    name_signature, path_to_save, 
    width = 6, height = 3, ...
){
    
    coefs_surv <- exp(coef(fit))
    confint_surv <- exp(confint(fit))
    base_data <- tibble(
        coefficient = names_coefficients[[cohort]][names(coefs_surv)],
        mean = coefs_surv, 
        lower = confint_surv[, 1],
        upper = confint_surv[, 2],
        HR = round(coefs_surv, 2) %>% as.character
    ) %>% dplyr::mutate(
        dplyr::across(
            where(is.numeric), 
            round, 
            2
        )
    )
    
    header <- tibble(
        coefficient = c("Coefficient"),
        HR = c("HR"),
        summary = TRUE
    )
    
    
    output_df <- dplyr::bind_rows(
        header,
        base_data
    )
    
    forest_plot <- output_df %>% forestplot::forestplot(
        labeltext = c(coefficient, HR), 
        hrzl_lines = list(
            "2" = gpar(col = "#444444", columns = 1:3, lwd = 1)
        ),
        is.summary = summary,
        xlog = TRUE, 
        col = fpColors(
            box = viridis::viridis(2)[1],
            line = "black",
            summary = viridis::viridis(3)[2]
        ),
        title = paste0(
            type_survival, " ",
            toupper(cohort), ", ", patients,
            "\n",
            "#Patients: ", fit$n,
            ", ",
            "#Events: ", fit$nevent
        ),
        boxsize = 0.2,
        graph.pos = 2,
        ...
    )
    
    filename <- paste0(
        path_to_save, "/", c("png", "pdf"), "/",
        cohort, "_", type_survival,
        "_", name_signature
    )
    
    
    # note that pdf and png have different arguments. the default
    # unit for pdf is inches and for png is pixels. for png when
    # specifying the unit, if inches, then specify also resolution. 
    pdf(file = paste0(filename[2], ".pdf"), width = width, height = height)
    plot(forest_plot)
    dev.off()
     
    png(
        filename = paste0(filename[1], ".png"), 
        width = width, 
        height = height,
        units = "in", 
        res = 300
    )
    plot(forest_plot)
    dev.off()
    
    forest_plot
}


#' Prepare the data for normalization
#'
#' @description 
#' `prepare_data_normalization` checks if the data structure
#'     is right and if it has the appropriate slots 
#' 
#' @details 
#' The function check if the object is a summarized experiment object
#' first, so the normalization procedures can be saved later on. It also
#' checks if the assay_to_use is available on the summarized experiment
#' object, without this it cannot proceed with the normalization. Lastly
#' it proceeds to subset the dataframe to only contain the necessary genes
#' to calculate the ranking obtained in the procedure described
#' in the paper.
#' 
#' @param sum_exp Summarized experiment object. It should contain a "rank"
#'     slot that will be used to calculate the average rank of the stable 
#'     genes.
#' @param assay_to_use A string. Which assay to use when calculating the
#'     average expression of the stable genes.
#' @param stable_genes A vector. Stable genes.
#' @param stable_genes A vector. Most variable genes determined by some 
#'     procedure.
#' @param verbose. An integer. If verbose equals to 1, then it prints the
#'     number of stable genes, total genes in the dataframe and total number
#'     of samples
#' @examples
#' \dontrun{
#' prepare_data_normalization(
#'     tcga,
#'     "logFPKM_TMM",
#'     c("GAPDH"),
#'     c("ESR1", "GREB1")
#' )
#' }
prepare_data_normalization <- function(
    sum_exp, 
    assay_to_use, 
    stable_genes,
    most_variable_genes,
    verbose = 1
){

    sum_exp_classes <- c(
        "SummarizedExperiment", 
        "RangedSummarizedExperiment"
    )
    
    if( !(is(sum_exp, "SummarizedExperiment")) ){
        stop(
            paste0(
                "Object is not from class(es) ", 
                paste(sum_exp_classes, collapse = ", "),
                ". Please check your data structure."
            )
        )
    }
    
    if( !(assay_to_use %in% SummarizedExperiment::assayNames(sum_exp)) ){
        stop(
            paste0(
                "Assay not available, check if you specified the right ", 
                "name or you included your assay in the data structure."
            )
        )
    }
    
    # subselect the necessary genes for downstream procedures
    sum_exp <- sum_exp[
        intersect(rownames(sum_exp), c(stable_genes, most_variable_genes)),
    ]
    
    if (verbose == 1){
        cat(paste0(
            "Total number of stable genes: ", 
            length(intersect(rownames(sum_exp), stable_genes)), "\n",
            "Total number of genes: ", nrow(sum_exp), "\n",
            "Number of samples: ", ncol(sum_exp), "\n"
        ))
    }
    
    # calculate the ranking that will be used in the other functions
    assay(sum_exp, "rank") <- singscore::rankGenes(
        as.matrix(assay(sum_exp, assay_to_use))
    )
    
    sum_exp
}


#' Calculate the average ranking and average expression of stable genes
#'
#' @param sum_exp Summarized experiment object. It should contain a "rank"
#'     slot that will be used to calculate the average rank of the stable 
#'     genes.
#' @param assay_to_use A string. Which assay to use when calculating the
#'     average expression of the stable genes.
#' @param stable_genes A vector. Stable genes.
#' @examples
#' \dontrun{
#' avg_ranking(tcga, "logFPKM_TMM", c("GAPDH"))
#' }
avg_ranking <- function(sum_exp, assay_to_use, stable_genes){

    available_stable_genes <- intersect(rownames(sum_exp), stable_genes)
    
    rank_stable_genes <- assay(sum_exp, "rank")[available_stable_genes, ]
    expr_stable_genes <- assay(sum_exp, assay_to_use)[available_stable_genes, ]
    
    
    values <- list(
        avg_ranking = colMeans(rank_stable_genes)/nrow(sum_exp),
        avg_expression = colMeans(expr_stable_genes)
    )
    
    values$sd_ranking <- sd(values$avg_ranking)
    values$sd_expression <- sd(values$avg_expression)
    
    values
}

#' Calculate the average expression and ranking based on stable genes
#'
#' @inheritParams avg_ranking
#' @param avg_ranking_sds A list. Output from the function avg_ranking
#' @examples
#' \dontrun{
#' calculate_norm_ranks(tcga, "logFPKM_TMM", output_avg_ranking)
#' }
calculate_norm_ranks <- function(sum_exp, assay_to_use, avg_ranking_sds){
        
    # for each patient calculate the division for each gene in each
    # element of average expression obtained from the avg_ranking function
    df <- sapply(
        1:ncol(sum_exp),
        function(i, sum_exp, avg_ranking_sds){
            expression_levels <- assay(sum_exp, assay_to_use)[,i]
            expression_levels/avg_ranking_sds$avg_expression[i]
        },
        sum_exp = sum_exp,
        avg_ranking_sds = avg_ranking_sds
    ) %>% 
        `colnames<-`(colnames(sum_exp)) %>%
        `rownames<-`(rownames(sum_exp))
    
    # we now add the normalized average expression to the original summarized
    # experiment so it can be used to other analysis
    assay(sum_exp, "avg_expression") <- df
    
    # notice here we are converting the average normalized ranking back to 
    # average rank, since we will divide the rank of the gene and not
    # the normalized rank. Also all this operations are sums, so we 
    # can go back and forth here.
    df <- sapply(
        1:ncol(sum_exp),
        function(i, sum_exp, avg_ranking_sds){
            rank_genes <- assay(sum_exp, "rank")[,i]
            rank_genes/(avg_ranking_sds$avg_ranking[i]*nrow(sum_exp))
        },
        sum_exp = sum_exp,
        avg_ranking_sds = avg_ranking_sds
    ) %>% 
        `colnames<-`(colnames(sum_exp)) %>%
        `rownames<-`(rownames(sum_exp))
    
    # we now add the normalized average ranking to the original summarized
    # experiment so it can be used to other analysis
    assay(sum_exp, "avg_ranking") <- df
    
    sum_exp$avg_ranking <- avg_ranking_sds$avg_ranking
    sum_exp$avg_expression <- avg_ranking_sds$avg_expression
    
    sum_exp
    
}

#' Wrapper for prepare_data_normalization avg_ranking and calculate_norm_ranks
#'
#' @inheritParams prepare_data_normalization
#' @examples
#' \dontrun{
#' get_final_ranking_values(
#'     tcga,
#'     "logFPKM_TMM",
#'     c("GAPDH"),
#'     c("ESR1", "GREB1")
#' )
#' }
get_final_ranking_values <- function(
    sum_exp, 
    assay_to_use,
    stable_genes,
    most_variable_genes,
    verbose = 1
){
    
    sum_exp <- prepare_data_normalization(
        sum_exp = sum_exp,
        assay_to_use = assay_to_use,
        stable_genes = stable_genes,
        most_variable_genes = most_variable_genes,
        verbose = verbose
    )
    
    avg_ranking_sds <- avg_ranking(
        sum_exp = sum_exp, 
        assay_to_use = assay_to_use, 
        stable_genes = stable_genes
    )
    
    sum_exp <- calculate_norm_ranks(
        sum_exp = sum_exp,
        avg_ranking_sds = avg_ranking_sds,
        assay_to_use = assay_to_use
    )
    
    sum_exp
    
}

#' Extract gene expression levels and add to colData
#'
#' @param sum_exp A summarized experiment object
#' @param genes A character vector. A vector containing the name of genes
#'     of interest to extract the values
#' @param assay_to_use A string. Specify from which assay slot to extract
#'     the gene values
#' @return a dataframe
#' @examples
#' \dontrun{
#' get_gene_col_data(
#'     tcga,
#'     c("ESR1", "GREB1"),
#'     "logFPKM_TMM"
#' )
#' }
get_gene_col_data <- function(
    sum_exp,
    genes,
    assay_to_use
){
    common_genes <- intersect(rownames(sum_exp), genes)
    
    if(length(common_genes) != length(genes)) {
        print("Not all genes are available")
    }
    
    colData(sum_exp)[, common_genes] <- 
        assay(sum_exp, assay_to_use)[common_genes, ] %>%
        as.matrix %>% t
    
    colData(sum_exp) %>% data.frame

}

#' Get the PCA projection coordinates for a dataset
#'
#' @param sum_exp A summarized experiment object
#' @param pca_fit A PCA fit returned from PCATools. This object is the one
#'     trained using the samples from TCGA and METABRIC.
#' @return A dataframe. The PC coordinates are returned for the new
#'     dataset
#' @examples
#' \dontrun{
#' get_pca_coordinates(
#'     scanb,
#'     pca_fit
#' )
#' }
get_pca_coordinates <- function(
    sum_exp,
    pca_fit,
    which_exp = "avg_ranking"
){
    
    # first check if avg_ranking is there if it was
    # selected
    if (which_exp == "avg_ranking"){
        if( !("avg_ranking" %in% SummarizedExperiment::assayNames(sum_exp)) ){
            stop(
                paste0(
                    "avg_ranking not available in the object, make sure",
                    " the normalization was done."
                )
            )
        }    
    }
    
    
    loadings_pca <- pca_fit$loadings
    
    # we now add 0 to average ranking for the genes that are not
    # available in the summarized experiment 
    genes_for_pca <- rownames(loadings_pca)
    assay_matrix <- assay(sum_exp, which_exp) %>% as.matrix
    genes_not_available <- setdiff(genes_for_pca, rownames(assay_matrix))
    if (length(genes_not_available) > 0){
        assay_matrix <- rbind(
            assay_matrix,
            matrix(
                0, 
                nrow = length(genes_not_available), 
                ncol = ncol(sum_exp),
                dimnames = list(
                    genes_not_available,
                    colnames(assay_matrix)
                )
            )
        )
    }
    
    # calculate the pc coordinates
    assay_matrix <- assay_matrix[genes_for_pca, ]
    t(assay_matrix) %*% (loadings_pca %>% as.matrix)
    
}

#' Plot projection given the PCA coordinates
#'
#' @param df_pca A dataframe. This is a data frame containing the principal 
#'     components given in the get_pca_coordinates and any other 
#'     metadata related to the samples
#' @param color A string. ggplot2 aesthetic parameter
#' @param x A string. ggplot2 aesthetic parameter
#' @param y A string. ggplot2 aesthetic parameter
#' @param base_size An integer. ggplot2 aesthetic parameter
#' @param size An integer. ggplot2 aesthetic parameter
#' @param title A string. ggplot2 aesthetic parameter
#' @return A ggplot2 plot
#' @examples
#' \dontrun{
#' plot_pca_coordinates(
#'     df_pca, x = "PC1", y = "PC2", color = "cohort", base_size = 20
#' )
#' }
plot_pca_coordinates <- function(
    df_pca,
    x = "PC1",
    y = "PC2",
    base_size = 20,
    size = 2,
    title = "",
    color = NULL
){
  
    ggplot2::ggplot(
        df_pca,
        aes_string(x = x, y = y, colour = color)
    ) + 
        ggplot2::geom_point(size = size) +
        ggplot2::labs(title = title) + 
        ggplot2::theme_bw(base_size = base_size)
    
}


#' Plot projection given the PCA coordinates
#'
#' @param seed_nb An integer.
#' @param merged_col_data A dataframe. Contains the name of the samples from
#'     the cohorts
#' @param datasets_normalized A list of Summarized experiment. A summarized experiment 
#'     with normalized data from all cohorts
#' @param which_cohorts_trainings A vector of strings. Specify which
#'     cohorts to use when training the PCA.
#' @return A list with the new PCA fit and the embeddings from
#'     all samples.
get_new_pca <- function(
    seed_nb, 
    merged_col_data, 
    datasets_normalized,
    which_cohorts_trainings
){
 
    print(paste("Seed:", seed_nb))
    set.seed(seed_nb)
    samples_for_training <- merged_col_data %>% 
        data.frame %>%
        dplyr::filter(cohort %in% which_cohorts_training) %>%
        dplyr::pull(sample_name) %>%
        sample(., size = 1000)
    
    training_set <- lapply(
        datasets_normalized[which_cohorts_training], 
        function(sum_exp, i, genes_for_pca) 
            assay(sum_exp[genes_for_pca, ], i = i) %>% 
            data.frame(check.names = FALSE), 
        i = "avg_ranking",
        genes_for_pca = genes_for_pca
    ) %>% 
        dplyr::bind_cols() %>%
        .[, samples_for_training]

    pca_fit <- PCAtools::pca(
        training_set,
        metadata = dplyr::bind_rows(
            lapply(
                datasets_normalized[which_cohorts_training],
                function(df){
                    colData(df) %>% data.frame %>%
                        dplyr::filter(sample_name %in% colnames(training_set))
                }
            ),
            .id = "cohort"
        ) %>% .[colnames(training_set), ]
    )

    datasets_pca_coordinates <- lapply(
        datasets_normalized,
        get_pca_coordinates,
        pca_fit = pca_fit
    )
    
    list(
        train = pca_fit,
        test = datasets_pca_coordinates
    )
}


#' Calculate fuzziness score
#'
#' @param genes_to_remove A character vector. Genes that are not available
#'     in the sample
#' @param pca_fit A PCAtools object. Output of pca from PCAtools.
#' @return A vector with the fuzziness score for each component.
get_fuzziness_score <- function(
    genes_to_remove,
    pca_fit,
    which_pcs = 2:3
){
    
    
    mapply(
        function(genes, pc, pca_fit){
            sum(abs(pca_fit$loadings[genes, pc]))
        },
        genes = genes_to_remove %>% as.data.frame,
        pc = names(pca_fit$loadings)[which_pcs],
        MoreArgs = list(pca_fit = pca_fit)
    )
    
}
