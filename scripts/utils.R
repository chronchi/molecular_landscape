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
    cat(cohort, type_survival, name_signature, "\n")
    
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
        aes_string(x = x, y = y, color = color)
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
        ) %>% .[colnames(training_set), ],
        center = FALSE, 
        scale = FALSE
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
    which_pcs = 3:4
){
    
    
    mapply(
        function(genes, pc, pca_fit){
            sum(abs(pca_fit$loadings[genes, pc]))
        },
        genes = genes_to_remove %>% as.data.frame,
        pc = names(pca_fit$loadings)[which_pcs],
        MoreArgs = list(pca_fit = pca_fit)
    ) %>% `names<-`(paste0("PC", which_pcs))
    
}

#' Plot PCA of samples with random missing genes
#'
#' @param pca_random_genes_patients A list. PCA coordinates from 
#'     patients with missing genes
#' @param patient A string. Name of the patient
#' @param df_pca_coordinates A dataframe. Contains the principal
#'     components from TCGA, METABRIC and SCANB
#' @return A plot showing where points would be if genes were
#'     missing.
plot_pca_random_genes <- function(
    pca_random_genes_patients,
    patient,
    df_pca_coordinates
){
    df <- pca_random_genes_patients[[patient]]
    
    df %>%
    ggplot2::ggplot(
        aes(
            x = PC3, 
            y = PC4, 
            color = proportion, 
            shape = embedding, 
            size = embedding,
            alpha = embedding
        )
    ) +
    ggplot2::geom_line(
        data = df %>% dplyr::filter(embedding == "original") %>%
            dplyr::bind_rows(., data.frame("PC3" = c(0), "PC4" = c(0))),
        aes(x = PC3, y = PC4), 
        inherit.aes = FALSE,
        alpha = 0.4,
        linetype = "dashed"
    ) +
    ggplot2::geom_point(
        data = df_pca_coordinates %>% dplyr::filter(cohort == "tcga"), 
        aes(x = PC3, y = PC4),
        alpha = 0.2, 
        color = "gray",
        inherit.aes = FALSE
    ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_viridis_d() +
    ggplot2::scale_size_manual(
        values = c("original" = 5, "random" = 2)
    ) +
    ggplot2::scale_alpha_manual(values = c("original" = 1, "random" = 0.5)) +
    ggplot2::labs(
        title = patient
    ) +
    ggplot2::theme_bw()
}

#' Plot PCA embedding of the dataframe using new samples
#'
#' @param df_pca A dataframe. PCA coordinates from 
#'     patients, including the SCANB, METABRIC and TCGA
#'     cohorts.
#' @param color A string. Which column to use when coloring
#'     the dots.
#' @param name_cohort A string. Name of the cohort that 
#'     is being used in the embedding
#' @param x A string. X-axis PC component
#' @param y A string. Y-axis PC component
#' @param title A string. Title for the plot
#' @return A plot of the embedding for new samples with the 
#'     TCGA, SCANB and METABRIC in the background.
get_plot_new_samples <- function(
    df_pca, 
    color, 
    name_cohort = "smc",
    x = "PC2",
    y = "PC3",
    title = "SMC samples overlayed on the molecular landscape"
){
    
    df_pca %>% 
        dplyr::filter(pam50 %in% 
            c("basal", "her2", "lumb", "luma", "normal", "not_available")
        ) %>%
        ggplot2::ggplot(
            aes_string(x = x, y = y, color = color)
        ) +
        ggplot2::geom_point(aes(alpha = cohort, size = cohort)) +
        ggplot2::scale_size_manual(
            values = c(3, 1, 1, 1) %>% 
                `names<-`(c(
                    name_cohort, 
                    "tcga", 
                    "scanb", 
                    "metabric"
                ))
        ) +
        ggplot2::scale_alpha_manual(
            values = c(1, .1, .1, .1) %>% 
                `names<-`(c(
                    name_cohort, 
                    "tcga", 
                    "scanb", 
                    "metabric"
                ))
        ) +
        ggplot2::labs(
            alpha = "Cohort",
            size = "Cohort",
            title = title,
            subtitle = "All samples from TCGA, METABRIC and SCANB are plotted"
        ) + 
        ggplot2::theme_bw(base_size = 20)
}

#' Plot PCA embedding of the 3 big cohorts to serve as a base plot
#'
#' @param df_pca A dataframe. PCA coordinates from 
#'     patients, including the SCANB, METABRIC and TCGA
#'     cohorts.
#' @param x A string. PC for x-axis
#' @param y A string. PC for y-axis
#' @param color A string. Which column to use when coloring
#'     the dots.
#' @param size_dots An integer. Size of the points in the plot
#' @param alpha_val A number. Alpha value of the points
#' @param size_legend An integer. Size of the color legend
#' @param base_size An integer. Parameter for ggplot2::theme_bw
#' @return A plot of the TCGA, SCANB and METABRIC embedding
get_base_plot <- function(
    df_pca,
    x = "PC3",
    y = "PC4",
    color = "pam50",
    size_dots = 2,
    alpha_val = 0.1,
    size_legend = 4,
    base_size = 10
){
    
    df_pca %>%
        ggplot2::ggplot(aes_string(x = x, y = y, color = color)) + 
        ggplot2::geom_point(size = size_dots, alpha = alpha_val) +
        ggplot2::scale_alpha(guide = 'none') +
        ggplot2::labs(
            color = "PAM50"
        ) + 
        ggplot2::guides(
            colour = ggplot2::guide_legend(
                override.aes = list(size = size_legend, alpha = 1)
            )
        ) + 
        ggplot2::theme_bw(base_size = base_size)
}


#' Get list of samples that are in the radius of another sample
#'
#' @param df_pca A dataframe. PCA coordinates from 
#'     patients, including the SCANB, METABRIC and TCGA
#'     cohorts.
#' @param components_sample A vector containing the PC values from the 
#'     center sample.
#' @param sample_name A string. Name of the sample, it is used to remove it
#'     from the list of samples in the neighborhood.
#' @param radius A number. Size of the ball
#' @param which_components A character vector. Name of the components used to
#'     to calculate the distance. Length should match the components_sample.
#' @param column_patients A string. name of the column from df_pca where the
#'     patient_names are stored.
#' @return A list containing two vectors, one with the names from the 
#'     samples and another with the distance to the center
get_samples_neighborhood <- function(
    components_sample, 
    sample_name,
    radius,
    df_pca,
    which_components = c("PC3", "PC4"),
    column_patients = "sample_name"
){
    
    x <- df_pca[, which_components[1]]
    y <- df_pca[, which_components[2]]
    
    distance_to_center <- (
        (components_sample[1] - x)^2 + 
            (components_sample[2] - y)^2
    )
    
    samples_to_select <- distance_to_center < radius^2
    
    samples_in_neighborhood <- df_pca[samples_to_select, ] %>% 
        dplyr::pull(!!sym(column_patients)) %>%
        setdiff(., sample_name)
    
    list(
        samples = samples_in_neighborhood,
        distances = sqrt(distance_to_center[samples_to_select])
    )
    
}

#' Get formatted dataframe to calculate the average scores
#'
#' @param df_pca A dataframe. PCA coordinates from 
#'     patients, including the SCANB, METABRIC and TCGA
#'     cohorts.
#' @param samples_to_use A character vector. Vector with all samples in
#'     the neighborhood of a specific sample
#' @param scores_to_use A character vector. List with scores to be 
#'     returned.
#' @param column_patients A string. name of the column from df_pca where the
#'     patient_names are stored.
#' @return A long dataframe containing the scores and pathways along with
#'     with the patient names.
get_scores_for_average <- function(
    df_pca, 
    samples_to_use, 
    scores_to_use,
    column_patients = "sample_name"
){
    df_pca %>% 
        dplyr::filter(!!sym(column_patients) %in% samples_to_use) %>%
        dplyr::select(all_of(c(scores_to_use, column_patients))) %>%
        tidyr::pivot_longer(
            cols = all_of(scores_to_use),
            names_to = "pathway",
            values_to = "score"
        )
}

#' Get fitted models for all pathways selected in get_scores_for_average
#'
#' @param scores_for_average A dataframe. Output from get_scores_for_average
#' @return A list with the fitted models from rstanarm for each pathway
#'     individually
get_average_neighboorhood <- function(scores_for_average){
    
    sapply(
        unique(scores_for_average$pathway),
        function(pathway_name, scores_for_average){

            
            rstanarm::stan_glm(
                data = scores_for_average %>% 
                    dplyr::filter(pathway == pathway_name),
                formula = "score ~ 1", 
                chains = 4,
                prior_intercept = rstanarm::normal(
                    location = 0, 
                    scale = 1
                ),
                refresh = 0
            )   
            
        },
        scores_for_average = scores_for_average,
        USE.NAMES = TRUE,
        simplify = FALSE
    )
}



plots_estimates <- function(
    models, 
    sample_name, 
    scores_to_use, 
    scores_patient,
    base_size
){
    
    tidy_draws <- lapply(
        models, 
        function(x) x %>% tidybayes::spread_draws(`(Intercept)`)
    ) %>% dplyr::bind_rows(
        .id = "pathway"
    ) 
    
    tidy_draws$pathway <- scores_to_use[tidy_draws$pathway] %>% as.character
    
    tidy_draws %>% 
    ggplot2::ggplot(
        aes(x = `(Intercept)`)
    ) + 
        ggplot2::geom_vline(
            data = scores_patient,
            mapping = aes(xintercept = score),
            linetype = "dashed",
            color = "red"
        ) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
        tidybayes::stat_dotsinterval(
            quantiles = 100,
            size = 3
        ) + 
        ggplot2::facet_wrap(~pathway, ncol = 3) +
        ggplot2::labs(
            title = paste0(
                "Distribution of average scores\n",
                "in ", sample_name, " neighborhood"
            ),
            y = "",
            x = "Average score",
            caption = paste0(
                "Red dashed line is patient score\n",
                "Black dashed line is centered at 0"
            )
        ) + 
        ggplot2::theme_bw(base_size = base_size)
}


plots_estimates_without_patient <- function(
    models, 
    sample_name, 
    scores_to_use, 
    scores_patient,
    base_size,
    size_dots = 5,
    which_direction = "top_down",
    color_title = c("top_down" = "red", "left_right" = "black")
){
    
    # specify names for the title depending on the direction
    proper_path_name <- c(
        "top_down" = "Path 1 (right)",
        "left_right" = "Path 2 (left)"
    )
    
    # get draws from the intercept as they are considered the
    # averages
    tidy_draws <- lapply(
        models, 
        function(x) x %>% tidybayes::spread_draws(`(Intercept)`)
    ) %>% dplyr::bind_rows(
        .id = "pathway"
    ) 
    
    # change to better names when plotting
    tidy_draws$pathway <- scores_to_use[tidy_draws$pathway] %>% as.character
    
    tidy_draws %>% 
    ggplot2::ggplot(
        aes(x = `(Intercept)`)
    ) + 
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
        tidybayes::stat_dotsinterval(
            quantiles = 100,
            fill = color_title[which_direction],
            alpha = 0.5,
            show.legend = FALSE,
            size = size_dots
        ) + 
        ggplot2::facet_wrap(~pathway, ncol = 3) +
        ggplot2::labs(
            title = paste0(
                proper_path_name[which_direction], 
                "\nDistribution of average scores"
            ),
            y = "",
            x = "Average score",
            caption = ifelse(
                which_direction == "top_down",
                "",
                "Black dashed line is centered at 0"
            )
        ) + 
        ggplot2::theme_bw(base_size = base_size) +
        ggplot2::coord_cartesian(xlim = c(-0.6, 0.6)) +
        ggplot2::theme(
            plot.title = element_text(
                colour = color_title[which_direction]
            )
        )
}

get_patient_scores_distributions <- function(
    patient_name,
    df_pca,
    scores_to_use,
    radius = 0.3,
    components = c("PC3", "PC4"),
    column_patients = "sample_name",
    which_direction = NULL
){

    components_sample <- df_pca %>% 
        dplyr::filter(sample_name == patient_name) %>%
        dplyr::select(all_of(components)) %>%
        dplyr::slice(1) %>%
        unlist

    samples_and_distance_neighborhood <- get_samples_neighborhood(
        components_sample = components_sample, 
        sample_name = patient_name, 
        radius = radius, 
        df_pca = df_pca %>% 
            dplyr::filter(cohort %in% c("tcga", "metabric", "scanb")), 
        which_components = components, 
        column_patients = column_patients
    )


    scores_avg <- get_scores_for_average(
        df_pca = df_pca %>%
            dplyr::filter(cohort %in% c("tcga", "metabric", "scanb")), 
        samples_to_use = samples_and_distance_neighborhood$samples,
        scores_to_use = names(scores_to_use),
        column_patients = column_patients
    )
    
    
    average_scores_neighborhood <- get_average_neighboorhood(scores_avg)

    scores_patient <- df_pca %>% 
        dplyr::filter(sample_name == patient_name) %>%
        dplyr::select(all_of(names(scores_to_use))) %>%
        tidyr::pivot_longer(
            cols = all_of(names(scores_to_use)),
            names_to = "pathway",
            values_to = "score"
        ) %>% 
        dplyr::mutate(
            pathway = scores_to_use[pathway]
        )

    list(
        average_scores_neighborhood = average_scores_neighborhood,
        samples_distance = samples_and_distance_neighborhood,
        scores_patient = scores_patient
    )
}

get_patient_scores_distributions_all <- function(
    patient_name,
    df_pca,
    scores_to_use,
    radius = 1,
    components = c("PC3", "PC4"),
    column_patients = "sample_name",
    which_direction = NULL,
    base_size = 15
){

    components_sample <- df_pca %>% 
        dplyr::filter(sample_name == patient_name) %>%
        dplyr::select(all_of(components)) %>%
        dplyr::slice(1) %>%
        unlist

    samples_and_distance_neighborhood <- get_samples_neighborhood(
        components_sample = components_sample, 
        sample_name = patient_name, 
        radius = radius, 
        df_pca = df_pca %>% 
            dplyr::filter(cohort %in% c("tcga", "metabric", "scanb")), 
        which_components = components, 
        column_patients = column_patients
    )


    scores_avg <- get_scores_for_average(
        df_pca = df_pca %>%
            dplyr::filter(cohort %in% c("tcga", "metabric", "scanb")), 
        samples_to_use = samples_and_distance_neighborhood$samples,
        scores_to_use = names(scores_to_use),
        column_patients = column_patients
    )
    
    
    average_scores_neighborhood <- get_average_neighboorhood(scores_avg)

    scores_patient <- df_pca %>% 
        dplyr::filter(sample_name == patient_name) %>%
        dplyr::select(all_of(names(scores_to_use))) %>%
        tidyr::pivot_longer(
            cols = all_of(names(scores_to_use)),
            names_to = "pathway",
            values_to = "score"
        ) %>% 
        dplyr::mutate(
            pathway = scores_to_use[pathway]
        )

    p <- plots_estimates(
        average_scores_neighborhood,
        patient_name, 
        scores_to_use, 
        scores_patient,
        base_size
    )
    
    list(
        average_scores_neighborhood = average_scores_neighborhood,
        samples_distance = samples_and_distance_neighborhood,
        scores_patient = scores_patient,
        plot = p
    )
}

get_plot_patient_distribution <- function(
    avg_rstan_samples,
    scores_to_use,
    which_direction,
    radius = 1,
    base_size = 10,
    size_dots = 5
){

    average_scores_neighborhood <- avg_rstan_samples$average_scores_neighborhood
    scores_patient <- avg_rstan_samples$scores_patient
    
    if (is.null(which_direction)){
        plots_estimates(
            average_scores_neighborhood,
            scores_to_use, 
            base_size = base_size
        )   
    } else {
        plots_estimates_without_patient(
            average_scores_neighborhood, 
            scores_to_use = scores_to_use, 
            scores_patient = scores_patient, 
            base_size = base_size,
            which_direction = which_direction
        )
    }

}

get_closest_sample <- function(
    pca_values,
    df_pca,
    sample_name = "sample_name"
){
    
    df_pca %>% 
        dplyr::rowwise() %>%
        dplyr::mutate(
            distance = norm(as.matrix(c(PC3, PC4) - pca_values))
        ) %>%
        dplyr::ungroup() %>%
        dplyr::slice_min(order_by = distance, n = 1) %>%
        dplyr::pull(!!sym(sample_name))
     
}

plot_selected_samples <- function(
    df_pca, 
    pca_values, 
    scores_plots_movie,
    radius = 0.3,
    title_plot = "",
    base_size = 20,
    size_points = 3,
    size_line = 3
){


    lapply(
        1:length(scores_plots_movie$top_down),
        function(i, pca_values, scores_plots_movie, df_pca, ...){
            
            
            samples_to_highlight <- lapply(
                scores_plots_movie,
                function(x, i = i) x[[i]]$samples_distance$samples,
                i = i
            ) %>% unlist %>% unique
            
            ggplot2::ggplot(
                df_pca,
                aes(
                    x = PC3, 
                    y = PC4, 
                    color = pam50, 
                    shape = cohort, 
                    alpha = ifelse(
                        sample_name %in% samples_to_highlight, 
                        1, 
                        0.2
                    )
                )
            ) +
                ggplot2::geom_path(
                    inherit.aes = FALSE,
                    data = pca_values[["top_down"]] %>% data.frame,
                    aes(x = PC3, y = PC4),
                    color = "red", 
                    show.legend = FALSE, 
                    size = size_line
                ) + 
                ggplot2::geom_path(
                    inherit.aes = FALSE,
                    data = pca_values[["left_right"]] %>% data.frame,
                    aes(x = PC3, y = PC4),
                    color = "black",
                    show.legend = FALSE,
                    size = size_line
                ) + 
                ggplot2::geom_point(size = size_points) +
                ggforce::geom_circle(
                    data = data.frame(
                        x_center = pca_values[["top_down"]][i, 1], 
                        y_center = pca_values[["top_down"]][i, 2], 
                        radius = radius
                    ),
                    mapping = aes(
                        x0 = x_center, 
                        y0 = y_center, 
                        r = radius*1.1
                    ),
                    alpha = 0.7,
                    inherit.aes = FALSE,
                    color = "red"
                ) + 
                ggforce::geom_circle(
                    data = data.frame(
                        x_center = pca_values[["left_right"]][i, 1], 
                        y_center = pca_values[["left_right"]][i, 2], 
                        radius = radius
                    ),
                    mapping = aes(
                        x0 = x_center, 
                        y0 = y_center, 
                        r = radius*1.1
                    ),
                    alpha = 0.7,
                    inherit.aes = FALSE,
                    color = "black"
                ) + 
                ggplot2::scale_alpha(guide = 'none') +
                ggplot2::labs(title = title_plot) +
                ggplot2::theme_bw(base_size = base_size)
        },
        pca_values = pca_values, 
        scores_plots_movie = scores_plots_movie, 
        df_pca = df_pca,
        title_plot = title_plot,
        base_size = base_size,
        size_points = size_points,
        size_line = size_line
    )
}

