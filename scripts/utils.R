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
    
    filename <- paste0(
        path_to_save, "/",
        cohort, "_", type_survival,
        "_", name_signature, ".pdf"
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
    
    pdf(file = filename, width = width, height = height)
    plot(forest_plot)
    dev.off()
    
    forest_plot
}