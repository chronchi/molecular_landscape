# Author: Carlos Ronchi
# Date: 2021.07.22
# Email: carlos.ronchi@epfl.ch
# Description: A set of helper functions that are reused across the different
#   documents in order to decluter them.

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


plot_umap <- function(
    dataset, umap_projection, name_dataset, color_by,
    base_size = 10,
    point_size = 1,
    alpha_val = 0.8
){
    
    umap_embedding <- umap_projection %>% data.frame %>%
        `colnames<-`(c("umap1", "umap2")) %>%
        tibble::rownames_to_column(var = "sample_name") %>%
        dplyr::inner_join(
            .,
            colData(dataset) %>% data.frame,
            by = "sample_name"
        )
    
    ggplot2::ggplot(
        umap_embedding,
        aes_string(x = "umap1", y = "umap2", color = color_by)
    ) + 
        ggplot2::geom_point(size = point_size, alpha = alpha_val) + 
        ggplot2::labs(title = toupper(name_dataset)) +
        ggplot2::theme_bw(base_size = base_size)
        
}

forest_plot_fits <- function(
    fit, cohort, names_coefficients, patients, type_survival, ...
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
    
    output_df %>% forestplot::forestplot(
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
}