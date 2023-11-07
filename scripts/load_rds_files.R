db_files <- list(
    surv_analysis_estrogen = c(),
    pca_merging = c(
        datasets = "datasets_with_scores",
        stable_genes = "stable_genes",
        which_exp = "which_exp",
        datasets_normalized = "datasets_normalized",
        correlation_scores = "correlation_scores",
        genes_for_pca = "genes_for_pca",
        pca_fit = "pca_fit",
        datasets_pca_coordinates = "datasets_pca_coordinates",
        merged_col_data = "merged_col_data",
        which_cohorts_training = "which_cohorts_training",
        df_pca_coordinates = "df_pca_coordinates",
        new_pcas = "new_pcas",
        pca_fit_no_norm = "pca_fit_no_norm",
        jaccard_indices = "jaccard_indices"
    ),
    validation = c(
        df_pca_coordinates = "df_pca_coordinates",
        pca_fit = "pca_fit",
        datasets_normalized = "datasets_normalized",
        datasets = "datasets_with_scores",
        stable_genes = "stable_genes",
        genes_for_pca = "genes_for_pca",
        samples_to_use = "samples_to_use",
        pca_random_genes_patients = "pca_random_genes_patients",
        fuzziness_scores_random = "fuzziness_scores_random",
        smc_df_pca = "smc_df_pca",
        smc_normalized = "smc_normalized",
        pdx_df_pca = "pdx_df_pca",
        pdx_normalized = "pdx_normalized",
        normal_swiss_df_pca = "normal_swiss_df_pca",
        normal_swiss_normalized = "normal_swiss_normalized",
        pdx_tpm = "pdx_tpm",
        gene_sets = "gene_sets",
        plot_estimates_all = "plot_estimates_all",
        nb_samples_avg = "nb_samples_avg"
    ),
    scoring = c(
        df_pca_coordinates = "poetic_df_pca_with_scores",
        df_pca = "df_pca_coordinates",
        pca_fit = "pca_fit",
        datasets = "datasets_with_poetic",
        stable_genes = "stable_genes",
        poetic_df_pca = "poetic_df_pca_with_scores",
        poetic_normalized = "poetic_normalized",
        datasets_normalized = "datasets_normalized",
        pipeline_scores_plots = "pipeline_scores_plots",
        pipeline_scores_plots_er_only = "pipeline_scores_plots_er_only",
        pipeline_scores_plots_er_androgen = "pipeline_scores_plots_er_androgen",
        comparison_to_average_neighbors = "comparison_to_average_neighbors"
    ),
    trying = c(
        datasets = "datasets_with_poetic",
        stable_genes = "stable_genes",
        which_exp = "which_exp",
        which_cohorts_training = "which_cohorts_training",
        merged_col_data = "merged_col_data_trying",
        gene_sets = "gene_sets",
        genes_for_pca = "genes_for_pca",
        df_pca = "poetic_df_pca_with_scores",
        pca_fit = "pca_fit",
        datasets_normalized_og = "datasets_normalized",
        datasets_normalized = "datasets_normalized_all_genes",
        pca_fit_all_genes = "pca_fit_all_genes",
        pca_fit_all_genes_wo_outliers = "pca_fit_all_genes_wo_outliers",
        df_pca_coordinates = "df_pca_coordinates_all_genes",
        df_pca_coordinates_wo_outliers = "df_pca_coordinates_wo_outliers",
        df_pca_wo_pcs = "df_pca_wo_pcs",
        pca_removed_pcs = "pca_removed_pcs",
        df_pcs_regressed = "df_pcs_regressed",
        gsva_scores_embeddings_tcga = "gsva_scores_embeddings_tcga",
        gsva_scores_embeddings_scanb = "gsva_scores_embeddings_scanb",
        gsva_scores_embeddings_scanb_within_tcga = "gsva_scores_embeddings_scanb_within_tcga",
        gsva_scores_embeddings_scanb_within_tcga_scanb = "gsva_scores_embeddings_scanb_within_tcga_scanb",
        df_pca_scores = "df_pca_scores",
        pipeline_scores_plots_new_scores = "pipeline_scores_plots_new_scores",
        singscore_dfs = "singscore_dfs",
        gene_sets_gsea = "gene_sets_gsea",
        singscore_dfs_regressed = "singscore_dfs_regressed",
        gsva_one_pathway = "gsva_one_pathway",
        gsva_scores = "gsva_scores_trying",
        gsva_scores_embeddings_scanb_within_tcga_regressed = "gsva_scores_embeddings_scanb_within_tcga_regressed",
        datasets_nanostring = "datasets_with_scores_nanostring",
        scores_nanostring = "scores_nanostring_only",
        pdx_tpm = "pdx_tpm",
        pdx_regressed = "pdx_regressed",
        scores = "scores_regressed_only",
        ssgsea_dfs = "ssgsea_dfs"
    ),
    risk_score = c(
        datasets = "datasets_with_scores",
        stable_genes = "stable_genes",
        which_exp = "which_exp",
        which_cohorts_training = "which_cohorts_training",
        merged_col_data = "merged_col_data",
        gene_sets = "gene_sets",
        genes_for_pca = "genes_for_pca",
        df_pca = "poetic_df_pca_with_scores",
        pca_fit = "pca_fit",
        datasets_normalized_og = "datasets_normalized",
        poetic_normalized = "poetic_normalized",
        risk_score_model = "risk_score_model",
        coefs_risk_model = "coefs_risk_model",
        coefs_risk_model_pcs = "coefs_risk_model_pcs"
    )
)

# find where the files are in the results folder
path_to_files <- sapply(
    paste0("^", db_files[[name_document]], "\\.rds$"),
    list.files,
    path = "../results",
    full.names = TRUE,
    recursive = TRUE, 
    USE.NAMES = FALSE
)
#path_to_files <- paste0("../", path_to_files)
names(path_to_files) <- names(db_files[[name_document]])

# check if they have the same length, if not stop. it means there are more
# files than it should
if (length(path_to_files) != length(db_files[[name_document]])){
    stop(paste0(
        "The number of files are not matching. ",
        "There might be files with the same name"
    ))
}

# exclude entries that have character 0 as path
path_to_files_ <- path_to_files[path_to_files != "../character(0)"]

# it is not possible to use parallel, as the workers don't have 
# access to R main's environment. for this we need to load first using
# mcmapply and then use the function list2env to load to the global 
# environment. after that we can remove the list
nb_cores <- parallel::detectCores() - 1
files_loaded <- parallel::mclapply(
    path_to_files,
    function(x){
        readRDS(x)  
    },
    mc.cores = nb_cores
)
list2env(files_loaded, envir = .GlobalEnv)
rm(files_loaded)