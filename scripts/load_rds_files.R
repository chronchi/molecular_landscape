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
        datasets_pca_coordinates = "datasets_pca_coordinates"
    )
)

# find where the files are in the results folder
path_to_files <- sapply(
    paste0(db_files[[name_document]], "\\.rds$"),
    list.files,
    path = "..",
    recursive = TRUE, 
    USE.NAMES = FALSE
) 
path_to_files <- paste0("../", path_to_files)
names(path_to_files) <- names(db_files[[name_document]])

# check if they have the same length, if not stop. it means there are more
# files than it should
if (length(path_to_files) != length(db_files[[name_document]])){
    stop(paste0(
        "The number of files are not matching. ",
        "There might be files with the same name"
    ))
}

# it is not possible to use parallel, as the workers don't have 
# access to R main's environemnt. for this we need to load first using
# mcmapply and then use the function list2env to load to the global 
# environment. after that we can remove the list
nb_cores <- parallel::detectCores() - 1
files_loaded <- parallel::mclapply(
    path_to_files,
    readRDS,
    mc.cores = nb_cores
)
list2env(files_loaded, envir = .GlobalEnv)
rm(files_loaded)