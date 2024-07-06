# Import environment variables as global variables
rds_dir <- Sys.getenv("RDS_DIR")
image_dir <- Sys.getenv("IMAGE_DIR")
image_ext <- Sys.getenv("IMAGE_EXT")

####### SAVE LIST OF ITEMS #######

save_plots <- function(plots_list, folder_path, image_extension = ".png") {
  # Ensure the folder exists, if not, create it
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }
  
  # Helper function to flatten the list and collect plots and data frames
  # Images are also lists so simply unflattening the list does not work
  flatten_list <- function(lst) {
    flat_list <- list()
    
    for (name in names(lst)) {
      if (inherits(lst[[name]], "ggplot") || inherits(lst[[name]], "trellis")) {
        flat_list[[name]] <- lst[[name]]
      } else if (is.data.frame(lst[[name]])) {
        flat_list[[name]] <- lst[[name]]
      } else if (is.list(lst[[name]])) {
        flat_list <- c(flat_list, flatten_list(lst[[name]]))
      }
    }
    
    return(flat_list)
  }
  
  # Flatten the plots list
  flat_items <- flatten_list(plots_list)
  
  # Iterate over the flattened list of items and save each plot or data frame
  for (item_name in names(flat_items)) {
    # Construct the file path
    if (inherits(flat_items[[item_name]], "ggplot") || inherits(flat_items[[item_name]], "trellis")) {
      file_path <- file.path(folder_path, paste0(item_name, image_extension))
      # Save the plot to the specified file path
      ggsave(filename = file_path, plot = flat_items[[item_name]])
    } else if (is.data.frame(flat_items[[item_name]])) {
      file_path <- file.path(folder_path, paste0(item_name, ".csv"))
      # Save the data frame to the specified file path
      write.csv(flat_items[[item_name]], file_path, row.names = FALSE)
    }
  }
}

####### GET PATIENT NUM #######

get_patient_num <- function(patient_data) {
  # The double [[ are necessary to make it return an unnamed number
  patient_num <- patient_data$Patient.ID[[1]]
  return(patient_num)
}

####### GET PATIENT DIRECTORIES #######

get_patient_dir_img <- function(patient_num) {
  patient_dir_img <- here(image_dir, paste0("Patient_", patient_num, "_plots/"))
  return(patient_dir_img)
}

get_patient_dir_rds <- function(patient_num) {
  patient_dir_rds <- here(rds_dir, paste0("Patient_", patient_num, "_data/"))
  return(patient_dir_rds)
}

####### GENERATE HEATMAP #######

# Heatmap plot, define a function to dynamically calculate label size based on the number of features
calculate_label_size <- function(num_features) {
  base_size <- 400
  scaled_size <- base_size / (num_features)
  return(scaled_size)
}

# Define the size parameter (size of the cluster names on top) of the heatmap
# proportionally to the size of the smallest cluster
calculate_clusters_names_size <- function(smallest_cluster_proportion) {
  base_size <- 9  # Base size for the heatmap
  # scaling_factor <- 0.1
  return(base_size * exp(log(smallest_cluster_proportion)/3))
}

generate_dyn_text_heatmap <- function(
    patient_data,
    cluster_var,
    assay_name,
    cluster_name = NULL,
    color_lookup_table = NULL) {
  
  # If a human friendly name is not given, use the name of the column in the Seurat object
  if (is.null(cluster_name)) {
    cluster_name <- cluster_var
  }
  
  # Use a better palette by default
  if (is.null(color_lookup_table)) {
    color_lookup_table <- generate_colors_lookup_table(patient_data, cluster_var)
  }
  
  # Select the cluster as the identity
  Idents(patient_data) <- cluster_var
  
  most_significant_markers <- find_most_significant_markers(
    patient_data,
    cluster_var,
    assay_name)
  
  # Size of smallest cluster
  smallest_cluster_proportion <- min(get_clusters_proportions(patient_data))
  # Size of cluster label
  cluster_label_size <- calculate_clusters_names_size(smallest_cluster_proportion)
  
  # Calculate the number of features (rows)
  num_features <- length(most_significant_markers$gene)
  # Calculate the label size for the current figure
  features_label_size <- calculate_label_size(num_features)
  
  # Read patient number
  patient_num <- get_patient_num(patient_data)
  
  # Create the heatmap with the scaled text
  diff_expr_genes_heatmap <- DoHeatmap(
    patient_data,
    features = most_significant_markers$gene,
    assay = assay_name,
    label = TRUE,
    size = 2,
    group.colors = color_lookup_table,
  ) + theme(
    axis.text.y = element_text(size = features_label_size),
  ) + labs(
    title = paste("Patient", patient_num, cluster_name),
    subtitle = "Top 10 Differentially Expressed Genes per Cluster"
  )
  
  result_list <- list(
    diff_expr_genes_heatmap,
    most_significant_markers
  )
  
  # Dynamically assign the desired names
  names(result_list) <- c(
    paste("Patient", patient_num, cluster_var, "diff_expr_genes_heatmap", sep = "_"),
    paste("Patient", patient_num, cluster_var, "diff_expr_genes_heatmap_list", sep = "_")
  )
  
  # Return the named list
  return(result_list)
}

find_most_significant_markers <- function(
    patient_data,
    cluster_var,
    assay_name,
    log2fc_threshold = 1,
    genes_per_cluster = 10) {
  
  # Select the cluster as the identity
  Idents(patient_data) <- cluster_var
  
  # Find markers (differentially expressed genes) for each of the identity classes in the filtered dataset
  # If you have more than one assay it is necessary to specify the assay parameter
  markers_data <- FindAllMarkers(
    patient_data,
    assay = assay_name,
    only.pos = TRUE,
    random.seed = 5)
  
  # Filter markers to get the most significant ones per cluster
  most_significant_markers <- markers_data %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > log2fc_threshold) %>%
    slice_head(n = genes_per_cluster) %>%
    ungroup()
  
  return(most_significant_markers)
}


get_clusters_proportions <- function(patient_data) {
  
  # Calculate the size of each cluster
  cluster_sizes <- table(Idents(patient_data))
  total_cells <- sum(cluster_sizes)
  cluster_proportions <- cluster_sizes / total_cells
  
  return(cluster_proportions)
}

####### INSITUTYPE SEMISUPERVISED #######

generate_flightpath <- function(IST_object, cols, patient_num) {
  
  fp <- flightpath_plot(flightpath_result = NULL,
                        insitutype_result = IST_object,
                        col = cols[IST_object$clust],
                        showclusterconfidence = TRUE) +
    labs(title = paste("Patient", patient_num),
         subtitle = "InSituType Semisupervised Clustering")
  
  fp_name <- paste("Patient",  patient_num, "InSituType_semisup_clusters_flightpath_plot", sep = "_")
  return(setNames(list(fp), fp_name))
  }

run_IST_semisup_extract_data <- function(
    patient_data,
    assay) {
 
  # Get the data necessary to run InSituType
  patient_rna_data <- extract_patient_rna_data(patient_data, assay)
  patient_cohort <- get_patient_cohort(patient_data)
  patient_rna_counts <- get_seurat_layer_data(patient_data, assay, "counts")
  patient_avg_neg_probes <- get_avg_neg_probes(patient_data, assay)
  
  # Call the InSituType semisupervised function
  patient_semisup <- run_InSituType_semisupervised(
    patient_data = patient_rna_data,
    patient_cohort = patient_cohort,
    patient_rna_counts = patient_rna_counts,
    patient_avg_neg_probes = patient_avg_neg_probes)
  
  # Return the Seurat object with the added clusters
  return(patient_semisup)
}

run_InSituType_semisupervised <- function(
    patient_data,
    patient_cohort,
    patient_rna_counts,
    patient_avg_neg_probes,
    io_profiles = NULL) {
  
  if (is.null(io_profiles)) {
    # Get cell reference profile data from NanoString
    # Use this reference profile as it is the only one available from CosMx data, originally from:
    # https://raw.githubusercontent.com/Nanostring-Biostats/CosMx-Cell-Profiles/main/Human/IO/IO.profiles.csv
    # I have added a name to the first column
    io_profiles <- read_csv(
      file = here("Analysis", "metadata", "NanoString.CosMx.Human.IO.profiles.csv"),
      col_types = cols(
        `Gene` = col_character(),
        `B cell` = col_double(),
        `Dendritic cell` = col_double(),
        Endothelial = col_double(),
        Fibroblast = col_double(),
        Macrophage = col_double(),
        `Mast cell` = col_double(),
        Monocyte = col_double(),
        Neutrophil = col_double(),
        `NK cell` = col_double(),
        Plasma = col_double(),
        Plasmablast = col_double(),
        `Plasmacytoid dendritic cell` = col_double(),
        `T cell CD4` = col_double(),
        `T cell CD8` = col_double(),
        `T cell regulatory` = col_double()
      ))
    row_name = io_profiles$Gene
    io_profiles %<>% select(-Gene) %>% as.matrix
    rownames(io_profiles) = row_name
  }
  
  set.seed(6)
  clusts_num <- chooseClusterNumber(
    counts = patient_rna_counts,
    neg = patient_avg_neg_probes,
    fixed_profiles = io_profiles,
    n_clusts = 1:8)
  
  # Semi-supervised learning with insitutype and reference profiles
  # InSituType needs integers, if given floating point numbers it fails with misleading errors
  patient_semisup <- insitutype(
    x = patient_rna_counts,
    neg = patient_avg_neg_probes,
    cohort = patient_cohort,
    reference_profiles = io_profiles,
    
    # Enter your own per-cell background estimates here if you
    # have them; otherwise insitutype will use the negprobes to
    # estimate background for you.
    bg = NULL,
    # Group the cells the do not correspond to any type in the reference matrix
    n_clusts = clusts_num$best_clust_number,
    # reference_profiles = updatedprofiles$updated_profiles,
    # Update the reference profile based on the current data
    update_reference_profiles = FALSE,
    # choosing inadvisably low numbers to speed the vignette; using the defaults
    # in recommended.
    # This is the number of cells used in each phase, because of random sampling
    n_phase1 = 20,
    n_phase2 = 50,
    n_phase3 = 200,
    n_starts = 1,
    max_iters = 5
  )
  
  return(patient_semisup)
}

####### EXTRACT PATIENT DATA #######

extract_patient_data <- function(all_patients_data, patient_num) {
  
  # Extract patient data from global Seurat object
  patient_data <- subset(x = all_patients_data, subset = Patient.ID == patient_num)
  
  return(patient_data)
}

extract_patient_rna_data <- function(patient_data, assay, start_index = 1, end_index = 1000) {
  
  assay_to_subset <- patient_data[[assay]]
  
  # Check if the provided indices are within the valid range
  if (start_index < 1 || end_index > nrow(assay_to_subset) || start_index > end_index) {
    stop("Invalid range of indices provided.")
  }
  
  # Create Seurat object with only the specified range of RNA data
  features_to_include <- rownames(assay_to_subset)[start_index:end_index]
  rna_only_assay <- subset(x = assay_to_subset, features = features_to_include)
  patient_rna_only <- patient_data
  patient_rna_only[[assay]] <- NULL
  patient_rna_only[["RNA"]] <- rna_only_assay
  
  return(patient_rna_only)
}
  
get_avg_neg_probes <- function(patient_data, assay) {
  
  assay_to_subset <- patient_data[[assay]]
  
  # Extract the negative probes from the Seurat object
  patient_neg_probes <- GetAssayData(subset(
        assay_to_subset,
        features = row.names(assay_to_subset) %>%
  grep("Negative", ., value = TRUE))) %>%
  as.matrix() %>%
  t()
  # Calculate the average negative probes per cell
  patient_avg_neg_probes <- Matrix::rowMeans(patient_neg_probes)
  # Return a large numeric
  return(patient_avg_neg_probes)
}

move_counts_to_new_assay <- function(patient_data, assay_name, new_assay_name, pattern) {
  
  assay_to_subset <- patient_data[[assay_name]]
  
  features_to_extract <- row.names(assay_to_subset) %>%
    grep(pattern, ., value = TRUE)
  
  # Extract the data from the specified assay
  extracted_data <- subset(
    assay_to_subset,
    features = features_to_extract)
  
  # Change the key of the newly created assay
  extracted_data@key <- paste0(new_assay_name, "_")
  
  # Save the extracted data into a new assay
  patient_data[[new_assay_name]] <- extracted_data
  
  # Remove data from the original assay
  patient_data[[assay_name]] <- subset(
    assay_to_subset,
    features = setdiff(row.names(assay_to_subset), features_to_extract))
  
  return(patient_data)
}

get_seurat_layer_data <- function(patient_data, assay_name, layer_name) {
  
  # Extract the data from the specified assay and layer
  patient_assay_layer_data <- GetAssayData(
    patient_data,
    assay = assay_name,
    layer = layer_name) %>%
  as.matrix() %>%
  t()
  # Return a matrix of the data
  return(patient_assay_layer_data)
}

get_patient_cohort <- function(patient_data) {
  
  # Features to be used for the cohorting
  features <- c("Mean.PanCK", "Mean.CD45", "Mean.CD68")
  # Cohort of the patient
  patient_immunofluorescence <- patient_data@meta.data %>% select(all_of(features))
  # fastCohorting is stochastic, so set the seed for reproducibility
  set.seed(42);
  # "Gaussian_transform = TRUE" maps variables to gaussians in order to place dramatically different variables on the same scale
  patient_cohort <- fastCohorting(patient_immunofluorescence, gaussian_transform = TRUE, n_cohorts = 5)
  # check clusters and cohort numbers
  table(patient_cohort)
  
  return(patient_cohort)
}

####### PATIENT INFO #######

get_patient_info <- function(patient_data) {
  
  patient_num <- get_patient_num(patient_data)
  # Initialize list to store the cores, stamps, and fovs associated with the patient
  cores <- list()
  # Get cores for the current patient
  patient_cores <- sort(unique(patient_data@meta.data$core_serial))
  
  for (core in patient_cores) {
    stamps <- list()
    # Get stamps for the current core
    patient_core_stamps <- sort(unique(patient_data@meta.data$stamp[patient_data@meta.data$core_serial == core]))
    
    for (stamp in patient_core_stamps) {
      fovs <- table(patient_data@meta.data$fov
                    [patient_data@meta.data$core_serial == core & patient_data@meta.data$stamp == stamp])
      stamps[[stamp]] <- as.list(fovs)
    }
    cores[[core]] <- stamps
  }
  patient_info <- list(patient_num = patient_num, cores = cores)
  return(invisible(patient_info))
}

print_patient_info <- function(patient_data) {
    
  patient_info <- get_patient_info(patient_data)
  patient_num <- patient_info$patient_num
  # Print cores and fovs associated with the patient
  print(paste("FOVs and cell count associated with patient", patient_num))
  
  cores <- patient_info$cores
  for (core in names(cores)) {
    print(paste("----------", "CORE", core, "----------"))
    
    stamps <- cores[[core]]
    # Iterate over stamps using indices
    for (i in seq_along(stamps)) {
      print(paste("-----", "STAMP", i, "-----"))
      # Print FOVs from the stamps list
      fovs <- stamps[[i]]
      
      if (length(fovs) == 0) {
        cat("Missing data\n\n")
      } else {
        # Print FOVs in a table
        cat("FOV\tCell_count\n")
        for (fov_name in names(fovs)) {
          cat(fov_name, "\t", fovs[[fov_name]], "\n")
        }
        cat("\n")
      }
    }
  }
}

####### CLUSTER INFO #######

# Define function to extract cluster information
get_cluster_info <- function(seurat_obj, clustering_column) {
  
  # Check if clustering_column is valid
  if (!clustering_column %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", clustering_column, "not found in meta.data of Seurat object."))
  }
  
  seurat_obj@meta.data[[clustering_column]] <- as.character(seurat_obj@meta.data[[clustering_column]])
  
  # Get cluster information
  clusters <- unique(seurat_obj@meta.data[[clustering_column]])
  
  # Initialize list to store cluster names and cell counts
  cluster_info <- c()
  
  # Loop through clusters and get cell counts
  for (i in seq_along(clusters)) {
    cluster_name <- clusters[i]
    cell_count <- sum(seurat_obj@meta.data[[clustering_column]] == cluster_name)
    cluster_info[cluster_name] <- cell_count
  }
  
  # Sort cluster_info vector by cell_count in descending order
  cluster_info <- cluster_info[order(cluster_info, decreasing = TRUE)]
  
  return(cluster_info)
}

# Define function to create cluster summary dataframe with dynamic columns
create_cluster_summary <- function(patient_data, clustering_column) {
  
  # Check if clustering_column is valid
  if (!clustering_column %in% colnames(patient_data@meta.data)) {
    stop(paste("Column", clustering_column, "not found in meta.data of Seurat object."))
  }
  
  # Initialize a dataframe to store cluster summary
  cluster_summary <- data.frame(cluster_name = character(),
                                stringsAsFactors = FALSE)
  
  # Get unique clusters
  clusters <- as.vector(unique(patient_data[[clustering_column]][[clustering_column]]))
  
  # Combine two existing columns to create a new column
  dataa <- paste(patient_data$core_serial, patient_data$stamp, sep = "_")
  
  patient_data <- AddMetaData(patient_data, dataa, "combined_id")
  
  # Get unique combinations of core_serial and stamp in the current cluster
  unique_combinations <- sort(unique(dataa))
  
  
  # Loop through clusters
  for (cluster_name in clusters) {
    
    # Subset data for the current cluster
    cluster_data <- patient_data[ , patient_data[[clustering_column]] == cluster_name]
    
    # Get total cell count for the cluster
    cell_count <- nrow(cluster_data@meta.data)
    
    # Create a new row for the cluster in cluster_summary dataframe
    new_row <- data.frame(cluster_name = cluster_name,
                          total_cells = cell_count,
                          stringsAsFactors = FALSE)
    
    # Loop through unique combinations of core_serial and stamp
    for (core_stamp in unique_combinations) {
      
      # Count cells for the current core_stamp combination
      cell_count <- sum(cluster_data$combined_id == core_stamp)
      
      # Add a column for the current core_stamp combination
      new_row[[core_stamp]] <- cell_count
    }
    
    # Append new_row to cluster_summary dataframe
    cluster_summary <- rbind(cluster_summary, new_row)
  }
  
  # Reset row names of cluster_summary dataframe
  rownames(cluster_summary) <- NULL
  
  return(cluster_summary)
}


####### ANALYZE RNA #######

normalize_cluster_data <- function(patient_data, assay, patient_dims = 1:25, patient_res = 0.8) {
  
  # Set the assay to the one containing the RNA data
  DefaultAssay(patient_data) <- assay
  # Normalize the count data present in a given assay
  patient_data <- Seurat::NormalizeData(
    object = patient_data,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE)
  # Scales and centers features in the dataset
  patient_data <- Seurat::ScaleData(
    object = patient_data,
    model.use = "linear",
    verbose = FALSE)
  # Detect highly variable genes for the pca
  # Identifies features that are outliers on a 'mean variability plot'
  patient_data <- Seurat::FindVariableFeatures(
    object = patient_data,
    verbose = FALSE)
  # Run a PCA dimensionality reduction
  patient_data <- Seurat::RunPCA(
    object = patient_data,
    reduction.name = "pca_RNA",
    reduction.key = "PCRNA_",
    seed.use = 1,
    verbose = FALSE)
  # Computes the k.param nearest neighbors
  patient_data <- Seurat::FindNeighbors(
    object = patient_data,
    dims = patient_dims,
    reduction = "pca_RNA",
    verbose = FALSE)
  # Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm
  # Use the resolution parameter to fine tune the number of expected clusters
  patient_data <- Seurat::FindClusters(
    object = patient_data,
    resolution = patient_res,
    cluster.name = "RNA_clusters",
    graph.name = paste0(assay, "_snn"),
    random.seed = 1,
    verbose = FALSE)
  # Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique
  patient_data <- Seurat::RunUMAP(
    object = patient_data,
    n.neighbors = 30L,
    dims = patient_dims,
    repulsion.strength = 5,
    reduction = "pca_RNA",
    reduction.name = "umap_RNA",
    reduction.key = "UMAPRNA_",
    seed.use = 1,
    verbose = FALSE)
  
  return(patient_data)
}

####### COMPARE CLUSTERING METHODS #######

compare_clustering_methods <- function(patient_rna_data) {
  
  # Read patient number
  patient_num <- patient_rna_data$Patient.ID[1]
  
  # Create the contingency table (table used to study the correlation between the two variables)
  contingency_tab_clusters <- table(
    patient_rna_data$InSituType_semisup_clusters,
    patient_rna_data$RNA_clusters,
    dnn = c("InSituType", "Seurat"))
  
  # Convert the table to a data frame for ggplot2
  df_compare_clusters <- as.data.frame(as.table(contingency_tab_clusters))
  df_compare_clusters$Freq <- log10(df_compare_clusters$Freq + 10)
  
  # Plot using ggplot2 with geom_tile
  heatmap_seurat_vs_insitutype <- ggplot(df_compare_clusters, aes(InSituType, Seurat, fill = Freq)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c() +
    labs(fill = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = paste("Patient", patient_num),
      subtitle = "Seurat vs InSituType clusters"
    )
  
  return(heatmap_seurat_vs_insitutype)
}

####### ANALYZE PROTEINS #######

analyze_proteins <- function(patient_data) {
  
  metadata_proteins_columns <- c(
    "Mean.PanCK",
    "Max.PanCK",
    "Mean.CD68",
    "Max.CD68",
    "Mean.Membrane",
    "Max.Membrane",
    "Mean.CD45",
    "Max.CD45",
    "Mean.DAPI",
    "Max.DAPI")
  
  for (col in metadata_proteins_columns) {
    patient_data[[col]] <- as.numeric(as.factor(patient_data@meta.data[[col]]))
  }
  
  # Extract proteins features and replace underscores in column names
  proteins_features <- patient_data@meta.data[, metadata_proteins_columns]
  colnames(proteins_features) <- gsub("_", "-", colnames(proteins_features))
  proteins_matrix <- as.matrix(proteins_features)
  
  # Ensure row names (cell names) match the Seurat object cell names
  rownames(proteins_matrix) <- rownames(patient_data@meta.data)
  
  # Create a new assay with the proteins data
  proteins_assay <- CreateAssay5Object(counts = t(proteins_matrix))
  
  # Add the new assay to the Seurat object with the name "proteins"
  patient_data[["proteins"]] <- proteins_assay
  DefaultAssay(patient_data) <- "proteins"
  
  # Scale the raw metadata values
  patient_data <- ScaleData(
    layer = "counts",
    patient_data,
    assay = "proteins",
    features = metadata_proteins_columns,
    do.center = TRUE,
    do.scale = TRUE)
  
  # Run PCA and specify the number of PCs to compute
  # Only 10 features available, use 10 PCs
  n_pcs <- 9
  patient_data <- RunPCA(
    patient_data,
    assay = "proteins",
    features = metadata_proteins_columns,
    npcs = n_pcs,
    reduction.name = "pca_proteins",
    reduction.key = "PCPR_",
    seed.use = 2)
  
  # Use the available number of PCs for FindNeighbors and clustering
  patient_data <- FindNeighbors(patient_data, dims = 1:n_pcs, assay = "proteins", reduction = "pca_proteins")
  patient_data <- FindClusters(
    patient_data,
    resolution = 0.5,
    assay = "proteins",
    cluster.name = "protein_clusters",
    graph.name = "proteins_snn",
    random.seed = 2)
  
  # Run UMAP for visualization
  patient_data <- RunUMAP(
    patient_data,
    assay = "proteins",
    dims = 1:n_pcs,
    reduction = "pca_proteins",
    reduction.name = "umap_proteins",
    reduction.key = "UMAPPR_",
    seed.use = 2)
  
  return(patient_data)
  
}

####### GENERATE PLOTS #######

generate_feature_plot <- function(patient_data, reduction, features, max_cutoff = NA) {
  
  print(paste("Generate FeaturePlots from", reduction, "reduction of features", features))
  
  patient_num <- get_patient_num(patient_data)
  
  features_plots <- FeaturePlot(
    object = patient_data,
    features = features,
    reduction = reduction,
    max.cutoff = max_cutoff) +
    plot_annotation(
      title = 'Patient 1',
      subtitle = reduction,
    )  & NoLegend() & NoAxes()
    
  features_plot_name <- paste0("Patient_",  patient_num, "_featureplots_", reduction)
  return(setNames(list(features_plots), features_plot_name))
}

generate_elbow_plot<- function(patient_data, reduction, dims) {
  
  print(paste("Generate ElbowPlot from", reduction))
  
  patient_num <- get_patient_num(patient_data)
  
  elbow_plot <- ElbowPlot(patient_data, reduction = reduction, ndims = dims) +
    labs(title = paste("Patient", patient_num), subtitle = reduction)
  
  elbow_plot_name <- paste("Patient",  patient_num, "elbow_plot", reduction, sep = "_")
  return(setNames(list(elbow_plot), elbow_plot_name))
  
}

####### COLOR CLUSTERS #######

# Assign random colors to unknown clusters
generate_colors_lookup_table <- function(data, cluster_column_name, known_clusters_colors = NULL, color_palette = DiscretePalette(36, palette = "polychrome")) {
  
  # Remove the colors that are used by the known clusters
  usable_color_palette <- setdiff(color_palette, known_clusters_colors)
  
  # Random colors for the unknown clusters
  # The cluster 0 in Seurat for example is different for every patient
  usable_color_palette <- sample(usable_color_palette)
  
  unknown_clusters_colors <- c()
  next_color <- 1
  # Loop over the unknown clusters and assign each of them a random color
  for (cluster in setdiff(
    as.vector(as.character(unique(data[[cluster_column_name]])[[cluster_column_name]])),
    as.vector(names(known_clusters_colors)))) {
    unknown_clusters_colors <- c(unknown_clusters_colors, setNames(usable_color_palette[next_color], cluster))
    next_color <- next_color + 1
  }
  # Return the colors for the known and unknown clusters
  return(c(known_clusters_colors, unknown_clusters_colors))
}

####### PRINT CLUSTERING PLOTS #######

generate_clustering_plots <- function(
    patient_data,
    cluster_var,
    cluster_assay,
    cluster_reduction = "umap",
    create_heatmap = FALSE,
    cluster_name = NULL,
    color_lookup_table = NULL
    ) {
  
  # If a human friendly name is not given, use the name of the column in the Seurat object
  if (is.null(cluster_name)) {
    cluster_name <- cluster_var
  }
  
  # If no color lookup table is given create one
  if (is.null(color_lookup_table)) {
    color_lookup_table <- generate_colors_lookup_table(patient_data, cluster_var)
  }

  # Get name of the first image
  patient_image <- Images(patient_data)[1]
  
  # List to be returned with all the plots
  clustering_plots <- list()
  
  patient_num <- get_patient_num(patient_data)

  # Print some information about the clusters
  print(paste(cluster_name, "and number of cells in each of them associated with patient", patient_num))
  print(table(patient_data[[cluster_var]]))
  
  # Select the cluster as the identity
  Idents(patient_data) <- cluster_var
  # Plot the cells using their polygonal boundaries
  DefaultBoundary(patient_data[[patient_image]]) <- "segmentation"
  
  # Graphs the output of a dimensional reduction technique on a 2D scatter plot
  # Each point is a cell and it's positioned based on the cell embeddings determined by the reduction technique
  umap_clusters <- Seurat::DimPlot(
    patient_data,
    reduction = cluster_reduction,,
    group.by = cluster_var,
    label=TRUE,
    label.box=TRUE,
    repel=TRUE,
    cols = color_lookup_table) +
    labs(
      title = paste("Patient", patient_num),
      subtitle = cluster_name) +
    NoLegend()
  # Save plot to list
  clustering_plots[[paste("Patient",  patient_num, cluster_var, "umap", sep = "_")]] <- umap_clusters
  
  # Plot cells in their spatial context
  stamps_list <- list()
  for(curr_core in sort(unique(patient_data@meta.data$core_serial))) {
    for (curr_stamp in sort(unique(patient_data@meta.data$stamp[patient_data@meta.data$core_serial == curr_core]))) {
      
      # Subset data from current core and stamp
      core_stamp_subset <- subset(patient_data, subset = core_serial == curr_core & stamp == curr_stamp)
      
      stamp_plot <- ImageDimPlot(
        core_stamp_subset,
        fov = patient_image,
        # Set border color to 'NA' as 'white' masks all cells when zoomed out
        border.color = NA,
        flip_xy = FALSE,
        cols = color_lookup_table) + theme(
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8),
          legend.key.size = unit(0.5, 'lines'), # Adjust the size of the legend keys
          legend.spacing = unit(0.5, 'lines') # Adjust the spacing between legend items
        ) +
        labs(
          title = paste("Patient", patient_num, "Core", curr_core, ", Stamp", curr_stamp),
          subtitle = cluster_name
        )
      stamp_plot_name <- paste("Patient",  patient_num, cluster_var, "core", curr_core, "stamp", as.character(curr_stamp), sep = "_")
      clustering_plots[[stamp_plot_name]] <- stamp_plot
    }
  }
  
  # Create heatmap of most differentially expressed genes
  if (create_heatmap) {
    diff_expr_genes_heatmap <- generate_dyn_text_heatmap(
      patient_data,
      cluster_var,
      cluster_assay,
      color_lookup_table = color_lookup_table,
      cluster_name = cluster_name)
    # Save plot to list
    clustering_plots <- c(clustering_plots, diff_expr_genes_heatmap)
  }
  
  return(clustering_plots)
}

remove_clusters <- function(patient_data, cluster_col, cluster_vals) {
  
  # Identify cells that do not belong to the clusters specified in cluster_vals
  cells_to_keep <- !patient_data[[cluster_col, drop=TRUE]] %in% cluster_vals
  
  # Subset the Seurat object to keep only the desired cells
  patient_data_clusters_removed <- subset(patient_data, cells = Cells(patient_data)[cells_to_keep])
  
  return(patient_data_clusters_removed)
}

####### GENERATE PLOTS #######

# Know clusters that should have consistent colors
known_clusters_colors <- c(
  "B cell" = "#5A5156",
  "Dendritic cell" = "#E4E1E3",
  "Endothelial" = "#F6222E",
  "Fibroblast" = "#FE00FA",
  "Macrophage" = "#16FF32",
  "Mast cell" = "#3283FE",
  "Monocyte" = "#FEAF16",
  "Neutrophil" = "#B00068",
  "NK cell" = "#1CFFCE",
  "Plasma" = "#90AD1C",
  "Plasmablast" = "#2ED9FF",
  "Plasmacytoid dendritic cell" = "#DEA0FD",
  "T cell CD4" = "#AA0DFE",
  "T cell CD8" = "#F8A19F",
  "T cell regulatory" = "#325A9B",
  "Tumor" = "#C4451C",
  "T cell" = "#1C8356"
)

generate_proteins_plots <- function(patient_data, assay) {
  
  print("Generate proteins plots")
  
  patient_num <- get_patient_num(patient_data)
  
  # List with all the plots to be returned
  plot_list <- list()
  
  # Show the significance of every principal component of the PCA
  # It can be used to decide the number of dims of the FindNeighbors function
  elbow_plot_red = "pca_proteins"
  proteins_elbow_plot<- generate_elbow_plot(patient_data, elbow_plot_red, 9)
  plot_list <- c(plot_list, proteins_elbow_plot)
  
  proteins_features_plots <- generate_feature_plot(
    patient_data = patient_data,
    reduction = "umap_proteins",
    features = c("Mean.PanCK", "Mean.CD45", "Mean.CD68", "Mean.Membrane", "Mean.DAPI", "Area" ),
    max_cutoff = "q95")
  plot_list <- c(plot_list, proteins_features_plots)
  
  protein_cluster_var <- "protein_clusters"
  protein_color_lookup_table <- generate_colors_lookup_table(patient_data, protein_cluster_var, known_clusters_colors)
  protein_cluster <- generate_clustering_plots(
    patient_data,
    protein_cluster_var,
    cluster_assay = assay, 
    cluster_reduction = "umap_proteins",
    create_heatmap = FALSE,
    cluster_name = "Protein Clusters",
    color_lookup_table = protein_color_lookup_table)
  plot_list[[protein_cluster_var]] <- protein_cluster
  
  return(plot_list)
}

generate_rna_plots <- function(patient_data, assay, RNA_cluster_var) {
  
  print("Generate RNA plots")
  
  # List with all the plots to be returned
  plot_list <- list()
  
  patient_num <- get_patient_num(patient_data)
  
  # Show the significance of every principal component of the PCA
  # It can be used to decide the number of dims of the FindNeighbors function
  elbow_plot_red = "pca_RNA"
  # By default the RunPCA function uses 50 dimensions, plot all of them
  RNA_elbow_plot <- generate_elbow_plot(patient_data, elbow_plot_red, 50)
  plot_list <- c(plot_list, RNA_elbow_plot)
  
  RNA_features_plots <- generate_feature_plot(
    patient_data = patient_data,
    reduction = "umap_RNA",
    features = c("Mean.PanCK", "Mean.CD45", "Mean.CD68", "Mean.Membrane", "Mean.DAPI", "Area" ),
    max_cutoff = "q100")
  plot_list <- c(plot_list, RNA_features_plots)
  
  RNA_color_lookup_table <- generate_colors_lookup_table(patient_data, RNA_cluster_var, known_clusters_colors)
  RNA_cluster <- generate_clustering_plots(
    patient_data,
    RNA_cluster_var,
    cluster_assay = assay, 
    cluster_reduction = "umap_RNA",
    create_heatmap = TRUE,
    cluster_name = "RNA Clusters",
    color_lookup_table = RNA_color_lookup_table)
  plot_list[[RNA_cluster_var]] <- RNA_cluster
  
  return(plot_list)
}
 
generate_IST_plots <- function(patient_data, assay, IST_object) {
  
  print("Generate InSituType plots")
  
  # List with all the plots to be returned
  InSituType_clusters <- list()
  patient_num <- get_patient_num(patient_data)
  
  InSituType_cluster_var <- "InSituType_semisup_clusters"
  InSituType_color_lookup_table <- generate_colors_lookup_table(
    patient_data,
    InSituType_cluster_var,
    known_clusters_colors)
  
  InSituType_clusters <- c(
    InSituType_clusters,
    generate_flightpath(IST_object, InSituType_color_lookup_table, patient_num))
 
  InSituType_clusters <- c(InSituType_clusters, generate_clustering_plots(
    patient_data,
    InSituType_cluster_var,
    cluster_assay = assay, 
    cluster_reduction = "umap_RNA",
    create_heatmap = TRUE,
    cluster_name = "InSituType Semisupervised Clusters",
    color_lookup_table = InSituType_color_lookup_table))
  
  return(InSituType_clusters)
}

generate_comparison_plots <- function(patient_data) {
  
  print("Generate Comparison plots")
  
  # List with all the plots to be returned
  plot_list <- list()
  
  patient_num <- get_patient_num(patient_data)

  seurat_vs_insitutype_plot_name <- paste0("Patient_",  patient_num, "_heatmap_seurat_vs_insitutype")
  seurat_vs_insitutype_plot <- compare_clustering_methods(patient_data)
  plot_list[[seurat_vs_insitutype_plot_name]] <- seurat_vs_insitutype_plot
  
  return(plot_list)
}

####### ANALYZE PATIENT #######

analyze_patient <- function(all_patients_data, patient_num) {
  
  # List with all the plots to be returned
  plot_list <- list()
  
  # Create directories to save patient data if they do not exist
  patient_directories <- list(get_patient_dir_img(patient_num), get_patient_dir_rds(patient_num))
  for (dir in patient_directories) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }
  
  # Load patient data from RDS files if they exist, otherwise generate them
  patient_rna_rds <- here(
    rds_dir, 
    paste0("Patient_", patient_num, "_data/", "Patient_", patient_num, "_rna_data.rds"))
  patient_misc_rds <- here(
    rds_dir,
    paste0("Patient_", patient_num, "_data/", "Patient_", patient_num, "_misc_data.rds"))
  if ((!file.exists(patient_rna_rds)) || (!file.exists(patient_misc_rds))) {
    print("Generating patient data")
    
    # Extract patient data
    all_patient_data <- extract_patient_data(all_patients_data, patient_num)
    patient_rna_only <- all_patient_data[[1]]
    patient_cohort <- all_patient_data[[2]]
    patient_rna_counts <- all_patient_data[[3]]
    patient_avg_neg_probes <- all_patient_data[[4]]
    
    # Normalize, scale, cluster, ...
    patient_rna_only <- normalize_cluster_data(
      patient_rna_only,
      assay = "Nanostring",
      patient_dims = 1:25,
      patient_res = 0.8)
    
    # Run InSituType semisupervised clustering
    patient_semisup <- run_InSituType_semisupervised(
      patient_rna_only,
      patient_cohort,
      patient_rna_counts,
      patient_avg_neg_probes)
    # add phenotypes to the metadata for plotting
    patient_rna_only$InSituType_semisup_clusters <- patient_semisup$clust
    
    # Analyze protein data
    patient_rna_only <- analyze_proteins(patient_rna_only)

    # Save the data to RDS files
    saveRDS(patient_rna_only, file = patient_rna_rds)
    patient_misc_data <- list(patient_cohort, patient_rna_counts, patient_avg_neg_probes)
    saveRDS(patient_misc_data, file = patient_misc_rds)
    
  } else {
    print("Reading patient data from RDS files")
    
    patient_rna_only <- readRDS(patient_rna_rds)
    patient_misc_data <- readRDS(patient_misc_rds)
    patient_cohort <- patient_misc_data[[1]]
    patient_rna_counts <- patient_misc_data[[2]]
    patient_avg_neg_probes <- patient_misc_data[[3]]
  }
  
  # Print patient information
  print_patient_info(patient_rna_only)
  
  # Return/print all plots together, otherwise only the last one is shown
  plot_list <- generate_plots(patient_rna_only)
  return(list(patient_rna_only, plot_list))
}
