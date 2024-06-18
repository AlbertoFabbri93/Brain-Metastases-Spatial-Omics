# Import environment variables as global variables
rds_dir <- Sys.getenv("RDS_DIR")
image_dir <- Sys.getenv("IMAGE_DIR")
image_ext <- Sys.getenv("IMAGE_EXT")

####### SAVE PLOTS #######

save_plots <- function(plots_list, folder_path, file_extension = ".png") {
  # Ensure the folder exists, if not, create it
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }
  
  # Helper function to flatten the list and collect plots
  # Images are also lists so simply unflattening the list does not work
  flatten_list <- function(lst) {
    flat_list <- list()
    
    for (name in names(lst)) {
      if (inherits(lst[[name]], "ggplot") || inherits(lst[[name]], "trellis")) {
        flat_list[[name]] <- lst[[name]]
      } else if (is.list(lst[[name]])) {
        flat_list <- c(flat_list, flatten_list(lst[[name]]))
      }
    }
    
    return(flat_list)
  }
  
  # Flatten the plots list
  flat_plots <- flatten_list(plots_list)
  
  # Iterate over the flattened list of plots and save each plot
  for (plot_name in names(flat_plots)) {
    # Construct the file path
    file_path <- file.path(folder_path, paste0(plot_name, file_extension))
    
    # Save the plot to the specified file path
    # Using png() resulted in low quality images
    ggsave(filename = file_path, plot = flat_plots[[plot_name]])
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

  # Calculate the size of each cluster
  cluster_sizes <- table(Idents(patient_data))
  smallest_cluster_size <- min(cluster_sizes)
  total_cells <- sum(cluster_sizes)
  cluster_proportions <- cluster_sizes / total_cells
  smallest_cluster_proportion <- tail(cluster_proportions, n = 1)
  
  # Identify clusters above threshold
  clusters_above_threshold <- names(cluster_proportions[cluster_proportions >= 0.01])
  # Subset the Seurat object to include only cells from clusters above the threshold
  data_filtered <- subset(patient_data, idents = clusters_above_threshold)
  # Find markers (differentially expressed genes) for each of the identity classes in the filtered dataset
  # If you have more than one assay it is necessary to specify the assay parameter
  markers.data_filtered <- FindAllMarkers(data_filtered, assay = assay_name, only.pos = TRUE)

  # Calculate the size of each cluster
  cluster_sizes_filtered <- table(Idents(data_filtered))
  smallest_cluster_size_filtered <- min(cluster_sizes_filtered)
  total_cells_filtered <- sum(cluster_sizes_filtered)
  cluster_proportions_filtered <- cluster_sizes_filtered / total_cells_filtered
  smallest_cluster_proportion_filtered <- min(cluster_proportions_filtered)
  
  # Filter markers to get the most significant ones per cluster
  most_significant_markers <- markers.data_filtered %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  # Calculate the number of features (rows)
  num_features <- length(most_significant_markers$gene)
  
  # Calculate the label size for the current figure
  label_size <- calculate_label_size(num_features)
  
  # Read patient number
  patient_num <- get_patient_num(patient_data)
  
  # Create the heatmap with the scaled text
  diff_expr_genes_heatmap <- DoHeatmap(
    data_filtered,
    features = most_significant_markers$gene,
    assay = assay_name,
    label = TRUE,
    size = calculate_clusters_names_size(smallest_cluster_proportion_filtered),
    group.colors = color_lookup_table,
  ) + theme(
    axis.text.y = element_text(size = label_size),
  ) + labs(
    title = paste("Patient", patient_num, cluster_name),
    subtitle = "Top 10 Differentially Expressed Genes per Cluster"
  )
  
  return(diff_expr_genes_heatmap)
}

####### INSITUTYPE SEMISUPERVISED #######

runInSituTypeSemisupervised <- function(patient_data, ioprofiles, patient_cohort, patient_rna_counts, patient_avg_neg_probes) {
  
  # Semi-supervised learning with insitutype and reference profiles
  # InSituType needs integers, if given floating point numbers it fails with misleading errors
  patient_semisup <- insitutype(
    x = patient_rna_counts,
    neg = patient_avg_neg_probes,
    cohort = patient_cohort,
    reference_profiles = ioprofiles,
    
    # Enter your own per-cell background estimates here if you
    # have them; otherwise insitutype will use the negprobes to
    # estimate background for you.
    bg = NULL,
    # condensed to save time. n_clusts = 5:15 would be more optimal
    # Group the cells the do not correspond to any type in the reference matrix
    n_clusts = c(5),
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
  
  # Create Seurat object with only the RNA data
  patient_rna_only <- subset(x = patient_data, features = rownames(patient_data)[1:1000])
  
  # Cohort of the patient
  patient_immunofluorescence <- patient_data@meta.data %>% select("Mean.PanCK", "Mean.CD45", "Mean.CD68")
  # "Gaussian_transform = TRUE" maps variables to gaussians in order to place dramatically different variables on the same scale
  patient_cohort <- fastCohorting(patient_immunofluorescence, gaussian_transform = TRUE, n_cohorts = 5)
  # check clusters and cohort numbers
  table(patient_cohort)
  
  # Extract the count data from the Seurat object
  patient_rna_counts <- GetAssayData(
    patient_rna_only,
    layer = "counts") %>%
  as.matrix() %>%
  t()
  # Extract the negative probes from the Seurat object
  patient_neg_probes <- GetAssayData(
      subset(
        patient_data,
        features = row.names(GetAssayData(patient_data)) %>%
  grep("Negative", ., value = TRUE))) %>%
  as.matrix() %>%
  t()
  # Calculate the average negative probes per cell
  patient_avg_neg_probes <- Matrix::rowMeans(patient_neg_probes)
  
  return(list(patient_rna_only, patient_cohort, patient_rna_counts, patient_avg_neg_probes))
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
      fovs <- table(patient_data@meta.data$fov[patient_data@meta.data$core_serial == core & patient_data@meta.data$stamp == stamp])
      
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

library(Seurat)

# Define function to extract cluster information
get_cluster_info <- function(seurat_obj, clustering_column) {
  
  # Check if clustering_column is valid
  if (!clustering_column %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", clustering_column, "not found in meta.data of Seurat object."))
  }
  
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


####### ANALYZE RNA #######

normalize_cluster_data <- function(patient_data, assay, patient_dims, patient_res) {
  
  # Set the assay to the one containing the RNA data
  DefaultAssay(patient_data) <- assay
  # Normalize the count data present in a given assay
  patient_data <- Seurat::NormalizeData(patient_data)
  # Scales and centers features in the dataset
  patient_data <- Seurat::ScaleData(patient_data)
  # Detect highly variable genes for the pca
  # Identifies features that are outliers on a 'mean variability plot'
  patient_data <- Seurat::FindVariableFeatures(patient_data)
  # Run a PCA dimensionality reduction
  patient_data <- Seurat::RunPCA(
    patient_data,
    reduction.name = "pca_RNA",
    reduction.key = "PCRNA_",
    seed.use = 1)
  # Computes the k.param nearest neighbors
  patient_data <- Seurat::FindNeighbors(patient_data, dims = patient_dims, reduction = "pca_RNA")
  # Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm
  # Use the resolution parameter to fine tune the number of expected clusters
  patient_data <- Seurat::FindClusters(
    patient_data,
    resolution = patient_res,
    cluster.name = "RNA_clusters",
    graph.name = paste0(assay, "_snn"),
    random.seed = 1)
  # Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique
  patient_data <- Seurat::RunUMAP(
    patient_data,
    dims = patient_dims,
    repulsion.strength = 5,
    reduction = "pca_RNA",
    reduction.name = "umap_RNA",
    reduction.key = "UMAPRNA_",
    seed.use = 1)
  
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
  proteins_assay <- CreateAssayObject(counts = t(proteins_matrix))
  
  # Add the new assay to the Seurat object with the name "proteins"
  patient_data[["proteins"]] <- proteins_assay
  DefaultAssay(patient_data) <- "proteins"
  
  # Scale the raw metadata values
  features_to_scale <- colnames(proteins_matrix)
  patient_data <- ScaleData(patient_data, assay = "proteins", features = features_to_scale, do.center = TRUE, do.scale = TRUE)
  
  # Run PCA and specify the number of PCs to compute
  # Only 10 features available, use 10 PCs
  n_pcs <- 9
  patient_data <- RunPCA(
    patient_data,
    assay = "proteins",
    features = features_to_scale,
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

####### GENERATE FEATUREPLOTS DATA #######

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
    
  features_plot_name <- paste0("Patient_",  patient_num, "_", reduction)
  return(setNames(list(features_plots), features_plot_name))
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
    clustering_plots[[paste("Patient", patient_num, cluster_var, "diff_expr_genes_heatmap", sep = "_")]] <- diff_expr_genes_heatmap
  }
  
  return(clustering_plots)
}

####### ANALYZE PATIENT #######

analyze_patient <- function(all_patients_data, patient_num) {
  
  # List with all the plots to be printed
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
    
    # Get cell reference profile data from NanoString
    # Use this reference profile as it is the only one available from CosMx data, originally from:
    # https://raw.githubusercontent.com/Nanostring-Biostats/CosMx-Cell-Profiles/main/Human/IO/IO.profiles.csv
    ioprofiles <- read.csv(
      file = here("Analysis", "Cell_profile_matrices", "NanoString.CosMx.Human.IO.profiles.csv"),
      header = T,
      sep = ",",
      fill = T)
    rownames(ioprofiles) <- ioprofiles[, 1]
    ioprofiles <- ioprofiles[, -1] %>% as.matrix()
    
    # Run InSituType semisupervised clustering
    patient_semisup <- runInSituTypeSemisupervised(
      patient_rna_only,
      ioprofiles,
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
  
  ################## PRINT CLUSTERING PLOTS ##################
  
  # Know clusters that should have consistent colors
  known_clusters_colors <- c(
    "B.cell" = "#5A5156",
    "Dendritic.cell" = "#E4E1E3",
    "Endothelial" = "#F6222E",
    "Fibroblast" = "#FE00FA",
    "Macrophage" = "#16FF32",
    "Mast.cell" = "#3283FE",
    "Monocyte" = "#FEAF16",
    "Neutrophil" = "#B00068",
    "NK.cell" = "#1CFFCE",
    "Plasma" = "#90AD1C",
    "Plasmablast" = "#2ED9FF",
    "Plasmacytoid.dendritic.cell" = "#DEA0FD",
    "T.cell.CD4" = "#AA0DFE",
    "T.cell.CD8" = "#F8A19F",
    "T.cell.regulatory" = "#325A9B"
  )
  
  proteins_features_plots <- generate_feature_plot(
    patient_data = patient_rna_only,
    reduction = "umap_proteins",
    features = c("Mean.PanCK", "Mean.CD45", "Mean.CD68", "Mean.Membrane", "Mean.DAPI", "Area" ),
    max_cutoff = "q95")
  plot_list <- c(plot_list, proteins_features_plots)
  
  protein_cluster_var <- "protein_clusters"
  protein_color_lookup_table <- generate_colors_lookup_table(patient_rna_only, protein_cluster_var, known_clusters_colors)
  protein_cluster <- generate_clustering_plots(
    patient_rna_only,
    protein_cluster_var,
    "proteins", 
    cluster_reduction = "umap_proteins",
    create_heatmap = FALSE,
    cluster_name = "Protein Clusters",
    color_lookup_table = protein_color_lookup_table)
  plot_list[[protein_cluster_var]] <- protein_cluster
  
  # Show the significance of every principal component of the PCA
  # It can be used to decide the number of dims of the FindNeighbors function
  print("Generate Elbow plot")
  elbow_plot_red = "pca_RNA"
  elbow_plot_name <- paste("Patient",  patient_num, elbow_plot_red, "elbow_plot", sep = "_")
  # elbow_plot_rds <- paste0(patient_dir_rds_img, elbow_plot_name, ".rds")
  # if (!file.exists(elbow_plot_rds)) {
  elbow_plot <- ElbowPlot(patient_rna_only, reduction = elbow_plot_red, ndims = 50) +
    labs(title = paste("Patient", patient_num), subtitle = elbow_plot_red)
  #   saveRDS(elbow_plot, file = elbow_plot_rds)
  # } else {
  #   elbow_plot <- readRDS(elbow_plot_rds)
  # }
  plot_list[[elbow_plot_name]] <- elbow_plot
  
  RNA_features_plots <- generate_feature_plot(
    patient_data = patient_rna_only,
    reduction = "umap_RNA",
    features = c("Mean.PanCK", "Mean.CD45", "Mean.CD68", "Mean.Membrane", "Mean.DAPI", "Area" ),
    max_cutoff = "q100")
  plot_list <- c(plot_list, RNA_features_plots)
  
  RNA_cluster_var <- "RNA_clusters"
  RNA_color_lookup_table <- generate_colors_lookup_table(patient_rna_only, RNA_cluster_var, known_clusters_colors)
  RNA_cluster <- generate_clustering_plots(
    patient_rna_only,
    RNA_cluster_var,
    "Nanostring", 
    cluster_reduction = "umap_RNA",
    create_heatmap = TRUE,
    cluster_name = "RNA Clusters",
    color_lookup_table = RNA_color_lookup_table)
  plot_list[[RNA_cluster_var]] <- RNA_cluster
  
  InSituType_cluster_var <- "InSituType_semisup_clusters"
  InSituType_color_lookup_table <- generate_colors_lookup_table(patient_rna_only, InSituType_cluster_var, known_clusters_colors)
  InSituType_cluster <- generate_clustering_plots(
    patient_rna_only,
    InSituType_cluster_var,
    "Nanostring", 
    cluster_reduction = "umap_RNA",
    create_heatmap = TRUE,
    cluster_name = "InSituType Semisupervised Clusters",
    color_lookup_table = InSituType_color_lookup_table)
  plot_list[[InSituType_cluster_var]] <- InSituType_cluster
  
  # List to be returned with all the plots
  clustering_plots_list <- list(RNA_cluster, InSituType_cluster, protein_cluster)
  
  print("Generate Comparing Clusters plot")
  seurat_vs_insitutype_plot_name <- paste0("Patient_",  patient_num, "_heatmap_seurat_vs_insitutype")
  # seurat_vs_insitutype_plot_rds <- paste0(patient_dir_rds_img, seurat_vs_insitutype_plot_name, ".rds")
  # if (!file.exists(seurat_vs_insitutype_plot_rds)) {
    seurat_vs_insitutype_plot <- compare_clustering_methods(patient_rna_only)
  #   saveRDS(seurat_vs_insitutype_plot, file = seurat_vs_insitutype_plot_rds)
  # } else {
  #   seurat_vs_insitutype_plot <- readRDS(seurat_vs_insitutype_plot_rds)
  # }
  # Return/print all plots together, otherwise only the last one is shown
  plot_list[[seurat_vs_insitutype_plot_name]] <- seurat_vs_insitutype_plot
  return(list(patient_rna_only, plot_list))
}
