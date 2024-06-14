# Import environment variables as global variables
rds_dir <- Sys.getenv("RDS_DIR")
image_dir <- Sys.getenv("IMAGE_DIR")
image_ext <- Sys.getenv("IMAGE_EXT")

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

generate_dyn_text_heatmap <- function(patient_data, clusters_column_name, assay_name, group_colors = NULL) {
  
  # Select the cluster as the identity
  Idents(patient_data) <- clusters_column_name

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
  patient_num <- patient_data$Patient.ID[1]
  
  # Create the heatmap with the scaled text
  diff_expr_genes_heatmap <- DoHeatmap(
    data_filtered,
    features = most_significant_markers$gene,
    assay = assay_name,
    label = TRUE,
    size = calculate_clusters_names_size(smallest_cluster_proportion_filtered),
    group.colors = group_colors
  ) + theme(
    axis.text.y = element_text(size = label_size),
  ) + labs(
    title = paste("Patient", patient_num, clusters_column_name),
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

####### PRINT PATIENT DATA #######

get_patient_info <- function(patient_data) {
  
  patient_num <- patient_data$Patient.ID[1]
  
  # Inilitize list to store the cores, stamps and fovs associated with the patient
  cores <- list()
  # Print cores and fovs associated with the patient
  print(paste("FOVs and cell count associated with patient", patient_num))
  
  # Get cores for the current patient
  patient_cores <- sort(unique(patient_data@meta.data$core_serial))
  for(core in patient_cores) {
    
    print(paste("----------", "CORE", core, "----------"))
    stamps <- list()
    
    # Get stamps for the current core
    patient_core_stamps <- sort(unique(patient_data@meta.data$stamp[patient_data@meta.data$core_serial == core]))
    for(stamp in patient_core_stamps) {
      
      print(paste("-----", "STAMP", stamp, "-----"))
      fovs <- table(patient_data@meta.data$fov[patient_data@meta.data$core_serial == core &  patient_data@meta.data$stamp == stamp])
      print(fovs)
      
      stamps[[stamp]] <- as.list(fovs)
    }
    cores[[core]] <- stamps
  }
  return(invisible(cores))
}

####### ANALYZE RNA #######

normalize_cluster_data <- function(patient_rna_only, patient_dims, patient_res) {
  
  # Normalize the count data present in a given assay
  patient_rna_only <- Seurat::NormalizeData(patient_rna_only, assay = "Nanostring")
  # Scales and centers features in the dataset
  patient_rna_only <- Seurat::ScaleData(patient_rna_only)
  # Detect highly variable genes for the pca
  # Identifies features that are outliers on a 'mean variability plot'
  patient_rna_only <- Seurat::FindVariableFeatures(patient_rna_only)
  # Run a PCA dimensionality reduction
  patient_rna_only <- Seurat::RunPCA(patient_rna_only, seed.use = 1)
  # Computes the k.param nearest neighbors
  patient_rna_only <- Seurat::FindNeighbors(patient_rna_only, dims = patient_dims)
  # Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm
  # Use the resolution parameter to fine tune the number of expected clusters
  patient_rna_only <- Seurat::FindClusters(
    patient_rna_only,
    resolution = patient_res,
    random.seed = 1,
    cluster.name = "RNA_clusters")
  # Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique
  patient_rna_only <- Seurat::RunUMAP(patient_rna_only, dims = patient_dims, repulsion.strength = 5, seed.use = 1)
  
  return(patient_rna_only)
  
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
    reduction.key = "PCPR_")
  
  # Use the available number of PCs for FindNeighbors and clustering
  patient_data <- FindNeighbors(patient_data, dims = 1:n_pcs, assay = "proteins", reduction = "pca_proteins")
  patient_data <- FindClusters(
    patient_data,
    resolution = 0.5,
    assay = "proteins",
    cluster.name="protein_clusters",
    graph.name = "proteins_snn")
  
  # Run UMAP for visualization
  patient_data <- RunUMAP(
    patient_data,
    assay = "proteins",
    dims = 1:n_pcs,
    reduction.name = "umap_proteins",
    reduction.key = "UMAPPR_")
  
  return(patient_data)
  
}

####### PRINT PROTEINS DATA #######

print_proteins_data <- function(patient_data, patient_num, patient_dir_img, patient_dir_rds_img) {
  
  print("Generate FeaturePlot from protein data")
  protein_data_feature_plots <- paste0("Patient_",  patient_num, "_protein_feature_plots")
  protein_data_feature_plots_rds <- paste0(patient_dir_rds_img, "Patient_",  protein_data_feature_plots, ".rds")
  if (!file.exists(protein_data_feature_plots_rds)) {
    
    protein_plots <- FeaturePlot(
      object = patient_data,
      features = c("Mean.PanCK", "Mean.CD45", "Mean.CD68", "Mean.Membrane", "Mean.DAPI", "Area" ),
      reduction = "umap_proteins",
      max.cutoff = "q95") +
      plot_annotation(
        title = 'Patient 1',
        subtitle = 'UMAP from protein data',
      )  & NoLegend() & NoAxes()
    
    saveRDS(protein_plots, file = protein_data_feature_plots_rds)
    protein_plots_image <- paste0(patient_dir_img, protein_data_feature_plots, image_ext)
    ggsave(filename = protein_plots_image, plot = protein_plots)
  } else {
    protein_plots <- readRDS(protein_data_feature_plots_rds)
  }
  
  return(list(protein_plots))
}

####### ANALYZE PATIENT #######

analyze_patient <- function(all_patients_data, patient_num) {
  
  # Create directories to save patient data if they do not exist
  patient_dir_img <- paste0(image_dir, "Patient_", patient_num, "_plots/")
  patient_dir_rds <- paste0(rds_dir, "Patient_", patient_num, "_data/")
  patient_dir_rds_img <- paste0(rds_dir, "Patient_", patient_num, "_data/", "Img_data/")
  patient_directories <- list(patient_dir_img, patient_dir_rds, patient_dir_rds_img)
  for (dir in patient_directories) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }
  
  # Load patient data from RDS files if they exist, otherwise generate them
  patient_rna_rds <- paste0(rds_dir, "Patient_", patient_num, "_data/", "Patient_", patient_num, "_rna_data.rds")
  patient_misc_rds <- paste0(rds_dir, "Patient_", patient_num, "_data/", "Patient_", patient_num, "_misc_data.rds")
  if ((!file.exists(patient_rna_rds)) || (!file.exists(patient_misc_rds))) {
    print("Generating patient data")
    
    # Extract patient data
    all_patient_data <- extract_patient_data(all_patients_data, patient_num)
    patient_rna_only <- all_patient_data[[1]]
    patient_cohort <- all_patient_data[[2]]
    patient_rna_counts <- all_patient_data[[3]]
    patient_avg_neg_probes <- all_patient_data[[4]]
    
    # Normalize, scale, cluster, ...
    patient_rna_only <- normalize_cluster_data(patient_rna_only, patient_dims = 1:25, patient_res = 0.8)
    
    # Get cell reference profile data from NanoString
    # Use this reference profile as it is the only one available from CosMx data, originally from:
    # https://raw.githubusercontent.com/Nanostring-Biostats/CosMx-Cell-Profiles/main/Human/IO/IO.profiles.csv
    ioprofiles <- read.csv("./Cell_profile_matrices/NanoString.CosMx.Human.IO.profiles.csv", header = T, sep = ",", fill = T)
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
  get_patient_info(patient_rna_only)
  
  # Show the significance of every principal component of the PCA
  # It can be used to decide the number of dims of the FindNeighbors function
  print("Generate Elbow plot")
  elbow_plot_name <- paste0("Patient_",  patient_num, "_elbow_plot")
  elbow_plot_rds <- paste0(patient_dir_rds_img, elbow_plot_name, ".rds")
  if (!file.exists(elbow_plot_rds)) {
    elbow_plot <- ElbowPlot(patient_rna_only, ndims = 50) + labs(title = paste("Patient", patient_num))
    saveRDS(elbow_plot, file = elbow_plot_rds)
    elbow_plot_image <- paste0(patient_dir_img, elbow_plot_name, image_ext)
    ggsave(filename = elbow_plot_image,plot = elbow_plot)
  } else {
    elbow_plot <- readRDS(elbow_plot_rds)
  }
  
  # Plot Mean PanCK
  print("Generate Mean PanCK plot")
  mean_panck_plot_name <- paste0("Patient_",  patient_num, "_panCK")
  mean_panck_plot_rds <- paste0(patient_dir_rds_img, mean_panck_plot_name, ".rds")
  if (!file.exists(mean_panck_plot_rds)) {
    mean_panck_plot <- FeaturePlot(object = patient_rna_only, features = "Mean.PanCK", min.cutoff = 2000) +
      labs(title = paste("Patient", patient_num), subtitle = "Mean PanCK")
    saveRDS(mean_panck_plot, file = mean_panck_plot_rds)
    mean_panck_plot_image <- paste0(patient_dir_img, mean_panck_plot_name, image_ext)
    ggsave(filename = mean_panck_plot_image, plot = mean_panck_plot)
  } else {
    mean_panck_plot <- readRDS(mean_panck_plot_rds)
  }
  
  # Plot KRT17
  print("Generate KRT17 plot")
  krt17_plot_name <- paste0("Patient_",  patient_num, "_krt17")
  krt17_plot_rds <- paste0(patient_dir_rds_img, krt17_plot_name, ".rds")
  if (!file.exists(krt17_plot_rds)) {
    KRT17_plot <- FeaturePlot(
      object = patient_rna_only,
      features = "KRT17",
      cols = c("white", "red")
    ) + 
      labs(
        title = paste("Patient", patient_num),
        subtitle = "KRT17"
      )
    saveRDS(KRT17_plot, file = krt17_plot_rds)
    KRT17_plot_image <- paste0(patient_dir_img, krt17_plot_name, image_ext)
    ggsave(filename = KRT17_plot_image, plot = KRT17_plot)
  } else {
    KRT17_plot <- readRDS(krt17_plot_rds)
  }
  
  protein_data_plots <- print_proteins_data(patient_rna_only, patient_num, patient_dir_img, patient_dir_rds_img)
  
  ################## PRINT CLUSTERING PLOTS ##################
  
  # Get name of the first image
  patient_image <- Images(patient_rna_only)[1]
  
  # Create a list of the clusters
  RNA_cluster <- list(name = "RNA Clusters", var = "RNA_clusters", assay = "Nanostring", reduction = "umap", heatmap = TRUE)
  InSituType_cluster <- list(name = "InSituType Semisupervised Clusters", var = "InSituType_semisup_clusters", assay = "Nanostring", reduction = "umap", heatmap = TRUE)
  protein_cluster <- list(name = "Protein Clusters", var = "protein_clusters", assay = "proteins", reduction = "umap_proteins", heatmap = FALSE)
  cell_clusters <- list(RNA_cluster, InSituType_cluster, protein_cluster)
  
  # Color of the clusters
  clusters_colors <- DiscretePalette(20, palette = "polychrome", shuffle = TRUE)
  # Print the color palette
  # pie(rep(1, length(clusters_colors)), col = clusters_colors , main="")
  
  # Save all the plots in a list to return them all together
  clustering_plots_list <- list()
  
  for (cell_cluster in cell_clusters) {
    
    print(paste(cell_cluster$name, "and number of cells in each of them associated with patient", patient_num))
    print(table(patient_rna_only[[cell_cluster$var]]))
    
    # Select the cluster as the identity
    Idents(patient_rna_only) <- cell_cluster$var
    # Plot the cells using their polygonal boundaries
    DefaultBoundary(patient_rna_only[[patient_image]]) <- "segmentation"
    
    if (cell_cluster$heatmap) {
      diff_expr_genes_heatmap <- generate_dyn_text_heatmap(patient_rna_only, cell_cluster$var, cell_cluster$assay, clusters_colors)
      
      # Save plot to list
      clustering_plots_list[[paste(cell_cluster$var, "heatmap", sep = "_")]] <- diff_expr_genes_heatmap
      
      # Save the heatmap to an image
      ggsave(
        filename = paste0(patient_dir_img, "Patient_",  patient_num, "_", cell_cluster$var, "_diff_expr_genes_heatmap", image_ext),
        plot = diff_expr_genes_heatmap
      )
    }
    
    # Graphs the output of a dimensional reduction technique on a 2D scatter plot
    # Each point is a cell and it's positioned based on the cell embeddings determined by the reduction technique
    umap_clusters <- Seurat::DimPlot(
      patient_rna_only,
      reduction = cell_cluster$reduction,,
      group.by = cell_cluster$var,
      label=TRUE,
      label.box=TRUE,
      repel=TRUE,
      cols = clusters_colors) +
      labs(
        title = paste("Patient", patient_num),
        subtitle = cell_cluster$name) +
      NoLegend()
    # Save plot to list
    clustering_plots_list[[paste(cell_cluster$var, "umap", sep = "_")]] <- umap_clusters
    # Save plot to image file
    ggsave(
      filename = paste0(patient_dir_img, "Patient_",  patient_num, "_", cell_cluster$var,"_umap", image_ext),
      plot = umap_clusters
    )
    
    # Plot cells in their spatial context
    stamps_list <- list()
    for(curr_core in sort(unique(patient_rna_only@meta.data$core_serial))) {
      for (curr_stamp in sort(unique(patient_rna_only@meta.data$stamp[patient_rna_only@meta.data$core_serial == curr_core]))) {
        
        # Subset data from current core and stamp
        core_stamp_subset <- subset(patient_rna_only, subset = core_serial == curr_core & stamp == curr_stamp)
        
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
            subtitle = cell_cluster$name
          )
        clustering_plots_list[[paste(cell_cluster$var, curr_core, as.character(curr_stamp), sep = "_")]] <- stamp_plot
        ggsave(
          filename = paste0(patient_dir_img, "Patient_",  patient_num, "_", cell_cluster$var, "_core_", curr_core, "_stamp_", curr_stamp, image_ext),
          plot = stamp_plot)
      }
    }
  }
  
  print("Generate Comparing Clusters plot")
  seurat_vs_insitutype_plot_name <- paste0("Patient_",  patient_num, "_heatmap_seurat_vs_insitutype")
  seurat_vs_insitutype_plot_rds <- paste0(patient_dir_rds_img, seurat_vs_insitutype_plot_name, ".rds")
  if (!file.exists(seurat_vs_insitutype_plot_rds)) {
    seurat_vs_insitutype_plot <- compare_clustering_methods(patient_rna_only)
    saveRDS(seurat_vs_insitutype_plot, file = seurat_vs_insitutype_plot_rds)
    seurat_vs_insitutype_plot_image <- paste0(patient_dir_img, seurat_vs_insitutype_plot_name, image_ext)
    ggsave(
      filename = seurat_vs_insitutype_plot_image,
      plot = seurat_vs_insitutype_plot
    )
  } else {
    seurat_vs_insitutype_plot <- readRDS(seurat_vs_insitutype_plot_rds)
  }
  
  # Return all plots together, otherwise only the last one is shown
  plot_list <- c(list(elbow_plot, mean_panck_plot, KRT17_plot, protein_data_plots, seurat_vs_insitutype_plot), clustering_plots_list)
  print(plot_list)
  
  return(patient_rna_only)
}