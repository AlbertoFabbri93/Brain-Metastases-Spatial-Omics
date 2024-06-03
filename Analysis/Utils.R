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

generateDynTextHeatmap <- function(data, clusters_column_name) {
  
  # Select the cluster as the identity
  Idents(data) <- clusters_column_name
  
  # Calculate the size of each cluster
  cluster_sizes <- table(Idents(data))
  smallest_cluster_size <- min(cluster_sizes)
  total_cells <- sum(cluster_sizes)
  cluster_proportions <- cluster_sizes / total_cells
  smallest_cluster_proportion <- tail(cluster_proportions, n = 1)
  
  # Identify clusters above threshold
  clusters_above_threshold <- names(cluster_proportions[cluster_proportions >= 0.01])
  # Subset the Seurat object to include only cells from clusters above the threshold
  data_filtered <- subset(data, idents = clusters_above_threshold)
  # Find markers (differentially expressed genes) for each of the identity classes in the filtered dataset
  markers.data_filtered <- FindAllMarkers(data_filtered, only.pos = TRUE)
  
  # Calculate the size of each cluster
  cluster_sizes_filtered <- table(Idents(data_filtered))
  smallest_cluster_size_filtered <- min(cluster_sizes_filtered)
  total_cells_filtered <- sum(cluster_sizes_filtered)
  cluster_proportions_filtered <- cluster_sizes_filtered / total_cells_filtered
  smallest_cluster_proportion_filtered <- min(cluster_proportions_filtered)
  
  # Filter markers to get the most significant ones per cluster
  most_significant_markers <- markers.data_filtered %>%
    group_by(cluster) %>%
    filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  # Calculate the number of features (rows)
  num_features <- length(most_significant_markers$gene)
  
  # Calculate the label size for the current figure
  label_size <- calculate_label_size(num_features)
  
  # Read patient number
  patient_num <- data$Patient.ID[1]
  
  # Create the heatmap with the scaled text
  diff_expr_genes_heatmap <- DoHeatmap(
    data_filtered,
    features = most_significant_markers$gene,
    assay = "Nanostring",
    label = TRUE,
    size = calculate_clusters_names_size(smallest_cluster_proportion_filtered),
  ) + theme(
    axis.text.y = element_text(size = label_size),
  ) + labs(
    title = paste("Patient", patient_num, clusters_column_name),
    subtitle = "Top 10 Differentially Expressed Genes per Cluster"
  )
  
  return(diff_expr_genes_heatmap)
}

####### INSITUTYPE SEMISUPERVISED #######

runInSituTypeSemisupervised <- function(patient_data, ioprofiles) {

  # Cohort of all the patients
  patient_immunofluorescence <- patient_data@meta.data %>% select("Mean.PanCK", "Mean.CD45", "Mean.CD68")
  # "Gaussian_transform = TRUE" maps variables to gaussians in order to place dramatically different variables on the same scale
  patient_cohort <- fastCohorting(patient_immunofluorescence, gaussian_transform = TRUE, n_cohorts = 5)
  # check clusters and cohort numbers
  table(patient_cohort)
  
  # Extract the count data from the Seurat object
  patient_rna.counts <- GetAssayData(subset(patient_data, features = row.names(GetAssayData(patient_data))[1:1000]), layer = "counts") %>%
    as.matrix() %>%
    t()
  # Extract the negative probes from the Seurat object
  patient_neg.probes <- GetAssayData(subset(patient_data, features = row.names(GetAssayData(patient_data)) %>% grep("Negative", ., value = TRUE))) %>%
    as.matrix() %>%
    t()
  # Calculate the average negative probes per cell
  patient_avg.neg.probes <- Matrix::rowMeans(patient_neg.probes)
  
  # Semi-supervised learning with insitutype and reference profiles
  # InSituType needs integers, if given floating point numbers it fails with misleading errors
  patient_semisup <- insitutype(
    x = patient_rna.counts,
    neg = patient_avg.neg.probes,
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