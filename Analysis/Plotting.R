###### COLORS ######

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
  "T cell" = "#1C8356",
  "Unknown" = "#f5ee59"
)

# Create a vector with a color for every cluster
# It can be used to create a palette which includes known clusters (e.g. B cell, T cell, etc.) and unknown clusters (e.g. 0, 1, ..., a, b, etc.)
generate_colors_lookup_table <- function(
    data,
    cluster_column_name,
    known_clusters_colors = NULL,
    color_palette = DiscretePalette(36, palette = "polychrome")) {
  
  # Remove the colors that are used by the known clusters
  usable_color_palette <- setdiff(color_palette, known_clusters_colors)
  
  # Random colors for the unknown clusters
  # The cluster 0 in Seurat for example is different for every patient
  usable_color_palette <- sample(usable_color_palette)
  
  # "data" must be either a data frame or a tibble
  # If the object is of type Seurat, extract the meta data
  if ("Seurat" %in% class(data)) {
    data <- data@meta.data
  }
  
  # Get the unique clusters in the data (works with both base R data frames and tibbles)
  clusters <- unique(dplyr::pull(data, cluster_column_name))
  
  # Remove the known clusters from the list of clusters as they already have a color assigned
  unknown_clusters <- setdiff(clusters, names(known_clusters_colors))
  
  # Assign a random colors to the unknown clusters
  unknown_clusters_colors <- setNames(usable_color_palette[1:length(unknown_clusters)], unknown_clusters)
  
  # Return the colors for the known and unknown clusters
  return(c(known_clusters_colors, unknown_clusters_colors))
}

###### GENERATE SINGLE PLOTS ######

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

generate_elbow_plot <- function(patient_data, reduction, dims) {
  
  print(paste("Generate ElbowPlot from", reduction))
  
  patient_num <- get_patient_num(patient_data)
  
  elbow_plot <- ElbowPlot(patient_data, reduction = reduction, ndims = dims) +
    labs(title = paste("Patient", patient_num), subtitle = reduction)
  
  elbow_plot_name <- paste("Patient",  patient_num, "elbow_plot", reduction, sep = "_")
  return(setNames(list(elbow_plot), elbow_plot_name))
}

generate_umap <- function(patient_data, cluster_var, cluster_reduction, cluster_name = NULL, color_table = NULL) {
  
  if (is.null(cluster_name)) {
    cluster_name <- cluster_var
  }
  
  patient_num <- get_patient_num(patient_data)
  
  # Graphs the output of a dimensional reduction technique on a 2D scatter plot
  # Each point is a cell and it's positioned based on the cell embeddings determined by the reduction technique
  umap_clusters <- Seurat::DimPlot(
    object = patient_data,
    reduction = cluster_reduction,,
    group.by = cluster_var,
    label=TRUE,
    label.box=TRUE,
    repel=TRUE,
    cols = color_table) +
    labs(
      title = paste("Patient", patient_num),
      subtitle = cluster_name) +
    NoLegend()
  
  # Return the plot inside a list with a name
  umap_clusters_plot_name <- paste("Patient",  patient_num, cluster_var, "umap", sep = "_")
  return(setNames(list(umap_clusters), umap_clusters_plot_name))
}