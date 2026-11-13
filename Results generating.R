# ==============================================================================
# CLUSTERING STEP 1: THE ELBOW METHOD FOR OPTIMAL K
# ==============================================================================

# Isolate just the first four master shape components from your data frame
pc_matrix <- pcscores_all_df[, c("Comp1", "Comp2", "Comp3", "Comp4")]

# Calculate the total within-groups sum of squares for 1 to 5 clusters
wss <- (nrow(pc_matrix) - 1) * sum(apply(pc_matrix, 2, var))
for (i in 2:5) {
  wss[i] <- sum(kmeans(pc_matrix, centers = i, nstart = 25)$withinss)
}

# Plot the Elbow Curve
plot(1:5, wss, type = "b", xlab = "Number of Clusters (K)",
     ylab = "Within Groups Sum of Squares",
     main = "Elbow Method for Optimal Cluster Selection",
     pch = 19, col = "darkblue")

# ==============================================================================
# CLUSTERING STEP 2: RUN K-MEANS ALGORITHM
# ==============================================================================

# Run K-means with K=2 clusters
set.seed(123) # Locks randomness
kmeans_results <- kmeans(pc_matrix, centers = 2, nstart = 25)

# Save the mathematical partition cluster right into your data frame
pcscores_all_df$KMeans_Cluster <- as.factor(kmeans_results$cluster)

# Look at the results! See which specimens fell into Cluster 1 vs Cluster 2
print(pcscores_all_df[, c("specimens", "category", "group", "KMeans_Cluster")])

# ==============================================================================
# REBUILDING THE DATA FRAME WITH SCORES AND METADATA
# ==============================================================================
library(dplyr)

# 1. Start fresh by turning your raw PCA scores matrix into a data frame
pcscores_all_df <- as.data.frame(cranium_pca_relaxed$x)

# 2. Add your metadata columns back in using your actual metadata object
# (Make sure to change 'my_metadata_df' to the exact name of your spreadsheet object!)
pcscores_all_df <- pcscores_all_df %>% 
  mutate(
    specimens = rownames(pcscores_all_df),          # Pulls your "aml_1", "aml_2" labels
    group     = gdf$group,                          # Pulls your island populations (Cuba vs Ven)
    category  = gdf$category,                       # Pulls modification status
    size      = gdf$Csize                           # Pulls centroid sizes
  )

# 3. NOW run the K-means clustering algorithm on just the numeric components
pc_matrix <- pcscores_all_df[, c("Comp1", "Comp2", "Comp3", "Comp4")]
set.seed(123)
kmeans_results <- kmeans(pc_matrix, centers = 2, nstart = 25)

# 4. Save the new cluster assignment as "KMeans_Cluster"
pcscores_all_df$KMeans_Cluster <- as.factor(kmeans_results$cluster)

# ==============================================================================
# VERIFY THE RESULTS
# ==============================================================================
# Now this print line will work perfectly because all columns exist!
print(pcscores_all_df[, c("specimens", "category", "group", "KMeans_Cluster")])


#######
# ==============================================================================
# RUNNING K-MEANS DIRECTLY ON YOUR PC SCORES
# ==============================================================================

# 1. Grab just the primary shape axes (PC1 to PC4) directly from your PCA object
pc_matrix <- cranium_pca_relaxed$x[, 1:4]

# 2. Run the K-means algorithm (Setting K=2 to test Modified vs Unmodified separation)
set.seed(123) # Locks randomness
kmeans_results <- kmeans(pc_matrix, centers = 2, nstart = 25)

# 3. Build a completely fresh, simple summary data frame to view the results
clustering_summary <- data.frame(
  Specimen = rownames(pc_matrix),               # Pulls your "aml_" names cleanly
  Cluster  = as.factor(kmeans_results$cluster)  # The mathematical group assigned by R
)

# ==============================================================================
# PRINT YOUR RESULTS
# ==============================================================================
print(clustering_summary)

#####Generate the Clustered Morphospace Plot#####
library(ggplot2)
library(ggrepel)

# 1. Combine your PCA scores, your metadata, and the new clusters into one plotting dataframe
plot_df <- data.frame(
  PC1 = cranium_pca_relaxed$x[, 1],
  PC2 = cranium_pca_relaxed$x[, 2],
  Specimen = rownames(cranium_pca_relaxed$x),
  Cluster = clustering_summary$Cluster
)

# 2. Build the plot
ggplot(plot_df, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 4, alpha = 0.8) +
  # Draw automatic statistical boundaries around the two clusters
  stat_ellipse(aes(fill = Cluster), geom = "polygon", alpha = 0.1, level = 0.95) +
  # Label each dot with its "aml_" ID so you can see who is who
  geom_text_repel(aes(label = Specimen), size = 3.5, color = "black", max.overlaps = 20) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(
    title = "Unsupervised K-Means Clustering on Cranial Morphospace",
    subtitle = "Partitions calculated purely from PC1-PC4 shape scores",
    x = "Principal Component 1",
    y = "Principal Component 2",
    color = "Calculated Cluster",
    fill = "Calculated Cluster"
  ) +
  theme(plot.title = element_text(face = "bold"))

# 3. Save the high-resolution image directly to your folder
ggsave("Output/5-PCA/KMeans_Cranial_Clusters.png", width = 8, height = 6, dpi = 300)

# ==============================================================================
# Revised Cranial Morphospace Configuration
# ==============================================================================

library(ggplot2)
library(ggrepel)

# 1. Reset the scores data frame fresh from the relaxed PCA object
pcscores_all_df <- as.data.frame(cranium_pca_relaxed$x)

# 2. Assign the specimen IDs directly from the row names ("aml_1", "aml_2", etc.)
pcscores_all_df$specimens <- rownames(pcscores_all_df)

# 3. Pull directly from your raw spreadsheet (REPLACE WITH YOUR EXACT EXCEL COLUMN HEADERS!)
# If your spreadsheet columns are named 'Island' and 'Status', change these accordingly:
pcscores_all_df$group     <- my_metadata_df$YOUR_EXCEL_ISLAND_COLUMN     # e.g., Cuba vs Venezuela
pcscores_all_df$category  <- my_metadata_df$YOUR_EXCEL_MODIFICATION_COLUMN # e.g., Modified vs Unmodified

# 4. Grab Centroid Size straight from the PCA output
pcscores_all_df$size      <- cranium_pca_relaxed$Csize

# 5. Calculate your exact axis variance percentages dynamically
pc1_percent <- round((cranium_pca_relaxed$sdev[1]^2 / sum(cranium_pca_relaxed$sdev^2)) * 100, digits = 2)
pc2_percent <- round((cranium_pca_relaxed$sdev[2]^2 / sum(cranium_pca_relaxed$sdev^2)) * 100, digits = 2)

# 6. Generate the publication plot
PCA_human_ggplot <- ggplot(pcscores_all_df, aes(x = Comp1, y = Comp2, label = specimens)) +
  # Color points by population group and shape them by modification status
  geom_point(aes(color = group, shape = category), size = 4, alpha = 0.8) +
  # Prevent labels from stacking
  geom_text_repel(colour = "black", size = 3, max.overlaps = 40) +
  scale_color_brewer(name = "Population Series", palette = "Set1") +
  theme_bw() + 
  xlab(paste0("PC 1 (", pc1_percent, "%) - Deformation Axis")) +
  ylab(paste0("PC 2 (", pc2_percent, "%) - Population Axis")) +
  theme(
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 12)
  ) +
  labs(title = "Cranial Morphospace Configuration")

# 7. Print the plot to your screen
PCA_human_ggplot
#####

# ==============================================================================
# PROCRUSTES ANOVA: TESTING GEOGRAPHY AND MODIFICATION EFFECT
# ==============================================================================
library(geomorph)

# 1. Cleanly map the variables into a fresh, isolated testing list
anova_data_list <- list(
  ShapeCoords = Y.gpa_relaxed$coords,               # Your 3D relaxed landmarks
  Population   = metadata_subset$Geographic_origin,     # e.g., Cuba vs Venezuela
  Modification = metadata_subset$Type # e.g., Modified vs Unmodified
)

# 2. Run the Procrustes ANOVA using the exact names in your isolated list
# We use 9,999 permutations for highly robust p-values for Gizeh and Agnese
shape_anova <- procD.lm(
  ShapeCoords ~ Population * Modification, 
  data = anova_data_list, 
  iter = 9999, 
  RRPP = TRUE
)

# ==============================================================================
# DISPLAY THE TABLE FOR YOUR POWERPOINT SLIDE
# ==============================================================================
summary(shape_anova)

# ==============================================================================
# AUTOMATING 40 ANATOMICAL SNAPSHOTS
# ==============================================================================
library(rgl)

# Define your standard anatomical orientations (Theta, Phi, Fov)
orientations <- list(
  frontal = c(0, 0, 30),
  lateral_left = c(90, 0, 30),
  lateral_right = c(-90, 0, 30),
  upper   = c(0, 90, 30),
  bottom  = c(0, -90, 30),
  back    = c(180, 0, 30)
)

# Example for PC1 min (Repeat this structure for your PC shapes)
# 1. Plot your warped mesh using your standard script line:
# plotRefToTarget(cranium_pca2$shapes$shapes.comp1$min, mean_shape, mesh = template_mesh)
# ==============================================================================
# RETRIEVING OR RECALCULATING THE MEAN SHAPE
# ==============================================================================

# Option A: Pull the consensus coordinates directly out of your PCA object
mean_shape <- cranium_pca2$consensus

# Option B: If Option A is empty, recalculate it directly from your aligned array
# (Uncomment the line below if option A doesn't solve it)
# mean_shape <- mshape(Y.gpa_relaxed$coords)

# ==============================================================================
# RE-RUN TASK 1 RE-ORIENTED PLOT
# ==============================================================================
# This will now find 'mean_shape' perfectly!
plotRefToTarget(cranium_pca2$shapes$shapes.comp1$min, mean_shape, mesh = template_mesh)

# 2. Run this loop to automatically save every view perfectly straight
for(view_name in names(orientations)) {
  # Snap camera to angle
  view3d(theta = orientations[[view_name]][1], 
         phi   = orientations[[view_name]][2], 
         fov   = orientations[[view_name]][3])
  
  # Save snapshot to your Output folder
  rgl.snapshot(paste0("Output/5-PCA/PC1min_relaxed_", view_name, ".png"))
}
clear3d()

