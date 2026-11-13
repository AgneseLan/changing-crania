## ---- Setup ----

library(Rvcg)
library(SURGE)
library(rgl)
library(Morpho)
library(tidyverse)
library(geomorph)

source("D:/Ada_MSproject/R project/custom-functions.R")

load("D:/Ada_MSproject/R project/ply/slidedlms.cranium.R")

patch_coords <- read.table(
  "D:/Ada_MSproject/R project/patching/aml_2_patch.pts"
) |>
  dplyr::select(-1) |>
  as.matrix()

template_mesh <- vcgImport(
  "D:/Ada_MSproject/R project/patching/aml_2.ply"
)

curve_table <- readr::read_csv(
  "D:/Ada_MSproject/R project/Modified_curve_def - copia.csv",
  show_col_types = FALSE
)

my_curves <- create_curve_info(curve_table, n_fixed = 35)

template_specimen <- slid.lms$out[,,2]

points <- my_curves$Fixed

curves <- unlist(my_curves$Curve.in)

curvein <- my_curves$Curve.in

shade3d(template_mesh, col="white")

spheres3d(template_specimen[points,], radius=.1, col="red")

spheres3d(template_specimen[curves,], radius=.05, col="blue")

atlas <- createAtlas(
  mesh = template_mesh,
  landmarks = template_specimen,
  patch = patch_coords,
  corrCurves = curvein,
  patchCurves = NULL,
  keep.fix = points
)

plotAtlas(atlas)

save(atlas,file="D:/Ada_MSproject/R project/atlas.R")

load("D:/Ada_MSproject/R project/atlas.R")

plypath <- "E:/Portatil copiar/Work/Articulos/deformed skulls/landmarks/Landmarks/pts - copia/ply/"

patched_specimens <- placePatch(
  atlas = atlas,
  dat.array = slid.lms$out,
  inflate = 1,
  path = plypath,
  fileext = ".ply"
)

checkLM.mod(
  patched_specimens,
  path = plypath,
  pt.size = NULL,
  pt.size.scale = 0.001,
  suffix = ".ply",
  render = "s",
  begin = 1
)

patchedLMs <- patched_specimens

save(
  patchedLMs,
  file="E:/Portatil copiar/Work/Articulos/deformed skulls/landmarks/Landmarks/pts - copia/patchedLMs.Rdata"
)

########################
# PCA
########################

load("E:/Portatil copiar/Work/Articulos/deformed skulls/landmarks/Landmarks/pts - copia/patchedLMs.Rdata")

patchedLMs2 <- patchedLMs[,,-3]

Y.gpa <- geomorph::gpagen(patchedLMs)

cranium_pca <- geomorph::gm.prcomp(Y.gpa$coords)

PCscores <- cranium_pca$x

summary(cranium_pca)

plot(cranium_pca, axis1 = 1, axis2 = 2)

pcscores <- tibble::as_tibble(PCscores[,1:6], rownames="Filename")

pcscores_plot <- pcscores %>%
  dplyr::select(PC1, PC2, PC3, PC4, Filename)

full_specimen_data <- readr::read_csv(
  "E:/Portatil copiar/Work/Articulos/deformed skulls/landmarks/Landmarks/pts - copia/c-full-specimen-data.csv",
  show_col_types = FALSE
)

pcscores_plot <- dplyr::left_join(
  pcscores_plot,
  full_specimen_data,
  by = "Filename"
)

ggplot(
  pcscores_plot,
  aes(x = PC1, y = PC2, shape = Filename, colour = Type)
) +
  geom_point(size = 2.5) +
  scale_shape_manual(values=c(3,4,8,7,17,18,20,21,22,23,24,25)) +
  labs(
    x="Principal Component 1",
    y="Principal Component 2",
    title="Cranial Morphospace"
  ) +
  theme_bw()