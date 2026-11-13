#Code for reading in curves, resampling, and sliding
#Make sure to run the SLIDER3D_2 function first
#Code adapted from Ryan Felice

rm(list=ls())
library(tidyverse)
library(Morpho)
library(geomorph)
library(Rvcg)
library(paleomorph)
library(EMMLi)
library(qgraph)
library(ape)
library(geiger)
library(abind)
library("devtools")
library(SURGE)
library(RColorBrewer)
library(magick)
library(readr)
#devtools::install_github("rnfelice/SURGE")
library(SURGE)
library(rgl)

# set working directory
wd <- "D:/OneDrive/my research project/R analyses/changing-crania/Data"
setwd(wd)

#source(paste(wd, "Modified crania/custom-functions.R", sep = ""))

#import table defining curves
curve_table <- readr::read_csv("./Modified_curve_def - copia.csv", show_col_types = FALSE)

#identify the folder where your pts files are
ptsfolder <- paste(wd,"/pts", sep = "")

#import the pts file names
ptslist <- dir(ptsfolder, pattern='\\.pts$', recursive=FALSE)

my_curves <- create_curve_info(curve_table, n_fixed = 35)

setwd(ptsfolder)

subsampled.lm <- import_chkpt_data(
  ptslist,
  my_curves,
  subsampl = TRUE,
  verbose = TRUE
)

# #if you have any missing points, Checkpoint will set them to 9999
subsampled.lm[subsampled.lm == 9999] <- NA

# #correct rare infinite values
subsampled.lm[subsampled.lm == Inf] <- NA
subsampled.lm[subsampled.lm == -Inf] <- NA

spheres3d(subsampled.lm[,,2], r = 1) 
------------------------------------------- #####
#SET WORKING DIRECTORY TO ASCII PLY
setwd("D:/Ada_MSproject/R project/ply")

checkLM.mod(
  subsampled.lm,
  path="D:/Ada_MSproject/R project/ply",
  pt.size=NULL,
  pt.size.scale=0.1,
  suffix=".ply",
  render="s",
  begin=36
)

newpts <- subsampled.lm

#Create missing list
# misslist <- createMissingList(dim(newpts)[3])
# 
# for (j in seq_len(dim(newpts)[3])) {
#   misslist[[j]] <- which(is.na(newpts[,1,j]))
# }
# 
# newpts2 <- fixLMtps(newpts)

#check that pts and ply names match
plyfolder <- "D:/Ada_MSproject/R project/ply"



ptslist2 <- gsub("\\.ply.pts$", "", ptslist)

plylist <- dir(plyfolder, pattern='\\.ply$', recursive=FALSE)
plylist2 <- gsub("\\.ply$", "", plylist)

setdiff(plylist2, ptslist2)

#newpts2 <- newpts2$out

spheres3d(newpts[,,2], r = 1)

#################################################
#        RUN THE SLIDER 3D_2 CODE
#################################################

###!!SET WD to ply from console!! -->

dimnames(newpts$out)[[3]]<-gsub(".ply”,"",dimnames(newpts$out)[[3]])

#Slide points on surface
slided4.all <- slider3d_2(newpts, SMvector= my_curves$Sliding.LMs,
                          outlines = my_curves$Curves, 
                          #copy ply folder path to ensure it works - set as working directory from console to print
                          sur.path = "D:/Ada_MSproject/R project/ply",
                          sur.name = NULL, 
                          meshlist = paste("",dimnames(newpts)[[3]],"",sep=""), ignore = NULL,
                          sur.type = "ply", tol = 1e-10, deselect = FALSE, inc.check = FALSE,
                          recursive = TRUE, iterations = 3, initproc = TRUE,
                          pairedLM = 0, mc.cores = 1, bending=TRUE,
                          fixRepro = F, stepsize=0.2)

# try slider3d instead? And it worked...
# Iteration 1: squared distance between means: 1.4567704158347e-06 
# Iteration 2: squared distance between means: 1.12019231726991e-06 
# Iteration 3: squared distance between means: 6.81562338380338e-07 
# Meaning: 

setwd("D:/Ada_MSproject/R project/")

slided4.all <- slider3d_2(
  newpts,
  SMvector = my_curves$Sliding.LMs,
  outlines = my_curves$Curves,
  sur.path = "./ply",
  sur.name = NULL,
  meshlist = dir(plyfolder, pattern='\\.ply$', recursive=FALSE),
  ignore = NULL,
  sur.type = "ply",
  tol = 1e-10,
  deselect = FALSE,
  inc.check = FALSE,
  recursive = TRUE,
  iterations = 3,
  initproc = TRUE,
  pairedLM = 0,
  mc.cores = 1,
  bending = TRUE,
  fixRepro = FALSE,
  stepsize = 0.2
)

dimnames(slided4.all[["dataslide"]])[3] <- dimnames(newpts2)[3]

#re-estimate missing post sliding
slided4.all$dataslide[which(is.na(newpts))] <- NA

#Fix missing landmarks
slid.lms <- fixLMtps(slided4.all$dataslide)

slid.lms.neurocranium <- slided4.all

save(slid.lms, file='slidedlms.cranium.R')

# export landmark data
write.csv(
  slid.lms$out,
  "D:/Ada_MSproject/R project/resampled__landmarks.csv",
  row.names = FALSE
)

save(slid.lms, file="./resampled__landmarks.R")

#---------------#
# APPLY PATCHES #
#---------------#
