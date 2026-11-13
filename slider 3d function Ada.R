
# Process for reading in landmarks and resampling curves

# 1.      Read in curves
# 2.      Resample - using the curve.csv (separate code for archs/mysts/odonts)
# 3.      Make a list of the missing data (using the createMissinglist function)
# 4.      Use fixLMtps to produce dummy variables
# 5.      SLIDE
# 6.      Then fixLMtps again to better estimate dummy variables
# 7.      Add absent bone landmarks for a couple of specimens
# 8.      Make a fake midline with the rostrum landmarks - different for odonts/archs/mysts
# 9.      MIRROR - different code for archs, mysts, and odonts
# 10.    Remove fake midline landmarks + additional midlines that aren't needed
# 11.    Analyse (edited) 

#This makes nice 3D outputs of the files for publications as well (such as the M.Mirus)
#script for converting binary to ASCII .ply files
#rm(list=ls())

library(Rvcg)
library(rgl)
#set working directory with ascii or binary ply files
#put the .ply files you want to change (for reading into R) in the .ply binary folder (not the .ply folder) 
wd = ("D:/Ada_MSproject/R project/ply - Copy")
setwd(wd)
#set output folder
#The extra forward slash at the end tells R that that is the output folder
outputfolder<-("D:/Ada_MSproject/R project/ply ASCII")
meshlist<-dir(pattern='.ply',recursive=F)
# meshlist<-dir(path = paste("/meshes/neurocranium"), pattern = '.ply', recursive = F)

for (i in 1:length(meshlist)){
  x<-vcgImport(file=meshlist[i]) #import meshes
  
  vcgPlyWrite(x,filename=paste(outputfolder,meshlist[i],sep=""),binary=FALSE)#export
}


#SLIDER3D_2 function 

calcGamma <- function(Gamma0,Lsubk3,U,dims,stepsize=1)
{
  
  U <- as(U,"sparseMatrix")
  tUL <- crossprod(U,Lsubk3)
  ULU <- forceSymmetric(tUL%*%U)
  B <- tUL%*%Gamma0
  B <- as(B,"sparseMatrix")
  T <- solve(ULU,B)
  ULUT <- U%*%T
  Gamma0 <- Gamma0-stepsize*ULUT
  Gamma0 <- matrix(Gamma0,length(Gamma0)/dims,dims)
  return(Gamma0)
}
calcProcDGamma <- function(U,Gamma0,mshape,dims,stepsize=1) {
  Tpart <- tcrossprod(U)
  mshape <- as.vector(mshape)
  tmpdiff <- Gamma0-mshape
  slided <- Gamma0-stepsize*(Tpart%*%tmpdiff)
  slided <- matrix(slided,length(slided)/dims,dims)
  return(slided)
}
calcTang_U_s2 <- function(datamatrix,normalmatrix=NULL,SMvector,outlines=NULL,surface=NULL,free=NULL,deselect=FALSE)
{
  
  dims <- dim(datamatrix)[2]
  k <- dim(datamatrix)[1]
  
  if (deselect==TRUE)
    SMvector <- c(1:k)[-SMvector]
  
  SMvector <- unique(SMvector)
  m <- length(SMvector)
  if ( !is.null(free)) {
    udims <- c(dims*k,m*3)
    tanvec <- matrix(0,k,dims*3)
    #U <- matrix(0,dims*k,m*3)
    type <- 2
  } else if(!is.null(surface)) {
    udims <- c(dims*k,m*2)
    tanvec <- matrix(0,k,dims*2)
    #U <- matrix(0,dims*k,m*2)
    type <- 1
  } else {
    udims <- c(dims*k,m)
    tanvec <- matrix(0,k,dims)
    #U <- matrix(0,dims*k,m)
    type <- 0
  }
  Gamma0 <- c(datamatrix)
  
  if (is.null(outlines) == FALSE) {             
    if (is.list(outlines)==FALSE) {
      outlines <- list(outlines)
    }
    for ( j in 1:length(outlines)) {
      lt <- length(outlines[[j]])
      temp <- outlines[[j]]
      
      ### procedure for open curves ####            
      if (outlines[[j]][1]!= outlines[[j]][lt]) {
        for (i in 1:lt) {
          if (temp[i]%in%SMvector==TRUE && i!=1 && i!=lt) {
            tanvec[temp[i],1:3] <- (datamatrix[temp[i-1],]-datamatrix[temp[i+1],])/sqrt(sum((datamatrix[temp[i-1],]-datamatrix[temp[i+1],])^2))
          } else if (temp[i]%in%SMvector==TRUE && i==1) {
            tanvec[temp[i],1:3] <- (datamatrix[temp[i],]-datamatrix[temp[i+1],])/sqrt(sum((datamatrix[temp[i],]-datamatrix[temp[i+1],])^2))
          } else if (temp[i]%in%SMvector==TRUE && i==lt) {
            tanvec[temp[i],1:3] <- (datamatrix[temp[i-1],]-datamatrix[temp[i],])/sqrt(sum((datamatrix[temp[i-1],]-datamatrix[temp[i],])^2))
          }
        } 
      }
      ### procedure for closed curves ####
      else if (outlines[[j]][1]== outlines[[j]][lt]) {
        for (i in 1:(lt-1)) {
          if (temp[i]%in%SMvector==TRUE && i!=1 && i!=lt) {
            tanvec[temp[i],1:3] <- (datamatrix[temp[i-1],]-datamatrix[temp[i+1],])/sqrt(sum((datamatrix[temp[i-1],]-datamatrix[temp[i+1],])^2))
          } else if (temp[i]%in%SMvector==TRUE && i==1) {
            tanvec[temp[i],1:3] <- (datamatrix[temp[lt-1],]-datamatrix[temp[i+1],])/sqrt(sum((datamatrix[temp[lt-1],]-datamatrix[temp[i+1],])^2))
          }
        } 
      }
    }
  }
  ### procedure for surfaces ###
  if (is.null (surface) ==F) {  
    lt <- length(surface)
    temp <- surface
    for (i in 1:lt) {
      tanp <- tangentPlane(normalmatrix[temp[i],])
      tanvec[temp[i],1:6] <- c(tanp$y,tanp$z)                
    }
  }
  #### end surfaces ####
  
  ### procedure for free sliding points ####
  
  if (!is.null(free)) {
    lf <- length(free)
    tmp <- free
    for (i in 1:lf)
      tanvec[tmp[i],] <- c(1,0,0,0,1,0,0,0,1)
  }
  ### end free sliding ##
  gc()  
  SMsort <- sort(SMvector)
  xinfo <- .Call("tweakU",tanvec,m, type,SMsort)
  U <- sparseMatrix(i=xinfo$rows,j=xinfo$cols+1, x=xinfo$x,dims=udims)
  #U <- weights*U
  
  
  outOnly <- outlines
  surfOnly <- surface
  freeOnly <- free
  if (!is.null(surface))
    surfOnly <- unique(sort(surface))
  if (!is.null(free))
    freeOnly <- unique(sort(free))
  if(!is.null(outlines))
    outOnly <- unique(sort(unlist(outlines)))
  
  allsurf <- c(outOnly,surfOnly,freeOnly)
  
  if (length(which(!SMvector %in% allsurf)))
    stop("all semi-landmarks must to be tagged as outlines or surfaces")
  ## remove fix columns
  
  Ured0 <- as(U[,1:m],"sparseMatrix")
  
  if (!is.null(surface) || !is.null(free)) {
    Ured1 <- as(U[,(m+1):(2*m)],"sparseMatrix")
    smsurffree <- which(! SMvector %in% c(surfOnly,freeOnly))
    if (length(smsurffree))
      Ured1 <- Ured1[,-smsurffree]
    Ured0 <- cbind(Ured0,Ured1)
    
  }
  if (!is.null(free)) {
    Ured1 <- as(U[,(2*m+1):(3*m)],"sparseMatrix")
    smsurffree <- which(! SMvector %in% c(freeOnly))
    if (length(smsurffree))
      Ured1 <- Ured1[,-smsurffree]
    Ured0 <- cbind(Ured0,Ured1)
    
  }
  
  return(list(SMvector=SMvector,Gamma0=Gamma0,U=Ured0))             
}     
#dat.array<-newpts2$out[,,c(260:270)]
#SMvector= my_curves$Sliding.LMs
#
#outlines = my_curves$Curves
#sur.path = "./Raw_Data/ply"
#sur.name = NULL
#meshlist = paste("./Raw_Data/ply/",dimnames(newpts2$out[,,c(260:270)])[[3]],".ply",sep="")
#ignore = NULL
#
#sur.type = "ply"
#tol = 1e-10
#deselect = FALSE
#inc.check = FALSE
#
#recursive = TRUE
#iterations = 3
#initproc = TRUE
#pairedLM = 0
#mc.cores = 1
#bending=TRUE
#fixRepro = FALSE
#stepsize=0.2
#
#missingList=misslist[c(260:270)]
slider3d_2<-function (dat.array, SMvector, outlines = NULL, surp = NULL, 
                      sur.path = NULL, sur.name = NULL, meshlist = NULL, ignore = NULL, 
                      sur.type = "ply", tol = 1e-05, deselect = FALSE, inc.check = TRUE, 
                      recursive = TRUE, iterations = 0, initproc = TRUE, fullGPA = FALSE, 
                      pairedLM = 0, bending = TRUE, stepsize = ifelse(bending, 
                                                                      1, 0.5), mc.cores = parallel::detectCores(), fixRepro = TRUE, 
                      missingList = NULL, use.lm = NULL, silent = FALSE) 
{
  require(Matrix)
  require(foreach)
  require(doParallel)
  
  
  if (.Platform$OS.type == "windows" && mc.cores > 1) {
    cl <- makeCluster(mc.cores)            
    doParallel::registerDoParallel(cl=cl)
  } else if (mc.cores > 1) {
    doParallel::registerDoParallel(cores = mc.cores)
  } else
    foreach::registerDoSEQ()
  
  if (iterations == 0)
    iterations <- 1e10
  
  if (is.null(outlines) && is.null(surp))   
    stop("nothing to slide")
  fixLM <- integer(0)
  n <- dim(dat.array)[3]
  k <- dim(dat.array)[1]
  m <- dim(dat.array)[2]
  nomesh <- FALSE
  if (length(sur.name) && is.null(sur.path))
    sur.path <- ""
  if (is.null(meshlist) && is.null(sur.path))
    nomesh <- TRUE
  if (pairedLM[1]!=0 && is.vector(pairedLM))# check if there are only 2 symmetric lms
    pairedLM <- t(as.matrix(pairedLM))
  if(!is.null(missingList))
    if(length(missingList) != n)
      stop(paste0("missingList must be of length", n," - same as samplesize"))
  ### update indexing for after ignored landmarks are removed ###   
  if (!is.null(ignore)) {
    li <- length(ignore)
    lm.old <- c(1:k)[-ignore]
    mat.ptr <- matrix(c(1:(k-li),lm.old),k-li,2)
    ptr <- function(xo)  ### define pointer function for indexing
    {
      if (length(which(ignore %in% xo))!= 0)
        xo <- xo[-which(xo %in% ignore)]
      for (i in 1:(k-li))
        xo[which(xo==mat.ptr[i,2])] <- mat.ptr[i,1]
      return(xo)
    }
    if (!is.null(missingList))
      missingList <- lapply(missingList,ptr)
    if (!is.null(outlines)) ### update outline indices
      outlines <- lapply(outlines,ptr)
    if (!is.null(surp))     ### update surface indices
      surp <- ptr(surp)
    
    if (!is.null(SMvector)) ### of fixed/sliding definition
      SMvector <- ptr(SMvector)
    
    if (pairedLM[1]!=0){    ### update paired landmarks indices
      count <- 0
      del <- NULL
      for (i in 1:dim(pairedLM)[1]) {   
        if (length(which(ignore %in% pairedLM[i,]))!=0) {
          count <- count+1
          del[count] <- i
        }
      }
      pairedLM <- pairedLM[-del,]
      if (is.vector(pairedLM))
        pairedLM <- t(as.matrix(pairedLM))
      
      if (dim(pairedLM)[1]==0) {
        pairedLM <- 0
      } else {
        pairedLM <- apply(pairedLM,2,ptr)
        if (is.vector(pairedLM))
          pairedLM <- t(as.matrix(pairedLM))
      }
    }
    dat.array <- dat.array[-ignore,,]
    k <- dim(dat.array)[1]
  }
  
  vn.array <- dat.array
  data.orig <- dat.array
  if (deselect)
    fixLM <- SMvector
  else if (length(SMvector) < k)
    fixLM <- c(1:k)[-SMvector]
  else
    fixRepro <- TRUE
  
  weights <- NULL
  if (!is.null(use.lm)) {
    weights <- rep(0,dim(dat.array)[1])
    weights[use.lm] <- 1
  }
  if(length(sur.name)==0 && !is.null(sur.path)) {
    sur.name <- dimnames(dat.array)[[3]]
    sur.name <- paste(sur.path,"/",sur.name,".",sur.type,sep="")
    
  }
  p1 <- 10^12
  
  ## check for existing surfaces (if needed)
  if (length(sur.name)) {
    if (length(sur.name) != dim(dat.array)[3])
      stop("length of sur.name does not match number of landmark configurations")
    checkfiles <- file.exists(sur.name)
    if (prod(checkfiles) == 0) {
      warning(paste0("missing mesh files: ", sur.name[checkfiles == 0]))
      stop("Some mesh files do not exist, please check your data and/or correct array naming")
    }
  }
  if (length(meshlist)) {
    if (length(meshlist) != dim(dat.array)[3])
      stop("length of meshlist does not match number of landmark configurations")
  }
  
  ini <- rotonto(dat.array[,,1],dat.array[,,2],signref=FALSE) # create mean between first tow configs to avoid singular BE Matrix
  mshape <- (ini$Y+ini$X)/2
  if (!silent && !nomesh)
    cat(paste("Points will be initially projected onto surfaces","\n","-------------------------------------------","\n"))
  ## parallel function in case meshlist != NULL
  parfunmeshlist <- function(i,data) {
    if (!is.list(data))
      tmpdata <- data[,,i]
    else
      tmpdata <- data[[i]]
    out <- projRead(tmpdata,meshlist[[i]])
    if (!is.null(missingList))
      if(length(missingList[[i]]))
        out$vb[1:3,missingList[[i]]] <- t(tmpdata[missingList[[i]],])
    return(out)
  }
  parfunmeshfile <- function(i, data) {
    if (!is.list(data))
      tmpdata <- data[,,i]
    else
      tmpdata <- data[[i]]
    
    out <- projRead(tmpdata,sur.name[i])
    if (!is.null(missingList))
      if(length(missingList[[i]]))
        out$vb[1:3,missingList[[i]]] <- t(tmpdata[missingList[[i]],])
    
    return(out)
    
  }
  parfunnomesh <- function(i, data) {
    if (!is.list(data))
      tmpdata <- data[,,i]
    else
      tmpdata <- data[[i]]
    
    out <- vcgUpdateNormals(tmpdata,silent=TRUE)
    return(out)
  }
  
  if (is.null(meshlist) && !nomesh) {
    repro <- foreach(i=1:n, .inorder=TRUE,.errorhandling="pass",.packages=c("Morpho","Rvcg")) %dopar% parfunmeshfile(i,dat.array)
  } else if (!nomesh) {
    repro <- foreach(i=1:n, .inorder=TRUE,.errorhandling="pass",.packages=c("Morpho","Rvcg")) %dopar% parfunmeshlist(i,dat.array)
    
  } else {
    repro <- foreach(i=1:n, .inorder=TRUE,.errorhandling="pass",.packages=c("Morpho","Rvcg")) %dopar% parfunnomesh(i,dat.array)
    message("no surfaces specified - surface is approximated from point cloud")
  }    
  for (ii in 1:n){
    if(class(repro[[ii]]) !=  "mesh3d")
      writeLines(paste0("ERROR: Specimen ", dimnames(dat.array)[[3]][[ii]], " failed projection of points on to surface."))
  }
  for (j in 1:n) {
    reprotmp <- repro[[j]]         
    dat.array[,,j] <- t(reprotmp$vb[1:3,])
    vn.array[,,j] <- t(reprotmp$normals[1:3,])
    
  }
  
  
  if (!fixRepro)# use original positions for fix landmarks
    dat.array[fixLM,,] <- data.orig[fixLM,,]
  
  if (!silent && !nomesh)
    cat(paste("\n","-------------------------------------------","\n"),"Projection finished","\n","-------------------------------------------","\n")
  
  if (initproc==TRUE) { # perform proc fit before sliding
    if (!silent)
      cat("Inital procrustes fit ...")    
    procini <- ProcGPA(dat.array,scale=fullGPA,silent=silent)
    mshape <- procini$mshape
  }
  dataslide <- dat.array
  
  if (pairedLM[1]!=0) {# create symmetric mean to get rid of assymetry along outlines/surfaces after first relaxation
    Mir <- diag(c(-1,1,1))
    A <- mshape
    Amir <- mshape%*%Mir
    Amir[c(pairedLM),] <- Amir[c(pairedLM[,2:1]),]
    symproc <- rotonto(A,Amir)
    mshape <- (symproc$X+symproc$Y)/2
  }
  if (!silent)
    cat(paste("Start sliding...","\n","-------------------------------------------","\n"))
  gc(verbose=F)
  ## calculation for a defined max. number of iterations
  count <- 1
  while (p1>tol && count <= iterations) {
    dataslide_old <- dataslide
    mshape_old <- mshape
    if (!silent)
      cat(paste("Iteration",count,sep=" "),"..\n")  # reports which Iteration is calculated
    if (recursive==TRUE)    # slided Semilandmarks are used in next iteration step
      dat.array <- dataslide
    if (bending)
      L <- CreateL(mshape,output="Lsubk3")
    else
      fixRepro=TRUE
    a.list <- as.list(1:n)
    slido <- function(j)                 
    {
      free <- NULL
      if (!is.null(missingList))
        if(length(missingList[[j]])) {
          if (!is.null(outlines)) {
            notoutlines <- which(!missingList[[j]] %in% unlist(outlines))
            if (length(notoutlines))
              free <- missingList[[j]][notoutlines]
          } else {
            free <- missingList[[j]]
          }
        }
      tmpdata <- dat.array[,,j]
      tmpvn <- vn.array[,,j]
      if (!bending) {
        rot <- rotonto(mshape,tmpdata,reflection=FALSE,scale=TRUE,weights=weights,centerweight=TRUE)
        tmpdata <- rot$yrot
        tmpvn <- tmpvn%*%rot$gamm
      }
      U <- calcTang_U_s2(tmpdata,tmpvn,SMvector=SMvector,outlines=outlines,surface=surp,deselect=deselect,free=free)
      if (bending) {
        dataslido <- calcGamma(U$Gamma0,L$Lsubk3,U$U,dims=m,stepsize=stepsize)
      } else {
        dataslido <- calcProcDGamma(U$U,U$Gamma0,mshape,dims=m,stepsize=stepsize)
        dataslido <- rotreverse(dataslido,rot)
      }
      return(dataslido)
    }
    
    a.list <- foreach(i=1:n, .inorder=TRUE,.errorhandling="pass",.export=c("calcGamma","calcTang_U_s2"),.packages=c("Morpho","Rvcg")) %dopar% slido(i)
    
    for (ii in 1:n){
      if(!any(class(a.list[[ii]]) ==  c("matrix")))
        writeLines(paste0("ERROR: sliding failed for specimen ", dimnames(dat.array)[[3]][[ii]]))
    }
    ###projection onto surface
    
    if (is.null(meshlist) && !nomesh) {
      repro <- foreach(i=1:n, .inorder=TRUE,.errorhandling="pass",.packages=c("Morpho","Rvcg")) %dopar% parfunmeshfile(i,a.list)  
    } else if (!nomesh) {
      repro <- foreach(i=1:n, .inorder=TRUE,.errorhandling="pass",.packages=c("Morpho","Rvcg")) %dopar% parfunmeshlist(i,a.list)
    } else {
      repro <- foreach(i=1:n, .inorder=TRUE,.errorhandling="pass",.packages=c("Morpho","Rvcg")) %dopar% parfunnomesh(i,a.list)
    }
    
    for (ii in 1:n){
      if(class(repro[[ii]]) !=  "mesh3d")
        writeLines(paste0("ERROR: Specimen ", dimnames(dat.array)[[3]][[ii]], " failed projection of points on to surface."))
    }
    
    for (j in 1:n) {
      reprotmp <- repro[[j]]         
      dataslide[,,j] <- t(reprotmp$vb[1:3,])
      vn.array[,,j] <- t(reprotmp$normals[1:3,])
    }
    
    if (!fixRepro)# use original positions for fix landmarks
      dataslide[fixLM,,] <- data.orig[fixLM,,]
    
    if (!silent)
      cat("estimating sample mean shape...")              
    proc <- ProcGPA(dataslide,scale=fullGPA,silent=silent)
    mshape <- proc$mshape
    if (pairedLM[1]!=0) {# create symmetric mean to get rid of assymetry along outline after first relaxation
      Mir <- diag(c(-1,1,1))
      A <- mshape
      Amir <- mshape%*%Mir
      Amir[c(pairedLM),] <- Amir[c(pairedLM[,2:1]),]
      symproc <- rotonto(A,Amir)
      mshape <- (symproc$X+symproc$Y)/2
    }     
    p1_old <- p1
    testproc <- rotonto(mshape_old,mshape)               
    p1 <- sum(diag(crossprod((testproc$X/cSize(testproc$X))-(testproc$Y/cSize(testproc$Y)))))
    
    ### check for increasing convergence criterion ###      
    if (inc.check) {
      if (p1 > p1_old) {
        dataslide <- dataslide_old
        if (!silent)
          cat(paste("Distance between means starts increasing: value is ",p1, ".\n Result from last iteration step will be used. \n"))
        p1 <- 0
      } else {
        if (!silent)
          cat(paste("squared distance between means:",p1,sep=" "),"\n","-------------------------------------------","\n")
        count <- count+1         
      }
    } else {
      if (!silent)
        cat(paste("squared distance between means:",p1,sep=" "),"\n","-------------------------------------------","\n")
      count <- count+1         
    }
    gc(verbose = FALSE)
  }
  gc(verbose = FALSE)
  out <- list(dataslide=dataslide,vn.array=vn.array)
  class(out) <- "slider3d"
  out$sliderinfo <- list(fixLM=fixLM,outlineLM=outlines,surfaceLM=surp)
  if (.Platform$OS.type == "windows" && mc.cores > 1)
    stopCluster(cl)
  return(out)
}

# >>> Modified checkLM function of "Morpho" package <<< #
# Scales the landmark size to centroid size if:
# @pt.size = NULL
# @pt.size.scale: number to scale centroid size to sphere radius (0.0001 works well)
checkLM.mod <- function (dat.array, path = NULL, prefix = "", suffix = ".ply", 
                         col = "white", pt.size = NULL, pt.size.scale=NULL, alpha = 1, begin = 1,
                         render = c("w", "s"), point = c("s", "p"), add = FALSE, Rdata = FALSE, 
                         atlas = NULL, text.lm = FALSE) 
{
  k <- NULL
  marked <- NULL
  j <- 1
  if (!Rdata) 
    load <- file2mesh
  outid <- NULL
  point <- point[1]
  radius <- pt.size
  if (is.null(radius)) {
    if (point == "s") 
      
      # ! modified section start ! #
      
      # radius <- (cSize(dat.array[, , 1])/sqrt(nrow(dat.array[,,1]))) * (1/30) # original code in checkLM
      radius <- NULL
    for(i in 1:dim(dat.array)[3]) {
      radius[i] <- cSize(dat.array[,,i])*pt.size.scale
    }
    # ! modified section end ! #
  }
  else { radius <- 10
  }
  size <- radius
  render <- render[1]
  arr <- FALSE
  point <- point[1]
  if (point == "s") {
    rendpoint <- spheres3d
  }
  else if (point == "p") {
    rendpoint <- points3d
  }
  #else {
  #stop("argument \\"point\\" must be \\"s\\" for spheres or \\"p\\" for points")
  #}
  dimDat <- dim(dat.array)
  if (length(dimDat) == 3) {
    n <- dim(dat.array)[3]
    name <- dimnames(dat.array)[[3]]
    arr <- TRUE
  }
  else if (is.list(dat.array)) {
    n <- length(dat.array)
    name <- names(dat.array)
  }
  else {
    stop("data must be 3-dimensional array or a list")
  }
  i <- begin
  if (render == "w") {
    back <- front <- "lines"
    rend <- wire3d
  }
  else {
    back <- front <- "filled"
  }
  if (!add || rgl.cur() == 0) 
    open3d()
  if (!is.null(atlas)) {
    k <- dim(atlas$landmarks)[1]
  }
  meshnames <- paste(path, prefix, name, suffix, sep = "")
  while (i <= n) {
    rgl.bringtotop()
    tmp.name <- meshnames[i]
    if (arr) 
      landmarks <- dat.array[, , i]
    else landmarks <- dat.array[[i]]
    if (is.null(atlas)) {
      outid <- rendpoint(landmarks, radius = radius[i], size = size[i])
      if (text.lm) 
        outid <- c(outid, text3d(landmarks, texts = paste(1:dim(landmarks)[1], 
                                                          sep = ""), cex = 1, adj = c(1, 1.5)))
      if (!is.null(path)) {
        if (!Rdata) {
          tmpmesh <- file2mesh(tmp.name)
        }
        else {
          input <- load(tmp.name)
          tmp.name <- gsub(path, "", tmp.name)
          tmpmesh <- get(input)
        }
        outid <- c(outid, shade3d(tmpmesh, col = col, 
                                  alpha = alpha, back = back, front = front))
        rm(tmpmesh)
        if (Rdata) 
          rm(list = input)
        gc()
      }
    }
    else {
      atlas.tmp <- atlas
      atlas.tmp$mesh <- NULL
      atlas.tmp$landmarks <- landmarks[1:k, ]
      atlas.tmp$patch <- landmarks[-c(1:k), ]
      if (!is.null(path)) {
        if (!Rdata) {
          atlas.tmp$mesh <- file2mesh(tmp.name)
        }
        else {
          input <- load(tmp.name)
          tmp.name <- gsub(path, "", tmp.name)
          atlas.tmp$mesh <- get(input)
        }
      }
      outid <- plotAtlas(atlas.tmp, add = TRUE, alpha = alpha, 
                         pt.size = radius[i], render = render, point = point, 
                         meshcol = col, legend = FALSE)
    }
    answer <- readline(paste("viewing #", i, "(return=next | m=mark current | s=stop viewing)\\n"))
    if (answer == "m") {
      marked[j] <- i
      j <- j + 1
    }
    else if (answer == "s") {
      i <- n + 1
    }
    else i <- i + 1
    rgl.pop(id = outid)
  }
  invisible(marked)
}

