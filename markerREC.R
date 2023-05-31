
#' Extract expression matrix from Seurat object or matrix
#'
#' This function extracts the expression matrix from either a Seurat object or a matrix object. If a Seurat object is provided, the function can optionally specify the assay and slot to be used for extracting the expression matrix. If a matrix is provided, the function will create a Seurat object with the provided matrix and default metadata fields. The function can also specify the minimum number of cells and features required to keep genes in the expression matrix.
#'
#' @name assay2mtx
#' @author Qiong Zhang
#' @param object A Seurat object or matrix containing expression data
#' @param assay A character string specifying the assay name to extract expression data from the Seurat object
#' @param slot A character string specifying the slot name to extract expression data from the Seurat object
#' @param min.cells An integer specifying the minimum number of cells required to keep genes in the expression matrix
#' @param min.features An integer specifying the minimum number of features required to keep genes in the expression matrix
#' @return A matrix containing the expression data, with gene names as row names and cell barcodes as column names
#' @examples
#' expmtx <- assay2mtx(object = seurat_object, assay = "RNA", slot = "counts", min.cells = 10, min.features = 5)
#' expmtx <-assay2mtx(object = expression_matrix, min.cells = 10, min.features = 5)
#'
assay2mtx <- function(object, assay = NULL, slot = "counts", min.cells = 3, min.features = 1, update = F) {
  if (any(methods::is(object) %in% "Seurat")) {
    assay <- if (is.null(assay)) Seurat::DefaultAssay(object) else assay
    object <- if (isTRUE(update)) Seurat::UpdateSeuratObject(object) else object
  } else {
    object <- Seurat::CreateSeuratObject(counts = object, project = "GSES", assay = "RNA", min.cells = min.cells, min.features = min.features)
    assay <- "RNA"
  }
  Seurat::GetAssayData(object, assay = assay, slot = slot)
}

#' Calculate gini specific score for markers
#'
#' @param x vector of gene
#' @name ginirev
#' @author Qiong Zhang
#' @return feature-specific gini scores.
#'
ginirev <- function(x, rev = T){
  .x <- sort(x)
  .G <- sum(.x * 1L:length(.x))
  .G <- 2 * .G/sum(.x) - (length(.x) + 1L)
  .G <- .G/length(.x)
  if (rev) .G <- 1 - .G
  return(.G)
}

#' Transform expression matrix to binary type based on give mode or customize value
#'
#' @param mat vector
#' @param pattern transform pattern
#' @name expMat2Bin
#' @return binary matrix
#' @author Qiong Zhang
#'
expMat2Bin <- function(mat, pattern = "min", cutoff = NULL){
  .stabe2b <- function(mat, pattern, cutoff = NULL){
    if (pattern == "rowMean"){
      .x.m <- Matrix::rowMeans(mat)
      return( list(mat = mat > .x.m[rownames(mat)], cutoff = cutoff))
    }
    .x <- mat@x
    if (is.null(cutoff)){
      .cutoff <- switch(pattern,
                        "min" = min(.x),
                        "max" = max(.x),
                        "mean" = mean(.x),
                        "median" = median(.x),
                        stop("Error: pattern need to be one of 'min', 'max', 'median', 'mean' and 'adpmed'. ")
      )
    }else{
      .cutoff <- as.numeric(cutff)
    }
    .index <- which(.x > .cutoff)
    .dp <- diff(mat@p)
    .col.index <- (rep(seq_along(.dp), .dp))[.index]
    .row.index <- (mat@i+1)[.index]
    .mat.new <- Matrix::sparseMatrix(.row.index[1], .col.index[1], dims = mat@Dim)
    .mat.new[cbind(.row.index, .col.index)] <- 1
    .mat.new <- as(.mat.new, "dgCMatrix")
    rownames(.mat.new) <- rownames(mat)
    return(list(mat = .mat.new, cutoff = .cutoff))
  }
  .adpmed <- function(mat){
    .m.multi <- Mclust(Matrix::rowSums(mat), 1:30, modelNames=c("V"), verbose = FALSE)
    if (.m.multi$G == 1) return(.stabe2b(mat, type = "median") )
    .mat.new <- list()
    .cutff <- rep(0, .m.multi$G)
    for (i in seq_along(1:.m.multi$G)) {
      .x.clu <- mat[which(.m.multi$classification == i), ]
      .cutff[i] <- median(.x.clu@x)
      .index <- which(.x.clu@x >= .cutff[i])
      .dp <- diff(.x.clu@p)
      .col.index <- (rep(seq_along(.dp), .dp))[.index]
      .row.index <- (.x.clu@i+1)[.index]
      .mat.new[[i]] <- Matrix::sparseMatrix(.row.index[1], .col.index[1], dims=.x.clu@Dim)
      .mat.new[[i]][cbind(.row.index, .col.index)] <- 1
      .mat.new[[i]] <- as(.mat.new[[i]], "dgCMatrix")
      rownames(.mat.new[[i]]) <- rownames(.x.clu)
    }
    .mat.new <- do.call(rbind, .mat.new)
    return(list(mat = .mat.new, cutoff = .cutoff))
  }
  if (pattern == "adpmed"){
    .exp2bin <- .adpmed(mat)
  }else {
    .exp2bin <- .stabe2b(mat = mat, pattern = pattern, cutoff = cutoff)
  }
  return(.exp2bin)
}

#' rescale a vector between 0 (lowest) and 1 (highest)
#'
#' @param x vector
#' @name scaleZO
#' @return rescale vector with value between 0 and 1
#' @author Qiong Zhang
#'
scaleZO <- function(x) {
  (x-min(x))/(max(x)-min(x))
}

#' calculate prior probabilities for clusters
#'
#' @param x vector
#' @name calProb4C
#' @return probabilities of clusters
#' @author Qiong Zhang
#'
calProb4C <- function(x) {
  table(x) / length(x)
}

#' calculate prior probabilities for genes
#'
#' @param mat expression matrix
#' @name calProb4G
#' @return probabilities for genes
#' @author Qiong Zhang
#'
calProb4G <- function(mat) {
  .prob4G <- as.vector( (Matrix::rowSums(mat)) / ncol(mat) )
  names(.prob4G) <- rownames(mat)
  return(.prob4G)
}

#' calculate probabilities genes in clusters
#'
#' @param mat expression matrix
#' @param cluProb
#' @param cluID
#' @param cpu
#' @name calProbGIC
#' @return probabilities of genes in clusters
#' @author Qiong Zhang
#'
#get gene probability conditional on cluster
calProbGIC <- function(mat, cluProb, cluID, cpu = 1) {
  .clusterID <- names(cluProb)
  if (cpu > 1) {
    .GIC <- do.call(cbind, mclapply(1:length(.clusterID), function(i) Matrix::rowMeans(mat[,which(cluID == .clusterID[i])]), mc.cores = cpu))
  } else {
    .GIC <- do.call(cbind, lapply(1:length(.clusterID), function(i) Matrix::rowMeans(mat[,which(cluID == .clusterID[i])])))
  }
  .GIC <- as(.GIC, "dgCMatrix")
  colnames(.GIC) <- .clusterID
  return(.GIC)
}

#' calculate probabilities of genes in given condition
#'
#' @param mat expression matrix
#' @param cluProb
#' @param cluID
#' @param cpu
#' @name calProbCIG
#' @return probabilities of genes in given condition
#' @author Qiong Zhang
#'
calProbCIG <- function(binmat, probG, ProbC) {

  .prob <- t(apply(log(binmat), 1, function(x) x + log(ProbC)))
  .prob <- exp(apply(.prob, 2, function(x) x - log(probG)))
  colnames(.prob) <- names(ProbC)
  return(as(.prob, "dgCMatrix"))
}

#' calculate the feature score of genes
#'
#' @param ginc gene expression probability
#' @param binmat geneXcell binary matrix
#' @param scale scale the scores between groups
#' @name calGeneFea
#' @return feature score of genes
#' @author Qiong Zhang
#'
calGeneFea <- function(ginc, binmat, scale = F) {

  if (isTRUE(scale)){
    .genefeature <- apply(ginc * binmat, 2, scaleZO)
  } else {
    .genefeature <- exp(log(ginc) + log(binmat))
  }
  colnames(.genefeature) <- colnames(binmat)
  return(.genefeature)
}

#' Calculate the entropy of probability vector
#'
#' @param x probability vector
#' @name probEnt
#' @return the entropy of probability vector
#' @author Qiong Zhang
#'
probEnt <- function(x) {
  .et <- 0
  if (any(x > 0 )) {
    .t <- x[which(x > 0)]
    .et <- -(sum(.t * (log(.t))))
  }
  return(.et)
}

#' Calculate the shanno entropy for vector
#'
#' @param x vector
#' @name shannoEnt
#' @return the entropy value
#' @author Qiong Zhang
#'
shannoEnt <- function(x){
  .p <- x/sum(x,na.rm=TRUE)
  .entropy <- -sum(.p * log(.p), na.rm=TRUE)
  return(.entropy)
}

#' Calculate mutual information for two vectors
#'
#' @param a vector a
#' @param b vector b
#' @name calMutInfo
#' @return mutual information between two vectors
#' @author Qiong Zhang
#'
calMutInfo <- function(a, b) {
  .al <- length(a)
  .bl <- length(b)
  if (.al != .bl) stop("Inequal elements between the two vectors ")
  .a <- probEnt(calProb4C(a))
  .b <- probEnt(calProb4C(b))
  .ab <- probEnt(table(a,b) / .bl)
  .a + .b - .ab
}

#' Calculate log-likelihood similarity for two vectors
#'
#' @param a vector a
#' @param b vector b
#' @name calLLS
#' @return log-likelihood similarity between two vectors
#' @author Qiong Zhang
#'
calLLRSim <- function(a, b, n) {
  .a <- sum(a)
  .b <- sum(b)
  .aANDb <- sum(a & b)
  .aSPE <- .a - .aANDb
  .bSPE <- .b - .aANDb
  .abNOT <- sum(!(a + b))
  .entropy.row <- shannoEnt(c(.b, n - .b ))
  .entropy.col <- shannoEnt(c(.a, n - .a ))
  .entropy.mat <- shannoEnt(c(.aANDb, .bSPE, .aSPE, .abNOT))
  .llRSim <- 2*n*(.entropy.row + .entropy.col - .entropy.mat)
  return(.llRSim)
}

#' Calculate TFIDF for clusters
#'
#' @param x  a
#' @param b vector b
#' @name calTFIDF
#' @return calculate TFIDF for clusters
#' @author Qiong Zhang
#'
calTFIDF <- function(condMat) {
  .idf <- log(ncol(condMat)) - log(Matrix::rowSums(expMat2Bin(condMat, pattern = "min")$mat))
  return(exp(log(condMat) + log(.idf)))
}


#' QC for clusterIDs
#'
#' @param mat  gene expression matrix
#' @param clusterIDs a vector contains clusterID for cells
#' @name clusidQC
#' @return validated clusterIDs
#' @author Qiong Zhang
#'
clusidQC <- function(mat, clusterIDs){

  if ( any(table(clusterIDs) < 2)) stop("Several clusters comprised of only one cell. Please check clusterIDs")
  if (length(levels(clusterIDs)) < 2) stop("Only one cluster exists in clusterIDs!")
  .ncell <- ncol(mat)
  .cluIDs.len <- length(clusterIDs)
  if (.ncell != .cluIDs.len) stop("Inequal cell number and clusterIDs.")
  return(droplevels(clusterIDs))
}

#' Construct the pesudo binary matrix
#'
#' @param clusterIDs a vector contains clusterID for cells
#' @param value fill in value
#' @name buildPesudoBinMatrix
#' @return pesudo binmatrix matrix
#' @author Qiong Zhang
#'
buildPesudoBinMatrix <- function(clusterIDs, combine = FALSE, value = 3){
  .clusters.levels <- levels(clusterIDs)
  clusterIDs <- as.numeric(clusterIDs)
  .cluster.uniq <- unique(clusterIDs)
  .cluster.name <- .clusters.levels[.cluster.uniq]
  .pesu.lst <- as.list(.cluster.uniq)
  if (isTRUE(combine)){
    .defvalue <- length(.cluster.uniq) / 2
    if (value > .defvalue) value <- .defvalue
    for(i in 2:value){
      .tmp.combination <- as.data.frame(combn(.cluster.uniq, i))
      .pesu.lst <- c(.pesu.lst, as.list(.tmp.combination))
      .cluster.name <- c(.cluster.name, unlist(lapply(.tmp.combination, function(x) paste0(.clusters.levels[x], collapse = ".AND."))))
    }
  }
  .pesu.mat <- do.call(rbind, lapply(.pesu.lst, function(x) {
    .idx <- clusterIDs %in% x
    clusterIDs[.idx] <- 1
    clusterIDs[!.idx] <- 0
    return(clusterIDs)
  }))
  rownames(.pesu.mat) <- .cluster.name
  .pesu.mat <- as(.pesu.mat, "dgCMatrix")
  return(.pesu.mat)
}

#' Calculate marker specificity
#'
#' @param mat expression matrix
#' @param clusterIDs cluster IDs for each cell.
#' @param binmethod Method used to binarize the matrix
#' @param expR The proportion of feature expressed cells in a cluster.
#' @param cpu How many cores will be used for detection of markers
#' @param marker The number of markers returned.
#' @name LLRMarker
#' @author Qiong Zhang
#' @return data.frame with marker features of each cluster.
#'
LLRMarker <- function(mat, clusterIDs, binmethod = "min", expR = 0.1, cpu = 1, marker = 20) {

  clusterIDs <- clusidQC(mat, clusterIDs)
  .clusIDs <- as.factor(clusterIDs)
  .clusterIDs.int <- as.integer(clusterIDs)
  .clusterIDs.lev <- levels(clusterIDs)
  .mat <- as(mat, "dgCMatrix")
  .mat.bin <- expMat2Bin(.mat, pattern = binmethod)
  .clu.prob <- calProb4C(.clusterIDs.int)
  names(.clu.prob) <- .clusterIDs.lev
  .binmat <- calProbGIC(.mat.bin$mat, .clu.prob, clusterIDs, cpu = cpu)
  .binmat.f.l <- apply(.binmat, 1, function(x) any(x > expR))
  .binmat.f <- .binmat[.binmat.f.l, ]
  .gene.prob <- calProb4G(.mat.bin$mat[.binmat.f.l, ])
  .ginc <- calProbCIG(.binmat.f, .gene.prob, .clu.prob)
  .gfeat <- calGeneFea(.ginc, .binmat.f)
  .w <- replicate(ncol(.gfeat), Matrix::rowSums(.gfeat))
  .gfeat.n <- .gfeat^3 / .w ^ 2
  .gfeat.s <- sqrt(.gfeat / (.w - .gfeat))
  .clu.n <- ncol(.gfeat.s)
  .marker.n <- do.call(cbind, lapply(1:.clu.n, function(x, marker) {
    names(sort(.gfeat.n[, x], decreasing = T)[1:marker])}, marker))
  .marker.s <- do.call(cbind, lapply(1:.clu.n, function(x, marker) {
    .x.s <- names(sort(.gfeat.s[, x], decreasing = T)[1:marker])
    .k.b <- names(sort(.binmat[.x.s, x], decreasing = T))
    }, marker))
  .col.name <- colnames(.gfeat.s)
  colnames(.marker.s) <- .col.name
  colnames(.marker.n) <- .col.name
  return(list( NormM = .marker.n, SpecM = .marker.s, pct = .binmat))
}

#' Recommend marker genes for expression matrix or Seurat object
#'
#' @param obj Seurat object or expression matrix
#' @param clusterIDs The cluster ID for each cell or the column name which contains cluster ID information in Seurat object.
#' @param binmethod Method used to binarize the matrix
#' @param expR The proportion of feature expressed cells in a cluster.
#' @param cpu How many cores will be used for detection of markers
#' @param marker The number of markers returned.
#' @name markerREC
#' @author Qiong Zhang
#' @return data.frame of marker features in each cluster.
#'
markerREC <- function(obj, assay = "RNA", slot = "data", binmethod = "min", expR = 0.1, cpu = 1, marker = 20, clusterIDs = NULL) {
  require("Seurat", quietly = T)
  .mat <- assay2mtx(object = obj, assay = assay, slot = slot, min.cells = 3, min.features = 1)
  if (is.null(clusterIDs)) {
    clusterIDs <- Seurat::Idents(obj)
  } else{
    if ( length(clusterIDs) != ncol(.mat)) clusterIDs <- obj@meta.data[[clusterIDs]]
  }
  .markers <- LLRMarker(mat = .mat, clusterIDs = clusterIDs, binmethod = binmethod, expR = expR, cpu = cpu, marker = marker)
  return(.markers)
}
