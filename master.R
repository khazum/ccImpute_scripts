suppressPackageStartupMessages({
  library(mclust)
  library(Rtsne)
  library(SummarizedExperiment)
  library(SC3)
  library(SingleCellExperiment)
  library(stats)
  library(Rcpp)
  library(cluster)
  library(DrImpute)
  library(foreach)
  library(doParallel)
  library(Rmagic)
  library(scImpute)
  
  sourceCpp("~/ccImpute/cpp/wCorr_m.cpp")
  sourceCpp("~/ccImpute/cpp/solver.cpp")
})

ccImpute <- function(clust_alg, X, X_log, labels, num_clusters) {
  threshold <- 0
  if(clust_alg == "kmeans"){
    threshold <- .65
  }else if(clust_alg == "hclust"){
    threshold <- .50
  }else if(clust_alg == "dcem-star"){
    threshold <- .95
  }
  
  distances <- list()
  names <- c()
  names <- c(names, "Spearman")
  distances <-list(w_cor_dist_m(X_log, rowVars(X_log)))
  names(distances) <- names
  
  sce <- SingleCellExperiment(
    assays = list(
      counts = as.matrix(X),
      logcounts = X_log
    ),
    colData = labels
  )
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- sc3_prepare(sce, gene_filter = FALSE)
  
  metadata(sce)$sc3$distances <- distances
  sce <- sc3_calc_transfs(sce)
  sce <- sc3_kmeans(sce, num_clusters, clust_alg)
  sce <- sc3_calc_consens(sce)
  
  # Get consensus matrix from the SC3
  cm <- eval(parse(text=paste("metadata(sce)$sc3$consensus$'", toString(num_clusters), "'$consensus", sep="")))
  
  # Remove diagonal entries
  cm <- cm - diag(nrow(cm))
  
  cm[cm < threshold] <- 0
  cm2 <- t(apply((cm), 2,  # Normalize the entries to get weighted average
                 function(i) i/sum(i)))
  
  # Replace NA values with 0
  cm2[is.na(cm2)] <- 0
  
  xlog_t = t(X_log)
  x_imp <- xlog_t
  
  t2 = x_imp == 0
  x_t_vote <- x_imp
  
  x_t_vote[t2] <- -1
  x_t_vote[x_imp > 0] <- 1
  
  #compute votes - majority wins - negative means an actual 0, otherwise it is some positive value
  votes <- matrix(0L, nrow = nrow(x_imp), ncol = ncol(x_imp))
  votes[t2] <- (cm2 %*% x_t_vote)[t2]
  
  t3 <- votes > 0
  x_imp <- solve_dropouts(cm2, x_imp, which(t3, arr.ind = TRUE))
  
  return(t(x_imp))
}

impute <- function(method, X, X_log, labels, num_clusters, dataset) {
  start_time <- Sys.time()
  X_imp <- NULL
  
  if(method == "ccImpute_kmeans" || method == "ccImpute_greedy" || method == "ccImpute_dcem-star" || method == "ccImpute_hclust"){
    X_imp <- ccImpute(strsplit(method, '_')[[1]][2], X, X_log, labels, num_clusters)
  }else if(method == "none"){
    X_imp <- X_log
  }else if(method == "drimpute"){
    X_imp <- DrImpute(X_log)
  }else if(method == "magic"){
    X_imp <- t(as.matrix(magic(t(X_log), genes="all_genes")))
  }else if(method == "scimpute"){
    write.csv(x=X, paste("~/ccImpute/datasets/temp_sci_", dataset, ".csv", sep=""))
    
    scimpute(# full path to raw count matrix
      count_path = paste("~/ccImpute/datasets/temp_sci_", dataset, ".csv", sep=""),
      infile = "csv",           # format of input file
      outfile = "csv",          # format of output file
      out_dir = paste("~/ccImpute/datasets/temp_sci_", dataset, "/", sep=""), # full path to output directory
      labeled = FALSE,          # cell type labels not available
      drop_thre = 0.5,          # threshold set on dropout probability
      Kcluster = num_clusters,  # 2 cell subpopulations
      ncores = 15)              # number of cores used in parallel co
    end_time <- Sys.time()
    
    X2 = read.csv(paste("~/ccImpute/datasets/temp_sci_", dataset, "/scimpute_count.csv", sep=""),row.names = 1, header=TRUE)
    
    librarySizes <- colSums(X2)
    X_norm <- t(t(X2)/librarySizes)*1000000
    
    X_imp <- log2(X_norm + 1)
    
    unlink(paste("~/ccImpute/datasets/temp_sci_", dataset, "/", sep=""), recursive = TRUE)
    unlink(paste("~/ccImpute/datasets/temp_sci_", dataset, ".csv", sep=""))
  }else if(method == "dca"){
    write.csv(x=X, paste("~/ccImpute/datasets/temp_dca_", dataset, ".csv", sep=""))
    system(paste("dca ~/ccImpute/datasets/temp_dca_", dataset, ".csv", " ~/ccImpute/datasets/temp_dca_", dataset, "/", " --threads 15", sep=""))
    X2 <- read.csv(paste("~/ccImpute/datasets/temp_dca_", dataset, "/mean.tsv", sep=""), sep="\t", row.names = 1, header=TRUE)
    librarySizes <- colSums(X2)
    X_norm <- t(t(X2)/librarySizes)*1000000
    X_imp <- log2(X_norm + 1)
    unlink(paste("~/ccImpute/datasets/temp_dca_", dataset, "/", sep=""), recursive = TRUE)
    unlink(paste("~/ccImpute/datasets/temp_dca_", dataset, ".csv", sep=""))
  }else if(method == "deepimpute"){
    write.csv(x=X, paste("~/ccImpute/datasets/temp_", dataset, ".csv", sep=""))
    system(paste("deepImpute ~/ccImpute/datasets/temp_", dataset, ".csv", " -o ~/ccImpute/datasets/temp2_", dataset, ".csv --cores 15", sep=""))
    X2 <- read.csv(paste("~/ccImpute/datasets/temp2_", dataset, ".csv", sep=""), row.names = 1, header=TRUE)
    librarySizes <- colSums(X2)
    X_norm <- t(t(X2)/librarySizes)*1000000
    X_imp <- log2(X_norm + 1)
    unlink(paste("~/ccImpute/datasets/temp1_", dataset, ".csv", sep=""))
    unlink(paste("~/ccImpute/datasets/temp2_", dataset, ".csv", sep=""))
  }
  
  end_time <- Sys.time()
  
  xlog_t <- t(X_imp)
  
  p <- 30
  
  cells <- ncol(X)
  
  if(cells > 1000){
      print("Reducing rank")
      pca_red <- prcomp(as.matrix(xlog_t), rank. = 1000)$x
      restarts <- 50
      
  }
  else{
      pca_red <- prcomp(as.matrix(xlog_t))$x
      if (ncol(t(as.matrix(xlog_t))) <= p*2){
          p <- 9
      }
      restarts <- 1000
  }
  
  tsne_red <- Rtsne(as.matrix(xlog_t), perplexity = p, check_duplicates = FALSE)$Y
  
  cl <- parallel::makeCluster(6, outfile = "")
  doParallel::registerDoParallel(cl, cores = 6)
  
  # calculate distances in parallel
  results <- foreach::foreach(i = 1:6, .combine=c, .packages = c("mclust","Rtsne", "cluster", "caret", "e1071")) %dopar% {
      if(i==1){
          adjustedRandIndex(kmeans(
              pca_red,
              centers = num_clusters,
              iter.max = 1e+09,
              nstart = restarts
          )$cluster,
          labels)
      }else if(i==2){
          adjustedRandIndex(kmeans(
              tsne_red,
              centers = num_clusters,
              iter.max = 1e+09,
              nstart = restarts
          )$cluster,
          labels)
      }else if(i==3){
          xlog_t_do <- t(X_log)
          xlog_t_do[xlog_t_do==0] <- xlog_t[xlog_t_do==0]
          dist <- as.matrix(stats::dist(xlog_t_do, method = "euclidean", p=2))
          int_labels <- as.numeric((as.factor(labels)))
          silh <- silhouette(int_labels, dist)
          as.numeric(summary(silh)['avg.width'])
      }
      else if(i==4){
          factor_labels4 <- as.factor(labels)
          folds4 = caret::createFolds(factor_labels4, k = 10)
          
          cv4 = lapply(folds4, function(x) { # start of function
              train_labels4 <- factor_labels4[-x]
              classifier4 = e1071::svm(x = pca_red[-x,1:2], y = train_labels4, type = "C-classification", kernel = "linear")
              if(length(x) == 1){
                  x = rep(x,2)
              }
              y_pred4 = predict(classifier4, as.matrix(pca_red[x,1:2]))
              test_labels4 <- factor_labels4[x]
              cm4 = table(test_labels4, y_pred4)
              accuracy4 = sum(diag(cm4)) / sum(cm4)
              return(accuracy4)
          })
          round(mean(as.numeric(cv4)), digits = 4)
      }
      else if(i==5){
          factor_labels5 <- as.factor(labels)
          folds5 = caret::createFolds(factor_labels5, k = 10)
          cv5 = lapply(folds5, function(x) { # start of function
              train_labels5 <- factor_labels5[-x]
              classifier5 = e1071::svm(x = pca_red[-x,1:4], y = train_labels5, type = "C-classification", kernel = "linear")
              if(length(x) == 1){
                  x = rep(x,2)
              }
              y_pred5 = predict(classifier5, as.matrix(pca_red[x,1:4]))
              test_labels5 <- factor_labels5[x]
              cm5 = table(test_labels5, y_pred5)
              accuracy5 = sum(diag(cm5)) / sum(cm5)
              return(accuracy5)
          })
          round(mean(as.numeric(cv5)), digits = 4)
      }
      else if(i==6){
          factor_labels6 <- as.factor(labels)
          folds6 = caret::createFolds(factor_labels6, k = 10)
          cv6 = lapply(folds6, function(x) { # start of function
              train_labels6 <- factor_labels6[-x]
              classifier6 = e1071::svm(x = tsne_red[-x,], y = train_labels6, type = "C-classification", kernel = "linear")
              if(length(x) == 1){
                  x = rep(x,2)
              }
              y_pred6 = predict(classifier6, as.matrix(tsne_red[x,]))
              test_labels6 <- factor_labels6[x]
              cm6 = table(test_labels6, y_pred6)
              accuracy6 = sum(diag(cm6)) / sum(cm6)
              return(accuracy6)
          })
          round(mean(as.numeric(cv6)), digits = 4)
      }
  }
  
  # stop local cluster
  parallel::stopCluster(cl)
  
  xlog_t[xlog_t<0.176]<-0
  prop_zeroes_removed <- (sum(X_log==0)-sum(xlog_t==0))/(nrow(X_log) * ncol(X_log))
  
  return(c(results[1], results[2], results[4], results[5], results[6], difftime(end_time, start_time, units="secs") , prop_zeroes_removed, results[3]))
}

driver <- function(method, dataset, repeats){
  sce <- readRDS(file = paste("~/ccImpute/datasets/", dataset, ".rds", sep=""))
  
  X <- assays(sce)$counts
  X_log <- assays(sce)$logcounts
  
  print(paste(method,dataset,repeats, "||", "Genes(rows):", nrow(X), "Cells(cols):", ncol(X), "||", sep=" "))
  
  labels<-if(is.null(colData(sce)$cell_type2)) colData(sce)$cell_type1 else colData(sce)$cell_type2
  row_sums <- rowSums(X[,-1])
  X <- X[row_sums>0,] # remove genes that are not expressed at all
  X_log <- X_log[row_sums>0,] # remove genes that are not expressed at all
  
  num_clusters = length(unique(labels))
  
  data_aris <- replicate(repeats, impute(method, X, X_log, labels, num_clusters, dataset))
  
  means <- rowMeans(data_aris)
  stdevs <- rowSds(data_aris)
  
  print(c(method, "Clustering results: ", dataset))
  print("c1, c2, time, prop_zeroes_rm, silh_pca_avr")
  print(means)
  print(stdevs)
  fileConn<-eval(parse(text=paste('file("~/ccImpute/results/', method, "_", dataset, '_', repeats, '_', Sys.time(), '")', sep="")))
  writeLines(c("c1,c2,time, prop_zeroes removed, silh_pca_avr",paste(method,dataset,repeats, "||", "Genes(rows):", nrow(X), "Cells(cols):", ncol(X), "||", sep=" "), means, stdevs), fileConn)
  close(fileConn)
}

args = commandArgs(trailingOnly=TRUE)

driver(args[1],args[2],strtoi(args[3], base=10L))
