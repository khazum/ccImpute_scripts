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
    library(devtools)
    library(e1071)
    library(caret)
})

devtools::install_github("khazum/ccImpute_exp")

impute <- function(method, X, X_log, labels, num_clusters, dataset) {
    start_time <- Sys.time()
    X_imp <- NULL
    
    if(method == "ccImpute"){
        X_imp <- ccImpute::impute(X_log, k = num_clusters, nCores=15)
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
        pca_red <- prcomp(as.matrix(xlog_t), rank. = 1000, center = TRUE, scale = TRUE)$x
        restarts <- 50
        
    }
    else{
        pca_red <- prcomp(as.matrix(xlog_t), center = TRUE, scale = TRUE)$x
        if (ncol(t(as.matrix(xlog_t))) <= p*2){
            p <- 9
        }
        restarts <- 1000
    }
    
    cl <- parallel::makeCluster(6, outfile = "")
    doParallel::registerDoParallel(cl, cores = 5)
    
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
            tsne_red <- Rtsne(as.matrix(xlog_t), perplexity = p, check_duplicates = FALSE)$Y
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
            factor_labels <- as.factor(labels)
            folds = caret::createFolds(factor_labels, k = 10)
            cv = lapply(folds, function(x) { # start of function
                classifier = e1071::svm(x = pca_red[-x,1:2], y = factor_labels[-x], type = "C-classification", kernel = "linear")
                y_pred = predict(classifier, pca_red[x,1:2])
                cm = table(factor_labels[x], y_pred)
                accuracy = sum(diag(cm)) / sum(cm)
                return(accuracy)
            })
            mean(as.numeric(cv))
        }
        else if(i==5){
            factor_labels <- as.factor(labels)
            folds = caret::createFolds(factor_labels, k = 10)
            cv = lapply(folds, function(x) { # start of function
                classifier = e1071::svm(x = pca_red[-x,1:4], y = factor_labels[-x], type = "C-classification", kernel = "linear")
                y_pred = predict(classifier, pca_red[x,1:4])
                cm = table(factor_labels[x], y_pred)
                accuracy = sum(diag(cm)) / sum(cm)
                return(accuracy)
            })
            mean(as.numeric(cv))
        }
        else if(i==6){
            tsne_red <- Rtsne(as.matrix(xlog_t), perplexity = p, check_duplicates = FALSE)$Y
            factor_labels <- as.factor(labels)
            folds = caret::createFolds(factor_labels, k = 10)
            cv = lapply(folds, function(x) { # start of function
                classifier = e1071::svm(x = tsne_red[-x,], y = factor_labels[-x], type = "C-classification", kernel = "linear")
                y_pred = predict(classifier, tsne_red[x,])
                cm = table(factor_labels[x], y_pred)
                accuracy = sum(diag(cm)) / sum(cm)
                return(accuracy)
            })
            mean(as.numeric(cv))
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
# args = c("ccImpute", "usoskin", 1)
# driver(args[1],args[2],strtoi(args[3], base=10L))

