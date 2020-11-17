# Functions to run Gaussian mixture modeling on SUPPA results

import_suppa_TCGA <- function(fpath){
  
  require(data.table)
  
  psi <- data.table::fread(file = fpath, sep = "\t", showProgress = TRUE,
                           data.table = FALSE  )
  rownames(psi) <- psi$V1
  psi <- psi[ , -1]
  colnames(psi) <- gsub("\\.", "-", x = colnames(psi))
  
  return(psi)
  
}

import_suppa_GTEX <- function(fpath){
  
  require(data.table)
  
  psi <- data.table::fread(file = fpath, sep = "\t", showProgress = TRUE,
                           data.table = FALSE  )
  rownames(psi) <- psi$V1
  psi <- psi[ , -1]
  colnames(psi) <- gsub("\\.", "-", x = colnames(psi))
  
  return(psi)
  
}

annotate_GTEX_barcode <- function(barcodes){
  
  gtex <- read.csv("./TCGA_annotation/gtex_annotation.csv", header = T,
                   stringsAsFactors = F)
  idx <- match(barcodes, gtex$Run)
  return(gtex[idx,])
}

annotate_TCGA_barcode <- function(barcodes){
  
  data_folder = "./TCGA_annotation/"
  
  #Annotate barcodes with sample type and tissues  
  tss.dict <- read.csv(paste0(data_folder, 
                              "/TCGA_TissueSourceSite_Dictionary.csv"), header = T,
                       stringsAsFactors = F)
  
  sampleType.dict <- read.csv(paste0(data_folder, 
                                     "/TCGA_SampleType_Dictionary.csv"), header = T,
                              stringsAsFactors = F, 
                              colClasses = c("character","character", "character"))
  
  studyAbbrev.dict <- read.csv(paste0(data_folder, 
                                      "/TCGA_StudyAbbreviations_Dictionary.csv"), header = T,
                               stringsAsFactors = F)
  
  project.dict <- read.csv(paste0(data_folder, 
                                  "/TCGA_Projects_Dictionary.csv"), header = T,
                           stringsAsFactors = F)
  
  bc_annot <- lapply(seq_along(barcodes), function(i){
    
    tokens <- strsplit(barcodes[i], "-")
    TSS <- tokens[[1]][2]
    participant <- tokens[[1]][3]
    sampleType <- substr(tokens[[1]][4], 1, 2)
    vial <- substr(tokens[[1]][4], 3, 3)
    
    TSS_desc <- tss.dict$Study.Name[match(TSS, tss.dict$TSS.Code)]
    Study_abbrev <- studyAbbrev.dict$Study.Abbreviation[match(TSS_desc,studyAbbrev.dict$Study.Name)]
    sample_desc <- sampleType.dict$Definition[match(sampleType, sampleType.dict$Code)]
    sample_shortLetterCode <- sampleType.dict$Short.Letter.Code[match(sampleType, sampleType.dict$Code)]
    primarySite <- project.dict$Primary.Site[match(Study_abbrev, project.dict$Project)]
    
    return(data.frame(
      Barcode = barcodes[i],
      TSS = TSS,
      TSS_desc = TSS_desc,
      Study_abbrev = Study_abbrev,
      participant = participant,
      sampleType = sampleType,
      vial = vial,
      sample_desc = sample_desc,
      sample_shortLetterCode = sample_shortLetterCode,
      primarySite = primarySite
    ))
    
  })
  
  return(do.call(rbind, bc_annot))
  
}

mixtureClusterAnalysis <- function(analysis_name, PSI_tumor, PSI_control,
                                       events, min_cluster, min_tumor, min_controls,
                                       clust_dir){
  
  require(doParallel)
  require(mclust)
  require(dplyr)
  require(parallel)
  
  if(dir.exists(clust_dir)){
    unlink(clust_dir, force = T, recursive = T)
  }
  dir.create(clust_dir)
  
  #Remove events based on coverage and number of samples
  flags <- data.frame(Junction.id = events$event_id,
                      Remove = rep(FALSE, nrow(events)))
  
  # res <- mclapply(seq_along(events$event_id), function(i){
  #   
  #   remove <- FALSE
  #   se.id <- events$event_id[i]
  #   gene <- events$gene_symbol[i]
  #   
  #   idx.case <- which(!is.na(PSI_tumor[se.id, ]))
  #   idx.control <- which(!is.na(PSI_control[se.id, ]))
  #   
  #   exp.case <- PSI_tumor[se.id, idx.case, drop = FALSE]
  #   exp.control <- PSI_control[se.id, idx.control, drop = FALSE]
  #   
  #   #samples event is detected
  #   if(length(exp.case) < min_tumor || length(exp.control) < min_controls){
  #     remove <- TRUE
  #   }
  #   
  #   if(sd(exp.case)==0 || sd(exp.control)==0){
  #     remove <- TRUE
  #   }
  #   
  #   return(data.frame(Junction.id = se.id, Remove = remove))
  #   
  # }, mc.cores = 16)
  # 
  # flags <- do.call(rbind, res)
  # head(flags)
  # table(flags$Remove)
  
  for(i in 1:nrow(events)){

    se.id <- events$event_id[i]
    gene <- events$gene_symbol[i]

    idx.case <- which(!is.na(PSI_tumor[se.id, ]))
    idx.control <- which(!is.na(PSI_control[se.id, ]))

    exp.case <- PSI_tumor[se.id, idx.case, drop = FALSE]
    exp.control <- PSI_control[se.id, idx.control, drop = FALSE]

    #samples event is detected
    if(length(exp.case) < min_tumor || length(exp.control) < min_controls){
      flags$Remove[i] <- TRUE
      next
    }

    if(sd(exp.case)==0 || sd(exp.control)==0){
      flags$Remove[i] <- TRUE
      next
    }

  }
  
  # Events for GMM clustering
  id.keep <- flags$Junction.id[!flags$Remove]
  events_filt <- events %>% filter(event_id %in% id.keep)
  
  cat("INFO: Events to be clustered: ", nrow(events_filt), "\n")
  
  #Run mixture in parallel
  # number of models to fit
  nmodels <- nrow(events_filt)
  
  # we first register a backend specifying the number
  # of cores that we wish to use
  registerDoParallel(cores = 16)
  
  foreach(i=1:nmodels, .combine = rbind) %dopar% {
    
    se.id <- events_filt$event_id[i]
    gene <- events_filt$gene_symbol[i]
    
    cat("Processing ", i, " ", se.id, "\n")
    
    idx.case <- which(!is.na(PSI_tumor[se.id, ]))
    idx.control <- which(!is.na(PSI_control[se.id, ]))
    
    exp.case <- PSI_tumor[se.id, idx.case, drop = FALSE]
    exp.control <- PSI_control[se.id, idx.control, drop = FALSE]
    
    df.exp <- data.frame(Junction = se.id,
                         Gene = gene,
                         PSI = c(as.numeric(exp.case), as.numeric(exp.control)),
                         Source = c(rep("Tumor", length(exp.case)),
                                    rep("Normal", length(exp.control))),
                         Sample = c(colnames(exp.case), colnames(exp.control)),
                         stringsAsFactors = FALSE)
    
    #First fit
    an.error.occured <- FALSE
    tryCatch( {mod.mixed <- mclust::densityMclust(df.exp$PSI, G=1:3) },
              error = function(e) {an.error.occured <<- TRUE} )
    
    if(!an.error.occured){
      
      nClust <- mod.mixed$G
      clust.assigned <- levels(factor(mod.mixed$classification))
      nAssignClust <- length(clust.assigned)
      
      if(length(clust.assigned)<nClust){ #empty clusters
        
        newClass <- 1:length(clust.assigned)
        newAssigClass <- rep(NA,nrow(df.exp))
        for(c in 1:length(newClass)){
          newAssigClass[ mod.mixed$classification == clust.assigned[c] ] = c
        }
        df.exp$Classification <- factor(newAssigClass)
        #table(newAssigClass)
        #table(mod.mixed$classification)
        
      }else{
        df.exp$Classification <- factor(mod.mixed$classification)
      }
      
      save(df.exp, file = paste0(clust_dir, analysis_name, "_",
                                 se.id, ".RData"))
      
    }
    
    
    
    
  }#foreach model fiting
  
} #

# build a data frame with cluster composition (% tumor and normal) accross
# all models
computeClusterComposition <- function(analysis_name, 
                                      clust_dir = "./clust_output", events){
  
  models <- list.files(clust_dir, pattern = ".RData")
  aux <- gsub(paste0(analysis_name, "_"), "", models)
  junction_ids <- gsub(paste0(".RData"), "", aux)
  
  main.counts <- mclapply(seq_along(models), function(i){
    
    jannot <- dplyr::filter(events, event_id == junction_ids[i])
    # Load df.exp data frame
    # Cluster classification for a given model (i.e junction) 
    load(file = file.path(clust_dir, models[i]))
    nClust <- max(as.vector(df.exp$Classification))
    
    #Compute clusters composition in each model (i.e. junction) 
    cluster.counts <- lapply(1:nClust, function(j){
      
      res <- df.exp %>% filter(Classification==j)
      counts <- table(res$Source)
      perc <- as.numeric(counts)/sum(as.numeric(counts))
      perc <- round(perc*100,digits=1)
      
      return(data.frame(Junction = jannot$event_id, Gene = jannot$gene_symbol,
                        Source = factor(names(counts), levels = c("Tumor", "Normal")),
                        Counts = as.numeric(counts),
                        Percent = perc,
                        Cluster = paste0("C", j),
                        stringsAsFactors = FALSE)
      )
      
      
    })
    # rbind clusters for the model
    return(do.call(rbind, cluster.counts))
    
  }, mc.cores = 16)
  
  #rbind for all models
  return(do.call(rbind, main.counts))
}

# build a data frame with cluster mean PSI (% tumor and normal) accross
# all models
#' @export
computeClusterMeanPSI <- function(analysis_name, 
                                  clust_dir = "./clust_output", events){
  
  models <- list.files(clust_dir, pattern = ".RData")
  aux <- gsub(paste0(analysis_name, "_"), "", models)
  junction_ids <- gsub(paste0(".RData"), "", aux)
  
  main.counts <- mclapply(seq_along(models), function(i){
    
    jannot <- dplyr::filter(events, event_id == junction_ids[i])
    # Load df.exp data frame
    # Cluster classification for a given model (i.e junction) 
    load(file = file.path(clust_dir, models[i]))
    nClust <- max(as.vector(df.exp$Classification))
    
    #Compute clusters composition in each model (i.e. junction) 
    cluster.counts <- lapply(1:nClust, function(j){
      
      res <- df.exp %>% filter(Classification==j)
      meanPSI <- mean(res$PSI)
      sdPSI <- sd(res$PSI)
      
      return(data.frame(Junction = jannot$event_id, Gene = jannot$gene_symbol,
                        Cluster = j, meanPSI = meanPSI,sdPSI = sdPSI, 
                        stringsAsFactors = FALSE))
      
    })
    # rbind clusters for the model
    return(do.call(rbind, cluster.counts))
    
  }, mc.cores = 16)
  
  #rbind for all models
  return(do.call(rbind, main.counts))
}


# Filter clusters based on selected parameters
selectClusters <- function(clust, params, clust_dir, analysis_name){
  
  require(dplyr)
  require(doParallel)
  
  res <- clust$composition %>% filter(Source=="Tumor", 
                                      Percent > params$tumor_purity, 
                                      Counts > params$min_size)
  res$Cluster <- as.numeric(gsub("C", "", res$Cluster))
  
  
  main.deltaPSI <- mclapply(seq(1, nrow(res)), function(i){
    
    # Load df.exp data frame
    # Cluster classification for a given model (i.e junction) 
    load(file = file.path(clust_dir, 
                          paste0(analysis_name, "_", res$Junction[i], ".RData")))
    
    tumor.cluster.exp <- df.exp %>% dplyr::filter(Junction==res$Junction[i],
                                           Classification==res$Cluster[i],
                                           Source == "Tumor")
    
    normal.exp <- df.exp %>% dplyr::filter(Junction==res$Junction[i],
                                    Source == "Normal")
    
    if (nrow(tumor.cluster.exp) > 10 || nrow(normal.exp) > 10) {
      deltaPSI <- mean(tumor.cluster.exp$PSI)-mean(normal.exp$PSI)
      test <- wilcox.test(tumor.cluster.exp$PSI, normal.exp$PSI)
      cluster.test <- data.frame(Junction = res$Junction[i], 
                                 Cluster = res$Cluster[i],
                                 meanPSI_T = mean(tumor.cluster.exp$PSI),
                                 meanPSI_N = mean(normal.exp$PSI),
                                 deltaPSI_TvsN = deltaPSI,
                                 absDeltaPSI_TvsN = abs(deltaPSI),
                                 Wilcox_pval = test$p.value,
                                 stringsAsFactors = FALSE)
      
    }else{
      cluster.test <- data.frame(Junction = res$Junction[i], 
                                 Cluster = res$Cluster[i],
                                 meanPSI_T = NaN,
                                 meanPSI_N = NaN,
                                 deltaPSI_TvsN = NaN,
                                 absDeltaPSI_TvsN = NaN,
                                 Wilcox_pval = NaN,
                                 stringsAsFactors = FALSE)
    }
    return(cluster.test)
    
    
  }, mc.cores = 16 )
  
  main.deltaPSI <- do.call(rbind, main.deltaPSI)
  
  res2 <- main.deltaPSI %>% filter(absDeltaPSI_TvsN>params$deltaPSI, 
                                   Wilcox_pval < params$pvalue)
  ts_clusters <- inner_join(res, res2, by = c("Junction", "Cluster") )
  
  return(ts_clusters)
  
}

# build a data frame with cluster mean PSI (% tumor and normal) accross
# all models
gatherClusterExp <- function(analysis_name, 
                             clust_dir = "./clust_output", events, ts_clusters){
  
  models <- list.files(clust_dir, pattern = ".RData")
  aux <- gsub(paste0(analysis_name, "_"), "", models)
  junction_ids <- gsub(paste0(".RData"), "", aux)
  
  idx_keep <- which(junction_ids %in% ts_clusters$Junction)
  models <- models[idx_keep]
  junction_selected <- junction_ids[idx_keep]
  
  list.dfexp <- mclapply(seq_along(models), function(i){
    
    jannot <- dplyr::filter(events, event_id == junction_selected[i])
    # Load df.exp data frame
    # Cluster classification for a given model (i.e junction) 
    load(file = file.path(clust_dir, models[i]))
    
    return(df.exp)
    
  }, mc.cores = 16)
  
  
  #rbind for all models
  return(do.call(rbind, list.dfexp))
}


nonParallel_mixtureClusterAnalysis <- function(analysis_name, PSI_tumor, PSI_control,
                                   events, min_cluster, min_tumor, min_controls,
                                   clust_dir){
  
  require(mclust)
  require(dplyr)
  
  if(dir.exists(clust_dir)){
    unlink(clust_dir, force = T, recursive = T)
  }
  dir.create(clust_dir)
  
  #Run mixture in parallel
  # number of models to fit
  nmodels <- nrow(events)
  
  
  for (i in 1:nmodels) {
    
    se.id <- events$event_id[i]
    gene <- events$gene_symbol[i]
    
    cat("Processing ", i, " ", se.id, "\n")
    
    idx.case <- which(!is.na(PSI_tumor[se.id, ]))
    idx.control <- which(!is.na(PSI_control[se.id, ]))
    
    exp.case <- PSI_tumor[se.id, idx.case, drop = FALSE]
    exp.control <- PSI_control[se.id, idx.control, drop = FALSE]
    
    #samples event is detected
    if(length(exp.case) < min_tumor || length(exp.control) < min_controls){
      next
    }
    
    if(sd(exp.case)==0 || sd(exp.control)==0){
      next
    }
    
    df.exp <- data.frame(Junction = se.id,
                         Gene = gene,
                         PSI = c(as.numeric(exp.case), as.numeric(exp.control)),
                         Source = c(rep("Tumor", length(exp.case)),
                                    rep("Normal", length(exp.control))),
                         Sample = c(colnames(exp.case), colnames(exp.control)),
                         stringsAsFactors = FALSE)
    
    #First fit
    an.error.occured <- FALSE
    tryCatch( {mod.mixed <- mclust::densityMclust(df.exp$PSI, G=1:3) },
              error = function(e) {an.error.occured <<- TRUE} )
    
    if(an.error.occured){
      next
    }
    
    nClust <- mod.mixed$G
    clust.assigned <- levels(factor(mod.mixed$classification))
    nAssignClust <- length(clust.assigned)
    
    if(length(clust.assigned)<nClust){ #empty clusters
      
      newClass <- 1:length(clust.assigned)
      newAssigClass <- rep(NA,nrow(df.exp))
      for(c in 1:length(newClass)){
        newAssigClass[ mod.mixed$classification == clust.assigned[c] ] = c
      }
      df.exp$Classification <- factor(newAssigClass)
      #table(newAssigClass)
      #table(mod.mixed$classification)
      
    }else{
      df.exp$Classification <- factor(mod.mixed$classification)
    }
    
    save(df.exp, file = paste0(clust_dir, analysis_name, "_",
                               se.id, ".RData"))
    
  }
  
} # function nonParallel_mixtureClusteringAnalysis

enrichment_subytpes <- function(clust, 
                                clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                                ts_clusters){
  
  # brca.dat RNA-seq gene level 
  # clinical data of tumors only
  load(clinDataPath)
  clin.data <- as.data.frame(SummarizedExperiment::colData(data))
  clin.data <- clin.data %>% dplyr::filter(shortLetterCode == "TP")
  
  res <- mclapply(seq_along(ts_clusters$Junction), function(i){
    
    cluster.exp <- clust$cluster_exp %>% 
      filter(Junction == ts_clusters$Junction[i],
             Classification == ts_clusters$Cluster[i],
             Source == "Tumor")
    
    
    #Test ER+ status
    #clin.filt <- clin.data %>% dplyr::filter(subtype_ER.Status %in% c("Positive", "Negative"))   
    clin.filt <- clin.data
    clin.pheno <- clin.filt %>% dplyr::filter(subtype_ER.Status == "Positive")
    
    m = nrow(clin.pheno) #has the phenotype in population
    n = nrow(clin.filt) - nrow(clin.pheno) #do not have the phenotype in population
    k = sum(cluster.exp$Sample %in% clin.filt$barcode) #effective cluster size
    x = sum(cluster.exp$Sample %in% clin.pheno$barcode) #effective cluster size with phenotype
    
    test_ERpos <- -1*log10(phyper(x, m, n, k, lower.tail = FALSE, log.p = FALSE))
    expected = m/(m+n)
    seen = x/k
    ratio_cluster_ERpos = paste0(x,"/",k)
    ratio_pop_ERpos = paste0(m, "/", m+n)
    
    #Test HER2+ status
    clin.filt <- clin.data 
    clin.pheno <- clin.filt %>% dplyr::filter(subtype_HER2.Final.Status == "Positive")
    
    m = nrow(clin.pheno) #has the phenotype in population
    n = nrow(clin.filt) - nrow(clin.pheno) #do not have the phenotype in population
    k = sum(cluster.exp$Sample %in% clin.filt$barcode) #effective cluster size
    x = sum(cluster.exp$Sample %in% clin.pheno$barcode) #effective cluster size with phenotype
    
    test_HER2pos <- -1*log10(phyper(x, m, n, k, lower.tail = FALSE, log.p = FALSE))
    expected = m/(m+n)
    seen = x/k
    ratio_cluster_HER2pos = paste0(x,"/",k)
    ratio_pop_HER2pos = paste0(m, "/", m+n)
    
    #Test PR+ status
    clin.filt <- clin.data 
    clin.pheno <- clin.filt %>% dplyr::filter(subtype_PR.Status == "Positive")
    
    m = nrow(clin.pheno) #has the phenotype in population
    n = nrow(clin.filt) - nrow(clin.pheno) #do not have the phenotype in population
    k = sum(cluster.exp$Sample %in% clin.filt$barcode) #effective cluster size
    x = sum(cluster.exp$Sample %in% clin.pheno$barcode) #effective cluster size with phenotype
    
    test_PR2pos <- -1*log10(phyper(x, m, n, k, lower.tail = FALSE, log.p = FALSE))
    expected = m/(m+n)
    seen = x/k
    ratio_cluster_PR2pos = paste0(x,"/",k)
    ratio_pop_PR2pos = paste0(m, "/", m+n)
    
    return(data.frame(Junction = ts_clusters$Junction[i], 
                      Gene = ts_clusters$Gene[i],
                      Cluster = ts_clusters$Cluster[i],
                      ERpos = test_ERpos,
                      ratio_cluster_ERpos = ratio_cluster_ERpos,
                      ratio_pop_ERpos = ratio_pop_ERpos,
                      HER2pos = test_HER2pos,
                      ratio_cluster_HER2pos = ratio_cluster_HER2pos,
                      ratio_pop_HER2pos = ratio_cluster_HER2pos,
                      PRpos = test_PR2pos,
                      ratio_cluster_PR2pos = ratio_cluster_PR2pos,
                      ratio_pop_PR2pos = ratio_pop_PR2pos))
    
    
  }, mc.cores = detectCores())
  
  return(do.call(rbind, res))

} #function enrichment_subtypes

# Enrichment analysis using clinical data
#' @export
enrichment_pam50 <- function(clust, 
                             clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                             ts_clusters){
  
  # brca.dat RNA-seq gene level 
  # clinical data of tumors only
  load(clinDataPath)
  clin.data <- as.data.frame(SummarizedExperiment::colData(data))
  clin.data <- clin.data %>% dplyr::filter(shortLetterCode == "TP")
  
  res <- mclapply(seq_along(ts_clusters$Junction), function(i){
    
    cluster.exp <- clust$cluster_exp %>% 
      filter(Junction == ts_clusters$Junction[i],
             Classification == ts_clusters$Cluster[i],
             Source == "Tumor")
    
    
    #Test PAM50 "Basal-like"
    clin.filt <- clin.data %>% dplyr::filter(subtype_PAM50.mRNA %in% 
                                               levels(clin.data$subtype_PAM50.mRNA))
    clin.pheno <- clin.filt %>% dplyr::filter(subtype_PAM50.mRNA == "Basal-like")
    
    m = nrow(clin.pheno) #has the phenotype in population
    n = nrow(clin.filt) - nrow(clin.pheno) #do not have the phenotype in population
    k = sum(cluster.exp$Sample %in% clin.filt$barcode) #effective cluster size
    x = sum(cluster.exp$Sample %in% clin.pheno$barcode) #effective cluster size with phenotype
    
    test_PAM50_BL <- -1*log10(phyper(x, m, n, k, lower.tail = FALSE, log.p = FALSE))
    expected = m/(m+n)
    seen = x/k
    ratio_cluster_PAM50_BL = paste0(x,"/",k)
    ratio_pop_PAM50_BL = paste0(m, "/", m+n)
    
    #Test PAM50 "Luminal A"
    clin.filt <- clin.data %>% dplyr::filter(subtype_PAM50.mRNA %in% 
                                               levels(clin.data$subtype_PAM50.mRNA))
    clin.pheno <- clin.filt %>% dplyr::filter(subtype_PAM50.mRNA == "Luminal A")
    
    m = nrow(clin.pheno) #has the phenotype in population
    n = nrow(clin.filt) - nrow(clin.pheno) #do not have the phenotype in population
    k = sum(cluster.exp$Sample %in% clin.filt$barcode) #effective cluster size
    x = sum(cluster.exp$Sample %in% clin.pheno$barcode) #effective cluster size with phenotype
    
    test_PAM50_LumA <- -1*log10(phyper(x, m, n, k, lower.tail = FALSE, log.p = FALSE))
    expected = m/(m+n)
    seen = x/k
    ratio_cluster_PAM50_LumA = paste0(x,"/",k)
    ratio_pop_PAM50_LumA = paste0(m, "/", m+n)
    
    #Test PAM50 "Luminal B"
    clin.filt <- clin.data %>% dplyr::filter(subtype_PAM50.mRNA %in% 
                                               levels(clin.data$subtype_PAM50.mRNA))
    clin.pheno <- clin.filt %>% dplyr::filter(subtype_PAM50.mRNA == "Luminal B")
    
    m = nrow(clin.pheno) #has the phenotype in population
    n = nrow(clin.filt) - nrow(clin.pheno) #do not have the phenotype in population
    k = sum(cluster.exp$Sample %in% clin.filt$barcode) #effective cluster size
    x = sum(cluster.exp$Sample %in% clin.pheno$barcode) #effective cluster size with phenotype
    
    test_PAM50_LumB <- -1*log10(phyper(x, m, n, k, lower.tail = FALSE, log.p = FALSE))
    expected = m/(m+n)
    seen = x/k
    ratio_cluster_PAM50_LumB = paste0(x,"/",k)
    ratio_pop_PAM50_LumB = paste0(m, "/", m+n)
    
    #Test PAM50 "HER2-enriched"
    clin.filt <- clin.data %>% dplyr::filter(subtype_PAM50.mRNA %in% 
                                               levels(clin.data$subtype_PAM50.mRNA))
    clin.pheno <- clin.filt %>% dplyr::filter(subtype_PAM50.mRNA == "HER2-enriched")
    
    m = nrow(clin.pheno) #has the phenotype in population
    n = nrow(clin.filt) - nrow(clin.pheno) #do not have the phenotype in population
    k = sum(cluster.exp$Sample %in% clin.filt$barcode) #effective cluster size
    x = sum(cluster.exp$Sample %in% clin.pheno$barcode) #effective cluster size with phenotype
    
    test_PAM50_HER2 <- -1*log10(phyper(x, m, n, k, lower.tail = FALSE, log.p = FALSE))
    expected = m/(m+n)
    seen = x/k
    ratio_cluster_PAM50_HER2 = paste0(x,"/",k)
    ratio_pop_PAM50_HER2 = paste0(m, "/", m+n)
    
    return(data.frame(Junction = ts_clusters$Junction[i], 
                      Gene = ts_clusters$Gene[i],
                      Cluster = ts_clusters$Cluster[i],
                      Basal_Like = test_PAM50_BL,
                      ratio_cluster_PAM50_BL = ratio_cluster_PAM50_BL,
                      ratio_pop_PAM50_BL = ratio_pop_PAM50_BL,
                      LuminalA = test_PAM50_LumA,
                      ratio_cluster_PAM50_LumA = ratio_cluster_PAM50_LumA,
                      ratio_pop_PAM50_LumA = ratio_cluster_PAM50_LumA,
                      LuminalB = test_PAM50_LumB,
                      ratio_cluster_PAM50_LumB = ratio_cluster_PAM50_LumB,
                      ratio_pop_PAM50_LumB = ratio_pop_PAM50_LumB,
                      HER2 = test_PAM50_HER2,
                      ratio_cluster_PAM50_HER2 = ratio_cluster_PAM50_HER2,
                      ratio_pop_PAM50_HER2 = ratio_pop_PAM50_HER2))
    
    
  }, mc.cores = detectCores())
  
  return(do.call(rbind, res))

} #function enrichment_pam50

# Enrichment analysis RNA-seq data
enrichment_RBP_expression <- function(clust, 
                                      clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq_FPKM.rda",
                                      ts_clusters,
                                      log_fc_threshold = log2(1.5),
                                      pval_threshold = 0.05){
  
  require(dplyr)
  require(SummarizedExperiment)
  
  rbp_file <- system.file("extdata", file.path("RBPs","RBP_shortList.txt"), 
                          package = "ts3")
  
  rbp <- read.table(rbp_file, header = TRUE, quote = NULL, sep = "\t")  
  
  # brca.dat RNA-seq gene level 
  # clinical data and FPKM matrix
  load(clinDataPath)
  clin.data <- as.data.frame(SummarizedExperiment::colData(data))
  clin.data <- clin.data %>% dplyr::filter(shortLetterCode == "TP")
  exp.matrix <- SummarizedExperiment::assay(data) #FPKM matrix
  exp.matrix.tp <- exp.matrix[ , colnames(exp.matrix) %in% clin.data$barcode ]
  
  res <- mclapply(seq_along(ts_clusters$Junction), function(i){
  
  #res <- lapply(seq_along(ts_clusters$Junction), function(i){  
    
    cat(ts_clusters$Junction[i], " i=", i, "\n")
    cluster.exp <- clust$cluster_exp %>% 
      dplyr::filter(Junction == ts_clusters$Junction[i],
             Classification == ts_clusters$Cluster[i],
             Source == "Tumor")
    
    cluster.fpkm <- exp.matrix.tp[, colnames(exp.matrix.tp) %in% cluster.exp$Sample]
    remain.fpkm <- exp.matrix.tp[, !colnames(exp.matrix.tp) %in% cluster.exp$Sample]
    
    rbp <- rbp[rbp$ensembl_gene_id %in% rownames(cluster.fpkm), ]
    
    rbp_enrich <- lapply(seq_along(rbp$ensembl_gene_id), function(j){
      
      cat(ts_clusters$Junction[i], " i=", i, "j=", j,  "\n")
      
      mu_cluster <- log2(mean(as.numeric(cluster.fpkm[rbp$ensembl_gene_id[j], ])))
      mu_rem <- log2(mean(as.numeric(remain.fpkm[rbp$ensembl_gene_id[j], ])))
      fc <- mu_cluster - mu_rem
      test <- wilcox.test(as.numeric(cluster.fpkm[rbp$ensembl_gene_id[j], ]),
                          as.numeric(remain.fpkm[rbp$ensembl_gene_id[j], ]))
      
      return(data.frame(Junction = ts_clusters$Junction[i],
                        Gene = ts_clusters$Gene[i],
                        Cluster = ts_clusters$Cluster[i],
                        ensembl_gene_id = rbp$ensembl_gene_id[j],
                        symbol = rbp$symbol[j],
                        log_fc = fc,
                        abs_log_fc = abs(fc),
                        Wilcox_pval = test$p.value))
      
    })
    
    rbp_enrich_df <- do.call(rbind, rbp_enrich)
    return(rbp_enrich_df)  
    
  }, mc.cores = detectCores())
  # })
  
  all_res <- do.call(rbind, res)
  
  all_res_sig <- all_res %>% dplyr::filter(abs_log_fc > log_fc_threshold, 
                                           Wilcox_pval < pval_threshold)
  
} #function enrichment_RBP_expression


reportClusterAnalysis <- function(clust, analysis_name, clusters){
  
  require(ggplot2)
  require(dplyr)
  
  folder <- paste0(analysis_name, "_plots")
  if(dir.exists(folder)){
    unlink(folder, force = T, recursive = T)
  }
  dir.create(folder)
  
  nplot <- nrow(clusters)
  
  for(i in 1:nplot) {
    
    se.id <- clusters$Junction[i]
    gene <- as.character(clusters$Gene[i])
    
    df.exp <- clust$cluster_exp %>% filter(Junction==se.id)
    df.counts <- clust$composition %>% filter(Junction==se.id)
    
    selected.clusters <- clusters %>% filter(Junction==se.id)
    selected.clusters <- gsub("C", "", as.vector(unique(selected.clusters$Cluster)))
    
    tcga_samples <- grepl("^TCGA-", df.exp$Sample, perl = TRUE)
    gtex_samples <- grepl("^SRR", df.exp$Sample, perl = TRUE)
    
    # Pheno TCGA
    pheno <- annotate_TCGA_barcode(as.character(df.exp$Sample[tcga_samples]))
    #as.character(head(df.exp$Sample,100))==as.character(head(pheno$Barcode,100))
    #identical(as.character(pheno$Barcode),as.character(df.exp$Sample))
    pheno$Type <- "TCGA-Tumor"
    pheno$Type [which(pheno$sample_shortLetterCode=="NT")] = "TCGA-Normal"
    pheno$sample_id <- pheno$Barcode
    head(pheno)
    table(pheno$Type)
    
    # Pheno GTEX
    pheno_gtex <- annotate_GTEX_barcode(as.character(df.exp$Sample[gtex_samples]))
    head(pheno_gtex)
    pheno_gtex$Type <- paste0("", pheno_gtex$histological_type)
    pheno_gtex$sample_id <- pheno_gtex$Run
    
    all_pheno <- rbind(pheno %>% dplyr::select(sample_id, Type),
                       pheno_gtex %>% dplyr::select(sample_id, Type))
    
    table(all_pheno$Type)
    
    df.exp <- inner_join(df.exp, all_pheno, by = c("Sample" = "sample_id"))
    table(df.exp$Type)
    
    df.exp$Source2 <- df.exp$Source
    df.exp$Source2[gtex_samples] <- "GTEx"
    df.exp$Source2[df.exp$Type == "TCGA-Normal" ] <- "TCGA-Normal"
    df.exp$Source2[df.exp$Type == "TCGA-Tumor" ] <- "TCGA-Tumor"
    table(df.exp$Source2)
    
    tissue_counts <- df.exp %>% group_by(Type) %>% dplyr::count() %>% as.data.frame()
    
    all_types <- c("TCGA-Tumor", "TCGA-Normal", "Breast", "Adipose Tissue", "Blood",
                   "Blood Vessel",  "Brain", "Colon", "Esophagus",
                   "Heart", "Liver", "Lung", "Pancreas", "Skin")
    
    # Add tissues for which the splice event was not detected
    na_tissues <- all_types[!all_types %in% tissue_counts$Type]
    if(length(na_tissues) > 0){
      na_df <- data.frame(Type = na_tissues, n = rep(0, length(na_tissues)))
      tissue_counts <- rbind(tissue_counts, na_df)
    }
    
    # plotResults(folder,se.id, gene, df.exp, df.counts, tissue_counts,
    #             selected.clusters, clust$rbp_enrich, clust$subtype_enrich,
    #             clust$pam50_enrich) 
    
    # plotResults_simple(folder,se.id, gene, df.exp, df.counts, tissue_counts,
    #             selected.clusters) 

    plotResults_simple_V2(folder,se.id, gene, df.exp, df.counts, tissue_counts,
                       selected.clusters)
    
  }
  
  
}

plotResults <- function(path, se.id, gene, df.exp, df.counts, tissue_counts,
                        selected.clusters, rbp_enrich,
                        subtype_enrich, pam50_enrich){

  require(ggplot2)
  require(cowplot)
  require(ggpubr)
  require(ggsci)
  require(RColorBrewer)
  require(dplyr)
  
  df.exp$Source <- factor(df.exp$Source, levels = c("Tumor", "Normal"))
  palette_npg <- pal_npg("nrc")(6)
  
  #Percentage tumor versus normal in clusters
  perplot <- ggplot(df.counts, aes(x=Cluster, y=Percent, fill = Source)) +
    geom_bar(stat = "identity", width = 0.7) +
    theme(axis.text.x = element_text(size=12, angle = 0, hjust = 0.5),
          axis.text.y = element_text(size=12),
          axis.title.x = element_text(face="plain", colour="black", size=12),
          axis.title.y = element_text(face="plain", colour="black", size=12),
          legend.title=element_blank(),
          legend.text = element_text(face="plain", colour="black", size=12)) +
    #scale_fill_brewer(palette="Dark2") +
    scale_fill_manual(values = palette_npg) +
    ylab("Source (%)") +
    xlab("PSI Clusters")
  
  #Absolute counts tumor versus normal in clusters
  countplot <- ggplot(df.counts, aes(x=Source, y=Counts, fill = Source, label=Counts)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(size = 3, position = position_stack(vjust = 0.5)) +
    theme_classic() +
    #geom_text(size = 3, position = position_dodge(0.5)) +
    theme(axis.text.x = element_text(size=12, angle = 0, hjust = 0.5),
          axis.text.y = element_text(size=8),
          axis.title.x = element_text(face="plain", colour="black", size=12),
          axis.title.y = element_text(face="plain", colour="black", size=12),
          legend.title=element_blank(),
          legend.text = element_text(face="plain", colour="black", size=12)) +
    #scale_fill_brewer(palette="Dark2") +
    scale_fill_manual(values = palette_npg) +
    ylab("Counts") +
    xlab("") +
    facet_grid(Cluster ~ .)
  
  tissues <- length(unique(tissue_counts$Type))
  colors <- colorRampPalette( brewer.pal(11, "BrBG") )(tissues)
  
  res.normal <- tissue_counts %>% filter(!Type %in% "Tumor-TCGA")
  res.tumor <- tissue_counts %>% filter(Type == "Tumor-TCGA")
  sumCounts <- c(sum(res.normal$n),sum(res.tumor$n))
  
  #Absolute counts by tissue type
  tissueplot <- ggplot(tissue_counts, aes(x=Type, y=n, fill = Type, label=n)) +
    geom_bar(stat = "identity", width = 0.7) +
    #geom_text(size = 3, position = position_stack(vjust = 0.5)) +
    #geom_text(size = 3, position = position_dodge(0.5)) +
    theme(axis.text.x = element_text(size=12, angle = 0, hjust = 0.5),
          axis.text.y = element_text(size=8),
          axis.title.x = element_text(face="plain", colour="black", size=12),
          axis.title.y = element_text(face="plain", colour="black", size=12),
          legend.title=element_blank(),
          legend.text = element_text(face="plain", colour="black", size=12)) +
    scale_fill_manual(values = colors) +
    ylab("Counts") +
    xlab("") +
    expand_limits(y=c(0, max(tissue_counts$n)+150)) +
    coord_flip() +
    annotate(geom = "text", x=c(1:nrow(tissue_counts)), y=tissue_counts$n+70,
            label = tissue_counts$n,
            size =3, hjust=0.5, angle=0)
  
  
  
  
  density.p <- ggdensity(df.exp, x = "PSI",
                         fill = "Classification",
                         color = "Classification",
                         #add = "mean",
                         font.label = list(size = 12, color = "black",
                                           style = "plain"),
                         rug = TRUE,
                         ylab = "Density",
                         xlab = "PSI",
                         palette = "Dark2")
  
  
  
  bp_theme <- theme_classic2()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=12),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),
      legend.title=element_blank(),
      legend.text = element_text(face="plain", colour="black", size=12)
      #panel.border = element_blank(),
      #panel.grid=element_blank(),
      #axis.ticks = element_blank(),
      #plot.title=element_text(size=14, face="bold"),
      #legend.title=element_blank(),
      #legend.text = element_text(face="bold", colour="black", size=20)
    )
  my_comparisons <- list( c("Tumor","Normal"))  #boxplot plotting
  bp<-ggboxplot(df.exp, x = "Source", y = "PSI",
                color = "Source",
                palette = palette_npg,
                #fill = "Status",
                xlab = NULL,
                #title = "PSI level - TCGA subsets",
                add = "jitter",
                orientation = "vertical",
                ggtheme = bp_theme,
                facet.by = "Classification")
  # stat_compare_means(comparisons = my_comparisons,
  #                    method = "wilcox.test",
  #                    label = "p.signif") # Add pairwise comparisons p-value
  
  
  #Boxplot selected clusters versus normal
  cluster.exp <- df.exp %>% filter(Source == "Tumor", Classification %in% selected.clusters)
  cluster.exp$Category <- paste0("C", cluster.exp$Classification)
  for (i in 1:nrow(cluster.exp)){
    count <- df.counts %>% filter(Cluster == cluster.exp$Category[i],
                                  Source == "Tumor")
    cluster.exp$Category[i] <- paste0(cluster.exp$Category[i],
                                      " (n=", count$Counts, ")")
  }
  
  normal.exp <- df.exp %>% filter(Source == "Normal")
  normal.exp$Category <- paste0(normal.exp$Type)
  for (i in 1:nrow(normal.exp)){
    count <- tissue_counts %>% filter(Type == normal.exp$Type[i])
    normal.exp$Category[i] <- paste0(normal.exp$Category[i],
                                     " (n=", count$n, ")")
  }
  
  
  bp_theme2 <- theme_classic2()+
    theme(
      axis.title.x = element_text(size=12),
      axis.title.y = element_text(size=12),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=10),
      legend.title=element_blank(),
      legend.text = element_text(face="plain", colour="black", size=12)
    )
  bp.selectedClusters <-ggboxplot(rbind(cluster.exp, normal.exp),
                                  x = "Category", y = "PSI",
                                  color = "Source",
                                  palette = "npg",
                                  #fill = "Status",
                                  #title = "PSI level - TCGA subsets",
                                  add = "jitter",
                                  orientation = "horizontal",
                                  ggtheme = bp_theme2,
                                  repel = T)
  
  
  rbp_enrich_plot <- plot_RBP_enrich(se.id, gene, selected.clusters[1], rbp_enrich)
  # subtype_enrich_plot <- plot_subytpe_enrich(se.id, gene, selected.clusters[1], 
  #                                            subtype_enrich, pam50_enrich)
  
  rbp_enrich_plot <- rbp_enrich_plot + bp_theme2
  
  # Survival analysis
  kmplot <- plot_KM_survival(se.id, gene, df.exp, selected.clusters[1], 
                             clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda")
  
  plot.file <- paste0("./", path, "/",
                      se.id, "_", gene, ".pdf")
  
  page <- cowplot::plot_grid(density.p, perplot, 
                               countplot, bp, bp.selectedClusters, 
                               tissueplot, rbp_enrich_plot, kmplot$plot, 
                               kmplot$table,
                               #labels = c("A", "B", "C","D", "E", "F"),
                               ncol = 2, 
                               nrow = 5,
                               label_size = 4)
    save_plot(plot.file, page,
              ncol = 2, # we're saving a grid plot of 2 columns
              nrow = 4,
              # each individual subplot should have an aspect ratio of 1.3
              base_aspect_ratio = 1.3
    )
  
  
  
  
}

plotResults_simple <- function(path, se.id, gene, df.exp, df.counts, tissue_counts,
                        selected.clusters){
  
  require(ggplot2)
  require(cowplot)
  require(ggpubr)
  require(ggsci)
  require(RColorBrewer)
  require(dplyr)
  require(ggridges)
  
  df.exp$Source <- factor(df.exp$Source, levels = c("Tumor", "Normal"))
  df.exp$Source2 <- factor(df.exp$Source2, levels = c("TCGA-Tumor", "TCGA-Normal",
                                                     "GTEx"))
  palette_npg <- pal_npg("nrc")(6)
  
  tissue_counts$Type[tissue_counts$Type == "Adjacent-TCGA"] <- "TCGA-Adjacent"
  
  #Percentage tumor versus normal in clusters
  perplot <- ggplot(df.counts, aes(x=Cluster, y=Percent, fill = Source)) +
    geom_bar(stat = "identity", width = 0.7) +
    theme(axis.text.x = element_text(size=12, angle = 0, hjust = 0.5),
          axis.text.y = element_text(size=12),
          axis.title.x = element_text(face="plain", colour="black", size=12),
          axis.title.y = element_text(face="plain", colour="black", size=12),
          legend.title=element_blank(),
          legend.text = element_text(face="plain", colour="black", size=12)) +
    #scale_fill_brewer(palette="Dark2") +
    scale_fill_manual(values = palette_npg) +
    ylab("Source (%)") +
    xlab("PSI Clusters")
  
  tissues <- length(unique(tissue_counts$Type))
  colors <- colorRampPalette( brewer.pal(11, "BrBG") )(tissues)
  
  res.normal <- tissue_counts %>% filter(!Type %in% "TCGA-Tumor")
  res.tumor <- tissue_counts %>% filter(Type == "TCGA-Tumor")
  sumCounts <- c(sum(res.normal$n),sum(res.tumor$n))
  
  density.p <- ggdensity(df.exp, x = "PSI",
                         fill = "Classification",
                         color = "Classification",
                         #add = "mean",
                         font.label = list(size = 12, color = "black",
                                           style = "plain"),
                         rug = TRUE,
                         ylab = "Density",
                         xlab = "PSI",
                         palette = "Dark2")
  
  
  
  bp_theme <- theme_classic2()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=12),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),
      legend.title=element_blank(),
      legend.text = element_text(face="plain", colour="black", size=12)
      #panel.border = element_blank(),
      #panel.grid=element_blank(),
      #axis.ticks = element_blank(),
      #plot.title=element_text(size=14, face="bold"),
      #legend.title=element_blank(),
      #legend.text = element_text(face="bold", colour="black", size=20)
    )
  my_comparisons <- list( c("Tumor","Normal"))  #boxplot plotting
  
  
  #Boxplot selected clusters versus normal
  cluster.exp <- df.exp %>% filter(Source == "Tumor", Classification %in% selected.clusters)
  cluster.exp$Category <- paste0("C", cluster.exp$Classification)
  for (i in 1:nrow(cluster.exp)){
    count <- df.counts %>% filter(Cluster == cluster.exp$Category[i],
                                  Source == "Tumor")
    cluster.exp$Category[i] <- paste0("TCGA-Tumor ", cluster.exp$Category[i],
                                      " (n=", count$Counts, ")")
  }
  
  table(cluster.exp$Category)
  
  normal.exp <- df.exp %>% filter(Source == "Normal")
  normal.exp$Type[normal.exp$Type == "Adjacent-TCGA"] <- "TCGA-Adjacent"
  normal.exp$Category <- paste0(normal.exp$Type)
  
  for (i in 1:nrow(normal.exp)){
    count <- tissue_counts %>% filter(Type == normal.exp$Type[i])
    normal.exp$Category[i] <- paste0(normal.exp$Category[i],
                                     " (n=", count$n, ")")
  }
  
  bp_theme2 <- theme_classic2()+
    theme(
      axis.title.x = element_text(size=12),
      axis.title.y = element_text(size=12),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=10),
      legend.title=element_blank(),
      legend.text = element_text(face="plain", colour="black", size=12)
    )
  
  tplot <- tissue_counts %>% filter(n > 10) %>% pull(Type)
  
  bp_jitter_df <- rbind(cluster.exp, 
                        normal.exp %>% filter(Type %in% tplot))
  bp_jitter_df$Category <- factor(bp_jitter_df$Category,
                                  levels = sort(unique(bp_jitter_df$Category)))
  bp.selectedClusters <-ggboxplot(bp_jitter_df,
                                  x = "Category", y = "PSI",
                                  add = "jitter",
                                  outlier.shape = NA,
                                  color = "Source",
                                  palette = "npg",
                                  #fill = "Status",
                                  #title = "PSI level - TCGA subsets",
                                  orientation = "horizontal",
                                  ggtheme = bp_theme2,
                                  repel = T,
                                  xlab = "Tissue")
  
  
  bp_df <- rbind(cluster.exp, 
                 normal.exp %>% filter(Type %in% tplot)) 
  
  bp_noOutlier <- ggplot(bp_df, 
         aes(x=Category, y=PSI, 
                        fill = Source, 
                        color = Source)) +
    #geom_violin(trim = FALSE, alpha = 1, na.rm = TRUE) +
    #geom_boxplot(outlier.shape = NA, na.rm = TRUE, color = "black", fill = "white", width = 0.2) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, na.rm = TRUE, width = 0.5) +
    stat_summary(fun.y=median, geom="point", size=1, color="black", na.rm = TRUE) +
    #scale_y_log10(labels=trans_format('log10', math_format(10^.x))) +
    expand_limits(x=c(-0.2, 1.2)) +
    scale_y_continuous(breaks = seq(0,1, .25), labels = c(0, 0.25, 0.5, 0.75, 1)) +
    theme_bw() +
    theme(axis.text.x = element_text(size=12),
          #axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.x = element_text(face="plain", colour="black", size=12),
          axis.title.y = element_text(face="plain", colour="black", size=12),
          #legend.title=element_blank(),
          panel.background=element_blank(),
          #panel.border=element_blank(),
          panel.grid.major=element_blank(),
          legend.text = element_text(face="plain", colour="black", size=12),
          legend.title = element_text(face="plain", colour="black", size=12),
          legend.position = "right",
          plot.title = element_text(face="plain", colour="black", size=14, hjust = 0)
    ) +
    scale_fill_manual(values = palette_npg) +
    scale_color_manual(values = palette_npg) +
    #ylab(sprintf("%s\n%s", "Minimum junction coverage", "Total RNA-seq counts" )) +
    ylab("PSI") +
    xlab("Tissues") +
    coord_flip()
  
  bp_df2 <- 
    bp_df %>%
    group_by(Category) %>%
    mutate(outlier = PSI > median(PSI) + IQR(PSI) * 1.5 | PSI < median(PSI) - IQR(PSI) * 1.5) %>%
    ungroup
  
  ridges_plot <- ggplot(bp_df2, 
         aes(y=Category, x=PSI, 
             fill = Source, 
             color = Source)) +
    geom_density_ridges(alpha = 0.5, scale = 1, size = 0.8, rel_min_height = 0.001) +
    #theme_ridges()
    #theme_ridges() +
    theme_bw() +
    theme(axis.text.x = element_text(size=12),
          #axis.ticks.y = element_blank(),
          #axis.text.y = element_blank(),
          axis.title.x = element_text(face="plain", colour="black", size=12),
          axis.title.y = element_text(face="plain", colour="black", size=12),
          #legend.title=element_blank(),
          #panel.background=element_blank(),
          #panel.border=element_blank(),
          #panel.grid.major=element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(face="plain", colour="black", size=12),
          legend.title = element_text(face="plain", colour="black", size=12),
          legend.position = "right",
          strip.background = element_rect(color="white", fill = "white"),
          strip.text = element_text(size=12)
    ) +
    scale_x_continuous(limits = c(0, 1.2), breaks = seq(0,1, .25), labels = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_x_continuous(breaks = seq(0,1, .25), labels = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_fill_manual(values = palette_npg) +
    scale_color_manual(values = palette_npg) +
    ylab("Tissues") +
    xlab("PSI") +
    labs(fill = "", color = "", 
         title = "") 
  
    
  
  # Survival analysis
  kmplot <- plot_KM_survival(se.id, gene, df.exp, selected.clusters[1], 
                             clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda")
  
  # Multivariate Coxph Cluster status, Age, Tumor Stage
  forestplot <- plot_Coxph(se.id, gene, df.exp, selected.clusters[1], 
                           clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda")
  
  plot.file <- paste0("./", path, "/",
                      se.id, "_", gene, ".pdf")
  
  page <- cowplot::plot_grid(density.p, 
                             perplot, 
                             #bp.selectedClusters,
                             bp_noOutlier,
                             ridges_plot,
                             kmplot$plot,
                             kmplot$table,
                             forestplot,
                             #labels = c("A", "B", "C","D", "E", "F"),
                             rel_widths = c(1, 1, 1, 1, 1, 1, 2 ),
                             ncol = 2, 
                             nrow = 4,
                             label_size = 4)
  save_plot(plot.file, page,
            ncol = 2, # we're saving a grid plot of 2 columns
            nrow = 4,
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1.3
  )
  
  
  
  
}

plotResults_simple_V2 <- function(path, se.id, gene, df.exp, df.counts, tissue_counts,
                               selected.clusters){
  
  require(ggplot2)
  require(cowplot)
  require(ggpubr)
  require(ggsci)
  require(RColorBrewer)
  require(dplyr)
  require(ggridges)
  
  
  
  # *********PSI ranges of subpopulations - combined
  PSI_density <- ggplot(df.exp, 
                            aes(PSI, 
                            fill = Classification, 
                            color = Classification)) +
    geom_density(alpha = 0.9) +
    #geom_density_ridges(alpha = 0.5, scale = 1, size = 0.8, rel_min_height = 0.001) +
    #theme_ridges()
    #theme_ridges() +
    theme_bw() +
    theme(axis.text.x = element_text(size=12),
          #axis.ticks.y = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.x = element_text(face="plain", colour="black", size=12),
          axis.title.y = element_text(face="plain", colour="black", size=12),
          #legend.title=element_blank(),
          #panel.background=element_blank(),
          #panel.border=element_blank(),
          #panel.grid.major=element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(face="plain", colour="black", size=12),
          legend.title = element_text(face="plain", colour="black", size=12),
          legend.position = "bottom",
          strip.background = element_rect(color="white", fill = "white"),
          strip.text = element_text(size=12)
    ) +
    scale_colour_grey() +
    scale_fill_grey() +
    labs(fill = "Subpopulation", color = "Subpopulation", x = "PSI", y = "Density")
    # scale_x_continuous(limits = c(0, 1), breaks = seq(0,1, .25), labels = c(0, 0.25, 0.5, 0.75, 1))
  
  PSI_density
  
  # *********PSI ranges of subpopulations - one subpoulation per row
  PSI_ridges <- ggplot(df.exp, 
                        aes(y = Classification, x = PSI, 
                            fill = Classification, 
                            color = Classification)) +
    #geom_density(alpha = 0.9) +
    geom_density_ridges(alpha = 0.9, scale = 1, size = 0.8, rel_min_height = 0.001) +
    #theme_ridges()
    #theme_ridges() +
    scale_x_continuous(breaks = seq(0,1, .25), limits = c(-0.1, 1.1 )) +
    #expand_limits(x = c(0,1)) +
    theme_bw() +
    theme(axis.text.x = element_text(size=12),
          #axis.ticks.y = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.x = element_text(face="plain", colour="black", size=12),
          axis.title.y = element_text(face="plain", colour="black", size=12),
          #legend.title=element_blank(),
          #panel.background=element_blank(),
          #panel.border=element_blank(),
          #panel.grid.major=element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(face="plain", colour="black", size=12),
          legend.title = element_text(face="plain", colour="black", size=12),
          legend.position = "bottom",
          strip.background = element_rect(color="white", fill = "white"),
          strip.text = element_text(size=12)
    ) +
    scale_colour_grey() +
    scale_fill_grey() +
    labs(fill = "Subpopulation", color = "Subpopulation", x = "PSI", y = "Subpopulations")
  
  
  PSI_ridges
  
  
  #Boxplot selected clusters versus normal
  df.exp$Source <- factor(df.exp$Source, levels = c("Tumor", "Normal"))
  df.exp$Source2 <- factor(df.exp$Source2, levels = c("TCGA-Tumor", "TCGA-Normal",
                                                      "GTEx"))
  palette_npg <- pal_npg("nrc")(6)
  
  cluster.exp <- df.exp %>% filter(Source == "Tumor")
  cluster.exp$Category <- paste0("S", cluster.exp$Classification)
  
  head(cluster.exp)
  cluster.exp$n_detected <- rep(NA, nrow(cluster.exp))
  table(cluster.exp$Category)
  
  
  for (i in 1:nrow(cluster.exp)){
    count <- df.counts %>% filter(Cluster == paste0("C", cluster.exp$Classification[i]),
                                  Source == "Tumor")
    cluster.exp$Category2[i] <- paste0("TCGA-Tumor ", cluster.exp$Category[i],
                                      " (n=", count$Counts, ")")
    cluster.exp$Category3[i] <- paste0("TCGA-Tumor ", cluster.exp$Category[i])
    cluster.exp$n_detected[i] <- as.numeric(count$Counts)
  }
  
  table(cluster.exp$Category)
  table(cluster.exp$Category2)
  table(cluster.exp$Category3)
  
  normal.exp <- df.exp %>% filter(Source == "Normal")
  normal.exp$Category <- paste0(normal.exp$Type)
  normal.exp$n_detected <- rep(NA, nrow(normal.exp))
  
  for (i in 1:nrow(normal.exp)){
    count <- tissue_counts %>% filter(Type == normal.exp$Type[i])
    normal.exp$Category2[i] <- paste0(normal.exp$Category[i],
                                     " (n=", count$n, ")")
    normal.exp$Category3[i] <- paste0(normal.exp$Category[i])
    normal.exp$n_detected[i] <- as.numeric(count$n)
  }
  head(normal.exp, 20)
  
  table(normal.exp$Category)
  table(normal.exp$Category2)
  table(normal.exp$Category3)
  
  
  
  #tplot <- tissue_counts %>% filter(n > 10) %>% pull(Type)
  # bp_df <- rbind(cluster.exp, 
  #                normal.exp %>% filter(Type %in% tplot)) 
  
  bp_df <- rbind(cluster.exp, normal.exp)
  table(bp_df$Category3)
  
  bp_df2 <-
    bp_df %>%
    group_by(Category) %>%
    mutate(outlier = PSI > median(PSI) + IQR(PSI) * 1.5 | PSI < median(PSI) - IQR(PSI) * 1.5) %>%
    ungroup
  
  order_categories <- c(unique(cluster.exp$Category3), 
                        "TCGA-Normal", "Breast", "Adipose Tissue", "Blood",
                         "Blood Vessel",  "Brain", "Colon", "Esophagus",
                         "Heart", "Liver", "Lung", "Pancreas", "Skin")
  
  bp_df$Category3 <- factor(bp_df$Category3,
                            levels = rev(order_categories))
  
  sample_size <- bp_df %>% select(Category3, n_detected, Source2) %>% distinct()
  absent_tissues <- tissue_counts %>% filter(n == 0) %>% pull(Type)
  
  if(length(absent_tissues) > 0){
    sample_size2 <- rbind(sample_size,
                          data.frame(Category3 = absent_tissues,
                                     n_detected = 0, Source2 = "GTEx")
    )  
  }else{
    sample_size2 = sample_size
  }
  
  
  sample_size2$Category3 <- factor(sample_size2$Category3,
                                   levels = rev(order_categories))
 
  
  bp_noOutlier <- ggplot(bp_df, 
                         aes(x=Category3, y=PSI, 
                             fill = Source2, 
                             color = Source2)) +
    #geom_violin(trim = FALSE, alpha = 1, na.rm = TRUE) +
    #geom_boxplot(outlier.shape = NA, na.rm = TRUE, color = "black", fill = "white", width = 0.2) +
    geom_boxplot(outlier.shape = NA, alpha = 1, na.rm = TRUE, width = 0.5, fill = "white") +
    stat_summary(fun.y=median, geom="point", size=1, color="black", na.rm = TRUE) +
    #scale_y_log10(labels=trans_format('log10', math_format(10^.x))) +
    expand_limits(y=c(0.0, 1)) +
    scale_y_continuous(breaks = seq(0,1, .25), labels = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_x_discrete(drop = FALSE) +
    # geom_text(size = 4, label = paste0("n=(", bp_df$n_detected, ")"), 
    #           aes(y = 1.1), show.legend = FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(size=12),
          #axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.x = element_text(face="plain", colour="black", size=12),
          axis.title.y = element_text(face="plain", colour="black", size=12),
          #legend.title=element_blank(),
          panel.background=element_blank(),
          #panel.border=element_blank(),
          panel.grid.major=element_blank(),
          legend.text = element_text(face="plain", colour="black", size=12),
          legend.title = element_text(face="plain", colour="black", size=12),
          legend.position = "bottom",
          plot.title = element_text(face="plain", colour="black", size=14, hjust = 0)
    ) +
    scale_fill_manual(values = palette_npg) +
    scale_color_manual(values = palette_npg) +
    labs(y = "PSI", x = "Tissues", color = "", fill = "") +
    coord_flip()
  
  bp_noOutlier2 <- bp_noOutlier +
  geom_text(data = sample_size2, 
            size = 4, label = paste0("n=(", sample_size2$n_detected, ")"), 
             aes(y = 1.1), show.legend = FALSE)
  
  
 bp_noOutlier2 
  
 
 # Boxplot horizontal
 bp_df$Category3 <- factor(bp_df$Category3,
                           levels = order_categories)
 
 order_categories2 <- c(unique(cluster.exp$Category2), 
                       "TCGA-Normal", "Breast", "Adipose Tissue", "Blood",
                       "Blood Vessel",  "Brain", "Colon", "Esophagus",
                       "Heart", "Liver", "Lung", "Pancreas", "Skin")
 
 normal_labels <- unique(normal.exp$Category2)
 idx_tcga_normal <- grep("^TCGA-Normal", normal_labels)
 idx_breast_gtex <- grep("^Breast", normal_labels)
 normal_labels_other <- sort(normal_labels[-1*c(idx_breast_gtex, idx_tcga_normal)])
 
 bp_df$Category2 <- factor(bp_df$Category2,
                           levels = c(sort(unique(cluster.exp$Category2)),
                                      normal_labels[idx_tcga_normal],
                                      normal_labels[idx_breast_gtex],
                                      normal_labels_other)
                           )
 bp_horizontal <- ggplot(bp_df, 
                        aes(x=Category2, y=PSI, 
                            fill = Source2, 
                            color = Source2)) +
   #geom_violin(trim = FALSE, alpha = 1, na.rm = TRUE) +
   #geom_boxplot(outlier.shape = NA, na.rm = TRUE, color = "black", fill = "white", width = 0.2) +
   geom_boxplot(outlier.shape = NA, alpha = 1, na.rm = TRUE, width = 0.5, fill = "white") +
   stat_summary(fun.y=median, geom="point", size=1, color="black", na.rm = TRUE) +
   #scale_y_log10(labels=trans_format('log10', math_format(10^.x))) +
   expand_limits(y=c(0.0, 1)) +
   scale_y_continuous(breaks = seq(0,1, .25), labels = c(0, 0.25, 0.5, 0.75, 1)) +
   scale_x_discrete(drop = FALSE) +
   # geom_text(size = 4, label = paste0("n=(", bp_df$n_detected, ")"), 
   #           aes(y = 1.1), show.legend = FALSE) +
   theme_bw() +
   theme(axis.text.x = element_text(size=12, angle = -45, hjust = 0),
         #axis.ticks.x = element_blank(),
         axis.text.y = element_text(size=12),
         axis.title.x = element_text(face="plain", colour="black", size=12),
         axis.title.y = element_text(face="plain", colour="black", size=12),
         #legend.title=element_blank(),
         panel.background=element_blank(),
         #panel.border=element_blank(),
         panel.grid.major=element_blank(),
         legend.text = element_text(face="plain", colour="black", size=12),
         legend.title = element_text(face="plain", colour="black", size=12),
         legend.position = "bottom",
         plot.title = element_text(face="plain", colour="black", size=14, hjust = 0)
   ) +
   scale_fill_manual(values = palette_npg) +
   scale_color_manual(values = palette_npg) +
   labs(y = "PSI", x = "Tissues", color = "", fill = "")
 

 bp_horizontal
  
  
  # ridges_plot <- ggplot(bp_df2, 
  #                       aes(y=Category, x=PSI, 
  #                           fill = Source2, 
  #                           color = Source2)) +
  #   geom_density_ridges(alpha = 0.5, scale = 1, size = 0.8, rel_min_height = 0.001) +
  #   #theme_ridges()
  #   #theme_ridges() +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(size=12),
  #         #axis.ticks.y = element_blank(),
  #         #axis.text.y = element_blank(),
  #         axis.title.x = element_text(face="plain", colour="black", size=12),
  #         axis.title.y = element_text(face="plain", colour="black", size=12),
  #         #legend.title=element_blank(),
  #         #panel.background=element_blank(),
  #         #panel.border=element_blank(),
  #         #panel.grid.major=element_blank(),
  #         panel.grid.minor = element_blank(),
  #         legend.text = element_text(face="plain", colour="black", size=12),
  #         legend.title = element_text(face="plain", colour="black", size=12),
  #         legend.position = "bottom",
  #         strip.background = element_rect(color="white", fill = "white"),
  #         strip.text = element_text(size=12)
  #   ) +
  #   scale_x_continuous(limits = c(0, 1), breaks = seq(0,1, .25), labels = c(0, 0.25, 0.5, 0.75, 1)) +
  #   #scale_x_continuous(breaks = seq(0,1, .25), labels = c(0, 0.25, 0.5, 0.75, 1)) +
  #   scale_fill_manual(values = palette_npg) +
  #   scale_color_manual(values = palette_npg) +
  #   ylab("Tissues") +
  #   xlab("PSI") +
  #   labs(fill = "", color = "", 
  #        title = "") 
  # 
  # ridges_plot
  
  # Survival analysis
  kmplot <- plot_KM_survival_by_subpop(se.id, gene, df.exp, selected.clusters[1], 
                             clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda")
  
  
  
  plot.file <- paste0("./", path, "/",
                      se.id, "_", gene, ".pdf")
  
  plot.file1 <- paste0("./", path, "/",
                      se.id, "_", gene, "_part1.pdf")
  plot.file2 <- paste0("./", path, "/",
                       se.id, "_", gene, "_part2.pdf")
  
  page1 <- cowplot::plot_grid(PSI_density, 
                             PSI_ridges, 
                             kmplot$plot,
                             kmplot$table,
                             #labels = c("A", "B", "C","D", "E", "F"),
                             #rel_widths = c(0.5, 0.5, 1, 1, 1, 1),
                             #base_width = 6,
                             ncol = 2, 
                             nrow = 2,
                             label_size = 20)
  page1
  save_plot(plot.file1, 
            page1,
            ncol = 2, # we're saving a grid plot of 2 columns
            nrow = 2,
            base_width = 5,
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1.3
  )
  
  # page <- cowplot::plot_grid(PSI_density, 
  #                            PSI_ridges, 
  #                            bp_noOutlier,
  #                            ridges_plot,
  #                            kmplot$plot,
  #                            kmplot$table,
  #                            labels = c("A", "B", "C","D", "E", "F"),
  #                            rel_widths = c(0.5, 0.5, 1, 1, 1, 1),
  #                            #base_width = 6,
  #                            ncol = 1, 
  #                            nrow = 6,
  #                            label_size = 20)
  
  page2 <- cowplot::plot_grid(bp_noOutlier2,
                             bp_horizontal,
                             labels = c("A", "B", "C","D", "E", "F"),
                             #rel_widths = c(0.5, 0.5, 1, 1, 1, 1),
                             #base_width = 6,
                             ncol = 1, 
                             nrow = 2,
                             label_size = 20)
  
  page2
  save_plot(plot.file2, 
            page2,
            ncol = 1, # we're saving a grid plot of 2 columns
            nrow = 2,
            base_width = 6,
            # each individual subplot should have an aspect ratio of 1.3
            base_aspect_ratio = 1.3
  )
  
  
  
  
}

plot_Coxph <- function(se.id, gene, df.exp, selected.cluster,
                    clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda"){
  
  require(survival)
  require(parallel)
  require(ggpubr)
  require(survminer)
  require(ggsci)
  require(dplyr)
  
  # clinical data
  load(clinDataPath)
  clin.data <- SummarizedExperiment::colData(data)
  
  cluster.exp <- df.exp %>% 
    filter(Junction == se.id,
           Classification == selected.cluster,
           Source == "Tumor")
  
  # Cluster co-variate
  cluster.status <- rep(0, nrow(clin.data))
  cluster.status[clin.data$barcode %in% cluster.exp$Sample] <- 1
  clin.data$cluster_status <- cluster.status
  
  OS_ind <- rep(0, nrow(clin.data))
  OS_ind[clin.data$vital_status == "dead"] <- 1
  clin.data$OS_ind <- OS_ind
  
  # Age co-variate
  clin.data$age_years <- round(clin.data$age_at_diagnosis/365.25)
  clin.data$age_categ <- ifelse(clin.data$age_years > 60, ">60","<= 60")
  
  # Tumor Stage co-variate
  tumor_stage_num <- rep(NA, length(clin.data$tumor_stage))
  
  idx_stage_1 <- grep("stage i$|stage ia$|stage ib$", 
                      clin.data$tumor_stage, perl = TRUE)
  idx_stage_2 <- grep("stage ii$|stage iia$|stage iib$", 
                      clin.data$tumor_stage, perl = TRUE)
  idx_stage_3 <- grep("stage iii$|stage iiia$|stage iiib$", 
                      clin.data$tumor_stage, perl = TRUE)
  idx_stage_4 <- grep("stage iv$", 
                      clin.data$tumor_stage, perl = TRUE)
  idx_stage_5 <- grep("stage x$", 
                      clin.data$tumor_stage, perl = TRUE)
  
  tumor_stage_num[idx_stage_1] <- "stage_1"
  tumor_stage_num[idx_stage_2] <- "stage_2"
  tumor_stage_num[idx_stage_3] <- "stage_3"
  tumor_stage_num[idx_stage_4] <- "stage_4"
  
  clin.data$tumor_stage_num <- tumor_stage_num
  #table(clin.data$tumor_stage_num)
  
  # Cox model with Cluster, Age, Cancer Stages
  # The final step is to consider all the predictors to determine the effect of
  #each predictor on the hazard rate, while accounting for all other predictors.
  fit <- coxph( Surv(subtype_OS.Time, OS_ind) ~ cluster_status + age_years + tumor_stage_num, 
                data = clin.data)
  
  #summary(fit)
  return(ggforest(fit, data = clin.data, noDigits = 3, fontsize = 0.5))
  
}

plot_KM_survival <- function(se.id, gene, df.exp, selected.cluster,
                             clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                             grayscale = FALSE){
  
  require(survival)
  require(survminer)
  require(dplyr)
  # brca.dat RNA-seq gene level 
  # clinical data
  load(clinDataPath)
  clin.data <- SummarizedExperiment::colData(data)
  
  cluster.exp <- df.exp %>% 
    filter(Junction == se.id,
           Classification == selected.cluster,
           Source == "Tumor")
  
  test <- cluster.exp$Sample %in% clin.data$barcode
  
  cluster.status <- rep(0, nrow(clin.data))
  cluster.status[clin.data$barcode %in% cluster.exp$Sample] <- 1
  clin.data$cluster_status <- cluster.status
  
  OS_ind <- rep(0, nrow(clin.data))
  OS_ind[clin.data$vital_status == "dead"] <- 1
  clin.data$OS_ind <- OS_ind
  
  fit <- survival::survfit(Surv(subtype_OS.Time, OS_ind) ~ cluster_status, data = clin.data)
  
  # clin.data %>% as_data_frame() %>% dplyr::group_by(cluster_status, 
  #                   subtype_Metastasis.Coded) %>% count()
  # 
  # clin.data %>% as_data_frame()  %>% dplyr::group_by(cluster.status) %>% count()
  
  if(grayscale){
    palette_km <- c("grey", "black")
  }else{
    palette_km <- rev(c(pal_npg("nrc")(1), "grey"))
  }
  
  kmplot <- ggsurvplot(fit, data = clin.data,
                       conf.int = F, 
                       pval = T,
                       ggtheme = theme_survminer(),
                       xscale = "d_y",
                       break.time.by = 5*365.25,
                       risk.table = T,
                       risk.table.col = "strata", # Change risk table color by groups
                       palette = palette_km,
                       legend.title = "Stratification", 
                       legend.labs = c("0", "1"),
                       title = paste("Event ID=", se.id, 
                                     gene),
                       font.main = c(16, "bold", "black"),
                       font.legend = c(12, "plain", "black"),
                       xlab = "Time (Years)")
  return(kmplot)
  
}

plot_KM_survival_by_subpop <- function(se.id, gene, df.exp, selected.cluster,
                             clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda"){

  require(survival)
  require(parallel)
  require(ggpubr)
  require(survminer)
  require(ggsci)
  
  
  load(clinDataPath)
  clin.data <- SummarizedExperiment::colData(data)
  
  
  cluster.exp <- df.exp %>%
    filter(Junction == se.id,
           Classification == selected.cluster,
           Source == "Tumor")
  
  # Groups by subpopulation
  all.exp.cluster <- df.exp %>% 
    filter(Junction == se.id,
           Source == "Tumor") %>%
    dplyr::select(Sample, Classification) %>%
    dplyr::rename(barcode = Sample, subpopulation = Classification)
  
  head(all.exp.cluster)
  all.exp.cluster$subpopulation <- paste("S", all.exp.cluster$subpopulation,
                                         sep = "")
  table(all.exp.cluster$subpopulation)
  
  ts_cluster_number = as.numeric(selected.cluster)
  
  # Filter clin.data: Only patients with analyzed by GMM
  clin.data <- dplyr::inner_join(as.data.frame(clin.data), all.exp.cluster,
                                 by = "barcode")
  
  # Filter clin.data: Only patients with survival available
  clin.data <- clin.data %>% dplyr::filter(!is.na(subtype_OS.Time))
  
  # Filter clin.data : Only subpopulations with survival data for at least 20 individuals
  sub_keep <- names(table(clin.data$subpopulation))[table(clin.data$subpopulation) > 20]
  clin.data <- clin.data %>% filter(subpopulation %in% sub_keep)
  
  table(clin.data$subpopulation)
  
  # Count number of subpopulations with clinical data
  Npop <- length(unique(table(clin.data$subpopulation)))
  
  OS_ind <- rep(0, nrow(clin.data))
  OS_ind[clin.data$vital_status == "dead"] <- 1
  clin.data$OS_ind <- OS_ind
  
  
  an.error.occured <- FALSE
  tryCatch( {
    fit <- survival::survfit(Surv(subtype_OS.Time, OS_ind) ~ subpopulation, 
                             data = clin.data)
    
    t <- survdiff(Surv(subtype_OS.Time, OS_ind) ~ subpopulation, data = clin.data)
    
    # Pairwise comparisons between group levels with corrections for multiple testing
    fit_pair <- survminer::pairwise_survdiff(Surv(subtype_OS.Time, OS_ind) ~ subpopulation, 
                                             data = clin.data)
    
    # Obtain median survival for each subpopulation
    median_surv <- surv_median(fit)
    
    
  },
  error = function(e) {an.error.occured <<- TRUE} )
  
  p.val <- 1 - pchisq(t$chisq, length(t$n) - 1)
  
  
  Survival_pval_S1_vs_S2 <- NA 
  Survival_pval_S1_vs_S3 <- NA
  Survival_pval_S2_vs_S3 <- NA
  pairwise_sig = FALSE
  
  if(Npop == 3){
    Survival_pval_S1_vs_S2 <- fit_pair$p.value["S2", "S1"]
    Survival_pval_S1_vs_S3 <- fit_pair$p.value["S3", "S1"]
    Survival_pval_S2_vs_S3 <- fit_pair$p.value["S3", "S2"] 
    
    if(ts_cluster_number == 1){
      if(Survival_pval_S1_vs_S2 < 0.05 & Survival_pval_S1_vs_S3 < 0.05){
        pairwise_sig = TRUE
      }
    }
    
    if(ts_cluster_number == 2){
      if(Survival_pval_S2_vs_S3 < 0.05 & Survival_pval_S1_vs_S2 < 0.05){
        pairwise_sig = TRUE
      }
    }
    
    if(ts_cluster_number == 3){
      if(Survival_pval_S2_vs_S3 < 0.05 & Survival_pval_S1_vs_S3 < 0.05){
        pairwise_sig = TRUE
      }
    }
    
  }
  
  if(Npop == 2 ){
    Survival_pval_S1_vs_S2 <- fit_pair$p.value[1,1]
    
    if(Survival_pval_S1_vs_S2 < 0.05){
      pairwise_sig = TRUE
    }
  }
  
  
  # Significant events: Global p-value < 0.01 and adjusted pairwise comparison < 0.05
  Survival_overall_sig = (p.val < 0.01) & (pairwise_sig == TRUE)
  
  
  #Determine prognosis for subpopulation with differential splicing
  Survival_prognosis <- "NA"
  if(Survival_overall_sig){
    
    
    # Obtain median survival for each subpopulation
    median_surv <- surv_median(fit)
    
    labels <- median_surv$strata
    target_label <- paste0("subpopulation=S",ts_cluster_number)
    other_labels <- labels[!labels %in% target_label]
    
    target_median <- median_surv$median[median_surv$strata %in% target_label]
    other_median <- median_surv$median[median_surv$strata %in% other_labels]
    
    # If one of the groups has not yet dropped to 50% survival at the end of
    # the available data, you cannot compute a median survival 
    if (all(!is.na(c(target_median, other_median)))){
      
      if(all(target_median > other_median)){
        Survival_prognosis <- "Favorable"
      }
      
      if(all(target_median < other_median)){
        Survival_prognosis <- "Unfavorable"
      }
      
    }else{
      Survival_prognosis <- "Check_plot"
    }
    
  } #survival prognosis
  
  kmplot <- ggsurvplot(fit, data = clin.data,
                       conf.int = F,
                       pval = T,
                       xscale = "d_y",
                       break.time.by = 5*365.25,
                       ggtheme = theme_survminer(),
                       risk.table = T,
                       risk.table.col = "strata", # Change risk table color by groups
                       palette = rev(pal_npg("nrc")(5)),
                       legend.title = "",
                       #legend.labs = c("0", "1"),
                       title = sprintf("Event= %s\nGene= %s\nSubpopulation= S%s\nPrognosis=%s",
                                       se.id,
                                       gene,
                                       selected.cluster,
                                       Survival_prognosis),
                       font.main = c(12, "bold", "black"),
                       font.legend = c(12, "plain", "black"),
                       xlab = "Time (Years)")
  
  
  return(kmplot)
  
}

plot_RBP_enrich <- function(se.id, gene, selected.cluster, rbp_enrich){
  
  rbp_cluster <- rbp_enrich %>% filter(Junction == se.id,
                                       Gene == gene,
                                       Cluster %in% selected.cluster) %>%
    arrange(log_fc)
  #remove repetead RPBs - take the best
  idx.map <- match(unique(rbp_cluster$symbol), rbp_cluster$symbol)
  rbp_cluster <- rbp_cluster[idx.map, ]
  
  rbp_cluster$symbol <- factor(rbp_cluster$symbol,
                               levels = rbp_cluster$symbol)
  
  bp_theme <- theme_classic2()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=12),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),
      legend.title=element_blank(),
      legend.text = element_text(face="plain", colour="black", size=12)
      #panel.border = element_blank(),
      #panel.grid=element_blank(),
      #axis.ticks = element_blank(),
      #plot.title=element_text(size=14, face="bold"),
      #legend.title=element_blank(),
      #legend.text = element_text(face="bold", colour="black", size=20)
    )
  
  #Percentage tumor versus normal in clusters
  barplot <- ggplot(rbp_cluster, aes(x=symbol, y=log_fc, fill = log_fc)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_hline(yintercept=0) +
    geom_hline(yintercept=log2(1.5), linetype="dashed") +
    geom_hline(yintercept=-1*log2(1.5), linetype="dashed") +
    #geom_point() +
    scale_fill_gradient2(low="blue", mid = "white", high="red", limits = c(-2.5, 2.5),
                         breaks = seq(-2,2, by=1)) +
    #scale_colour_gradient2(low="blue", mid = "white", high="red") +
    #scale_fill_gradient2(low="blue", high="red") +
    #scale_fill_brewer(palette="RdYlBu") +
    #scale_fill_manual(values = palette_npg) +
    bp_theme +
    ylab("Log2 fold-change") +
    xlab("RBP") +
    expand_limits(y = c(-2, 2)) +
    coord_flip()
  
  return(barplot)  
  
}

# BRCA Subtype and PAM50 enrichment plot
#' @export
plot_subytpe_enrich <- function(se.id, gene, selected.cluster, subtype_enrich,
                                pam50_enrich){
  
  subtype_cluster <- subtype_enrich %>% filter(Junction == se.id,
                                               Gene == gene,
                                               Cluster %in% selected.cluster)
  pam50_cluster <- pam50_enrich %>% filter(Junction == se.id,
                                           Gene == gene,
                                           Cluster %in% selected.cluster)
  
  subtype_cluster <- subtype_cluster[, -1*grep("^ratio", colnames(subtype_cluster))]
  pam50_cluster <- pam50_cluster[, -1*grep("^ratio", colnames(pam50_cluster))]
  
  subytpe_plot <- reshape2::melt(subtype_cluster, id = c("Junction", "Gene", "Cluster"))
  subytpe_plot$Type <- "Receptor Status"
  pam50_plot <- reshape2::melt(pam50_cluster, id = c("Junction", "Gene", "Cluster"))
  pam50_plot$Type <- "PAM50"
  
  comb_plot <- rbind(subytpe_plot, pam50_plot)
  comb_plot <- comb_plot %>% dplyr::arrange(value)
  
  comb_plot$variable <- factor(comb_plot$variable,
                               levels = comb_plot$variable)
  
  bp_theme <- theme_classic2()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=12),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),
      legend.title=element_blank(),
      legend.text = element_text(face="plain", colour="black", size=12)
      #panel.border = element_blank(),
      #panel.grid=element_blank(),
      #axis.ticks = element_blank(),
      #plot.title=element_text(size=14, face="bold"),
      #legend.title=element_blank(),
      #legend.text = element_text(face="bold", colour="black", size=20)
    )
  
  #Percentage tumor versus normal in clusters
  barplot <- ggplot(comb_plot, aes(x=variable, y=value, fill = value)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_hline(yintercept=log10(0.05)*-1, linetype="dashed") +
    #geom_point() +
    scale_fill_gradient2(low="blue", mid = "white", high="red", limits = c(0, 4),
                         breaks = seq(0,4, by=1)) +
    #scale_colour_gradient2(low="blue", mid = "white", high="red") +
    #scale_fill_gradient2(low="blue", high="red") +
    #scale_fill_brewer(palette="RdYlBu") +
    #scale_fill_manual(values = palette_npg) +
    bp_theme +
    ylab("-Log10 P") +
    xlab("Subtype") +
    expand_limits(y = c(0, 5)) +
    coord_flip() +
    facet_grid(. ~ Type)
  
  return(barplot)  
  
}

BRCA_survivalAnalysis <- function(clust, 
                                  clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                                  ts_clusters, analysis_name){
  require(survival)
  require(parallel)
  require(ggpubr)
  require(survminer)
  require(ggsci)
  
  folder <- paste0(analysis_name, "_survival")
  if(dir.exists(folder)){
    unlink(folder, force = T, recursive = T)
  }
  dir.create(folder)
  
  # brca.dat RNA-seq gene level 
  # clinical data
  load(clinDataPath)
  clin.data <- SummarizedExperiment::colData(data)
  
  res <- mclapply(seq_along(ts_clusters$Junction), function(i){
    
    cluster.exp <- clust$cluster_exp %>% 
      filter(Junction == ts_clusters$Junction[i],
             Classification == ts_clusters$Cluster[i],
             Source == "Tumor")
    
    #test <- cluster.exp$Sample %in% clin.data$barcode
    
    cluster.status <- rep(0, nrow(clin.data))
    cluster.status[clin.data$barcode %in% cluster.exp$Sample] <- 1
    clin.data$cluster_status <- cluster.status
    
    OS_ind <- rep(0, nrow(clin.data))
    OS_ind[clin.data$vital_status == "dead"] <- 1
    clin.data$OS_ind <- OS_ind
    
    fit <- survival::survfit(Surv(subtype_OS.Time, OS_ind) ~ cluster_status, data = clin.data)
    
    an.error.occured <- FALSE
    tryCatch( {
      t <- survdiff(Surv(subtype_OS.Time, OS_ind) ~ cluster_status, data = clin.data)
    },
      error = function(e) {an.error.occured <<- TRUE} )
    #print(an.error.occured)
    if(an.error.occured){
      return(data.frame(Gene = ts_clusters$Gene[i], 
                        Junction = ts_clusters$Junction[i],
                        Cluster = ts_clusters$Cluster[i],
                        Survival_pval = NA))
    }
    
    p.val <- 1 - pchisq(t$chisq, length(t$n) - 1)
    
    if (p.val < 0.01){
      plot.file <- paste0(folder, "/Signficant_", ts_clusters$Gene[i], "_",
                          ts_clusters$Junction[i], "_Cluster",
                          ts_clusters$Cluster[i], ".pdf")
      kmplot <- ggsurvplot(fit, data = clin.data,
                           conf.int = F,
                           pval = T,
                           ggtheme = theme_survminer(),
                           risk.table = T,
                           risk.table.col = "strata", # Change risk table color by groups
                           palette = rev(pal_npg("nrc")(2)),
                           legend.title = "Cluster status",
                           legend.labs = c("0", "1"),
                           title = paste("Event ID=", ts_clusters$Junction[i],
                                         ts_clusters$Gene[i]),
                           font.main = c(16, "bold", "black"),
                           font.legend = c(12, "plain", "black"))
      
      # Save the plot
      pdf(file = plot.file)
      print(kmplot)
      dev.off()
    }

    
    
    return(data.frame(Gene = ts_clusters$Gene[i], 
                      Junction = ts_clusters$Junction[i],
                      Cluster = ts_clusters$Cluster[i],
                      Survival_pval = p.val,
                      Survival_N0 = fit$n[1],
                      Survival_N1 = fit$n[2]))
    
    
  }, mc.cores = detectCores())
  
  return(do.call(rbind, res))
  
}

BRCA_survivalAnalysis_by_subpop <- function(clust, 
                                  clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                                  ts_clusters, analysis_name){
  require(survival)
  require(parallel)
  require(ggpubr)
  require(survminer)
  require(ggsci)
  
  folder <- paste0(analysis_name, "_survival")
  if(dir.exists(folder)){
    unlink(folder, force = T, recursive = T)
  }
  dir.create(folder)
  
  # brca.dat RNA-seq gene level 
  # clinical data
  load(clinDataPath)
  
  
  res <- mclapply(seq_along(ts_clusters$Junction), function(i){
    
    clin.data <- SummarizedExperiment::colData(data)
    
    
    cluster.exp <- clust$cluster_exp %>%
      filter(Junction == ts_clusters$Junction[i],
             Classification == ts_clusters$Cluster[i],
             Source == "Tumor")

    #test <- cluster.exp$Sample %in% clin.data$barcode
    
    # cluster.status <- rep(0, nrow(clin.data))
    # cluster.status[clin.data$barcode %in% cluster.exp$Sample] <- 1
    # clin.data$cluster_status <- cluster.status
    
    # Groups by subpopulation
    all.exp.cluster <- clust$cluster_exp %>% 
      filter(Junction == ts_clusters$Junction[i],
             Source == "Tumor") %>%
      dplyr::select(Sample, Classification) %>%
      dplyr::rename(barcode = Sample, subpopulation = Classification)
    
    head(all.exp.cluster)
    all.exp.cluster$subpopulation <- paste("S", all.exp.cluster$subpopulation,
                                           sep = "")
    table(all.exp.cluster$subpopulation)
    
    ts_cluster_number = ts_clusters$Cluster[i]
    
    # Filter clin.data: Only patients with analyzed by GMM
    clin.data <- dplyr::inner_join(as.data.frame(clin.data), all.exp.cluster,
                                  by = "barcode")
    
    # Filter clin.data: Only patients with survival available
    clin.data <- clin.data %>% dplyr::filter(!is.na(subtype_OS.Time))
    
    # Filter clin.data : Only subpopulations with survival data for at least 20 individuals
    sub_keep <- names(table(clin.data$subpopulation))[table(clin.data$subpopulation) > 20]
    clin.data <- clin.data %>% filter(subpopulation %in% sub_keep)
    
    table(clin.data$subpopulation)
    
    # Count number of subpopulations with clinical data
    Npop <- length(unique(table(clin.data$subpopulation)))
    
    # Only perform test if more than one subpopulation
    if(Npop <= 1){
      return(data.frame(Gene = ts_clusters$Gene[i], 
                        Junction = ts_clusters$Junction[i],
                        Cluster = ts_clusters$Cluster[i],
                        Number_clusters = Npop,
                        Survival_global_pval = NA,
                        Survival_number_S1 = NA,
                        Survival_number_S2 = NA,
                        Survival_number_S3 = NA,
                        Survival_pval_S1_vs_S2 = NA,
                        Survival_pval_S1_vs_S3 = NA,
                        Survival_pval_S2_vs_S3 = NA,
                        pairwise_sig = NA,
                        Survival_overall_sig = NA,
                        Survival_prognosis = "NA"
      )
      )
    }
    
    
    OS_ind <- rep(0, nrow(clin.data))
    OS_ind[clin.data$vital_status == "dead"] <- 1
    clin.data$OS_ind <- OS_ind
    
    
    an.error.occured <- FALSE
    tryCatch( {
      fit <- survival::survfit(Surv(subtype_OS.Time, OS_ind) ~ subpopulation, 
                               data = clin.data)
      
      t <- survdiff(Surv(subtype_OS.Time, OS_ind) ~ subpopulation, data = clin.data)
      
      # Pairwise comparisons between group levels with corrections for multiple testing
      fit_pair <- survminer::pairwise_survdiff(Surv(subtype_OS.Time, OS_ind) ~ subpopulation, 
                                               data = clin.data)
      
      # Obtain median survival for each subpopulation
      median_surv <- surv_median(fit)
      
      
    },
    error = function(e) {an.error.occured <<- TRUE} )
    #print(an.error.occured)
    if(an.error.occured){
      return(data.frame(Gene = ts_clusters$Gene[i], 
                        Junction = ts_clusters$Junction[i],
                        Cluster = ts_clusters$Cluster[i],
                        Number_clusters = Npop,
                        Survival_global_pval = NA,
                        Survival_number_S1 = NA,
                        Survival_number_S2 = NA,
                        Survival_number_S3 = NA,
                        Survival_pval_S1_vs_S2 = NA,
                        Survival_pval_S1_vs_S3 = NA,
                        Survival_pval_S2_vs_S3 = NA,
                        pairwise_sig = NA,
                        Survival_overall_sig = Survival_overall_sig,
                        Survival_prognosis = "NA"
                        
                        )
             )
    }
    
    p.val <- 1 - pchisq(t$chisq, length(t$n) - 1)
    
    
    Survival_pval_S1_vs_S2 <- NA 
    Survival_pval_S1_vs_S3 <- NA
    Survival_pval_S2_vs_S3 <- NA
    pairwise_sig = FALSE
    
    if(Npop == 3){
      Survival_pval_S1_vs_S2 <- fit_pair$p.value["S2", "S1"]
      Survival_pval_S1_vs_S3 <- fit_pair$p.value["S3", "S1"]
      Survival_pval_S2_vs_S3 <- fit_pair$p.value["S3", "S2"] 
      
      if(ts_cluster_number == 1){
        if(Survival_pval_S1_vs_S2 < 0.05 & Survival_pval_S1_vs_S3 < 0.05){
          pairwise_sig = TRUE
        }
      }
      
      if(ts_cluster_number == 2){
        if(Survival_pval_S2_vs_S3 < 0.05 & Survival_pval_S1_vs_S2 < 0.05){
          pairwise_sig = TRUE
        }
      }
      
      if(ts_cluster_number == 3){
        if(Survival_pval_S2_vs_S3 < 0.05 & Survival_pval_S1_vs_S3 < 0.05){
          pairwise_sig = TRUE
        }
      }
      
    }
    
    if(Npop == 2 ){
      Survival_pval_S1_vs_S2 <- fit_pair$p.value[1,1]
      
      if(Survival_pval_S1_vs_S2 < 0.05){
        pairwise_sig = TRUE
      }
    }
    
    
    # Significant events: Global p-value < 0.01 and adjusted pairwise comparison < 0.05
    Survival_overall_sig = (p.val < 0.01) & (pairwise_sig == TRUE)
    
   
    #Determine prognosis for subpopulation with differential splicing
    Survival_prognosis <- "NA"
    if(Survival_overall_sig){
      
      
      # Obtain median survival for each subpopulation
      median_surv <- surv_median(fit)
      
      labels <- median_surv$strata
      target_label <- paste0("subpopulation=S",ts_cluster_number)
      other_labels <- labels[!labels %in% target_label]
      
      target_median <- median_surv$median[median_surv$strata %in% target_label]
      other_median <- median_surv$median[median_surv$strata %in% other_labels]
      
      # If one of the groups has not yet dropped to 50% survival at the end of
      # the available data, you cannot compute a median survival 
      if (all(!is.na(c(target_median, other_median)))){
        
        if(all(target_median > other_median)){
          Survival_prognosis <- "Favorable"
        }
        
        if(all(target_median < other_median)){
          Survival_prognosis <- "Unfavorable"
        }
        
      }else{
        Survival_prognosis <- "Check_plot"
      }
      
    } #survival prognosis
    
    if (Survival_overall_sig){
      plot.file <- paste0(folder, "/Signficant_", ts_clusters$Gene[i], "_",
                          ts_clusters$Junction[i], "_Cluster",
                          ts_clusters$Cluster[i], ".pdf")
      
      kmplot <- ggsurvplot(fit, data = clin.data,
                           conf.int = F,
                           pval = T,
                           xscale = "d_y",
                           break.time.by = 5*365.25,
                           ggtheme = theme_survminer(),
                           risk.table = T,
                           risk.table.col = "strata", # Change risk table color by groups
                           palette = rev(pal_npg("nrc")(5)),
                           legend.title = "",
                           #legend.labs = c("0", "1"),
                           title = sprintf("Event= %s\nGene= %s\nSubpopulation= S%s\nPrognosis=%s",
                                           ts_clusters$Junction[i],
                                           ts_clusters$Gene[i],
                                           ts_clusters$Cluster[i],
                                           Survival_prognosis),
                           font.main = c(12, "bold", "black"),
                           font.legend = c(12, "plain", "black"),
                           xlab = "Time (Years)")
      
      # Save the plot
      pdf(file = plot.file)
      print(kmplot)
      dev.off()
    }
    
    
    return(data.frame(Gene = ts_clusters$Gene[i], 
                      Junction = ts_clusters$Junction[i],
                      Cluster = ts_clusters$Cluster[i],
                      Number_clusters = Npop,
                      Survival_global_pval = p.val, #global Survdiff p-value comparing all classes
                      Survival_number_S1 = fit$n[1],
                      Survival_number_S2 = fit$n[2],
                      Survival_number_S3 = fit$n[3],
                      Survival_pval_S1_vs_S2 = Survival_pval_S1_vs_S2,
                      Survival_pval_S1_vs_S3 = Survival_pval_S1_vs_S3,
                      Survival_pval_S2_vs_S3 = Survival_pval_S2_vs_S3,
                      pairwise_sig = pairwise_sig,
                      Survival_overall_sig = Survival_overall_sig,
                      Survival_prognosis = Survival_prognosis
                      )
                      
           )
    
    
  }, mc.cores = detectCores())
  
  return(do.call(rbind, res))
  
}

computeCoxph_BRCA <- function(clust, 
                clinDataPath = "/Users/dveiga/GDCdata/brca_GeneExp_HTSeq.rda",
                              ts_clusters){
  require(survival)
  require(parallel)
  require(survminer)
  
  # clinical data
  load(clinDataPath)
  clin.data <- SummarizedExperiment::colData(data)
  
  # Create co-variates Age, Tumor Stage
  OS_ind <- rep(0, nrow(clin.data))
  OS_ind[clin.data$vital_status == "dead"] <- 1
  clin.data$OS_ind <- OS_ind
  
  # Age co-variate
  clin.data$age_years <- round(clin.data$age_at_diagnosis/365.25)
  clin.data$age_categ <- ifelse(clin.data$age_years > 60, ">60","<= 60")
  
  # Tumor Stage co-variate
  tumor_stage_num <- rep(NA, length(clin.data$tumor_stage))
  
  idx_stage_1 <- grep("stage i$|stage ia$|stage ib$", 
                      clin.data$tumor_stage, perl = TRUE)
  idx_stage_2 <- grep("stage ii$|stage iia$|stage iib$", 
                      clin.data$tumor_stage, perl = TRUE)
  idx_stage_3 <- grep("stage iii$|stage iiia$|stage iiib$", 
                      clin.data$tumor_stage, perl = TRUE)
  idx_stage_4 <- grep("stage iv$", 
                      clin.data$tumor_stage, perl = TRUE)
  idx_stage_5 <- grep("stage x$", 
                      clin.data$tumor_stage, perl = TRUE)
  
  tumor_stage_num[idx_stage_1] <- "stage_1"
  tumor_stage_num[idx_stage_2] <- "stage_2"
  tumor_stage_num[idx_stage_3] <- "stage_3"
  tumor_stage_num[idx_stage_4] <- "stage_4"
  
  clin.data$tumor_stage_num <- tumor_stage_num
  #table(clin.data$tumor_stage_num)
  
  res <- mclapply(seq_along(ts_clusters$Junction), function(i){
    
    cluster.exp <- clust$cluster_exp %>% 
      filter(Junction == ts_clusters$Junction[i],
             Classification == ts_clusters$Cluster[i],
             Source == "Tumor")
    
    clust$cluster_exp %>% 
      filter(Junction == ts_clusters$Junction[i],
             Source == "Tumor") %>%
      dplyr::group_by(Classification) %>%
      dplyr::summarise(mean_psi = mean(PSI))
    
    # Co-variate PSI
    all.exp <- clust$cluster_exp %>% 
      filter(Junction == ts_clusters$Junction[i],
             Source == "Tumor") %>%
      dplyr::select(Sample, PSI) %>%
      dplyr::rename(barcode = Sample)
    
    clin.data <- dplyr::left_join(as.data.frame(clin.data), all.exp, by = "barcode")
    clin.data$subpopulation <- as.numeric(clin.data$subpopulation)
    
    #Co-variate Cluster status
    
    # #numeric
    # cluster.status <- rep(0, nrow(clin.data))
    # cluster.status[clin.data$barcode %in% cluster.exp$Sample] <- 1
    # clin.data$cluster_status <- cluster.status
    
    # categorical by subpopulation
    all.exp.cluster <- clust$cluster_exp %>% 
      filter(Junction == ts_clusters$Junction[i],
             Source == "Tumor") %>%
      dplyr::select(Sample, Classification) %>%
      dplyr::rename(barcode = Sample, subpopulation = Classification)
    
    head(all.exp.cluster)
    table(all.exp.cluster$subpopulation)
    all.exp.cluster$subpopulation <- paste("S", all.exp.cluster$subpopulation,
                                            sep = "")
    clin.data <- dplyr::left_join(as.data.frame(clin.data), all.exp.cluster,
                                  by = "barcode")
    head(clin.data)
    
    # Cox model with Cluster, Age, Cancer Stages
    # The final step is to consider all the predictors to determine the effect of
    #each predictor on the hazard rate, while accounting for all other predictors.
    fit <- coxph( Surv(subtype_OS.Time, OS_ind) ~ subpopulation + age_years + tumor_stage_num, 
                  data = clin.data)
    
    # Coefficients matrix
    coeffs <- coef(summary(fit))
    p.val <- coeffs[1,5]
    HR <- round(coeffs[1,2], 2)
    ggforest(fit)
    
    # Cox model with PSI, Age, Cancer Stages
    # The final step is to consider all the predictors to determine the effect of
    #each predictor on the hazard rate, while accounting for all other predictors.
    fit2 <- coxph( Surv(subtype_OS.Time, OS_ind) ~ PSI + age_years + tumor_stage_num, 
                  data = clin.data)
    
    # Coefficients matrix
    coeffs_2 <- coef(summary(fit2))
    p.val2 <- coeffs_2[1,5]
    HR_2 <- round(coeffs_2[1,2], 2)
    ggforest(fit2)
    
    return(data.frame(Gene = ts_clusters$Gene[i], 
                      Junction = ts_clusters$Junction[i],
                      Cluster = ts_clusters$Cluster[i],
                      Coxph_pval = p.val,
                      Coxph_HazardRatio = HR))
    
    
  }, mc.cores = detectCores())
  
  return(do.call(rbind, res))
  
}


getPatientInClusters <- function(clust, ts_clusters){
  
  require(parallel)
  require(dplyr)
  
  res <- mclapply(seq_along(ts_clusters$Junction), function(i){
    
    cluster.exp <- clust$cluster_exp %>% 
      filter(Junction == ts_clusters$Junction[i],
             Classification == ts_clusters$Cluster[i],
             Source == "Tumor")
    
    PSI <- round(cluster.exp$PSI, 2)
    return(data.frame(Gene = ts_clusters$Gene[i], 
                      Junction = ts_clusters$Junction[i],
                      Cluster = ts_clusters$Cluster[i],
                      Tumor_ids = paste(cluster.exp$Sample, collapse = ","),
                      Tumor_PSI = paste(PSI, collapse = ","))
           )
    
  }, mc.cores = detectCores())
  
  return(do.call(rbind, res))
  
}

getControlsInClusters <- function(clust, ts_clusters){
  
  require(dplyr)
  require(parallel)
  require(tidyr)
  
  res <- mclapply(seq_along(ts_clusters$Junction), function(i){
    
    cat(ts_clusters$Junction[i], " i=", i, "\n")
    
    # Retrieve PSI and sample ids for all controls used during 
    # Gaussian mixture clustering
    df.exp <- clust$cluster_exp %>% 
      filter(Junction == ts_clusters$Junction[i],
             Source == "Normal")
    
    df.counts <- clust$composition %>% 
      filter(Junction == ts_clusters$Junction[i])
    
    tcga_samples <- grepl("^TCGA-", df.exp$Sample, perl = TRUE)
    gtex_samples <- grepl("^SRR", df.exp$Sample, perl = TRUE)                      
    
    n_samples_tcga <- sum(tcga_samples)
    n_samples_gtex <- sum(gtex_samples)
    
    # Pheno TCGA
    if(n_samples_tcga>0){
      pheno <- annotate_TCGA_barcode(as.character(df.exp$Sample[tcga_samples]))
      #as.character(head(df.exp$Sample,100))==as.character(head(pheno$Barcode,100))
      #identical(as.character(pheno$Barcode),as.character(df.exp$Sample))
      pheno$Type <- "Adjacent-TCGA"
      pheno$Type [which(pheno$sample_shortLetterCode=="NT")] = "Adjacent-TCGA"
      pheno$sample_id <- pheno$Barcode
      head(pheno)
      table(pheno$Type)  
    }else{
      pheno <- data.frame(sample_id = NA, Type = NA)
    }
    
    # Pheno GTEX
    if(n_samples_gtex > 0){
      pheno_gtex <- annotate_GTEX_barcode(as.character(df.exp$Sample[gtex_samples]))
      head(pheno_gtex)
      pheno_gtex$Type <- paste0("GTEx-", pheno_gtex$histological_type)
      pheno_gtex$sample_id <- pheno_gtex$Run
      
    }else{
      pheno_gtex <- data.frame(sample_id = NA, Type = NA)
    }
    
    all_pheno <- rbind(pheno %>% dplyr::select(sample_id, Type),
                       pheno_gtex %>% dplyr::select(sample_id, Type)) %>%
      dplyr::filter(!is.na(sample_id))
    

    #table(all_pheno$Type)
    
    df.exp <- inner_join(df.exp, all_pheno, 
                         by = c("Sample" = "sample_id"))
    #table(df.exp$Type)
    df.exp$PSI <- round(df.exp$PSI, 2)
    
    
    controls <- c("Adjacent-TCGA", "GTEx-Adipose Tissue", "GTEx-Blood",
                  "GTEx-Blood Vessel", "GTEx-Brain", "GTEx-Breast", "GTEx-Colon",
                  "GTEx-Esophagus", "GTEx-Heart", "GTEx-Liver", "GTEx-Lung",
                  "GTEx-Pancreas", "GTEx-Skin")
    controls_gtex <- c("GTEx-Adipose Tissue", "GTEx-Blood",
                  "GTEx-Blood Vessel", "GTEx-Brain", "GTEx-Breast", "GTEx-Colon",
                  "GTEx-Esophagus", "GTEx-Heart", "GTEx-Liver", "GTEx-Lung",
                  "GTEx-Pancreas", "GTEx-Skin")
    
    df.exp$Type <- factor(df.exp$Type, levels = controls)
    
    # Tissue counts
    # tissue_counts <- df.exp %>% 
    #   dplyr::group_by(Type, .drop = FALSE) %>%
    #   dplyr::summarise(number = n())
    
    # **** Adjacent-TCGA
    # # # PSI per sample
    # PSI_adjacent_TCGA = paste(df.exp %>% 
    #                             filter(df.exp$Type %in% "Adjacent-TCGA") %>%
    #                             pull(PSI), collapse = ",")
    # ids_adjacent_TCGA = paste(df.exp %>% 
    #                             filter(df.exp$Type %in% "Adjacent-TCGA") %>%
    #                             pull(Sample), collapse = ",")
    
    number_adjacent_TCGA <- df.exp %>% 
      filter(df.exp$Type %in% "Adjacent-TCGA") %>%
      nrow()
    
    mean_PSI_adjacent_TCGA <- df.exp %>% 
      filter(df.exp$Type %in% "Adjacent-TCGA") %>%
      pull(PSI) %>% mean()
    
    
    
    # **** GTEx-Breast
    
    # # PSI per sample
    # PSI_gtex_breast = paste(df.exp %>% 
    #                             filter(df.exp$Type %in% "GTEx-Breast") %>%
    #                             pull(PSI), collapse = ",")
    # ids_gtex_breast = paste(df.exp %>% 
    #                             filter(df.exp$Type %in% "GTEx-Breast") %>%
    #                             pull(Sample), collapse = ",")
    
    number_gtex_breast <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Breast") %>%
      count() %>% as.numeric()
    
    mean_PSI_gtex_breast <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Breast") %>%
      pull(PSI) %>% mean()
    
    # **** GTEx-Adipose Tissue
    number_gtex_adipose <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Adipose Tissue") %>%
      count() %>% as.numeric()
    mean_PSI_gtex_adipose <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Adipose Tissue") %>%
      pull(PSI) %>% mean()
    
    # **** GTEx-Blood
    number_gtex_blood <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Blood") %>%
      count() %>% as.numeric()
    mean_PSI_gtex_blood <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Blood") %>%
      pull(PSI) %>% mean()
    
    # **** GTEx-Blood Vessel
    number_gtex_blood_vessel <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Blood Vessel") %>%
      count() %>% as.numeric()
    mean_PSI_gtex_blood_vessel <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Blood Vessel") %>%
      pull(PSI) %>% mean()
    
    # **** "GTEx-Brain"
    number_gtex_brain <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Brain") %>%
      count() %>% as.numeric()
    mean_PSI_gtex_brain <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Brain") %>%
      pull(PSI) %>% mean()
    
    # **** "GTEx-Colon"
    number_gtex_colon <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Colon") %>%
      count() %>% as.numeric()
    mean_PSI_gtex_colon <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Colon") %>%
      pull(PSI) %>% mean()
    
    # **** "GTEx-Esophagus"
    number_gtex_esophagus <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Esophagus") %>%
      count() %>% as.numeric()
    mean_PSI_gtex_esophagus <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Esophagus") %>%
      pull(PSI) %>% mean()
    
    # **** "GTEx-Heart"
    number_gtex_heart <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Heart") %>%
      count() %>% as.numeric()
    mean_PSI_gtex_heart <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Heart") %>%
      pull(PSI) %>% mean()
    
    # **** "GTEx-Liver"
    number_gtex_liver <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Liver") %>%
      count() %>% as.numeric()
    mean_PSI_gtex_liver <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Liver") %>%
      pull(PSI) %>% mean()
    
    # **** "GTEx-Lung"
    number_gtex_lung <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Lung") %>%
      count() %>% as.numeric()
    mean_PSI_gtex_lung <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Lung") %>%
      pull(PSI) %>% mean()
    
    # **** "GTEx-Pancreas"
    number_gtex_pancreas <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Pancreas") %>%
      count() %>% as.numeric()
    mean_PSI_gtex_pancreas <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Pancreas") %>%
      pull(PSI) %>% mean()
    
    # **** "GTEx-Skin"
    number_gtex_skin <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Skin") %>%
      count() %>% as.numeric()
    mean_PSI_gtex_skin <- df.exp %>% 
      filter(df.exp$Type %in% "GTEx-Skin") %>%
      pull(PSI) %>% mean()
    
    # **** All GTEx
    number_gtex_all <- df.exp %>% 
      filter(df.exp$Type %in% controls_gtex) %>%
      count() %>% as.numeric()
    mean_PSI_gtex_all <- df.exp %>% 
      filter(df.exp$Type %in% controls_gtex) %>%
      pull(PSI) %>% mean()
    
    # **** All controls
    number_controls_all <- df.exp %>% 
      filter(df.exp$Type %in% controls) %>%
      count() %>% as.numeric()
    mean_PSI_controls <- df.exp %>% 
      filter(df.exp$Type %in% controls) %>%
      pull(PSI) %>% mean()
    
    return(data.frame(Gene = ts_clusters$Gene[i], 
                      Junction = ts_clusters$Junction[i],
                      Cluster = ts_clusters$Cluster[i],
                      
                      number_controls_all = number_controls_all,
                      mean_PSI_controls = mean_PSI_controls,
                      
                      # PSI_adjacent_TCGA = PSI_adjacent_TCGA,
                      # ids_adjacent_TCGA = ids_adjacent_TCGA,
                      
                      # PSI_gtex_breast = PSI_gtex_breast,
                      # ids_gtex_breast = ids_gtex_breast,
                      # number_gtex_control = n_samples_gtex,
                      
                      number_adjacent_TCGA = n_samples_tcga,
                      mean_PSI_adjacent_TCGA = mean_PSI_adjacent_TCGA,
                      
                      number_gtex_adipose = number_gtex_adipose,
                      mean_PSI_gtex_adipose = mean_PSI_gtex_adipose,
                      
                      number_gtex_blood = number_gtex_blood,
                      mean_PSI_gtex_blood = mean_PSI_gtex_blood,
                      
                      number_gtex_blood_vessel = number_gtex_blood_vessel,
                      mean_PSI_gtex_blood_vessel = mean_PSI_gtex_blood_vessel,
                      
                      number_gtex_brain = number_gtex_brain,
                      mean_PSI_gtex_brain = mean_PSI_gtex_brain,
                      
                      number_gtex_breast = number_gtex_breast,
                      mean_PSI_gtex_breast = mean_PSI_gtex_breast,
                      
                      number_gtex_colon = number_gtex_colon,
                      mean_PSI_gtex_colon = mean_PSI_gtex_colon,
                      
                      number_gtex_esophagus = number_gtex_esophagus,
                      mean_PSI_gtex_esophagus = mean_PSI_gtex_esophagus,
                      
                      number_gtex_heart = number_gtex_heart,
                      mean_PSI_gtex_heart = mean_PSI_gtex_heart,
                      
                      number_gtex_liver = number_gtex_liver,
                      mean_PSI_gtex_liver = mean_PSI_gtex_liver,
                      
                      number_gtex_lung = number_gtex_lung,
                      mean_PSI_gtex_lung = mean_PSI_gtex_lung,
                      
                      number_gtex_pancreas = number_gtex_pancreas,
                      mean_PSI_gtex_pancreas = mean_PSI_gtex_pancreas,
                      
                      number_gtex_skin = number_gtex_skin,
                      mean_PSI_gtex_skin = mean_PSI_gtex_skin,
                      
                      number_gtex_all = number_gtex_all,
                      mean_PSI_gtex_all = mean_PSI_gtex_all
            )
    )
    
   
    
  }, mc.cores = detectCores())
  
  return(do.call(rbind, res))
  
}

getSuppaTranscripts <- function(ioe_file){
  
  ioe <- read.table(ioe_file, header = TRUE, 
                    sep = "\t", stringsAsFactors = FALSE)
  
  # Create alternative2 transcripts for all events
  alternative2_tr <- rep(NA, nrow(ioe))
  
  for (i in 1:nrow(ioe)) {
    
    # alternative transcripts are the alternative1 transcripts (numerator of the PSI formula)
    # total transcripts are the denominator of the PSI formula
    # alternative 2 transcripts are the difference between alternative1 and total transcritps
    skip <- setdiff(strsplit(ioe$total_transcripts[i], ",")[[1]], 
                    strsplit(ioe$alternative_transcripts[i], ",")[[1]])
    
    alternative2_tr[i] <- paste0(skip, collapse = ",")
    
  }
  
  ioe$alternative2_transcripts <- alternative2_tr
  return(ioe)
  
}

getTranscriptomeOrigin <- function(clusters){
  
  require(dplyr)
  
  event_type <- rep("NA", nrow(clusters))
  for (i in 1:nrow(clusters)) {
    
    # clusters$alternative_transcripts[i]
    # clusters$alternative2_transcripts[i]
    # clusters$total_transcripts[i]
    if(grepl("ENST", clusters$alternative_transcripts[i]) &
       grepl("ENST", clusters$alternative2_transcripts[i]) ){
      
      event_type[i] <- "GENCODE"
    }else{
      event_type[i] <- "LR-seq"
    }
    
  }
  return(event_type)
  
}

getNMDprediction <- function(clusters){
  
  require(dplyr)
  require(stringr)
  require(rtracklayer)
  
  # PacBio NMD predictions
  load("./ORF_effect/df_PTC_effect.Rd")
  head(df_PTC_effect)
  
  # GENCODE NMD transcripts
  genc_gtf <- import.gff("~/tools/Gencode/gencode.v30.annotation.sorted.gff.gz")
  #genc_tr <- genc_gtf %>% as.data.frame() %>% filter(type == "transcript")
  nmd_ids <- genc_gtf %>% as.data.frame() %>%
    filter(transcript_type == "nonsense_mediated_decay") %>%
    pull(transcript_id) %>% base::unique()
  
  genc_NMD <- data.frame(PBid = unique(genc_gtf$transcript_id),
                         PTC_status = "no",
                         stringsAsFactors = FALSE)
  
  idx_nmd <- base::match(nmd_ids, genc_NMD$PBid)
  sum(is.na(idx_nmd))

  genc_NMD$PTC_status[idx_nmd] <- "yes"
  table(genc_NMD$PTC_status)
  
  #combine PacBio and Gencode predictions
  df_PTC_effect <- rbind(df_PTC_effect, genc_NMD)
  
  ptc1 <- rep("NA", nrow(clusters))
  ptc2 <- rep("NA", nrow(clusters))
  ptc1_summary <- rep("NA", nrow(clusters))
  ptc2_summary <- rep("NA", nrow(clusters))
  
  for (i in 1:nrow(clusters)) {
    
    alt1 <- str_split(clusters$alternative_transcripts[i], ",")[[1]]
    
    df1 <- data.frame(PBid = alt1, stringsAsFactors = FALSE)
    df1 <- left_join(df1, df_PTC_effect, by = c("PBid"))
    ptc1[i] <- paste0(df1$PTC_status, collapse = ",")
    
    if(grepl("yes", ptc1[i])){
      ptc1_summary[i] <- 1
    }else{
      ptc1_summary[i] <- 0
    }
    
    #****** alt 2
    alt2 <- str_split(clusters$alternative2_transcripts[i], ",")[[1]]
    
    df2 <- data.frame(PBid = alt2, stringsAsFactors = FALSE)
    df2 <- left_join(df2, df_PTC_effect, by = c("PBid"))
    ptc2[i] <- paste0(df2$PTC_status, collapse = ",")
    
    if(grepl("yes", ptc2[i])){
      ptc2_summary[i] <- 1
    }else{
      ptc2_summary[i] <- 0
    }
    
  } #for
  
  return(data.frame(ptc1 = ptc1, 
                    ptc1_summary = ptc1_summary,
                    ptc2 = ptc2,
                    ptc2_summary = ptc2_summary
        )
  )
  
}
