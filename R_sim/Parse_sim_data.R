######################################################
##  Parsing raw data from simulations
##    Only shown for example data
######################################################

#----------------------------------------------------------
#  Load packages
#----------------------------------------------------------

library(data.table)
library(parallel)
library(SNPRelate)
library(factoextra)


#----------------------------------------------------------
# Helper: unpack / clean up a tarball
#----------------------------------------------------------
out_folder <- "./"

unpack_sim <- function(file_gz,n_comp=0) {
  
  if(!file.exists("./tmp/")) dir.create("./tmp/")
  folder <- sub("\\.tar\\.gz$", "", file_gz)
  folder <- basename(folder)
  folder <- paste0("./tmp/","",folder)
  
  if(!file.exists(folder)){
    system(paste0("tar -xvzf ", file_gz, " --strip-components=",n_comp, " -C ./tmp/")) ## on a mac you need the strip components in order not to recreate the entire path again
  } 
  
  return(folder)
  
}

cleanup_sim <- function(base_name) {
  system(paste("rm -r", paste0("./tmp/",base_name)))
}

#----------------------------------------------------------
# Core parser
#----------------------------------------------------------

parse_raw_data <- function(file_gz,
                           outfolder,
                           subfolder,
                           min_maf  = 0.05,
                           side,
                           keep,
                           cores = 4) {
  
  
  # ---------------------------- #
  # 1) Unpack and and prepare data files
  # ---------------------------- #
  
  folder <- unpack_sim(file_gz,n_comp = 4) ## unpack data into "tmp" folder
  base_name <- basename(folder)
  
  files <- list.files(folder, recursive = TRUE, full.names = TRUE)
  
  if(length(files)==0){
    folder <- unpack_sim(file_gz,n_comp = 0)
    files <- list.files(folder, recursive = TRUE, full.names = TRUE)
  } 
  
  
  message("Working on ", base_name)
  
  
  ## Parse parameter info from folder name
  params <- data.table(t(strsplit(base_name, "_", fixed = TRUE)[[1]]))[, 1:4]
  setnames(params, c("chr", "V", "c", "rep"))
  
  ## get files
  
  ## Check that simulation produced output
  if (!any(grepl("out", files))) {
    message("Simulation did not work\n\n")
    cleanup_sim(folder)
    return(invisible(NULL))
  }
  
  # ---------------------------- #
  # 2) Sim summary & map
  # ---------------------------- #
  
  sim_file <- files[grep(paste0(base_name, ".txt"), files, fixed = TRUE)][1]
  sim_data <- fread(sim_file)
  
  map <- readRDS(file.path(folder, "rec_map.rds"))
  map[, marker := paste(Chr, bp, sep = ":")]
  map <- cbind(params, map)
  
  # Chr_9sp is one of the 10 recombination maps used for simulations. There are two of them, one that contains QTN and one that does not
  setnames(map, "chr", "Chr_9sp") 
  map[, Chr_9sp := as.numeric(gsub("chr", "", Chr_9sp))]
  map[, Chr_type := ifelse(Chr == 1, "QTN", "ntrl")]
  map <- map[, .(V, c, rep, Chr_9sp, Chr, bp, marker, cM, pos_nemo, Chr_type, type, allelic_values)]
  
  # ---------------------------- #
  # 3) Environment / selection surface
  # ---------------------------- #
  
  env_raw <- scan(file.path(folder, "env.txt"), what = character(), quiet = TRUE)
  env_raw <- strsplit(env_raw,c("}{"),fixed=TRUE)[[1]]
  env_vals <- as.numeric(gsub("}}","",gsub("{{","",env_raw,fixed = TRUE),fixed = TRUE))
  
  env <- data.table(expand.grid(x = 1:side, y = side:1))
  env[, pop := 1:(side * side)]
  env[, env := env_vals]
  
  # ---------------------------- #
  # 4) Parse neutral genotypes
  # ---------------------------- #
  
  ntrl_files   <- list.files(file.path(folder, "out", "ntrl"),   full.names = TRUE)
  
  hapl_ntrl <- fread(ntrl_files[1])
  
  env_ind <- env[match(hapl_ntrl$pop, pop)]
  env_ind[, indx := .I]
  
  hapl_ntrl <- hapl_ntrl[keep]
  ID_ntrl   <- hapl_ntrl$ID
  env_ind   <- env_ind[keep]
  
  # Keep only neutral loci (two haplotype cols per locus)
  n_ntrl   <- sum(map$type == "ntrl")
  hapl_ntrl <- as.matrix(hapl_ntrl[, 2:(2 * n_ntrl + 1)])
  
  nInd   <- nrow(hapl_ntrl)
  
  loci_ntrl <- split(seq_len(ncol(hapl_ntrl)), rep(seq_len(ncol(hapl_ntrl) / 2), each = 2))
  
  names(loci_ntrl) <- map[type == "ntrl", marker]
  
  # get maf
  maf_ntrl <- unlist(mclapply(loci_ntrl, function(l) {
    maf <- sum(hapl_ntrl[, l] - 1) / (nInd * 2)
    min(maf, 1 - maf)
  }, mc.cores = cores))
  
  map[type == "ntrl", maf := maf_ntrl]
  
  # get genotypes
  GTs_ntrl <- do.call(cbind, lapply(loci_ntrl[maf_ntrl > min_maf], function(l) {
    rowSums(hapl_ntrl[, l] - 1)
  }))
  
  # ---------------------------- #
  # 5) QTN genotypes & phenotype
  # ---------------------------- #
  quanti_files <- list.files(file.path(folder, "out", "quanti"), full.names = TRUE)
  quanti <- fread(quanti_files[1])
  quanti <- quanti[match(ID_ntrl, ID)]
  
  env_ind[, pheno := quanti[match(ID_ntrl, ID), P1]]
  
  ColNames <- grep("l", colnames(quanti), value = TRUE)
  tmp      <- do.call(rbind, strsplit(do.call(rbind, strsplit(ColNames, "t"))[, 2], "l"))
  
  trait <- as.numeric(tmp[, 1])
  locus <- as.numeric(gsub("[xy]", "", tmp[, 2]))
  
  hapl_QTN <- as.matrix(quanti[, ..ColNames])
  loci_QTN <- split(seq_len(ncol(hapl_QTN)), rep(seq_len(ncol(hapl_QTN) / 2), each = 2))
  names(loci_QTN) <- map[type == "QTN", marker]
  
  # Locus evolvability
  Va_l <- sapply(loci_QTN, function(l) var(rowSums(hapl_QTN[, l])))
  evolv <- Va_l / mean(env_ind$pheno)^2
  
  # Phase QTN genotypes relative to first haplotype
  hapl_QTN_bin <- do.call(cbind, lapply(loci_QTN, function(l) {
    gt <- hapl_QTN[, l, drop = FALSE]
    focal_allele <- gt[1, 1]
    ifelse(gt == focal_allele, 0, 1)
  }))
  
  # maf
  maf_QTN <- sapply(loci_QTN, function(l) {
    maf <- sum(hapl_QTN_bin[, l]) / (nInd * 2)
    min(maf, 1 - maf)
  })
  
  # genotypes
  GTs_QTN <- do.call(cbind, lapply(loci_QTN[maf_QTN > min_maf], function(l) {
    rowSums(hapl_QTN_bin[, l])
  }))
  
  # add to map
  map[type == "QTN", maf := maf_QTN]
  map[type == "QTN", ev  := evolv]
  
  # ---------------------------- #
  # 6) Combine ntrl and QTN data sets
  # ---------------------------- #
  
  # Filter map to loci with sufficient MAF
  map <- rbind(
    map[type == "QTN" & maf > min_maf],
    map[type == "ntrl" & maf > min_maf]
  )
  setorder(map, Chr_9sp, Chr_type, bp)
  
  GTs <- matrix(NA_real_, nrow = nrow(GTs_ntrl), ncol = nrow(map))
  colnames(GTs) <- map$marker
  
  GTs[, map$type == "QTN"]   <- GTs_QTN
  GTs[, map$type == "ntrl"]  <- as.matrix(GTs_ntrl)
  
  # Remove duplicated markers (if any)
  dup <- duplicated(map$marker)
  GTs <- GTs[, !dup, drop = FALSE]
  map <- map[!dup]
  
  # ---------------------------- #
  # 7) Save and clean up
  # ---------------------------- #
  out_file <- paste0(file.path(outfolder,
                               paste(params, collapse = "_")), ".rds")
  saveRDS(
    list(
      GTs      = GTs,
      map      = map,
      env      = env_ind,
      sim_data = sim_data
    ),
    out_file
  )
  
  message("removing temporary data\n\n")
  cleanup_sim(base_name)
  
  invisible(NULL)
}

#----------------------------------------------------------
# Parse data and save to new location
#----------------------------------------------------------

raw_folder <- "./Nemo/" ## contains only one example file
files_gz   <- list.files(raw_folder, full.names = TRUE)
files_gz   <- files_gz[grepl("\\.gz$", files_gz)]

keep <- readRDS("./data/keep_500.rds") ## we first keep only 500 samples, this will be reduced to 125 for the final analyses. These individuals are fixed across all simulations!

## run parser for all data sets
# file_gz <- files_gz[1]
lapply(files_gz, function(f) {
  parse_raw_data(
    file_gz   = f,
    outfolder = "./parsed_data",
    min_maf   = 0.05,
    side      = 48,
    keep      = keep,
    cores     = 4
  )
})

system("rm -r tmp")


