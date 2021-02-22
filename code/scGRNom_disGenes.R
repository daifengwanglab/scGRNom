scGRNom_disGenes <- function(df, gwas_snps, extension_bps){
  
  library(GenomicRanges)
  library(biomaRt)
  #library(ReactomePA)
  #library(clusterProfiler)
  #library(UpSetR)
  library(stringr)
  #library(dplyr)
  #library(EnsDb.Hsapiens.v86)
  #library(MotifDb)
  library(motifbreakR)
  #library(tidyverse)
  library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
  library(BSgenome.Hsapiens.UCSC.hg19)

  colnames(gwas_snps) = c('CHR', 'BP', 'SNP')
  
  tfbs_setup_results_mg = tfbs_setup(df, gwas_snps,extension_bps)
  mg_en_mb = tfbs_interruption(tfbs_setup_results_mg[[1]])
  mg_pr_mb = tfbs_interruption(tfbs_setup_results_mg[[2]])
  
  mg_broken_tfbs = ct_tfbs_interruption(df, mg_en_mb, mg_pr_mb)
  disease_genes = data.frame(unique(c(as.character(mg_broken_tfbs[[1]]$TG),as.character(mg_broken_tfbs[[2]]$TG))))
  colnames(disease_genes) = 'disease_genes'
  return(disease_genes)
  
}


snps.from.rsid.update <- function(rsid = NULL, dbSNP = NULL,
                                  search.genome = NULL) {
  if (is.null(rsid)) {
    stop("no RefSNP IDs have been found, please include RefSNP ID numbers")
  }
  if (!inherits(dbSNP, "SNPlocs")) {
    stop(paste0("dbSNP argument was not provided with a valid SNPlocs object.\n",
                "Please run availible.SNPs() to check for availible SNPlocs"))
  }
  if (!inherits(search.genome, "BSgenome")) {
    stop(paste0("search.genome argument was not provided with a valid BSgenome object.\n",
                "Run availible.genomes() and choose the appropriate BSgenome object"))
  }
  if (all(!grepl("rs", rsid))) {
    bad.names <- rsid[!grepl("rs", rsid)]
    stop(paste(paste(bad.names, collapse = " "), "are not rsids, perhaps you want to import your snps from a bed or vcf file with snps.from.file()?"))
  }
  rsid <- unique(rsid)
  rsid.grange <- as(snpsById(dbSNP, rsid, ifnotfound = "warning"), "GRanges")
  rsid.grange <- change.to.search.genome(rsid.grange, search.genome)
  rsid.grange <- GRanges(rsid.grange)
  rsid.refseq <- getSeq(search.genome, rsid.grange)
  rsid.grange$UCSC.reference <- as.character(rsid.refseq)
  rsid.grange <- sapply(split(rsid.grange, rsid.grange$RefSNP_id), function(snp) {
    alt.allele <- determine.allele.from.ambiguous(snp$alleles_as_ambig, snp$UCSC.reference)
    if (length(alt.allele) == 0){
      if (length(alt.allele) > 1L) {
        snp <- do.call("c", replicate(length(alt.allele), snp))
        snp$UCSC.alternate <- NA
        names(snp) <- paste(snp$RefSNP_id, alt.allele, sep = ":")
      } else {
        snp$UCSC.alternate <- NA
        names(snp) <- snp$RefSNP_id
      }
    }
    else{
      if (length(alt.allele) > 1L) {
        snp <- do.call("c", replicate(length(alt.allele), snp))
        snp$UCSC.alternate <- alt.allele
        names(snp) <- paste(snp$RefSNP_id, alt.allele, sep = ":")
      } else {
        snp$UCSC.alternate <- alt.allele
        names(snp) <- snp$RefSNP_id
      }
    }
    
    return(snp)
  })
  
  
  rsid.grange <- unlist(do.call("GRangesList", rsid.grange), use.names = FALSE)
  colnames(mcols(rsid.grange)) <- c("RefSNP_id", "alleles_as_ambig", "REF", "ALT")
  rsid.grange = rsid.grange[!is.na(rsid.grange$ALT),]
  
  
  rsid.grange$REF <- DNAStringSet(rsid.grange$REF)
  rsid.grange$ALT <- DNAStringSet(rsid.grange$ALT)
  # rsid.grange$alleles_as_ambig <- DNAStringSet(rsid.grange$alleles_as_ambig)
  rsid.grange$alleles_as_ambig <- NULL
  colnames(mcols(rsid.grange))[1] <- "SNP_id"
  attributes(rsid.grange)$genome.package <- attributes(search.genome)$pkgname
  return(rsid.grange)
}



change.to.search.genome <- function(granges.object, search.genome) {
  sequence <- seqlevels(granges.object)
  ## sequence is in UCSC format and we want NCBI style
  newStyle <- mapSeqlevels(sequence,seqlevelsStyle(search.genome))
  newStyle <- newStyle[complete.cases(newStyle)] # removing NA cases.
  ## rename the seqlevels
  granges.object <- renameSeqlevels(granges.object,newStyle)
  seqlevels(granges.object) <- seqlevelsInUse(granges.object)
  seqinfo(granges.object) <- keepSeqlevels(seqinfo(search.genome),
                                           value = seqlevelsInUse(granges.object))
  return(granges.object)
}





determine.allele.from.ambiguous <- function(ambiguous.allele, known.allele) {
  neucleotide.ambiguity.code <- list(Y = c("C", "T"), R = c("A", "G"), W = c("A", "T"),
                                     S = c("G", "C"), K = c("T", "G"), M = c("C", "A"),
                                     D = c("A", "G", "T"), V = c("A", "C", "G"),
                                     H = c("A", "C", "T"), B = c("C", "G", "T"),
                                     N = c("A", "C", "G", "T"))
  specnac <- neucleotide.ambiguity.code[[ambiguous.allele]]
  unknown.allele <- specnac[-grep(known.allele, specnac)]
  return(unknown.allele)
}




# tfbs_setup:
#   - Sets up overlapping GRanges for SNPs between disease GWAS and GRN
tfbs_setup = function(ct_grn, gwas,extension) {
  
  # Remove NAs
  gwas = na.omit(gwas)
  # Convert to GRanges
  gwas$CHR = as.character(gwas$CHR)
  gr = GRanges(seqnames = gwas$CHR,
               IRanges(start = gwas$BP, end = gwas$BP))
  # Overlap GWAS with cell type enhancer positions
  gwas_index = grn_overlap(ct_grn, gr,extension)
  # Overlap GWAS with cell type promoter positions
  
  gwas_snps_en = snps.from.rsid.update(rsid = as.character(gwas[unlist(gwas_index[[3]]), ]$SNP),
                                       dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37,
                                       search.genome = BSgenome.Hsapiens.UCSC.hg19)
  # Get SNPs from RSID for promoter
  gwas_snps_pr = snps.from.rsid.update(rsid = as.character(gwas[unlist(gwas_index[[4]]), ]$SNP),
                                       dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37,
                                       search.genome = BSgenome.Hsapiens.UCSC.hg19)
  total_setup = list(
    gwas_snps_en,
    gwas_snps_pr)
  return(total_setup)
}

# tfbs_interruption:
#   - Performs TFBS interruption analysis using motifbreakR
tfbs_interruption = function(rsid_list) {
  temp_motifbreakR_result = motifbreakR(
    snpList = rsid_list,
    filterp = TRUE,
    pwmList = hocomoco,
    threshold = 1e-4,
    method = "ic",
    bkg = c(
      A = 0.25,
      C = 0.25,
      G = 0.25,
      T = 0.25
    )
  )
  return(temp_motifbreakR_result)
}



# grn_setup:
#   - Extracts enhancer and promoter information from given cell type GRN
grn_setup = function(ct_grn) {
  # Convert to DF
  curr_grn_df = data.frame(ct_grn)
  # Extract enhancer/promoter chromosomes
  curr_grn_df$en_chrs = str_extract(curr_grn_df$enhancer, "(?i)(?<=chr)\\d+")
  curr_grn_df$pr_chrs = str_extract(curr_grn_df$promoter, "(?i)(?<=chr)\\d+")
  # Extract enhancer/promoter start and end positions
  curr_grn_df$en_start = as.numeric(str_extract(curr_grn_df$enhancer, "(?i)(\\d+){5}"))
  curr_grn_df$en_end = abs(as.numeric(str_extract(
    curr_grn_df$enhancer, "(?i)-(\\d+){5}"
  )))
  curr_grn_df$pr_start = as.numeric(str_extract(curr_grn_df$promoter, "(?i)(\\d+){5}"))
  curr_grn_df$pr_end =  abs(as.numeric(str_extract(
    curr_grn_df$promoter, "(?i)-(\\d+){5}"
  )))
  # Return created grn dataframe
  return(curr_grn_df)
}

# grn_overlap:
#   - Creates GRanges information for current cell type GRN and disease GWAS
grn_overlap = function(ct_grn, gwas_gr,extension) {
  # Perform setup
  curr_grn = grn_setup(ct_grn)
  # Create GRanges objects for (en)hancer and (pr)omoter
  curr_en_gr = GRanges(seqnames = curr_grn$en_chrs,
                       IRanges(start = curr_grn$en_start-extension,
                               end = curr_grn$en_end+extension))
  curr_pr_gr = GRanges(seqnames = curr_grn$pr_chrs,
                       IRanges(start = curr_grn$pr_start-extension,
                               end = curr_grn$pr_end+extension))
  # Overlap (en)hancer and (pr)omoter with given disease GWAS
  overlap_en = as.data.frame(
    findOverlaps(
      query = curr_en_gr,
      subject = gwas_gr,
      type = 'any',
      select = "all",
      ignore.strand = T
    )
  )
  overlap_pr = as.data.frame(
    findOverlaps(
      query = curr_pr_gr,
      subject = gwas_gr,
      type = 'any',
      select = "all",
      ignore.strand = T
    )
  )
  # Get unique results for each, and return both
  unique_ens = unique(overlap_en$queryHits)
  unique_prs = unique(overlap_pr$queryHits)
  unique_ens_gwas = unique(overlap_en$subjectHits)
  unique_prs_gwas = unique(overlap_pr$subjectHits)
  unique_ids = list(unique_ens, unique_prs, unique_ens_gwas, unique_prs_gwas)
  return(unique_ids)
}


# ct_tfbs_interruption:
#   - Pulls cell type disease genes affected by TFBS interruption sites
ct_tfbs_interruption = function(ct_grn, mb_results_en, mb_results_pr) {
  # Perform setup
  curr_grn = grn_setup(ct_grn)
  # Create GRanges objects for (en)hancer and (pr)omoter
  curr_en_gr = GRanges(seqnames = curr_grn$en_chrs,
                       IRanges(start = curr_grn$en_start,
                               end = curr_grn$en_end))
  curr_pr_gr = GRanges(seqnames = curr_grn$pr_chrs,
                       IRanges(start = curr_grn$pr_start,
                               end = curr_grn$pr_end))
  
  mb_results_en_df = data.frame(mb_results_en)
  mb_results_en_df$chrs = as.numeric(gsub("[^0-9.-]", "", mb_results_en_df$seqnames))
  #mb_results_en_df = mb_results_en_df[which(mb_results_en_df$effect == 'strong'),]
  
  
  
  
  mb_results_en_gr = GRanges(
    seqnames = mb_results_en_df$chrs,
    IRanges(start = mb_results_en_df$start,
            end = mb_results_en_df$end)
  )
  
  grn_mb_en_overlap = as.data.frame(
    findOverlaps(
      query = curr_en_gr,
      subject = mb_results_en_gr,
      type = 'any',
      select = "all",
      ignore.strand = T
    )
  )
  
  grn_mb_en_overlap_df = unique(curr_grn[unlist(unique(grn_mb_en_overlap)), ])
  
  mb_results_pr_df = data.frame(mb_results_pr)
  mb_results_pr_df$chrs = as.numeric(gsub("[^0-9.-]", "", mb_results_pr_df$seqnames))
  #mb_results_pr_df = mb_results_pr_df[which(mb_results_pr_df$effect == 'strong'),]
  
  mb_results_pr_gr = GRanges(
    seqnames = mb_results_pr_df$chrs,
    IRanges(start = mb_results_pr_df$start,
            end = mb_results_pr_df$end)
  )
  
  grn_mb_pr_overlap = as.data.frame(
    findOverlaps(
      query = curr_pr_gr,
      subject = mb_results_pr_gr,
      type = 'any',
      select = "all",
      ignore.strand = T
    )
  )
  
  grn_mb_pr_overlap_df = unique(curr_grn[unlist(unique(grn_mb_pr_overlap)), ])
  
  grn_mb_total = list(grn_mb_en_overlap_df, grn_mb_pr_overlap_df)
  
  return(grn_mb_total)
}


