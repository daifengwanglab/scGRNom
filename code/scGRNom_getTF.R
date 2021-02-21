scGRNom_getTF <- function(df, database = JASPAR2018::JASPAR2018, species_type = 9606, min_score = 0.9,
                        pwm_type = 'prob',num_cores = 2){
  
  opts <- list()
  opts[["species"]] <- species_type
  opts[["all_versions"]] <- TRUE
  PFMatrixList <- TFBSTools::getMatrixSet(database, opts)
  pwmlist <- TFBSTools::toPWM(PFMatrixList, type = pwm_type)
  TF_names <- TFBSTools::name(pwmlist)
  names(TF_names) = NULL
  
  TF_names_splited = sapply(TF_names,data.table::tstrsplit,'::|\\(var.2\\)|\\(var.3\\)')
  
  df$promoter_id <- paste(df$gene_chr,':',df$promoter_start,'-',df$promoter_end,sep = '')
  df$enhancer_id <- paste(df$enh_chr,':',df$enh_start,'-',df$enh_end,sep = '')
  df = data.table::data.table(df)
  
  df_p <- dplyr::distinct(df[,c('gene_chr','promoter_start','promoter_end','promoter_id')])
  df_e <- dplyr::distinct(df[,c('enh_chr','enh_start','enh_end','enhancer_id')])
  
  suppressWarnings( G1 <- GenomicRanges::GRanges(seqnames = df_p$gene_chr,
                                                 IRanges::IRanges(start=df_p$promoter_start,
                                                                  end=df_p$promoter_end),
                                                 seqlengths = GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:24]))
  G1 <- GenomicRanges::trim(G1)
  suppressWarnings( G2 <- GenomicRanges::GRanges(seqnames = df_e$enh_chr,
                                                 IRanges::IRanges(start=df_e$enh_start,
                                                                  end=df_e$enh_end),
                                                 seqlengths = GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:24]))
  G2 <- GenomicRanges::trim(G2)
  
  
  cl <- parallel::makeCluster(num_cores) # not overload your computer
  doParallel::registerDoParallel(cl)
  
  `%dopar%` <- foreach::`%dopar%`
  df_p$promoter_TF <- foreach::foreach(i = 1:nrow(df_p), .combine = rbind,
                                       .packages = c('data.table','motifmatchr')) %dopar% {
                                         peak <- G1[i]
                                         motif_ix <- matchMotifs(pwmlist, peak,
                                                                 genome = "hg19",
                                                                 out = "scores"
                                         )
                                         result <- motifScores(motif_ix)[1,]
                                         curr_TF <- unique(unlist(TF_names_splited[result > quantile(result,min_score)]))
                                         if(length(curr_TF) == 0){
                                           curr_TF <- NA
                                         }
                                         data.table(promoter_TF = list(curr_TF))
                                         
                                       }
  parallel::stopCluster(cl)
  
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  df_e$enhancer_TF <- foreach::foreach(i = 1:nrow(df_e), .combine = rbind,
                                       .packages = c('data.table','motifmatchr')) %dopar% {
                                         peak <- G2[i]
                                         motif_ix <- matchMotifs(pwmlist, peak,
                                                                 genome = "hg19",
                                                                 out = "scores")
                                         result <- motifScores(motif_ix)[1,]
                                         curr_TF <- unique(unlist(TF_names_splited[result > quantile(result,min_score)]))
                                         if(length(curr_TF) == 0){
                                           curr_TF <- NA
                                         }
                                         data.table(promoter_TF = list(curr_TF))
                                       }
  parallel::stopCluster(cl)
  
  
  df$promoter_TF <- df_p$promoter_TF[match(df$promoter_id, df_p$promoter_id)]
  df$enhancer_TF <- df_e$enhancer_TF[match(df$enhancer_id, df_e$enhancer_id)]
  
  df <- df[, c('gene','promoter_id','enhancer_id',
               'promoter_TF','enhancer_TF')]
  colnames(df) <- c('gene','promoter','enhancer',
                    'promoter_TF','enhancer_TF')
  return(df)


}
