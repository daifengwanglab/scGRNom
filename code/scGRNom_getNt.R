scGRNom_getNt <- function(df, gexpr, df_gene_id = 'hgnc_symbol', gexpr_gene_id = 'hgnc_symbol',
                          cutoff_by = 'quantile', cutoff_percentage = 0.9, cutoff_absolute = 0.1,scaleby = 'no',
                          train_ratio = 0.7, num_cores = 2,
                          mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                                                  dataset="hsapiens_gene_ensembl",
                                                  host="uswest.ensembl.org"), open_chrom = NULL, extension_bps = 0){
   library(stringr)
   if(!is.null(open_chrom)){
      df$en_chrs = str_extract(df$enhancer, "(?i)(?<=chr)\\d+")
      # Extract enhancer start and end positions
      df$en_start = as.numeric(str_extract(df$enhancer, "(?i)(\\d+){5}"))
      df$en_end = abs(as.numeric(str_extract(df$enhancer, "(?i)-(\\d+){5}")))
      df$en_start_intend = df$en_start -extension_bps
      df$en_end_intend  = df$en_end + extension_bps
      
      colnames(open_chrom) = c('chrs', 'Start', 'Stop')
      open_chrom$chrs = gsub('chr','', open_chrom$chrs)
      
      library(GenomicRanges)
      en_scGRN = GRanges(seqnames = df$en_chrs,
                         IRanges(start = df$en_start_intend,
                                 end = df$en_end_intend))
      
      en_openchrom = GRanges(seqnames = open_chrom$chrs,
                        IRanges(start = open_chrom$Start,
                                 end = open_chrom$Stop))
      
      overlap_en = as.data.frame(
        findOverlaps(
          query = en_scGRN,
          subject = en_openchrom,
          type = 'any',
          select = "all",
          ignore.strand = T
        )
      )
      
      unique_ens = unique(overlap_en$queryHits)
      df = df[unique_ens,]
    }
  
    cl <- parallel::makeCluster(num_cores) # not overload your computer
    doParallel::registerDoParallel(cl)
    `%dopar%` <- foreach::`%dopar%`
    
    df$TFs <- foreach::foreach(i = 1:nrow(df), .combine = rbind,
                           .packages = c('data.table')) %dopar% {
                             curr_TF <- unique(c(df$enhancer_TF[[i]], df$promoter_TF[[i]]))
                             curr_TF <- curr_TF[!is.na(curr_TF)]
                             if(length(curr_TF) == 0){
                               curr_TF <- NA
                             }
                             data.table::data.table(TFs = list(curr_TF))
                           }
    parallel::stopCluster(cl)
    
    df <- df[, c('gene','enhancer','promoter','enhancer_TF','promoter_TF','TFs')]
    df <- df[!is.na(df$TFs),]
    df$id <- seq.int(nrow(df))
    df_TF <- df[,c('TFs','id')][,.(TF = unlist(TFs)), by = id]
    df <- dplyr::left_join(df,df_TF,by = 'id')
    
    df$id <- NULL
    df$TFs <- NULL
    
    cl <- parallel::makeCluster(num_cores) # not overload your computer
    doParallel::registerDoParallel(cl)
    df$TFbs <- foreach::foreach(i = 1:nrow(df), .combine = rbind
    ) %dopar% {
    
    if(is.na(df$TF[i])){
    print("error")
    }
    
    if((df$TF[i] %in% df$enhancer_TF[[i]]) & (df$TF[i] %in% df$promoter_TF[[i]])){
    binding_site <- 'both'
    }else if(df$TF[i] %in% df$enhancer_TF[[i]]){
    binding_site <- 'enhancer'
    }else if(df$TF[i] %in% df$promoter_TF[[i]]){
    binding_site <- 'promoter'
    }
    binding_site
    }
    parallel::stopCluster(cl)
    
    df <- df[,c('gene','promoter','enhancer','TFbs','TF')]
    df <- df[df$TF != df$gene, ]
    
    
    ###### change TF names from hgnc_symbol to emsembl_id
    ###### delete the TF whose ensembl_id is NA
    
    if(gexpr_gene_id == 'hgnc_symbol'){
    if(df_gene_id != 'hgnc_symbol'){
    gene_names <- biomaRt::getBM(attributes = c("hgnc_symbol","ensembl_gene_id"), filters = "ensembl_gene_id",
                                 values = unique(df$gene), mart = mart)
    df$gene <- gene_names$hgnc_symbol[match(df$gene, gene_names$ensembl_gene_id)]
    df <- na.omit(df)
    df <- df[df$gene != '',]
    }
    }
    
    if(gexpr_gene_id == 'ensembl_gene_id'){
    if(df_gene_id == 'ensembl_gene_id'){
    gene_names <- biomaRt::getBM(attributes = c("hgnc_symbol","ensembl_gene_id"), filters = "hgnc_symbol",
                                 values = unique(df$TF), mart = mart)
    df$TF <- gene_names$ensembl_gene_id[match(df$TF, gene_names$hgnc_symbol)]
    df <- na.omit(df)
    }else{
    gene_names <- biomaRt::getBM(attributes = c("hgnc_symbol","ensembl_gene_id"), filters = "hgnc_symbol",
                                 values = c(unique(df$gene),unique(df$TF)), mart = mart)
    df$gene <- gene_names$ensembl_gene_id[match(df$gene, gene_names$hgnc_symbol)]
    df$TF <- gene_names$ensembl_gene_id[match(df$TF, gene_names$hgnc_symbol)]
    df <- na.omit(df)
    }
    }
  
    library(glmnet)
    tgs <- as.character(unique(df$gene))
    output_df = data.frame(TG = NULL,TF = NULL, coef = NULL, mse = NULL)
    set.seed(123)
    
    for (i in 1:length(tgs)){
    
    selgene <- tgs[i]
    
    if(selgene %in% rownames(gexpr)){
    selTFs <- unique(df[df$gene == selgene, 'TF'])
    train_cols <- sample(1:ncol(gexpr), round(train_ratio*ncol(gexpr)))
    if(sum(rownames(gexpr) %in% selTFs$TF) > 1 ){
      
      x.train <- t(gexpr[rownames(gexpr) %in% selTFs$TF,train_cols])
      x.test <- t(gexpr[rownames(gexpr) %in% selTFs$TF,-train_cols])
      y.train <- t(gexpr[selgene,train_cols])
      y.test <- t(gexpr[selgene,-train_cols])
      
      yfit <- vector('list',11)
      yhat <- vector('list',length(yfit))
      mse <- rep(Inf,length(yfit))
      
      if(sd(y.train) == 0){
        curr_df <- data.frame(TG = NULL, TF = NULL, coef = NULL, mse = NULL)
      }else{
        
        for (j in 0:10) {
          yfit[[j+1]] <- glmnet::cv.glmnet(x.train, y.train, type.measure="mse",
                                           alpha=j/10,family="gaussian",standardize= F,intercept=T)
          yhat[[j+1]] <- as.numeric(predict(yfit[[j+1]], s=yfit[[j+1]]$lambda.1se, newx=x.test))
          mse[j+1] <- mean((y.test - yhat[[j+1]])^2)
        }
        fitcoef <- coef(yfit[[which.min(mse)]], s = "lambda.min")
        
        TF_coef <- as.matrix(fitcoef)
        TF_coef <- TF_coef[2:nrow(TF_coef),]
        if(cutoff_by == 'quantile'){
          TF_coef <- TF_coef[abs(TF_coef) > quantile(abs(TF_coef),1 - cutoff_percentage)]
        }else if(cutoff_by == 'absolute'){
          TF_coef <- TF_coef[abs(TF_coef) > cutoff_absolute]
        }else{
          print('cutoff_by can only take absolute or quantile')
        }
        
        if(length(TF_coef) > 0){
          curr_df <- data.frame(TG = rep(selgene,length(TF_coef)), TF = names(TF_coef),
                                coef = unname(TF_coef), mse =min(mse),  stringsAsFactors = F)
        }else{
          curr_df = data.frame(TG = NULL,TF = NULL, coef = NULL, mse = NULL)
        }
        
      }
    }else{
      curr_df = data.frame(TG = NULL,TF = NULL, coef = NULL, mse = NULL)
    }
    
    }else{
      curr_df = data.frame(TG = NULL,TF = NULL, coef = NULL, mse = NULL)
    }
    
    output_df = rbind(output_df,curr_df)
    }
    
    if(nrow(output_df) > 0){
    df$id <- paste(df$gene,'-',df$TF,sep = '')
    output_df$id <- paste(output_df$TG,'-',output_df$TF,sep = '')
    output_df$TF <- NULL
    df <- dplyr::full_join(output_df,df,by = 'id')
    
    df <- na.omit(df)
    df <- df[,c('TG','TF','enhancer','promoter','TFbs','coef', 'mse')]
    }else{
    df <- data.frame(TG = NULL, TF = NULL, enhancer = NULL, promoter = NULL, TFbs = NULL,coef = NULL, mse = NULL)
    }

    return(df)
}
