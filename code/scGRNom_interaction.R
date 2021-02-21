scGRNom_interaction = function(hic_interaction, enhancers, ref_promoters = 'all',up_stream = 2500,
                               down_stream = 2500, link_type = 'within' ,target_genes='all',
                      gene_id_option = 'hgnc_symbol',
                      mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                                     dataset="hsapiens_gene_ensembl",
                                     host="uswest.ensembl.org")){


  ehs <- GenomicRanges::GRanges(seqnames = enhancers$chr, ranges = IRanges::IRanges(start = enhancers$start,
                                                              end = enhancers$end))
  names(ehs) <- paste("ENH", as.character(ehs), sep = "_")

  # get hi-tad data
  # I assume there are strictly 5 columns or 6 columns
  if(ncol(hic_interaction)==6){
    anchor.one <- GenomicRanges::GRanges(hic_interaction$chr1,
                                   IRanges::IRanges(start = hic_interaction$start1, end = hic_interaction$end1))
    anchor.two <- GenomicRanges::GRanges(hic_interaction$chr2,
                                   IRanges::IRanges(start = hic_interaction$start2, end = hic_interaction$end2))
    hic <- GenomicInteractions::GenomicInteractions(anchor.one, anchor.two)

  }else{
    anchor.one <- GenomicRanges::GRanges(hic_interaction$chr,
                                   IRanges::IRanges(start = hic_interaction$start1, end = hic_interaction$end1))
    anchor.two <- GenomicRanges::GRanges(hic_interaction$chr,
                                   IRanges::IRanges(start = hic_interaction$start2, end = hic_interaction$end2))
    hic <- GenomicInteractions::GenomicInteractions(anchor.one, anchor.two)
  }


  # annotation
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  hg19.genes <- GenomicFeatures::genes(txdb)
  hg19.transcripts <- GenomicFeatures::transcriptsBy(txdb, by="gene")
  hg19.transcripts <- hg19.transcripts[ names(hg19.transcripts)  %in% unlist(hg19.genes$gene_id) ]
  hg19_refseq_promoters <- GenomicFeatures::promoters(hg19.transcripts,upstream = up_stream,
                                     downstream = down_stream)

  seqnames_record <- GenomeInfoDb::seqnames(hg19_refseq_promoters)
  hg19_refseq_promoters <- unlist(hg19_refseq_promoters[S4Vectors::`%in%`(seqnames_record,
                                                        c('chrX','chrY',paste('chr',1:22,sep='')))])
  hg19_refseq_promoters <- unique(hg19_refseq_promoters)
  gene_names <- biomaRt::getBM(attributes = c("hgnc_symbol",'entrezgene_id','ensembl_gene_id'),
                      filters = "entrezgene_id",
                      values = names(hg19_refseq_promoters), mart = mart)
  hg19_refseq_promoters$geneSymbol <- gene_names$hgnc_symbol[match(names(hg19_refseq_promoters),
                                                                   gene_names$entrezgene_id)]
  hg19_refseq_promoters$ensembl_id <- gene_names$ensembl_gene_id[match(names(hg19_refseq_promoters),
                                                                       gene_names$entrezgene_id)]

  if(gene_id_option=="hgnc_symbol"){
    names(hg19_refseq_promoters) <- hg19_refseq_promoters$geneSymbol
    hg19_refseq_promoters <- hg19_refseq_promoters[(!is.na(names(hg19_refseq_promoters)) )&
                                                     (names(hg19_refseq_promoters)!='')]
  }else if(gene_id_option=="ensembl_gene_id"){
    names(hg19_refseq_promoters) <- hg19_refseq_promoters$ensembl_id
    hg19_refseq_promoters <-hg19_refseq_promoters[!is.na(names(hg19_refseq_promoters))]
  }

  # if target gene is not string 'all' then
  if(target_genes != 'all'){
    hg19_refseq_promoters <- hg19_refseq_promoters[ names(hg19_refseq_promoters)
                                                    %in% target_genes]
  }


  if(typeof(ref_promoters) == typeof('all')){
    annotation_promoters <- hg19_refseq_promoters
  }else{
    ref_promoters_granges <-GenomicRanges::GRanges(
                             seqnames = ref_promoters$chr,
                             ranges = IRanges::IRanges(start = ref_promoters$start, end = ref_promoters$end))
    overlap_granges <- IRanges::findOverlapPairs(ref_promoters_granges,hg19_refseq_promoters,
                                        type=link_type, ignore.strand=T)

    annotation_promoters <- overlap_granges@first
    names(annotation_promoters) <- names(overlap_granges@second)
    annotation_promoters <- unique(annotation_promoters)
  }


   annotation.features <- list(promoter = annotation_promoters, enhancer = ehs)
   GenomicInteractions::annotateInteractions(hic, annotation.features)
   interaction_index <- GenomicInteractions::isInteractionType(hic, "promoter", "enhancer")
   hic_subset <- hic[interaction_index]

   if(ncol(hic_interaction) == 6){
     df1 <- data.table::data.table(promoters = S4Vectors::mcols(GenomicInteractions::anchorOne(hic_subset))$promoter.id,
                       promoter_chr = hic_interaction$chr1[interaction_index],
                       enhs = S4Vectors::mcols(GenomicInteractions::anchorTwo(hic_subset))$enhancer.id,
                       enh_chr = hic_interaction$chr2[interaction_index])
     df1 <- df1[ (!is.na(df1$promoters)) & (!is.na(df1$enhs)) ,]

     df2 <- data.table::data.table(promoters = S4Vectors::mcols(GenomicInteractions::anchorTwo(hic_subset))$promoter.id,
                       promoter_chr = hic_interaction$chr2[interaction_index],
                       enhs = S4Vectors::mcols(GenomicInteractions::anchorOne(hic_subset))$enhancer.id,
                       enh_chr = hic_interaction$chr1[interaction_index]
     )
     df2 <- df2[ (!is.na(df2$promoters)) & (!is.na(df2$enhs)) ,]
   }else if( ncol(hic_interaction) == 5){
     df1 <- data.table::data.table(promoters = S4Vectors::mcols(GenomicInteractions::anchorOne(hic_subset))$promoter.id,
                       promoter_chr = hic_interaction$chr[interaction_index],
                       enhs = S4Vectors::mcols(GenomicInteractions::anchorTwo(hic_subset))$enhancer.id,
                       enh_chr = hic_interaction$chr[interaction_index]
     )
     df1 <- df1[ (!is.na(df1$promoters)) & (!is.na(df1$enhs)) ,]

     df2 <- data.table::data.table(promoters = S4Vectors::mcols(GenomicInteractions::anchorTwo(hic_subset))$promoter.id,
                       promoter_chr = hic_interaction$chr[interaction_index],
                       enhs = S4Vectors::mcols(GenomicInteractions::anchorOne(hic_subset))$enhancer.id,
                       enh_chr = hic_interaction$chr[interaction_index]
     )
     df2 <- df2[ (!is.na(df2$promoters)) & (!is.na(df2$enhs)) ,]
   }


  df <- data.table::rbindlist(list(df1,df2))
  df$id <- seq.int(nrow(df))
  df_enhancers <- df[,c('enhs','id')]

  df_enhancers <- df_enhancers[, .(enhancer = unlist(enhs)),by = id]
  df <- dplyr::left_join(df,df_enhancers,by = 'id')[,c('promoters','promoter_chr','enhancer','enh_chr')]
  df <- data.table::data.table(df)
  df$id <- seq.int(nrow(df))
  df_promoters <- df[,c('promoters','id')]
  df_promoters<- df_promoters[, .(promoter = unlist(promoters)),by = id]
  df <- dplyr::left_join(df,df_promoters,by = 'id')[,c('promoter','promoter_chr','enhancer','enh_chr')]
  df <- dplyr::distinct(df)
  sp_list  <- data.table::tstrsplit(df$enhancer,"_|:|-", keep = c(2,3,4))
  df$enh_start <- as.numeric(sp_list[[2]])
  df$enh_end <- as.numeric(sp_list[[3]])


  ref_df <- data.frame(promoter = names(annotation_promoters),
                        promoter_start = annotation_promoters@ranges@start,
                        promoter_end = annotation_promoters@ranges@start +
                        annotation_promoters@ranges@width, stringsAsFactors = F)


  final_df <- dplyr::full_join(ref_df,df,by = 'promoter')
  final_df <- na.omit(final_df)
  final_df <- final_df[,c('promoter','promoter_chr',
                          'promoter_start',
                          'promoter_end','enh_chr',
                         'enh_start','enh_end')]
  colnames(final_df) <- c('gene','gene_chr',
                       'promoter_start',
                       'promoter_end','enh_chr',
                      'enh_start','enh_end')

  return(final_df)
}
