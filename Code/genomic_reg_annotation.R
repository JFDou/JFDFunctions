# add genomic regulatory element annotation (promoters, 1-5kb upstream TSS, 5UTR, 3UTR, exon, intron)
# data: GRanges object, or something that can be turned into one (chr, start, stop columns)
# preloaded.anno.data: GRanges object built by "build_annotations" function of "annotatr" package,
#     if using function multiple times will save lot of time loading it once outside function
#     and passing it to the function

genomic_reg_annotation <- function(data, preloaded.anno.data=NULL){
  require(annotatr)
  require(org.Mm.eg.db)
  require(GenomicRanges)
  require(dplyr)
  
  ### load in genomic annotation
  if(is.null(preloaded.anno.data)){
    genomic.function <- build_annotations(genome='mm10', annotations="mm10_basicgenes")
  }else{
    genomic.function <- preloaded.anno.data
  }
  genomic.function.sub <- genomic.function[!is.na(genomic.function$type)]
  
  ### make data dataframe so its easier add to it later 
  data <- as.data.frame(data)

  ### annotate DMR candidates
  anno_data <- findOverlaps(makeGRangesFromDataFrame(data), genomic.function.sub)
  
  matches <- data.frame(dat_row = anno_data@from,
                        gene = genomic.function.sub[anno_data@to,]$symbol,
                        type = genomic.function.sub[anno_data@to,]$type) 

  ### collapse and summarize regulatory element entries
  matches$type <- gsub('mm10_genes_','',matches$type)
  matches <- matches[!is.na(matches$gene),]
  matches_sm <- matches %>% group_by(dat_row, gene) %>%
    summarize(type=unique(type) %>% paste0(collapse=','))
  matches_sm <- matches_sm %>% group_by(dat_row) %>%
    reframe(gene=paste0(gene, collapse=';'),
            type=paste0(type, collapse=';'))
  
  data[matches_sm$dat_row,'gene'] <- matches_sm$gene
  data[matches_sm$dat_row,'type'] <- matches_sm$type
  
  data$gene <- ifelse(is.na(data$gene), "", data$gene)
  data$type <- ifelse(is.na(data$type), "", data$type)
  
  return(data)
}