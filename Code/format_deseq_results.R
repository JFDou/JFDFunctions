# function to clean output from deseq results of RNA expression analysis
#
# deseq.res: results object from results or lfcShrink function

format_deseq_results <- function(deseq.res, ){
  # rownames into a gene column
  deseq.res$gene <- rownames(deseq.res)
  rownames(deseq.res) <- NULL
  
  # round numbers
  deseq.res$baseMean <- round(deseq.res$baseMean,2)
  deseq.res$log2FoldChange <- round(deseq.res$log2FoldChange,3)
  deseq.res$lfcSE <- round(deseq.res$lfcSE,3)
  deseq.res$pvalue <- signif(deseq.res$pvalue,3) %>% as.character()
  deseq.res$padj <- signif(deseq.res$padj,3) %>% as.character()
  
  # order columns
  deseq.res[,c('gene','baseMean','log2FoldChange','pvalue','padj')]
}
