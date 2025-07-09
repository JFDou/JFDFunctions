# function to make volcano plot of epigenome-wide association results
#   effect estimates on x-axis and
#   -log10 pvalues on y-axis
# xlim/ylim: x and y axis limits for plot, will auto calculate something if NULL
# pct_include: whether to include % hyper/hypo methylated CpGs meeting pct_p_thresh
# pct_p_thresh: p-value threshold for coloring points differently and for 
#    considering & hyper/hypo methylated CpGs percentage
# pct_nudge: how much to move % hyper/hypo methylated labels along x-axis further from 0
# title: title for plot
# cols: three colors with positions corresponding to:
#     1: points meeting pct_p_thresh and effect estimate > 0
#     2: points meeting pct_p_thresh and effect estimate < 0
#     3: points not meeting pct_p_thresh 
plot_volcano <- function(ss.hits,xlim=NULL,ylim=NULL, 
                         pct_include=F, pct_p_thresh=0.05, pct_nudge=0, 
                         title="", 
                         cols=c('darkslategrey','cadetblue4','cadetblue3')){
  
  #get betas and pvals from single sites results
  betas <- ss.hits$Effect.Estimate
  pvals <- -log(ss.hits$P.Value,10)
  
  #assign colors
  pcol <- ifelse(pvals < -log(pct_p_thresh,10), cols[1], ifelse(betas>0,cols[2],cols[3]))
  
  #get percent hyper/hypo for sites above threshold
  if(include_pct){
    total <- table(pcol)[cols[2]] + table(pcol)[cols[3]]
    pct <- c(round(table(pcol)[cols[2]]/total * 100,1),round(table(pcol)[cols[3]]/total * 100,1))
    pct <- c(paste0(pct[1],'%'),paste0(pct[2],'%'))
  }else{
    pct=c('','')
  }
  
  
  #auto calculate axis range if not defined
  if(is.null(xlim)){
    max.x <- max(abs(betas))
    xlim <- c(-max.x-2, max.x+2)
  }
  if(is.null(ylim)){
    ylim <- c(0,ifelse(max(pvals)<7, 7, max(pvals)+0.5))
  }
  
  #make plot
  df <- data.frame(betas=betas, pvals=pvals, col=pcol)
  ggplot(df, aes(x=betas, y=pvals, col=I(col))) +
    geom_point() + ggtitle(title) +
    geom_hline(yintercept=-log(0.05,10),color='blue') +
    geom_hline(yintercept=7,color='red') +
    xlim(xlim) + ylim(ylim) +
    xlab(paste0('Effect Estimate')) + ylab(expression('-Log'[10]*'(P-Value)')) +
    geom_text(aes(x=max(abs(betas))+nudge,y=6,label=pct[1],color=I(cols[2]))) + 
    geom_text(aes(x=-max(abs(betas))-nudge,y=6,label=pct[2],color=I(cols[3]))) +
    theme_bw() +
    theme(axis.text=element_text(size=12),
          legend.position = "none")
}
