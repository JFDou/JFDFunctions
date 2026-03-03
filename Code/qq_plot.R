# function to make qq-plot from limma topTable results of epigenome-wide association results
#
# sshits: dataframe of results summary of EWAS
# p_col: column name of sshits for p-values
# xlim/ylim: x and y axis limits for plot, will auto calculate something if NULL

qq_plot <- function(ss.hits, p_col='P.Value'){
  ss.hits$P.Value <- ss.hits[,p_col]
  ggplot(ss.hits, 
         aes(y=-log(P.Value,10), 
             x=-log(ppoints(nrow(ss.hits)),10))) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    theme_bw() +
    ggtitle(paste0("Lambda = ", round(inflation(ss.hits$P.Value),3))) +
    xlab("Expected") + ylab("Observed") +
    theme(text = element_text(size=14))
}
