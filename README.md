Random Functions:

# plot_volcano:
function to make volcano plot of epigenome-wide association results with
  effect estimates on x-axis and
  -log10 pvalues on y-axis
  
- Inputs:
  - sshits: dataframe of results summary of EWAS  
  - est_col: column name of sshits for effect estimates  
  - p_col: column name of sshits for p-values  
  - xlim/ylim: x and y axis limits for plot, will auto calculate something if NULL  
  - pct_include: whether to include % hyper/hypo methylated CpGs meeting pct_p_thresh  
  - pct_p_thresh: p-value threshold for coloring points differently and for considering & hyper/hypo methylated CpGs percentage  
  - pct_nudge: how much to move % hyper/hypo methylated labels along x-axis further from 0  
  - pct_height: what height on y-axis should percentage annotations be placed  
  - title: title for plot  
  - cols: vector of three colors with positions corresponding to:  
     - 1: points not meeting pct_p_thresh   
     - 2: points meeting pct_p_thresh and effect estimate > 0   
     - 3: points meeting pct_p_thresh and effect estimate < 0  
     - if pct_include=F, then only first two slots for above/below threshold are used  

# qq_plot:
function to make qq-plot from summary of epigenome-wide association results  
  
- Inputs:
  - sshits: dataframe of results summary of EWAS  
  - p_col: column name of sshits for p-values  

# genomic_reg_annotation:
function to  add genomic regulatory element annotation (promoters, 1-5kb upstream TSS, 5UTR, 3UTR, exon, intron)    
  
- Inputs:
  - data: GRanges object, or something that can be turned into one (chr, start, stop columns)  
  - preloaded.anno.data: GRanges object built by "build_annotations" function of "annotatr" package, if using function multiple times will save lot of time loading it once outside function and passing it to the function  
