


plot_manhattan <- function(data_manh,chr_cols=c("bisque","white"),spacer=5000000){
  
  data_manh <- data_manh[order(as.numeric(gsub("Chr","",data_manh[,Chr])))]
  data_manh[,Chr:=as.factor(Chr)]
  data_manh[,Chr:=factor(Chr,levels=data_manh[,levels(Chr)][order(as.numeric(gsub("Chr","",data_manh[,levels(Chr)])))])]
  
  don <- data_manh %>%
    
    # Compute chromosome size
    group_by(Chr) %>%
    summarise(chr_len=max(bp)+spacer) %>%
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(data_manh, ., by=c("Chr"="Chr")) %>%
    
    # Add a cumulative position of each SNP
    arrange(Chr, bp) %>%
    mutate( BPcum=bp+tot)
  
  axisdf = don %>%
    group_by(Chr) %>%
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  
  #ylim <- c(0,max(-log10(don$P))*1.1)
  
  chrs <- don[,unique(Chr)]
  n_chr <- length(chrs)
  
  rect_data <- rbindlist(lapply(chrs,function(chr){
    data.table(x1=min(don[Chr==chr,BPcum]),x2=max(don[Chr==chr,BPcum]),y1=0,y2=Inf)
  }))
  
  rect_data[,col:=rep(chr_cols,ceiling(n_chr/2))[1:(n_chr)]]
  list(don=don,axisdf=axisdf,ylim=ylim,rect_data=rect_data)
}

col_vector <- c("#7FC97F", "#FDC086", "#F0027F", "#BF5B17", "#1B9E77", "#7570B3", 
                "#E6AB02", "#A6761D", "#1F78B4", "#B2DF8A", "#A6CEE3", "#33A02C", 
                "#FB9A99", "#E31A1C", "#CAB2D6", "#B15928", "#FBB4AE", "#6A3D9A", 
                "#FDDAEC", "#B3E2CD", "#FDCDAC", "#A6761D", "#FFF2AE", "#E41A1C", 
                "#4DAF4A", "steelblue", "salmon", "firebrick", "grey30", "#FF7F00", 
                "#FFFF33", "#66C2A5", "#8DA0CB", "#E78AC3", "#FFD92F", "#BEBADA", 
                "#FB8072", "forestgreen", "#E6AB02", "#1F78B4", "#BC80BD")