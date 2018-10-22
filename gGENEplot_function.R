gGENEplot <- function(source.file, min.size = 0.001,
                      line.color = "gray", line.size = 2, line.overhang = 0.05,
                      intron.color = "lightblue", exon.color = "blue", 
                      mRNA.names = NULL, CDS.names = NULL, 
                      check.introns = T, names.order = NULL,
                      mRNA.order = NULL, CDS.order = NULL,
                      bar = F, bar.pos = "bottom", bar.color = "black",
                      bar.size = line.size/2, bar.length = 0){
  require(gsubfn)
  require(ggplot2)
  
  # variables
  y.pos <- 3 # to construct the plot
  
  # messages
  ERROR.min.size <- "min.size argument cannot be negative. Reset to default value."
  ERROR.overh <- "line.overhang argument cannot be negative. Reset to default value."
  ERROR.bar <- "bar.pos argument must be 'bottom' or 'top'. bar is not shown."
  
  #################
  # PARSE gb file #
  #################
  input <- readLines(source.file)
  mRNA.lines <- which(grepl("^[ ]{5}mRNA", input))
  CDS.lines <- which(grepl("^[ ]{5}CDS", input))
  
  df <- data.frame("gene" = NA, "isoform" = NA, "feature" = NA, 
                   "exon" = NA, "x" = NA, "y" = NA)
  
  # Get mRNA data from mRNA.ines
  index <- 0
  mRNA <- c("")
  for (i in mRNA.lines){
    index <- index + 1
    line.num <- i
    line <- input[line.num]
    temp <- strapplyc(line, "^[ ]{5}mRNA[ ]{12}join.(.*)", simplify = TRUE)
    if (is.list(temp)){
      index <- index - 1
      next
    }
    while (!(grepl(")$", temp))){
      line.num <- line.num + 1
      line <- input[line.num]
      temp <- paste0(temp, strapplyc(line, "^[ ]{21}(.*)", simplify = TRUE))
    }
    #  Get gene.name from mRNA block
    line <- input[line.num + 1]
    gene.name <- strapplyc(line, "^[ ]+/gene=\\\"(.*)\\\"", simplify = TRUE)
    mRNA[index] <- strapplyc(temp, "(.*).+", simplify = TRUE)
    mRNA.blocks <- strsplit(mRNA[index], ",")[[1]]
    exon.index <- 0
    for (j in 1:length(mRNA.blocks)){
      exon.index <- exon.index + 1
      init <- strsplit(mRNA.blocks[j], "[..]")[[1]][1]
      if (!(grepl("^[0-9]+", init))){
        init <- as.numeric(strapplyc(init, "([0-9]+)", simplify = TRUE))
      } else {
        init <- as.numeric(init)
      }
      end <- strsplit(mRNA.blocks[j], "[..]")[[1]][3]
      if (!(grepl("^[0-9]+", end))){
        end <- as.numeric(strapplyc(end, "([0-9]+)", simplify = TRUE))
      } else {
        end <- as.numeric(end)
      }
      df <- rbind(df, 
                  c(gene.name, index, "mRNA", exon.index, init, y.pos - 0.2), 
                  c(gene.name, index, "mRNA", exon.index, init, y.pos + 0.2),
                  c(gene.name, index, "mRNA", exon.index, end, y.pos + 0.2), 
                  c(gene.name, index, "mRNA", exon.index, end, y.pos - 0.2))
    }
  }
  
  # Get CDS data from CDS.ines
  index <- 0
  CDS <- c("")
  for (i in CDS.lines){
    index <- index + 1
    line.num <- i
    line <- input[line.num]
    temp <- strapplyc(line, "^[ ]{5}CDS[ ]{13}join.(.*)", simplify = TRUE)
    while (!(grepl(")$", temp))){
      line.num <- line.num + 1
      line <- input[line.num]
      temp <- paste0(temp, strapplyc(line, "^[ ]{21}(.*)", simplify = TRUE))
    }
    #  Get gene.name from CDS block
    line <- input[line.num + 1]
    gene.name <- strapplyc(line, "^[ ]+/gene=\\\"(.*)\\\"", simplify = TRUE)
    CDS[index] <- strapplyc(temp, "(.*).+", simplify = TRUE)
    CDS.blocks <- strsplit(CDS[index], ",")[[1]]
    exon.index <- 0
    for (j in 1:length(CDS.blocks)){
      exon.index <- exon.index + 1
      init <- strsplit(CDS.blocks[j], "[..]")[[1]][1]
      if (!(grepl("^[0-9]+", init))){
        init <- as.numeric(strapplyc(init, "([0-9]+)", simplify = TRUE))
      } else {
        init <- as.numeric(init)
      }
      end <- strsplit(CDS.blocks[j], "[..]")[[1]][3]
      if (!(grepl("^[0-9]+", end))){
        end <- as.numeric(strapplyc(end, "([0-9]+)", simplify = TRUE))
      } else {
        end <- as.numeric(end)
      }
      df <- rbind(df, 
                  c(gene.name, index, "CDS", exon.index, init, y.pos - 0.2), 
                  c(gene.name, index, "CDS", exon.index, init, y.pos + 0.2),
                  c(gene.name, index, "CDS", exon.index, end, y.pos + 0.2), 
                  c(gene.name, index, "CDS", exon.index, end, y.pos - 0.2))
    }
  }
  
  df <- df[which(!is.na(df$gene)),]
  df$exon <- as.numeric(df$exon)
  df$x <- as.numeric(df$x)
  df$y <- as.numeric(df$y)
  df$isoform <- as.numeric(df$isoform)
  
  ######################
  # Name gene isoforms #
  ######################
  if (!is.null(mRNA.names) & !is.null(CDS.names)) {
    df$isoform.name <- NA
    for (g in 1:max(df$isoform, na.rm = T)){
      df$isoform.name[which(df$isoform == g & 
                              df$feature == "mRNA")] <- mRNA.names[g]
      df$isoform.name[which(df$isoform == g & 
                              df$feature == "CDS")] <- CDS.names[g]
    }
    ################
    # EXON matcher #
    ################
    df$new.exon <- df$exon
    for (g in unique(df$isoform.name)){
      mRNA.exons <- df$exon[which(df$isoform.name == g & df$feature == "mRNA")]
      mRNA.exons <- mRNA.exons[seq(1, length(mRNA.exons), 4)]
      CDS.exons <- df$exon[which(df$isoform.name == g & df$feature == "CDS")]
      CDS.exons <- CDS.exons[seq(1, length(CDS.exons), 4)]
      mRNA.index <- 1
      for (iii in CDS.exons){
        temp <- df$x[which(df$isoform.name == g & 
                             df$exon == iii & df$feature == "CDS")]
        init.CDS <- temp[1]
        end.CDS <- temp[4]
        for (jjj in mRNA.exons[mRNA.index:length(mRNA.exons)]){
          temp2 <- df$x[which(df$isoform.name == g & 
                                df$exon == jjj & df$feature == "mRNA")]
          init.mRNA <- temp2[1]
          end.mRNA <- temp2[4]
          if (init.mRNA <= init.CDS & end.mRNA == end.CDS){
            df$new.exon[which(df$isoform.name == g & 
                                df$feature == "CDS" & df$exon == iii)] <- jjj
            break
          }
          if (init.mRNA == init.CDS & end.mRNA > end.CDS){
            df$new.exon[which(df$isoform.name == g & 
                                df$feature == "CDS" & df$exon == iii)] <- jjj
            break
          }
          mRNA.index <- mRNA.index + 1
        }
      }
    }
    df$exon <- df$new.exon
    df <- df[,-which(colnames(df) == "new.exon")]
    #################################################
    # correct exon/intron size to make them visible #
    #################################################
    if (min.size != 0){
      if (min.size < 0) {
        min.size <- 0.001
        warning(ERROR.min.size)
      }
      tot.size <- (max(df$x) - min(df$x)) + 1
      if (min.size >= 1) {
        size.min <- min.size
      } else {
        size.min <- floor(tot.size * min.size)
      }
      # exons
      for (g in unique(df$isoform.name)){
        exons <- unique(df$exon[which(df$isoform.name == g & 
                                        df$feature == "mRNA")])
        exons.CDS <- unique(df$exon[which(df$isoform.name == g & 
                                            df$feature == "CDS")])
        for (ex in exons){
          lines.mRNA <- which(df$isoform.name == g & df$feature == "mRNA" & 
                                df$exon == ex)
          block.size.mRNA <- (df$x[lines.mRNA[3]] - df$x[lines.mRNA[2]]) + 1
          if (ex %in% exons.CDS){
            lines.CDS <- which(df$isoform.name == g & df$feature == "CDS" & 
                                 df$exon == ex)
            block.size.CDS <- (df$x[lines.CDS[3]] - df$x[lines.CDS[2]]) + 1
            if (block.size.CDS < size.min){
              diff <- size.min - block.size.CDS
              # mRNA block
              df$x[lines.mRNA[1]] <- df$x[lines.mRNA[1]] - floor(diff / 2)
              df$x[lines.mRNA[2]] <- df$x[lines.mRNA[2]] - floor(diff / 2)
              df$x[lines.mRNA[3]] <- df$x[lines.mRNA[3]] + floor(diff / 2)
              df$x[lines.mRNA[4]] <- df$x[lines.mRNA[4]] + floor(diff / 2)
              # CDS block
              df$x[lines.CDS[1]] <- df$x[lines.CDS[1]] - floor(diff / 2)
              df$x[lines.CDS[2]] <- df$x[lines.CDS[2]] - floor(diff / 2)
              df$x[lines.CDS[3]] <- df$x[lines.CDS[3]] + floor(diff / 2)
              df$x[lines.CDS[4]] <- df$x[lines.CDS[4]] + floor(diff / 2)
            }
          } else if (block.size.mRNA < size.min){
            diff <- size.min - block.size.mRNA
            df$x[lines.mRNA[1]] <- df$x[lines.mRNA[1]] - floor(diff / 2)
            df$x[lines.mRNA[2]] <- df$x[lines.mRNA[2]] - floor(diff / 2)
            df$x[lines.mRNA[3]] <- df$x[lines.mRNA[3]] + floor(diff / 2)
            df$x[lines.mRNA[4]] <- df$x[lines.mRNA[4]] + floor(diff / 2)
          }
        }
      }
      # introns
      if (check.introns){
        for (g in unique(df$isoform.name)){
          exons <- unique(df$exon[which(df$isoform.name == g & 
                                          df$feature == "mRNA")])
          exons.CDS <- unique(df$exon[which(df$isoform.name == g & 
                                              df$feature == "CDS")])
          for (ii in seq(2, length(exons))){
            lines.mRNA <- which(df$isoform.name == g & df$feature == "mRNA" & 
                                  df$exon == exons[ii])
            intron.size <- (df$x[lines.mRNA[1]] - df$x[lines.mRNA[1]-1]) + 1
            if (intron.size < size.min){
              diff <- size.min - intron.size
              start <- df$x[lines.mRNA[1]]
              end <- df$x[lines.mRNA[4]]
              df$x[lines.mRNA[1]] <- df$x[lines.mRNA[1]] + diff
              df$x[lines.mRNA[2]] <- df$x[lines.mRNA[2]] + diff
              df$x[lines.mRNA[3]] <- df$x[lines.mRNA[3]] + diff
              df$x[lines.mRNA[4]] <- df$x[lines.mRNA[4]] + diff
              if (ex %in% exons.CDS){
                lines.CDS <- which(df$isoform.name == g & df$feature == "CDS" & 
                                     df$exon == exons[ii])
                df$x[lines.CDS[1]] <- df$x[lines.CDS[1]] + diff
                df$x[lines.CDS[2]] <- df$x[lines.CDS[2]] + diff
                df$x[lines.CDS[3]] <- df$x[lines.CDS[3]] + diff
                df$x[lines.CDS[4]] <- df$x[lines.CDS[4]] + diff
              }
            }
          }
        }
      }
    }
    ##############
    # Build plot #
    ##############
    order.plot <- unique(df$isoform.name)
    if (!is.null(names.order)){
      order.plot <- names.order
    }
    
    if (line.overhang < 0) {
      line.overhang <- 0.05
      warning(ERROR.overh)
    }
    if (line.overhang >= 1) {
      x.corr <- line.overhang
    } else {
      x.corr <- max(df$x, na.rm = T) * line.overhang
    }
    
    df.segment <- data.frame("x" = NA, "xend" = NA, "y" = NA)
    corr <- 0
    for (g in order.plot){
      y.corr <- (0.5 * corr)
      corr <- corr + 1
      y.isoform <- y.pos + y.corr
      temp <- df[which(df$isoform.name == g),]
      df.segment <- rbind(df.segment,
                          c(min(temp$x, na.rm = T) - x.corr,
                            max(temp$x, na.rm = T) + x.corr,
                            y.isoform))
    }
    
    df.segment <- df.segment[which(!is.na(df.segment$x)),]
    
    p <- ggplot(data = df.segment)
    p <- p + geom_segment(aes(x=x, xend=xend, y=y, yend=y), 
                          color = line.color, size = line.size)
    
    # Add exon blocks
    corr <- 0
    for (g in order.plot){
      y.corr <- (0.5 * corr)
      corr <- corr + 1
      # mRNA blocks
      data.mRNA <- df[which(df$feature == "mRNA" & df$isoform.name == g),]
      data.mRNA$y <- data.mRNA$y + y.corr
      p <- p + geom_polygon(data=data.mRNA, aes(x=x, y=y), fill = intron.color)
      # CDS blocks
      data.CDS <- df[which(df$feature == "CDS" & df$isoform.name == g),]
      data.CDS$y <- data.CDS$y + y.corr
      p <- p + geom_polygon(data=data.CDS, aes(x=x, y=y), fill = exon.color)
    }
    # ggplot settings
    p <- p + theme(line = element_blank(), axis.title = element_blank(),
                   axis.text = element_blank(), 
                   plot.title = element_text(hjust = 0.5),
                   panel.background = element_blank()) 
    p <- p + labs(title = gene.name) 
  } else {
    # correct exon size to make it visible
    if (min.size != 0){
      if (min.size < 0) {
        min.size <- 0.001
        warning(ERROR.min.size)
      }
      tot.size <- (max(df$x) - min(df$x)) + 1
      if (min.size >= 1) {
        size.min <- min.size
      } else {
        size.min <- floor(tot.size * min.size)
      }
      for (ii in seq(1, length(df$x), 4)){
        block.size <- (df$x[ii+2] - df$x[ii]) + 1
        if (block.size < size.min){
          diff <- size.min - block.size
          df$x[ii] <- df$x[ii] - floor(diff / 2)
          df$x[ii+1] <- df$x[ii+1] - floor(diff / 2)
          df$x[ii+2] <- df$x[ii+2] + floor(diff / 2)
          df$x[ii+3] <- df$x[ii+3] + floor(diff / 2)
        }
      }
    }
    # mRNA order
    if (is.null(mRNA.order)){
      mRNA.order <- c(1:max(df$isoform, na.rm = T))
    }
    # CDS order
    if (is.null(CDS.order)){
      CDS.order <- c(1:max(df$isoform, na.rm = T))
    }
    ##############
    # Build plot #
    ##############
    # "genome" line
    if (line.overhang < 0) {
      line.overhang <- 0.05
      warning(ERROR.overh)
    }
    if (line.overhang >= 1) {
      x.corr <- line.overhang
    } else {
      x.corr <- max(df$x, na.rm = T) * line.overhang
    }
    
    df.segment <- data.frame("x" = NA, "xend" = NA, "y" = NA)
    corr <- 0
    for (g in mRNA.order){
      y.corr <- (0.5 * corr)
      corr <- corr + 1
      y.isoform <- y.pos + y.corr
      data.mRNA <- df[which(df$isoform == g),]
      df.segment <- rbind(df.segment,
                          c(min(data.mRNA$x, na.rm = T) - x.corr,
                            max(data.mRNA$x, na.rm = T) + x.corr,
                            y.isoform))
    }
    
    df.segment <- df.segment[which(!is.na(df.segment$x)),]
    
    p <- ggplot(data = df.segment)
    p <- p + geom_segment(aes(x=x, xend=xend, y=y, yend=y), 
                          color = line.color, size = line.size)
    
    # Add mRNA blocks
    corr <- 0
    for (g in mRNA.order){
      y.corr <- (0.5 * corr)
      corr <- corr + 1
      data.mRNA <- df[which(df$feature == "mRNA" & df$isoform == g),]
      data.mRNA$y <- data.mRNA$y + y.corr
      p <- p + geom_polygon(data=data.mRNA, aes(x=x, y=y), fill = intron.color)
    }
    
    # Add CDS blocks
    corr <- 0
    for (g in CDS.order){
      y.corr <- (0.5 * corr)
      corr <- corr + 1
      data.CDS <- df[which(df$feature == "CDS" & df$isoform == g),]
      data.CDS$y <- data.CDS$y + y.corr
      p <- p + geom_polygon(data=data.CDS, aes(x=x, y=y), fill = exon.color)
    }
    
    # ggplot settings
    p <- p + theme(line = element_blank(), axis.title = element_blank(),
                   axis.text = element_blank(), 
                   plot.title = element_text(hjust = 0.5),
                   panel.background = element_blank()) 
    p <- p + labs(title = gene.name) 
  }
  
  #############
  # add bar #
  #############
  if (bar){
    x.bar <- min(df$x, na.rm = T) + max(df$x, na.rm = T) * 0.05
    exp <- floor(log10(max(df$x, na.rm = T))) - 1
    xend.bar <- x.bar + 10^exp
    if (bar.length > 0) {
      xend.bar <- x.bar + bar.length
    }
    label.bar <- xend.bar - x.bar 
    label.bar <- paste(label.bar, "bp")
    if (length(df.segment$x) - 2 < 10){
      y.bar.corr <- 0.4 / (10 - (length(df.segment$x) - 2))
    } else {
      y.bar.corr <- 0.2
    }
    if (bar.pos == "bottom"){
      y.bar <- 2.8 - y.bar.corr
      p <- p + geom_segment(aes(x = x.bar, xend = xend.bar, y = y.bar, 
                                yend = y.bar), 
                            color = bar.color, size = bar.size)
      p <- p + annotate("text", x = (x.bar+xend.bar)/2,
                        y = y.bar - y.bar.corr, label = label.bar)
    } else if (bar.pos == "top") {
      y.bar <- max(data.CDS$y, na.rm = T) + y.bar.corr
      p <- p + geom_segment(aes(x = x.bar, xend = xend.bar, y = y.bar,
                                yend = y.bar), 
                            color = bar.color, size = bar.size)
      p <- p + annotate("text", x = (x.bar+xend.bar)/2, 
                        y = y.bar + y.bar.corr, label = label.bar)
    } else {
      warning(ERROR.bar)
      return (p)
    }
  }
  return (p)
}
