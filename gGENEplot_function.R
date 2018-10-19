gGENEplot <- function(source.file, size.percent = 0.001, line.color = "gray", 
                      line.size = 2, line.overhang = 0.05, intron.color = "lightblue", exon.color = "blue",
                      ruler = F, ruler.pos = "bottom", ruler.color = "black", ruler.size = line.size/2, ruler.length = F,
                      mRNA.order = NULL, CDS.order = NULL){
  require(gsubfn)
  require(ggplot2)
  
  input <- readLines(source.file)
  initdebut <- which(substring(input, 1, 8) == "FEATURES") + 1
  gene.lines <- which(grepl("^[ ]{5}gene", input))
  mRNA.lines <- which(grepl("^[ ]{5}mRNA", input))
  CDS.lines <- which(grepl("^[ ]{5}CDS", input))
  
  #  Get gene names from gene.lines
  gene.name.list <- c()
  for (jj in gene.lines){
    line <- input[jj+1]
    gene.name <- strapplyc(line, "^[ ]+/gene=\\\"(.*)\\\"", simplify = TRUE)
    gene.name.list <- c(gene.name.list, gene.name)
  }
  
  df <- data.frame("gene" = NA, "isoform" = NA, "feature" = NA, "x" = NA, "y" = NA)
  y.pos <- 3
  
  index <- 0
  mRNA <- c("")
  for (i in mRNA.lines){
    # print(i)
    index <- index + 1
    line.num <- i
    line <- input[line.num]
    # print(line)
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
    # print(mRNA[index])
    mRNA.blocks <- strsplit(mRNA[index], ",")[[1]]
    for (j in 1:length(mRNA.blocks)){
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
                  c(gene.name, index, "mRNA", init, y.pos - 0.2), 
                  c(gene.name, index, "mRNA", init, y.pos + 0.2),
                  c(gene.name, index, "mRNA", end, y.pos + 0.2), 
                  c(gene.name, index, "mRNA", end, y.pos - 0.2))
    }
  }
  
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
    for (j in 1:length(CDS.blocks)){
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
                  c(gene.name, index, "CDS", init, y.pos - 0.2), 
                  c(gene.name, index, "CDS", init, y.pos + 0.2),
                  c(gene.name, index, "CDS", end, y.pos + 0.2), 
                  c(gene.name, index, "CDS", end, y.pos - 0.2))
    }
  }
  
  df <- df[which(!is.na(df$gene)),]
  df$x <- as.numeric(df$x)
  df$y <- as.numeric(df$y)
  df$isoform <- as.numeric(df$isoform)
  
  
  # correct block size to make it visible
  if (size.percent != 0){
    tot.size <- (max(df$x) - min(df$x)) + 1
    size.min <- floor(tot.size * size.percent)
    for (ii in seq(1, length(df$x), 4)){
      block.size <- (df$x[ii+2] - df$x[ii]) + 1
      if (block.size < size.min){
        df$x[ii+2] <- df$x[ii] + size.min
        df$x[ii+3] <- df$x[ii] + size.min
      }
    }
  }
  
  # segment.order <- c(1:max(df$isoform, na.rm = T))
  
  if (is.null(mRNA.order)){
    mRNA.order <- c(1:max(df$isoform, na.rm = T))
  } 
  
  if (is.null(CDS.order)){
    CDS.order <- c(1:max(df$isoform, na.rm = T))
  } 
  
  x.corr <- max(df$x, na.rm = T) * line.overhang
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
  p <- p + geom_segment(aes(x=x, xend=xend, y=y, yend=y), color = line.color, size = line.size)

  corr <- 0
  for (g in mRNA.order){
    y.corr <- (0.5 * corr)
    corr <- corr + 1
    data.mRNA <- df[which(df$feature == "mRNA" & df$isoform == g),]
    data.mRNA$y <- data.mRNA$y + y.corr
    p <- p + geom_polygon(data=data.mRNA, aes(x=x, y=y), fill = intron.color)
  }
  
  corr <- 0
  for (g in CDS.order){
    y.corr <- (0.5 * corr)
    corr <- corr + 1
    data.CDS <- df[which(df$feature == "CDS" & df$isoform == g),]
    data.CDS$y <- data.CDS$y + y.corr
    p <- p + geom_polygon(data=data.CDS, aes(x=x, y=y), fill = exon.color)
  }  
  
  p <- p + theme(line = element_blank(), axis.title = element_blank(), 
                 axis.text = element_blank(), plot.title = element_text(hjust = 0.5),
                 panel.background = element_blank()) + 
           labs(title = gene.name)
  
  # add ruler
  if (ruler){
    x.ruler <- max(df$x, na.rm = T) * 0.05
    if (!ruler.length){
      exp <- floor(log10(max(df$x, na.rm = T))) - 1
      xend.ruler <- x.ruler + 10^exp
    } else {
      xend.ruler <- x.ruler + ruler.length
    }
    label.ruler <- xend.ruler - x.ruler 
    label.ruler <- paste(label.ruler, "bp")
    # y.ruler.corr <- (max(data.CDS$y, na.rm = T) - 3)/25
    # y.ruler.corr <- 0.2
    if (length(df.segment$x) - 2 < 10){
      y.ruler.corr <- 0.4 / (10 - (length(df.segment$x) - 2))
    } else {
      y.ruler.corr <- 0.2
    }
    if (ruler.pos == "bottom"){
      y.ruler <- 2.8 - y.ruler.corr
      p <- p + geom_segment(aes(x=x.ruler, xend=xend.ruler, y=y.ruler, yend=y.ruler), 
                            color = ruler.color, size = ruler.size)
      p <- p + annotate("text", x = (x.ruler+xend.ruler)/2, y = y.ruler - y.ruler.corr,
                        label = label.ruler)
    } else if (ruler.pos == "top") {
      y.ruler <- max(data.CDS$y, na.rm = T) + y.ruler.corr
      p <- p + geom_segment(aes(x=x.ruler, xend=xend.ruler, y=y.ruler, yend=y.ruler), 
                            color = ruler.color, size = ruler.size)
      p <- p + annotate("text", x = (x.ruler+xend.ruler)/2, y = y.ruler + y.ruler.corr,
                        label = label.ruler)
    } else {
      warning("ruler.pos argument must be 'bottom' or 'top'. Ruler is not shown.")
      return (p)
    }
  }
  return (p)
}

test <- gGENEplot("zfGHRb.gb")

ggsave("zfGHRb_0.001.tiff", plot = test, width = 30, height = 3, units = "cm", dpi = 300)
ggsave("zfGHRb_0.001_ggsavetest.tiff", test)
ggsave("zfGHRb_0.001_ggsavetest2.tiff", test, width = 30, height = 3, units = "in")
ggsave("zfGHRb_0.000_ggsavetest.tiff", test)
ggsave("zfGHRb_0.001_xcorr.tiff", plot = test, width = 30, height = 3, units = "cm", dpi = 300)
ggsave("zfGHRb_0.001_xcorr_ratio.tiff", plot = test)
ggsave("zfGHRb_0.001_xcorr_30_5.tiff", plot = test, width = 30, height = 5, units = "cm", dpi = 300)

test <- gGENEplot("humanGHR.gb", line.size = 1.5)
ggsave("humanGHR_0.001_xcorr_ggsave_line1.5.tiff", plot = test)


gGENEplot("humanGHR.gb", ruler = T, mRNA.order = c(3, 2, 4:9, 1, 11, 10))
ggsave("humanGHR_ordered.tiff", width = 30, height = 15, units = "cm", dpi = 300)
gGENEplot("humanGHR.gb", ruler = T, mRNA.order = c(10, 11, 1, 9:4, 3, 2), CDS.order = c(11:1), line.overhang = 0.01)
ggsave("humanGHR_ordered_opposite_overhang0.01.tiff", width = 30, height = 15, units = "cm", dpi = 300)
gGENEplot("humanGHR.gb", ruler = T, CDS.order = c(9,2,1,3:8,11,10), line.overhang = 0.01)