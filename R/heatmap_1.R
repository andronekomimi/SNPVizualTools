suppressMessages(library(GenomicRanges))
suppressMessages(library(ggplot2))
suppressMessages(library(plyr))
suppressMessages(library(reshape2))

create_gr_from_df_4_heatmap_range = function(DF, current_range) {
  
  chrom_start_idx = 1
  left_start_idx = 2
  chrom_end_idx = 4
  right_start_idx = 5
  left_end_idx = 3
  right_end_idx = 6
  
  left_start = c()
  left_end = c()
  right_start = c()
  right_end = c()
  chrom_start = c()
  chrom_end =  c()
  
  for (i in (1:nrow(DF)))
  {
    chrom_start_value =  sub("chr","",DF[i,chrom_start_idx])
    chrom_end_value =  sub("chr","",DF[i,chrom_end_idx])
    
    left_start_value = as.numeric(DF[i,left_start_idx])
    left_end_value = as.numeric(DF[i,left_end_idx])
    
    right_start_value = as.numeric(DF[i,right_start_idx])
    right_end_value = as.numeric(DF[i,right_end_idx])
    
    if(chrom_start_value == as.numeric(seqlevels(current_range))) {
      if((left_start_value > (start(current_range) - 1) && left_start_value < (end(current_range) - 1)) || (right_end_value > (start(current_range) - 1) && right_end_value < (end(current_range) - 1)))
      {
        chrom_start = c(chrom_start, chrom_start_value)
        chrom_end = c(chrom_end, chrom_end_value)
        left_start = c(left_start, left_start_value)
        left_end = c(left_end, left_end_value)
        right_start = c(right_start, right_start_value)
        right_end = c(right_end, right_end_value)
      }
    }    
  }
  
  left_ranges = IRanges(left_start,left_end)
  right_ranges = IRanges(right_start,right_end)
  
  my_ranges = NULL
  
  if(length(left_ranges) != 0) {
    my_ranges <- GRanges(seqnames = chrom_start, ranges = left_ranges)
    right_range <- GRanges(seqnames = chrom_end, ranges = right_ranges)
    values(my_ranges)$to.gr <- right_range
  }
  
  my_ranges
}


transpose_range_col = function(current_range, x, out_limit, pas) {
  
  if (x > end(current_range)) {ret = out_limit}
  else {
    if (x < start(current_range)) { ret = 2}
    else {ret = ((x - start(current_range)) %/% pas) + 3}
  }
  ret
}

transpose_range_row = function(current_range, x, out_limit, pas) {
  
  if (x > end(current_range)) {ret = out_limit}
  else {
    if (x < start(current_range)) { ret = 1}
    else {ret = ((x - start(current_range)) %/% pas) + 2}
  }
  ret
}


fill_matrix = function(m, current_range, range, pas) {
  for (i in seq(start(range), end(range), by = pas)) {
    for (j in seq(start(range$to.gr), end(range$to.gr), by = pas)) {
     # print(paste0("i = ",i,", j = ",j))
      I = transpose_range_row(current_range, i, dim(m)[1], pas)
      J = transpose_range_col(current_range, j, dim(m)[1], pas)
      d = m[I,J] + 1
      m[I,J] = d
     # print(paste0("I = ",I,", J = ",J," = ",d))
    }
  }
  m
}

################################################################
color = rev(heat.colors(256))

current_chr <- "12"

#current_range <- IRanges(min(start(snps))-border, max(end(snps))+border)
current_range <- IRanges(27650000, 28400000)
current_range <- GRanges(seqnames = current_chr, ranges = current_range)

chiapet_k562_file1 <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_lane13_chr12"
df2 <- read.table(chiapet_k562_file1, header=FALSE, stringsAsFactors=FALSE)

chiapet_k562_file2 <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_lane24_chr12"
df3 <- read.table(chiapet_k562_file2, header=FALSE, stringsAsFactors=FALSE)

k562_df = rbind(df2,df3)

chiapet_mcf7_file1 <- "/home/nekomimi/Workspace/SNPVIZU/data/CHIAPET_MCF7_rep3_chr12"
df4 <- read.table(chiapet_mcf7_file1, header=FALSE, stringsAsFactors=FALSE)
data <- strsplit(gsub("\\.\\.", "-", df4[,4]), '[:,-]')
df4 =  matrix(unlist(data), ncol = 7, byrow = TRUE)

chiapet_mcf7_file2 <- "/home/nekomimi/Workspace/SNPVIZU/data/CHIAPET_MCF7_rep4_chr12"
df5 <- read.table(chiapet_mcf7_file2, header=FALSE, stringsAsFactors=FALSE)
data <- strsplit(gsub("\\.\\.", "-", df5[,4]), '[:,-]')
df5 =  matrix(unlist(data), ncol = 7, byrow = TRUE)

chiapet_mcf7_file3 <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_lane11_chr12"
df5a <- read.table(chiapet_mcf7_file3, header=FALSE, stringsAsFactors=FALSE)

chiapet_mcf7_file4 <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_lane23_chr12"
df5b <- read.table(chiapet_mcf7_file4, header=FALSE, stringsAsFactors=FALSE)

mcf7_df_1 = rbind(df4,df5)
mcf7_df_2 = rbind(df5a,df5b)

k562_gr = create_gr_from_df_4_heatmap_range(k562_df,current_range)
mcf7_gr_1 = create_gr_from_df_4_heatmap_range(mcf7_df_1,current_range)
mcf7_gr_2 = create_gr_from_df_4_heatmap_range(mcf7_df_2,current_range)

ranges = c(k562_gr, mcf7_gr_1, mcf7_gr_2)

range_width = width(current_range)
pas = 10000
m_dim = range_width %/% pas

# M = array(0, dim = c(m_dim+2,m_dim+2))
y = as.character(as.integer(seq.int(from = start(current_range), to = end(current_range), length.out = m_dim)))
DF = data.frame(matrix(0,ncol = m_dim+3, nrow = m_dim+2))

#row.names(DF) <- c("Start_out3",y,"out5")
colnames(DF) <- c("pos",(start(current_range) - 1000),y,(end(current_range) + 1000))
DF$pos =  c((start(current_range) - 1000),y,(end(current_range) + 1000))

for(i in seq_along(ranges)){
  print(paste0(i,"eme range"))
  DF =fill_matrix(DF,current_range, ranges[i], pas)
}


#heatmap(as.matrix(DF), Rowv = NA, Colv = "Rowv", col = color, cexRow = .4, cexCol = .4)
DF.m = melt(DF)
DF.m = ddply(DF.m, .(variable))
p = ggplot(DF.m, aes(variable, pos)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient(low = "white", high="#6FDE01" )+ theme(text = element_text(size=10), axis.text.x = element_text(angle=90, vjust=1))
p
