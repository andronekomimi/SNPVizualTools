
######################### ARGUMENTS ###############################
args <- commandArgs(trailingOnly = TRUE)
# 
# if(length(args) < 5) {
#   stop("Argument missing! Usage : scrip.R chrX start stop target_name path_to_snp_list [highlight_region fac. start:end]")
# }
# 
# 
# current_chr = args[1]
# current_start = as.numeric(args[2])
# current_stop = as.numeric(args[3])
# current_target = args[4]
# snps_list = args[5]

# 
# current_chr = "10"
# current_start = 80710000
# current_stop = 80992000
# current_target = "ZMIZ1"
# snps_list = "~/Workspace/COLLAB/ZMIZ1_snps"

current_chr = "17"
current_start = 53000000
current_stop = 53500000
current_target = "COX11"
snps_list = "~/Workspace/COLLAB/COX11_snps"

######################### LOAD LIBRARY ###############################

suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(ggbio))

######################### FUNCTIONS ###############################

produce_circos_whole_genome <- function(gr_left, gr_right, current_name, ideo, snps) {
  seqlevels(gr_left) <- seqlevels(ideo)
  seqlengths(gr_left) <- seqlengths(ideo)
  seqlevels(gr_right) <- seqlevels(ideo)
  seqlengths(gr_right) <- seqlengths(ideo)
  values(gr_left)$to.gr <- gr_right
  gr <- gr_left
  
  seqlevels(snps) <- seqlevels(ideo)
  seqlengths(snps) <- seqlengths(ideo)
  
  p <- ggplot() +
    layout_circle(ideo,
                  geom = "ideo",
                  radius = 36,
                  fill = "gray",
                  space.skip = 0,
                  trackWidth = 3)
  
  p <- p + layout_circle(ideo,
                   geom = "ideo",
                   radius = 39,
                   fill = "black",
                   space.skip = 0,
                   trackWidth = 4)
  
  p <- p +
    layout_circle(ideo,
                  geom = "scale",
                  size = 2,
                  radius = 36,
                  space.skip = 0,
                  trackWidth = 3)
  p <- p +
    layout_circle(ideo,
                  geom = "text",
                  aes(label = seqnames),
                  vjust = 0,
                  radius = 40,
                  space.skip = 0,
                  trackWidth = 7)
  
  p <- p +
    layout_circle(gr,
                  geom = "link",
                  linked.to = "to.gr",
                  color = "cyan4",
                  radius = 30,
                  space.skip = 0,
                  trackWidth = 0.5)
  
  
  p <- p + ggtitle(current_name)
  
  p <- p +
    layout_circle(snps,
                  geom = "point",
                  aes(y = y, color = name),
                  grid = TRUE,
                  radius = 32,
                  size = 2,
                  space.skip = 0,
                  trackWidth = 4)
  
  
  p
}


create_gr_from_df_4_circos_range = function(DF,label, current_range) {
  
  chrom_start_idx = 1
  left_start_idx = 2
  chrom_end_idx = 4
  right_start_idx = 5
  left_end_idx = 3
  right_end_idx = 6
  
  chrom_start = c()
  chrom_end =  c()
  left_start = c()
  left_end = c()
  right_start = c()
  right_end = c()
  base_width = c()
  color = c()
  
  for (i in (1:nrow(DF)))
  {
    chrom_start_value =  sub("chr","",DF[i,chrom_start_idx])
    chrom_end_value =  sub("chr","",DF[i,chrom_end_idx])
    
    left_start_value = as.numeric(DF[i,left_start_idx])
    left_end_value = as.numeric(DF[i,left_end_idx])
    
    right_start_value = as.numeric(DF[i,right_start_idx])
    right_end_value = as.numeric(DF[i,right_end_idx])
    
    if(chrom_start_value == as.numeric(seqlevels(current_range))) {
      if((left_start_value >= start(current_range) && left_start_value <= end(current_range)) || (left_end_value >= start(current_range) && left_end_value <= end(current_range)))
      {
        chrom_start = c(chrom_start, chrom_start_value)
        chrom_end = c(chrom_end, chrom_end_value)
        left_start = c(left_start, left_start_value)
        left_end = c(left_end, left_end_value)
        right_start = c(right_start, right_start_value)
        right_end = c(right_end, right_end_value)
        base_width = c(base_width , left_end_value - left_start_value)
        color = c(color,label)
      }
    }
    else
    {
      if(chrom_end_value == as.numeric(seqlevels(current_range))) {
        if((right_start_value >= start(current_range) && right_start_value <= end(current_range)) || (right_end_value >= start(current_range) && right_end_value <= end(current_range)))
        {
          chrom_start = c(chrom_start, chrom_start_value)
          chrom_end = c(chrom_end, chrom_end_value)
          left_start = c(left_start, left_start_value)
          left_end = c(left_end, left_end_value)
          right_start = c(right_start, right_start_value)
          right_end = c(right_end, right_end_value)
          base_width = c(base_width , left_end_value - left_start_value)
          color = c(color,label)
        }
      }
    }
    
  }
  
  out = list()
  left_ranges = IRanges(left_start,left_end)
  right_ranges = IRanges(right_start,right_end)
  
  if(length(left_ranges) != 0) {
    
    left_range <- GRanges(seqnames = chrom_start, ranges = left_ranges, base_width = base_width, color = color)
    right_range <- GRanges(seqnames = chrom_end, ranges = right_ranges, base_width = base_width, color = color)
    
    out = list(left_range, right_range)
  }
  
  out
}

################ LOAD FILES #######################

############################ SINGAPORE DATA ############################  
chiapet_k562_file1 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_lane13_chr",current_chr,"_inter")
df2 <- read.table(chiapet_k562_file1, header=FALSE, stringsAsFactors=FALSE)
data <- strsplit(gsub("\\.\\.", "-", df2[,4]), '[:,-]')
df2 =  matrix(unlist(data), ncol = 7, byrow = TRUE)

chiapet_k562_file2 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_lane24_chr",current_chr,"_inter")
df3 <- read.table(chiapet_k562_file2, header=FALSE, stringsAsFactors=FALSE)
data <- strsplit(gsub("\\.\\.", "-", df3[,4]), '[:,-]')
df3 =  matrix(unlist(data), ncol = 7, byrow = TRUE)


chiapet_mcf7_file1 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_lane11_chr",current_chr,"_inter")
df4 <- read.table(chiapet_mcf7_file1, header=FALSE, stringsAsFactors=FALSE)
data <- strsplit(gsub("\\.\\.", "-", df4[,4]), '[:,-]')
df4 =  matrix(unlist(data), ncol = 7, byrow = TRUE)

chiapet_mcf7_file2 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_lane23_chr",current_chr,"_inter")
df5 <- read.table(chiapet_mcf7_file2, header=FALSE, stringsAsFactors=FALSE)
data <- strsplit(gsub("\\.\\.", "-", df5[,4]), '[:,-]')
df5 =  matrix(unlist(data), ncol = 7, byrow = TRUE)

snps_df <- read.table(snps_list, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")



###################### PROCESSIND ########################

data("CRC", package = "biovizBase")

##################### SNPS ###############################################

snp_ids <- snps_df$snp_id
snps_ranges <- IRanges(snps_df$location, snps_df$location)
snps <- GRanges(seqnames = current_chr, ranges = snps_ranges, imp = snps_df$metadata, y = seq(from = 1, to = length(snp_ids)))
snps$name <- snp_ids

##################### CURRENT RANGE ###############################################

current_range <- IRanges(current_start,current_stop)
current_range <- GRanges(seqnames = current_chr, ranges = current_range)
current_range_length = end(current_range) - start(current_range)


k562_df = rbind(df2,df3)
mcf7_df = rbind(df4,df5)


range_k562 = create_gr_from_df_4_circos_range(k562_df,"K562 - ChIAPET", current_range)
range_mcf7 = create_gr_from_df_4_circos_range(mcf7_df,"MCF7 - ChIAPET", current_range)

# current range VS genome
if(length(range_k562) == 2) {
  if(length(seqlevels(range_k562[[1]])) == 1 && length(seqlevels(range_k562[[2]])) == 1 ) {
    warning("No interchromosomic interaction for K562")
  }
  else
  {
    chiapet_k562_range = produce_circos_whole_genome(range_k562[[1]], range_k562[[2]] , "ChIA-PET K562" , hg19sub , snps)
  }
}

if(length(range_mcf7) == 2) {
  if(length(seqlevels(range_mcf7[[1]])) == 1 && length(seqlevels(range_mcf7[[2]])) == 1 ) {
    warning("No interchromosomic interaction for MCF7")
  }
  else
  {
    chiapet_mcf7_range = produce_circos_whole_genome(range_mcf7[[1]], range_mcf7[[2]] , "ChIA-PET MCF7" , hg19sub , snps)
  }
  
}
