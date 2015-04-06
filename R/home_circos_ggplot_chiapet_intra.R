
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

produce_circos <- function(gr_left, gr_right, current_name, ideo, snps) {
  seqlevels(gr_left) <- seqlevels(ideo)
  seqlengths(gr_left) <- seqlengths(ideo)
  seqlevels(gr_right) <- seqlevels(ideo)
  seqlengths(gr_right) <- seqlengths(ideo)
  values(gr_left)$to.gr <- gr_right
  gr <- gr_left
  
  seqlevels(snps) <- seqlevels(ideo)
  seqlengths(snps) <- seqlengths(ideo)
  
  # preparation affichage coord
  coord_name = seq(from = start(current_range), to = end(current_range), by = 10000)
  coord_pos = coord_name - start(current_range)
  coord_range = IRanges(coord_pos, coord_pos)
  coord_name = GRanges(seqnames = as.character(seqlevels(current_range)), 
                       ranges = coord_range, label = coord_name)
  
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
  
  
  p <- p + layout_circle(coord_name, 
                         geom = "rect", 
                         vjust = 0, radius = 43, 
                         size = 0.1,
                         color = "black", 
                         trackWidth = 0.5)
  
  p <- p + layout_circle(coord_name, 
                         geom = "text", 
                         aes(label = label), 
                         vjust = 0, radius = 44, 
                         size = 1, 
                         trackWidth = 0.5)
  
  p <- p +
    layout_circle(ideo,
                  geom = "scale",
                  size = 2,
                  radius = 36,
                  space.skip = 0,
                  scale.type = "B",
                  trackWidth = 3)
  
  p <- p +
    layout_circle(gr,
                  geom = "link",
                  linked.to = "to.gr",
                  #aes(color = color, alpha = 0.1, size = base_width),
                  aes(linetype = out_range),
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
  out_range = c()
  
  
  for (i in (1:nrow(DF)))
  {
    chrom_start_value =  sub("chr","",DF[i,chrom_start_idx])
    chrom_end_value =  sub("chr","",DF[i,chrom_end_idx])
    
    left_start_value = as.numeric(DF[i,left_start_idx])
    left_end_value = as.numeric(DF[i,left_end_idx])
    
    right_start_value = as.numeric(DF[i,right_start_idx])
    right_end_value = as.numeric(DF[i,right_end_idx])
    
    
    if(chrom_start_value == as.numeric(seqlevels(current_range)) && chrom_end_value == as.numeric(seqlevels(current_range))) {
      
      
      if(left_end_value >= start(current_range) && right_start_value <= end(current_range))
      {
        chrom_start = c(chrom_start, chrom_start_value)
        chrom_end = c(chrom_end, chrom_end_value)
        left_start = c(left_start, transform_coord(left_start_value))
        left_end = c(left_end, transform_coord(left_end_value))
        right_start = c(right_start, transform_coord(right_start_value))
        right_end = c(right_end, transform_coord(right_end_value))
        base_width = c(base_width , left_end_value - left_start_value)
        color = c(color,label)
        out_range = c(out_range, FALSE)
      }
      else
      {
        if(chrom_start_value == as.numeric(seqlevels(current_range))) {
          if((left_start_value >= start(current_range) && left_start_value <= end(current_range)) || (left_end_value >= start(current_range) && left_end_value <= end(current_range)))
          {
            chrom_start = c(chrom_start, chrom_start_value)
            chrom_end = c(chrom_end, chrom_end_value)
            left_start = c(left_start, transform_coord(left_start_value))
            left_end = c(left_end, transform_coord(left_end_value))
            right_start = c(right_start, transform_coord(right_start_value))
            right_end = c(right_end, transform_coord(right_end_value))
            base_width = c(base_width , left_end_value - left_start_value)
            color = c(color,label)
            out_range = c(out_range, TRUE)
          }
        }
        else
        {
          if(chrom_end_value == as.numeric(seqlevels(current_range))) {
            if((right_start_value >= start(current_range) && right_start_value <= end(current_range)) || (right_end_value >= start(current_range) && right_end_value <= end(current_range)))
            {
              chrom_start = c(chrom_start, chrom_start_value)
              chrom_end = c(chrom_end, chrom_end_value)
              left_start = c(left_start, transform_coord(left_start_value))
              left_end = c(left_end, transform_coord(left_end_value))
              right_start = c(right_start, transform_coord(right_start_value))
              right_end = c(right_end, transform_coord(right_end_value))
              base_width = c(base_width , left_end_value - left_start_value)
              color = c(color,label)
              out_range = c(out_range, TRUE)
            }
          }
        }
      }
      
    }
  }
  
  out = list()
  left_ranges = IRanges(left_start,left_end)
  right_ranges = IRanges(right_start,right_end)
  
  if(length(left_ranges) != 0) {
    
    left_range <- GRanges(seqnames = chrom_start, ranges = left_ranges, base_width = base_width, color = color, out_range = out_range)
    right_range <- GRanges(seqnames = chrom_end, ranges = right_ranges, base_width = base_width, color = color, out_range = out_range)
    
    out = list(left_range, right_range)
  }
  
  out
}

transform_coord <- function(x) {
  x - start(current_range)
}


###################### PROCESSIND ########################
###################### LOAD FILES #######################

############################ SINGAPORE DATA ############################  
chiapet_k562_file1 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_lane13_chr", current_chr)
df1A <- read.table(chiapet_k562_file1, header=FALSE, stringsAsFactors=FALSE)

chiapet_k562_file2 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_lane24_chr", current_chr)
df2A <- read.table(chiapet_k562_file2, header=FALSE, stringsAsFactors=FALSE)

chiapet_mcf7_file3 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_lane11_chr", current_chr)
df3A <- read.table(chiapet_mcf7_file3, header=FALSE, stringsAsFactors=FALSE)

chiapet_mcf7_file4 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_lane23_chr", current_chr)
df4A <- read.table(chiapet_mcf7_file4, header=FALSE, stringsAsFactors=FALSE)

############################ STANFORD DATA ############################  

chiapet_k562_file3 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_rep1_chr", current_chr)
df1B <- read.table(chiapet_k562_file3, header=FALSE, stringsAsFactors=FALSE)
data <- strsplit(gsub("\\.\\.", "-", df1B[,4]), '[:,-]')
df1B =  matrix(unlist(data), ncol = 7, byrow = TRUE)

chiapet_k562_file4 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_rep2_chr", current_chr)
df2B <- read.table(chiapet_k562_file4, header=FALSE, stringsAsFactors=FALSE)
data <- strsplit(gsub("\\.\\.", "-", df2B[,4]), '[:,-]')
df2B =  matrix(unlist(data), ncol = 7, byrow = TRUE)

chiapet_mcf7_file1 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_rep3_chr", current_chr)
df3B <- read.table(chiapet_mcf7_file1, header=FALSE, stringsAsFactors=FALSE)
data <- strsplit(gsub("\\.\\.", "-", df3B[,4]), '[:,-]')
df3B =  matrix(unlist(data), ncol = 7, byrow = TRUE)

chiapet_mcf7_file2 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_rep4_chr", current_chr)
df4B <- read.table(chiapet_mcf7_file2, header=FALSE, stringsAsFactors=FALSE)
data <- strsplit(gsub("\\.\\.", "-", df4B[,4]), '[:,-]')
df4B =  matrix(unlist(data), ncol = 7, byrow = TRUE)


snps_df <- read.table(snps_list, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")

print("Files loaded")

data("CRC", package = "biovizBase")


##################### CURRENT RANGE ###############################################

current_range <- IRanges(current_start,current_stop)
current_range <- GRanges(seqnames = current_chr, ranges = current_range)
current_range_length = end(current_range) - start(current_range)


##################### SNPS ###############################################
snps_df <- snps_df[order(snps_df$location),]
snp_ids <- snps_df$snp_id
snps_ranges <- IRanges(transform_coord(snps_df$location), transform_coord(snps_df$location))
snps <- GRanges(seqnames = current_chr, ranges = snps_ranges, imp = snps_df$metadata, y = seq(from = 1, to = length(snp_ids)), pos = snps_df$location)
snps$name <- snp_ids


###################### CREATE GENOMIC RANGES ###################### 

k562A_df = rbind(df1A,df2A)
mcf7A_df = rbind(df3A,df4A)

k562A = create_gr_from_df_4_circos_range(k562A_df,"K562 - ChIAPET", current_range)
mcf7A = create_gr_from_df_4_circos_range(mcf7A_df,"MCF7 - ChIAPET", current_range)


k562B_df = rbind(df1B,df2B)
mcf7B_df = rbind(df3B,df4B)

k562B = create_gr_from_df_4_circos_range(k562B_df,"K562 - ChIAPET", current_range)
mcf7B = create_gr_from_df_4_circos_range(mcf7B_df,"MCF7 - ChIAPET", current_range)


##################### CREATE IDEO ############################

ideo_current_range = hg19sub[as.numeric(current_chr)]
seqlevels(ideo_current_range) = as.character(unique(seqnames(ideo_current_range)))
start(ideo_current_range) = transform_coord(start(current_range))
end(ideo_current_range) = transform_coord(end(current_range))
seqlengths(ideo_current_range) = width(current_range)

ideo_current_chr = hg19sub[as.numeric(current_chr)]
seqlevels(ideo_current_chr) = as.character(unique(seqnames(ideo_current_chr)))

border = 100000
ideo_big_range = hg19sub[as.numeric(current_chr)]
seqlevels(ideo_big_range) = as.character(unique(seqnames(ideo_big_range)))
start(ideo_big_range) = transform_coord(start(current_range)) - border
end(ideo_big_range) = transform_coord(end(current_range)) + border
seqlengths(ideo_big_range) = width(ideo_big_range)


#################### RUN CIRCOS ##############################
# if(length(k562A) == 2) {
chiapet_k562A_on_range = suppressWarnings(produce_circos(k562A[[1]], k562A[[2]] , "ChIA-PET K562 - Singapore (range level)" , ideo_current_range , snps))
  chiapet_k562A_on_range$plot$labels$colour = "SNP IDs"
#   chiapet_k562_on_chr = suppressWarnings(produce_circos(k562A[[1]], k562A[[2]] , "ChIA-PET K562 - Singapore (chromosome level)" , ideo_current_chr , snps))
#   chiapet_k562_on_chr$plot$labels$colour = "SNP IDs"
# }
# 
# if(length(mcf7A) == 2) {
chiapet_mcf7A_on_range =  suppressWarnings(produce_circos(mcf7A[[1]], mcf7A[[2]] , "ChIA-PET MCF7 - Singapore (range level)" , ideo_current_range , snps))
  chiapet_mcf7A_on_range$plot$labels$colour = "SNP IDs"
#   chiapet_mcf7_on_chr = suppressWarnings(produce_circos(mcf7A[[1]], mcf7A[[2]] , "ChIA-PET MCF7 - Singapore (chromosome level)" , ideo_current_chr , snps))
#   chiapet_mcf7_on_chr$plot$labels$colour = "SNP IDs"
# }
# 
# if(length(k562B) == 2) {
chiapet_k562B_on_range = suppressWarnings(produce_circos(k562B[[1]], k562B[[2]] , "ChIA-PET K562 - Stanford (range level)" , ideo_current_range , snps))
   chiapet_k562B_on_range$plot$labels$colour = "SNP IDs"
#   chiapet_k562_on_chr = suppressWarnings(produce_circos(k562B[[1]], k562B[[2]] , "ChIA-PET K562 - Stanford (chromosome level)" , ideo_current_chr , snps))
#   chiapet_k562_on_chr$plot$labels$colour = "SNP IDs"
# }
# 
# if(length(mcf7B) == 2) {
chiapet_mcf7B_on_range =  suppressWarnings(produce_circos(mcf7B[[1]], mcf7B[[2]] , "ChIA-PET MCF7 - Stanford (range level)" , ideo_current_range , snps))
  chiapet_mcf7B_on_range$plot$labels$colour = "SNP IDs"
#   chiapet_mcf7_on_chr = suppressWarnings(produce_circos(mcf7B[[1]], mcf7B[[2]] , "ChIA-PET MCF7 - Stanford (chromosome level)" , ideo_current_chr , snps))
#   chiapet_mcf7_on_chr$plot$labels$colour = "SNP IDs"
# }



