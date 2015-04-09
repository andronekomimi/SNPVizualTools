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

current_chr = "chr10"
current_start = 80710000
current_stop = 80992000
current_target = "ZMIZ1"
snps_list = "~/Workspace/COLLAB/ZMIZ1_snps"

# current_chr = "chr12"
# current_start = 27950000
# current_stop = 28735000
# current_target = "PTHLH"
# snps_list = "~/Workspace/COLLAB/PTHLH_snps"


# current_chr = "chr17"
# current_start = 53000000
# current_stop = 53500000
# current_target = "COX11"
# snps_list = "~/Workspace/COLLAB/COX11_snps"

highlight_region = args[6] # facultatif
if (is.na(highlight_region)) {
  highlight_region = "0:0"
}
######################### LOAD LIBRARY ###############################
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(ggbio))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressMessages(library(biovizBase))
suppressMessages(library(org.Hs.eg.db))

######################### FUNCTIONS ###############################

get_chiapet_arch <- function(rep) {
  
  g = ggplot(rep) +
    geom_arch() +
    theme_bw() +
    #aes(size=base_width, color=color) +
    #aes(size = alpha) +
    #aes(color=color, alpha = alpha) +
    aes(color=color) +
    scale_colour_manual(values = c("gray","promoter"="black")) +
    xlim(current_range) +
    theme(axis.ticks = element_blank(), axis.text.y = element_blank()) +
    ylab("") + guides(alpha=FALSE, color=FALSE)
  
  if (!is.null(highlight_region)) {
    
    highlight_region = strsplit(x = highlight_region, fixed = T, split = ":")[[1]]
    d = data.frame(x1=as.numeric(highlight_region[1]), x2=as.numeric(highlight_region[2]), y1=0, y2=10)
    g = g + geom_rect(d, mapping=aes(xmin=x1, xmax=x2,ymin=y1, ymax=y2), color="black", alpha = 0.4) + xlim(current_range)
  }
  
  g
  
}

create_gr_from_df = function(DF,start_idx,end_idx,label) {
  label = ""
  left = c()
  right = c()
  base_width = c()
  color = c()
  alpha = c()	
  
  for (i in (1:nrow(DF)))
  {
    start_value = as.numeric(DF[i,start_idx])
    end_value = as.numeric(DF[i,end_idx])
    
    if((start_value > (start(current_range) - 1) 
        && start_value < (end(current_range) - 1)) 
       || (end_value > (start(current_range) - 1) 
           && end_value < (end(current_range) - 1))) {
      
      left = c(left, start_value)
      right = c(right, end_value)
      base_width = c(base_width ,as.numeric(DF[i,start_idx+1])-start_value)
      color = c(color,label)
    }
  }
  
  ranges = IRanges(left,right)
  
  if(length(left) > 0) {
    g_ranges <- GRanges(seqnames = current_chr, ranges = ranges, size = base_width, color = color, alpha = rep(0.5,length(left)))
    print(paste0("----Returning a ",length(g_ranges),"-long Gr"))
    invisible(g_ranges)
  }
  else
  {
    invisible(GRanges())
  }
}

create_gr_from_list = function(frame_data,start_idx,end_idx,label) {
  left = c()
  right = c()
  base_width = c()
  color = c()
  alpha = c()	
  
  
  for (i in (1: length(frame_data))){
    start_value = as.numeric(as.character(frame_data[[i]][start_idx]))
    end_value = as.numeric(as.character(frame_data[[i]][end_idx]))
    if((start_value > (start(current_range) - 1) 
        && start_value < (end(current_range) - 1)) 
       || (end_value > (start(current_range) - 1) 
           && end_value < (end(current_range) - 1))) {
      
      left = c(left, start_value)
      right = c(right, end_value)
      base_width = c(base_width ,(as.numeric(as.character(frame_data[[i]][start_idx+1]))-start_value))
      color = c(color,label)
      # in the spec zone?
      if((start_value > (start(special_range) - 1) 
          && start_value < (end(special_range) - 1)) 
         || (end_value > (start(special_range) - 1) 
             && end_value < (end(special_range) - 1))) {				
        
        alpha = c(alpha, 1)
      }
      else
      {
        alpha = c(alpha, 0,5)
      }
    }
  }
  
  ranges = IRanges(left,right)
  g_ranges <- GRanges(seqnames = current_chr, ranges = ranges, base_width = base_width, color = color, alpha = alpha)
  print(paste0("Returning a ",length(g_ranges),"-long Gr"))
  invisible(g_ranges)
}

make_emphasis = function(range, range_degree) {
  
  for (i in(1 : length(range))) {
    
    range_row = range[i]
    start_value = start(range_row)
    end_value = end(range_row)
    # in the spec zone?
    if((start_value > (start(range_degree) - 1) 
        && start_value < (end(range_degree) - 1)) 
       || (end_value > (start(range_degree) - 1) 
           && end_value < (end(range_degree) - 1))) {
      range[i]$alpha = range_degree$alpha
      range[i]$color = range_degree$color
    }
  }
  
  range
}

gr_list_transformation = function(g_range, labels) {
  
  gr = c()
  
  for(i in(1:length(g_range))) {
    left =  c()
    right =  c()
    metadata =  c()
    alpha = c()
    elt = g_range[[i]]
    left = start(elt$left_pair)
    right = end(elt$right_pair)
    metadata = elt$left_pair$PET_count
    ranges = IRanges(left,right)
    gr <- c(gr,GRanges(seqnames = current_chr, ranges = ranges, base_width = rep(10,length(left)), color=rep(labels[i],length(left), alpha = rep(0.5,length(left)))))
  }
  
  gr
}


######################### PROCESSING ###############################
##################### FILE LOADING ###############################################

hic_hmec_file_1 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/HIC_HMEC_Arrowhead_domainlist_", current_chr)
df1a <- read.table(hic_hmec_file_1, header=FALSE, stringsAsFactors=FALSE)

hic_hmec_file_2 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/HIC_HMEC_HiCCUPS_looplist_", current_chr)
df1b <- read.table(hic_hmec_file_2, header=FALSE, stringsAsFactors=FALSE)

hic_k562_file_1 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/HIC_K562_Arrowhead_domainlist_", current_chr)
df6a <- read.table(hic_k562_file_1, header=FALSE, stringsAsFactors=FALSE)

hic_k562_file_2 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/HIC_K562_HiCCUPS_looplist_", current_chr)
df6b <- read.table(hic_k562_file_2, header=FALSE, stringsAsFactors=FALSE)


############################ SINGAPORE DATA ############################
chiapet_mcf7_file3 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_lane11_", current_chr)
df5a <- read.table(chiapet_mcf7_file3, header=FALSE, stringsAsFactors=FALSE)

chiapet_mcf7_file4 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_lane23_", current_chr)
df5b <- read.table(chiapet_mcf7_file4, header=FALSE, stringsAsFactors=FALSE)

chiapet_k562_file1 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_lane13_", current_chr)
df2 <- read.table(chiapet_k562_file1, header=FALSE, stringsAsFactors=FALSE)

chiapet_k562_file2 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_lane24_", current_chr)
df3 <- read.table(chiapet_k562_file2, header=FALSE, stringsAsFactors=FALSE)

############################ STANFORD DATA ############################ 

chiapet_k562_file3 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_rep1_", current_chr)
df3a <- read.table(chiapet_k562_file3, header=FALSE, stringsAsFactors=FALSE)
data <- strsplit(gsub("\\.\\.", "-", df3a[,4]), '[:,-]')
df3a =  matrix(unlist(data), ncol = 7, byrow = TRUE)

chiapet_k562_file4 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_rep2_", current_chr)
df3b <- read.table(chiapet_k562_file4, header=FALSE, stringsAsFactors=FALSE)
data <- strsplit(gsub("\\.\\.", "-", df3b[,4]), '[:,-]')
df3b =  matrix(unlist(data), ncol = 7, byrow = TRUE)

chiapet_mcf7_file1 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_rep3_", current_chr)
df4 <- read.table(chiapet_mcf7_file1, header=FALSE, stringsAsFactors=FALSE)
data <- strsplit(gsub("\\.\\.", "-", df4[,4]), '[:,-]')
df4 =  matrix(unlist(data), ncol = 7, byrow = TRUE)

chiapet_mcf7_file2 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_rep4_", current_chr)
df5 <- read.table(chiapet_mcf7_file2, header=FALSE, stringsAsFactors=FALSE)
data <- strsplit(gsub("\\.\\.", "-", df5[,4]), '[:,-]')
df5 =  matrix(unlist(data), ncol = 7, byrow = TRUE)

chiapet_mcf7_file3 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_ER_rep1_", current_chr)
df7 <- read.table(chiapet_mcf7_file3, header=FALSE, stringsAsFactors=FALSE)

chiapet_mcf7_file4 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_ER_rep2_", current_chr)
df8 <- read.table(chiapet_mcf7_file4, header=FALSE, stringsAsFactors=FALSE)

chiapet_mcf7_file5 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_ER_rep3_", current_chr)
df9 <- read.table(chiapet_mcf7_file5, header=FALSE, stringsAsFactors=FALSE)

chiapet_mcf7_file6 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_CTCF_rep1_", current_chr)
df10 <- read.table(chiapet_mcf7_file6, header=FALSE, stringsAsFactors=FALSE)

chiapet_mcf7_file7 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_CTCF_rep2_", current_chr)
df11 <- read.table(chiapet_mcf7_file7, header=FALSE, stringsAsFactors=FALSE)

chiapet_k562_file3 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_CTCF_rep1_", current_chr)
df12 <- read.table(chiapet_k562_file3, header=FALSE, stringsAsFactors=FALSE)


############################ END OF STANDFORD DATA #########################################

snps_df <- read.table(snps_list, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")


### sup material
chiapet_hct116_file1 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_HCT116_rep1_", current_chr)
df15 <- read.table(chiapet_hct116_file1, header=FALSE, stringsAsFactors=FALSE)
data <- strsplit(gsub("\\.\\.", "-", df15[,4]), '[:,-]')
df15 =  matrix(unlist(data), ncol = 7, byrow = TRUE)

chiapet_hela_file1 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_HELAS3_rep1_", current_chr)
df16 <- read.table(chiapet_hela_file1, header=FALSE, stringsAsFactors=FALSE)
data <- strsplit(gsub("\\.\\.", "-", df16[,4]), '[:,-]')
df16 =  matrix(unlist(data), ncol = 7, byrow = TRUE)

chiapet_nb4_file1 <- paste0("/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_NB4_rep1_", current_chr)
df17 <- read.table(chiapet_nb4_file1, header=FALSE, stringsAsFactors=FALSE)
data <- strsplit(gsub("\\.\\.", "-", df17[,4]), '[:,-]')
df17 =  matrix(unlist(data), ncol = 7, byrow = TRUE)


print("Files loaded")


##################### SNPS ###############################################

snp_ids <- snps_df$snp_id
snps_ranges <- IRanges(snps_df$location, snps_df$location)
snps <- GRanges(seqnames = current_chr, ranges = snps_ranges, imp = snps_df$metadata)
snps$name <- snp_ids

##################### CURRENT RANGE ###############################################

#current_range <- IRanges(53000000,53500000) # cox11
#current_range <- IRanges(27950000, 28735000) # pthlh
current_range <- IRanges(current_start,current_stop) # default
current_range <- GRanges(seqnames = current_chr, ranges = current_range)

# special_range <- IRanges(28111017,28127138) # pthlh + promoter
# special_range <- GRanges(seqnames = current_chr, ranges = special_range, alpha = 0.6, color = "PTHLH promoter")
# special_range <- IRanges(53044646,53047646) # cox11 promoter -1500/+1500
special_range <- IRanges(80827292,80830292) # zmiz1 promoter -1500/+1500 
special_range <- GRanges(seqnames = current_chr, ranges = special_range, alpha = 0.6, color = "promoter")

# special_range_2 <- IRanges(28284682,28287682) #-1500/+1500 ccdc91 fow
# special_range_2 <- GRanges(seqnames = current_chr, ranges = special_range_2, alpha = 0.6, color = "CCDC91 promoter")


##################### HiC HMEC ###################################

hic_hmec_1 <- create_gr_from_df(df1a,2,3,"HiC-arrowhead")
hic_hmec_2 <- create_gr_from_df(df1b,2,6,"HiC-hiccups")

hic_hmec = c(hic_hmec_1, hic_hmec_2)


if(length(hic_hmec) > 0) {
  hic_hmec = make_emphasis(hic_hmec, special_range)
  hic_hmec_track = get_chiapet_arch(hic_hmec)
  print("Track hic_hmec_track -> DONE")
} else { print("No hic_hmec_track") }

####################### HiC K562 #################################

hic_k562_1 <- create_gr_from_df(df6a,2,3,"HiC-arrowhead")
hic_k562_2 <- create_gr_from_df(df6b,2,6,"HiC-hiccups")

hic_k562 = c(hic_k562_1, hic_k562_2)

if(length(hic_k562) > 0) {
  hic_k562 = make_emphasis(hic_k562, special_range)
  hic_k562_track = get_chiapet_arch(hic_k562)
  print("Track hic_k562_track -> DONE")
} else { print("No hic_k562_track") }

###################### CHIAPET K562 ##################################
chiapet_k562_rep1 <- create_gr_from_df(df3a,2,6,"ChIA-Pet rep 1")
chiapet_k562_rep2 <- create_gr_from_df(df3b,2,6,"ChIA-Pet rep 2")

chiapet_k562_lane13 <- create_gr_from_df(df2,2,6,"ChIA-Pet lane 13")
chiapet_k562_lane24 <- create_gr_from_df(df3,2,6,"ChIA-Pet lane 24")

#chiapet_k562 = c(chiapet_k562_rep1,chiapet_k562_rep2,chiapet_k562_lane13, chiapet_k562_lane24)
chiapet_k562 = c(chiapet_k562_rep1,chiapet_k562_rep2)

if(length(chiapet_k562) > 0) {
  chiapet_k562 = make_emphasis(chiapet_k562, special_range)
  chiapet_k562_track = get_chiapet_arch(chiapet_k562)
  print("Track chiapet_k562_track -> DONE")
} else { print("No chiapet_k562_track") }

####################### CHIAPET MCF7 #################################


chiapet_mcf7_rep3 <- create_gr_from_df(df4,2,6,"ChIA-Pet rep 3")
chiapet_mcf7_rep4 <- create_gr_from_df(df5,2,6,"ChIA-Pet rep 4")

chiapet_mcf7_lane11 <- create_gr_from_df(df5a,2,6,"ChIA-Pet lane 11")
chiapet_mcf7_lane23 <- create_gr_from_df(df5b,2,6,"ChIA-Pet lane 23")


chiapet_mcf7 = c(chiapet_mcf7_rep3,chiapet_mcf7_rep4, chiapet_mcf7_lane23)

if(length(chiapet_mcf7) > 0) {
  chiapet_mcf7 = make_emphasis(chiapet_mcf7, special_range)
  chiapet_mcf7_track = get_chiapet_arch(chiapet_mcf7)
  print("Track chiapet_mcf7_track -> DONE")
} else { print("No chiapet_mcf7_track") }

######################### CHIAPET MCF7 ER ###############################

chiapet_mcf7_er_rep1_data <- strsplit(gsub("\\.\\.", "-", df7[,4]), '[:,-]')
m =  matrix(unlist(chiapet_mcf7_er_rep1_data), ncol = 7, byrow = TRUE)
chiapet_mcf7_er_rep1 <- create_gr_from_df(m,2,6,"ChIA-Pet ER")

chiapet_mcf7_er_rep2_data <- strsplit(gsub("\\.\\.", "-", df8[,4]), '[:,-]')
m =  matrix(unlist(chiapet_mcf7_er_rep2_data), ncol = 7, byrow = TRUE)
chiapet_mcf7_er_rep2 <- create_gr_from_df(m,2,6,"ChIA-Pet ER")

chiapet_mcf7_er_rep3_data <- strsplit(gsub("\\.\\.", "-", df9[,4]), '[:,-]')
m =  matrix(unlist(chiapet_mcf7_er_rep3_data), ncol = 7, byrow = TRUE)
chiapet_mcf7_er_rep3 <- create_gr_from_df(m,2,6,"ChIA-Pet ER")

chiapet_mcf7_er = c(chiapet_mcf7_er_rep1,chiapet_mcf7_er_rep2, chiapet_mcf7_er_rep3)

if(length(chiapet_mcf7_er) > 0) {
  chiapet_mcf7_er = make_emphasis(chiapet_mcf7_er, special_range)
  chiapet_mcf7_er_track = get_chiapet_arch(chiapet_mcf7_er)
  print("Track chiapet_mcf7_er_track -> DONE")
} else {print("No ChiA-Pet MCF7 ER")}

######################### CHIAPET MCF7 CTCF ###############################

chiapet_mcf7_ctcf_rep1_data <- strsplit(gsub("\\.\\.", "-", df10[,4]), '[:,-]')
m =  matrix(unlist(chiapet_mcf7_ctcf_rep1_data), ncol = 7, byrow = TRUE)
chiapet_mcf7_ctcf_rep1 <- create_gr_from_df(m,2,6,"ChIA-Pet CTCF")

chiapet_mcf7_ctcf_rep2_data <- strsplit(gsub("\\.\\.", "-", df11[,4]), '[:,-]')
m =  matrix(unlist(chiapet_mcf7_ctcf_rep2_data), ncol = 7, byrow = TRUE)
chiapet_mcf7_ctcf_rep2 <- create_gr_from_df(m,2,6,"ChIA-Pet CTCF")

chiapet_mcf7_ctcf = c(chiapet_mcf7_ctcf_rep1,chiapet_mcf7_ctcf_rep2)

if(length(chiapet_mcf7_ctcf) > 0) {
  chiapet_mcf7_ctcf = make_emphasis(chiapet_mcf7_ctcf, special_range)
  chiapet_mcf7_ctcf_track = get_chiapet_arch(chiapet_mcf7_ctcf)
  print("Track chiapet_mcf7_ctcf_track -> DONE")
} else {print("No ChiA-Pet MCF7 CTCF")}

######################### CHIAPET K562 CTCF ###############################

chiapet_k562_ctcf_rep1_data <- strsplit(gsub("\\.\\.", "-", df12[,4]), '[:,-]')
m =  matrix(unlist(chiapet_k562_ctcf_rep1_data), ncol = 7, byrow = TRUE)
chiapet_k562_ctcf_rep1 <- create_gr_from_df(m,2,6,"ChIA-Pet CTCF")

chiapet_k562_ctcf = c(chiapet_k562_ctcf_rep1)

if(length(chiapet_k562_ctcf) > 0) {
  chiapet_k562_ctcf = make_emphasis(chiapet_k562_ctcf, special_range)
  chiapet_k562_ctcf_track = get_chiapet_arch(chiapet_k562_ctcf)
  print("Track chiapet_k562_ctcf_track -> DONE")
} else {print("No ChiA-Pet K562 CTCF")}

####################### CHIAPET HCT116 ####################################

chiapet_hct116 <- create_gr_from_df(df15,2,6,"ChIA-Pet rep 1")

if(length(chiapet_hct116) > 0) {
  chiapet_hct116_track = get_chiapet_arch(chiapet_hct116)
  print("Track chiapet_hct116_track -> DONE")
} else { print("No chiapet_hct116_track") }

####################### CHIAPET HELA S4 ####################################

chiapet_hela <- create_gr_from_df(df16,2,6,"ChIA-Pet rep 1")

if(length(chiapet_hela) > 0) {
  chiapet_hela_track = get_chiapet_arch(chiapet_hela)
  print("Track chiapet_hela_track -> DONE")
} else { print("No chiapet_hela_track") }

####################### CHIAPET NB4 ####################################

chiapet_nb4 <- create_gr_from_df(df17,2,6,"ChIA-Pet rep 1")

if(length(chiapet_nb4) > 0) {
  chiapet_nb4_track = get_chiapet_arch(chiapet_nb4)
  print("Track chiapet_nb4_track -> DONE")
} else { print("No chiapet_nb4_track") }

######################### SNPS TRACK ###############################
snps_track <- autoplot(snps, aes(color=imp)) +
  geom_text(aes(x = start, y = 1, label=name, angle = 90, vjust=-1), size = 1, color = "blue") +
  theme_bw() +
  xlim(current_range) + guides(color= FALSE)


print("Track snps_track -> DONE")


######################### ANNOTATIONS TRACK ###############################
gr_txdb <- crunch(TxDb.Hsapiens.UCSC.hg19.knownGene, which = current_range)
colnames(values(gr_txdb))[4] <- "model"
gr_txdb$symbols <- select(org.Hs.eg.db,
                          keys = as.character(gr_txdb$gene_id),
                          column = "SYMBOL",
                          keytype = "ENTREZID")$SYMBOL
i <- which(gr_txdb$model == "gap")
gr_txdb <- gr_txdb[-i]
levels(gr_txdb) <- c("cds", "exon", "utr")
grl_txdb <- split(gr_txdb, gr_txdb$symbols)


genes <- autoplot(grl_txdb, aes(type = model)) + theme_bw() + xlim(current_range)

print("Track gene -> DONE")

######################### PRINT TRACKS ###############################
#print(tracks( "K562 ChIAPET" = chiapet_k562_track,"MCF7 ChIAPET" = chiapet_mcf7_track, "SNPS"=snps_track, "Annotation" = genes , label.text.cex = 0.7))
#print(tracks( "K562 HiC" = hic_k562_track, "K562 ChIAPET" = chiapet_k562_track,"HMEC HiC" = hic_hmec_track, "MCF7 ChIAPET" = chiapet_mcf7_track, "SNPS"=snps_track, "Annotation" = genes , label.text.cex = 0.7))
#print(tracks( "K562 HiC" = hic_k562_track, "K562 ChIAPET" = chiapet_k562_track, "K562 ChIAPET CTCF" =  chiapet_k562_ctcf_track, "HMEC HiC" = hic_hmec_track, "MCF7 ChIAPET" = chiapet_mcf7_track, "MCF7 ChIAPET CTCF" = chiapet_mcf7_ctcf_track , "MCF7 ChIAPET ER" = chiapet_mcf7_er_track, "SNPS"=snps_track, "Annotation" = genes , label.text.cex = 0.7))