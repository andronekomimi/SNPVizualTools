#' Download data for a current_chromosome
#' After setting current_chr variable : current_chr <- "chr1"
#' 
#' @return list elements containing all data available for a particular 
#' chromosome
#'
#' @import GenomicRanges
#' @examples
#' current_chr <- "chr12"
#' my.data <- loadChrData()
#'  
#' @export
loadChrData <- function() {
  
  can_run()
  
  ### chiapet 
  sgp_k562_lane13 <- readRDS(paste0("data/",current_chr,"_sgp_k562_lane13.Rda"))
  sgp_k562_lane24 <- readRDS(paste0("data/",current_chr,"_sgp_k562_lane24.Rda")) 
  sgp_mcf7_lane11 <- readRDS(paste0("data/",current_chr,"_sgp_mcf7_lane11.Rda")) 
  sgp_mcf7_lane23 <- readRDS(paste0("data/",current_chr,"_sgp_mcf7_lane23.Rda"))
  sgp_k562_lane13_inter <- readRDS(paste0("data/",current_chr,"_sgp_k562_lane13_inter.Rda"))
  sgp_k562_lane24_inter <- readRDS(paste0("data/",current_chr,"_sgp_k562_lane24_inter.Rda")) 
  sgp_mcf7_lane11_inter <- readRDS(paste0("data/",current_chr,"_sgp_mcf7_lane11_inter.Rda"))
  sgp_mcf7_lane23_inter <- readRDS(paste0("data/",current_chr,"_sgp_mcf7_lane23_inter.Rda"))
  sfd_k562_rep1 <- readRDS(paste0("data/",current_chr,"_sfd_k562_rep1.Rda"))  
  sfd_k562_rep2 <- readRDS(paste0("data/",current_chr,"_sfd_k562_rep2.Rda"))
  sfd_mcf7_rep3 <- readRDS(paste0("data/",current_chr,"_sfd_mcf7_rep3.Rda"))
  sfd_mcf7_rep4 <- readRDS(paste0("data/",current_chr,"_sfd_mcf7_rep4.Rda"))
  
  ### hic 
  hmec_file_1 <- readRDS(paste0("data/",current_chr,"_hmec_file_1.Rda"))
  hmec_file_2 <- readRDS(paste0("data/",current_chr,"_hmec_file_2.Rda"))
  k562_file_1 <- readRDS(paste0("data/",current_chr,"_k562_file_1.Rda"))
  k562_file_2 <- readRDS(paste0("data/",current_chr,"_k562_file_2.Rda"))
  
  ### lncrna
  lncrna <- readRDS(paste0("data/",current_chr,"_mitrans.Rda"))
  lncrna_expr <- readRDS(paste0("data/lncrna_expr.Rda"))
  
  my.data = list(sgp_k562_lane13 = sgp_k562_lane13,
                 sgp_k562_lane24 = sgp_k562_lane24,
                 sgp_mcf7_lane11 = sgp_mcf7_lane11,
                 sgp_mcf7_lane23 = sgp_mcf7_lane23,
                 sgp_k562_lane13_inter = sgp_k562_lane13_inter,
                 sgp_k562_lane24_inter = sgp_k562_lane24_inter,
                 sgp_mcf7_lane11_inter = sgp_mcf7_lane11_inter,
                 sgp_mcf7_lane23_inter = sgp_mcf7_lane23_inter,
                 sfd_k562_rep1 = sfd_k562_rep1,
                 sfd_k562_rep2 = sfd_k562_rep2,
                 sfd_mcf7_rep3 = sfd_mcf7_rep3,
                 sfd_mcf7_rep4 = sfd_mcf7_rep4,
                 hmec_file_1 = hmec_file_1,
                 hmec_file_2 = hmec_file_1,
                 k562_file_1 = k562_file_1,
                 k562_file_2 = k562_file_2,
                 lncrna = lncrna,
                 lncrna_expr = lncrna_expr)
  
  invisible(my.data)
  
}

#' Check data availability and extract relevant data under the form of list of
#' GenomicRanges
#' After setting the current chromosome and study range
#' 
#' 
#' @return list of 3 lists containing GenomicRanges representing data in the 
#' current study range : arch, circos, lncrna
#'
#' @examples
#' current_chr <- "chr12"
#' my.data <- loadChrData()
#' current_range <- setStudyRange(27950000, 28735000)
#' my.ranges <- getDataOverview()
#' 
#' @export
getDataOverview <- function() {
  can_run_2()
  
  my.ranges = list(
       arch = list(sgp_k562_lane13 = create_gr_from_df(my.data$sgp_k562_lane13,2,6,"K562 ChIA-Pet lane 13"),
                   sgp_k562_lane24 = create_gr_from_df(my.data$sgp_k562_lane24, 2,6,"K562 ChIA-Pet lane 24"),
                   sgp_mcf7_lane11 = create_gr_from_df(my.data$sgp_mcf7_lane11,2,6,"MCF7 ChIA-Pet lane 11"),
                   sgp_mcf7_lane23 = create_gr_from_df(my.data$sgp_mcf7_lane23,2,6,"MCF7 ChIA-Pet lane 23"),
                   sfd_k562_rep1 = create_gr_from_df(my.data$sfd_k562_rep1,2,6,"K562 ChIA-Pet rep 1"),
                   sfd_k562_rep2 = create_gr_from_df(my.data$sfd_k562_rep2,2,6,"K562 ChIA-Pet rep 2"),
                   sfd_mcf7_rep3 = create_gr_from_df(my.data$sfd_mcf7_rep3,2,6,"MCF7 ChIA-Pet rep 3"),
                   sfd_mcf7_rep4 = create_gr_from_df(my.data$sfd_mcf7_rep4,2,6,"MCF7 ChIA-Pet rep 4"),
                   hmec_file_1 = create_gr_from_df(my.data$hmec_file_1,2,3,"HMEC HiC-arrowhead"),
                   hmec_file_2 = create_gr_from_df(my.data$hmec_file_2,2,6,"HMEC HiC-hiccups"),
                   k562_file_1 = create_gr_from_df(my.data$k562_file_1,2,3,"K562 HiC-arrowhead"),
                   k562_file_2 = create_gr_from_df(my.data$k562_file_2,2,6,"K562 HiC-hiccups")),
     circos = list(sgp_k562_lane13_inter = create_gr_from_df_4_circos_range(my.data$sgp_k562_lane13_inter,"K562 ChIA-Pet lane 13 inter-chrom"),
                   sgp_k562_lane24_inter = create_gr_from_df_4_circos_range(my.data$sgp_k562_lane24_inter,"K562 ChIA-Pet lane 24 inter-chrom"),
                   sgp_mcf7_lane11_inter = create_gr_from_df_4_circos_range(my.data$sgp_mcf7_lane11_inter,"MCF7 ChIA-Pet lane 11 inter-chrom"),
                   sgp_mcf7_lane23_inter = create_gr_from_df_4_circos_range(my.data$sgp_mcf7_lane23_inter,"MCF7 ChIA-Pet lane 23 inter-chrom")),
     lncrna = create_gr_from_microTfile(my.data$lncrna,4,5,"lncrna"))
   
  invisible(my.ranges)
}


#' Prepare tracks to plot them
#' 
#' @param ranges_list list of GenomicRanges generated with the functions 
#' \code{getDataOverview}
#' @param highlight_range_list list of GenomicRange generated with the function 
#' \code{setHighLight}
#' 
#' @return list of 3 lists containing GenomicRanges representing data in the 
#' current study range : arch, circos, lncrna
#'
#' @examples
#' current_chr <- "chr12"
#' my.data <- loadChrData()
#' current_range <- setStudyRange(27950000, 28735000)
#' my.ranges <- getDataOverview()
#' specific_range1 <- setHighLight(28111017,28127138,"alpha")
#' specific_range2 <- setHighLight(28284682,28287682,"alpha")
#' specific_range3 <- setHighLight(28284682,28287682,"color")
#' 
#' arch <- drawArchs(my.ranges$arch, list(specific_range1,specific_range2,
#' specific_range3))
#' tracks(arch) + xlim(current_range)
#' 
#' @export
drawArchs <- function(ranges_list, highlight_range_list = NULL) {
  #highlight_range_list : start, stop, highlight method
  
  my.tracks = c()
  
  for(i in seq_along(ranges_list)) {
    range = ranges_list[[i]]
    if(length(range) > 0) {
      for(highlight_range in highlight_range_list) {
        range = make_emphasis(range, highlight_range)      
      }
      track_title = gsub(x = range[1]$color, pattern = " ", replacement = "\n")
      range_track = get_chiapet_arch(range,track_title)
      my.tracks = c(my.tracks, range_track)
    }   
  }
  
  invisible(my.tracks)
}

#' Enable to set the current study range
#' 
#' @param current_start integer start of the range
#' @param current_stop integer end of the range
#' 
#' @return a GenomicRange
#'
#' @examples
#' current_range <- setStudyRange(27950000, 28735000)
#' 
#' @export
setStudyRange <- function(current_start, current_stop) {
  can_run()
  current_range <- IRanges::IRanges(current_start, current_stop)
  current_range <- GenomicRanges::GRanges(seqnames = current_chr, ranges = current_range)
  current_range
}


#' Enable to set the zone to higlight
#' 
#' @param current_start integer start of the range
#' @param current_stop integer end of the range
#' @param method character defining the way to highlight a specific zone (alpha 
#' or color)
#' 
#' @return a GenomicRange
#'
#' @examples
#' specific_range <- setHighLight(28111017,28127138,"alpha")
#' 
#' @export
setHighLight <- function(current_start, current_stop, method) {
  can_run_3(current_start, current_stop, method)
  
  # method : alpha, color, size
  special_range <- IRanges(current_start, current_stop)
  special_range <- GRanges(seqnames = current_chr, ranges = special_range, 
                           method = method)
  special_range
}


#' Extract annotation informations from the annotation package org.Hs.eg.db
#' 
#' @param file_path character absolute path to the file containing snp 
#' informations
#' @param label character label of the track, default "Annotations" 
#' @return a ggplot track
#'
#' @examples
#' current_chr <- "chr12"
#' current_range <- setStudyRange(27950000, 28735000)
#' annot_track <- drawAnnotations("My Genes")
#' 
#' @export
drawAnnotations <- function(label = "Annotations") {
  can_run_4()
  
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
  
  
  genes <- autoplot(grl_txdb, aes(type = model)) + 
    theme_bw() + xlim(current_range) + ylab(label) +
    theme(axis.title.y = element_text(size = rel(0.5), angle = 90))
  
  genes
}

#' Extract SNP informations from a flat file and build the resulting track
#' 
#' @param file_path character absolute path to the file containing snp 
#' informations
#' @param label character label of the track, default "SNPs"
#' 
#' @return a ggplot track
#'
#' @examples
#' snps_list = "~/Workspace/COLLAB/PTHLH_snps"
#' current_chr <- "chr12"
#' current_range <- setStudyRange(27950000, 28735000)
#' snps_track <- drawSNPs(snps_list, "My SNPs")
#' 
#' @export
drawSNPs <- function(file_path, label = "SNPs") {
  can_run_4()
  
  snps_df <- read.table(snps_list, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")
  snp_ids <- snps_df$snp_id
  snps_ranges <- IRanges(snps_df$location, snps_df$location)
  snps <- GRanges(seqnames = current_chr, ranges = snps_ranges, imp = snps_df$metadata)
  snps$name <- snp_ids
  
  snps_track <- autoplot(snps, aes(color=imp)) +
    geom_text(aes(x = start, y = 1, label=name, angle = 90, vjust=-1), size = 1, color = "blue") +
    theme_bw() + ylab(label) + xlim(current_range) + guides(color= FALSE) + 
    theme(axis.title.y = element_text(size = rel(0.5), angle = 90))
  
  snps_track
}




#' Merge range in a single one
#' 
#' @param ranges vector of GenomicRanges generated with the functions 
#' \code{getDataOverview} you want to merge
#' @param label character label of the merged range
#' 
#' @return a GenomicRange 
#'
#' @examples
#' current_chr <- "chr12"
#' my.data <- loadChrData()
#' current_range <- setStudyRange(27950000, 28735000)
#' my.ranges <- getDataOverview()
#' k562_chiapet <- mergeRanges(c(my.ranges$arch$sgp_k562_lane13,
#' my.ranges$arch$sgp_k562_lane24), label = "K562 ChiA-PET SGP")
#' 
#' @export
mergeRanges <- function(ranges, label) {
  ranges$color = rep(x = label, times = length(ranges))
  ranges
}



can_run <- function() {
  if(!exists("current_chr"))
  {
    stop(call. = FALSE, "Please define the chromosome to study before running 
         this function(ex : current_chr = \"chr12\")")
  }
}

can_run_2 <- function() {
  if(!exists("current_range"))
  {
    stop(call. = FALSE, "Please define the range of the chromosome to study 
         (ex : current_range = setStudyRange(27950000, 28735000))")
  }
  
  if(!exists("my.data"))
  {
    stop(call. = FALSE, "Please load data before running this function
         ex : my.data = loadChrData()")
  }
}

can_run_3 <- function(current_start, current_stop, method) {
  if(current_start < start(current_range) || current_stop > end(current_range))
  {
    stop(call. = FALSE, "The highlight zone must be included in the study range")
  }
  
  if(!(method %in% c("color","alpha")))
  {
    stop(call. = FALSE, "The highlight method must be \'color\' or \'alpha\'")
  }
}

can_run_4 <- function() {
  if(!exists("current_range"))
  {
    stop(call. = FALSE, "Please define the range of the chromosome to study 
         (ex : current_range = setStudyRange(27950000, 28735000))")
  }
}
