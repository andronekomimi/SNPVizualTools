#' Download data for a current_chrosome
#' 
#' @param current_chr current_chrosome name (ex : chr1)
#' 
#' @return void
#'
#' @import GenomicRanges
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

drawArch <- function(ranges_list, highlight_range_list = NULL) {
  #highlight_range_list : start, stop, highlight method
  
  my.tracks = c()
  
  for(range in ranges_list) {
    if(length(range) > 0) {
      for(highlight_range in highlight_range_list) {
        range = make_emphasis(range, highlight_range)      
      }
      range_track = get_chiapet_arch(range)
      my.tracks = c(my.tracks, range_track)
    }   
  }
  
  invisible(my.tracks)
}

setStudyRange <- function(current_start, current_stop) {
  can_run()
  current_range <- IRanges::IRanges(current_start, current_stop)
  current_range <- GenomicRanges::GRanges(seqnames = current_chr, ranges = current_range)
  current_range
}

setHighLight <- function(current_start, current_stop, method) {
  # method : alpha, color, size
  special_range <- IRanges(current_start, current_stop)
  special_range <- GRanges(seqnames = current_chr, ranges = special_range, 
                           alpha = 0.6, color = "promoter")
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
