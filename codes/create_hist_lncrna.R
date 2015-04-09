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

current_chr = "chr12"
current_start = 27950000
current_stop = 28735000
current_target = "PTHLH"
snps_list = "~/Workspace/COLLAB/PTHLH_snps"



######################### LOAD LIBRARY ###############################
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(ggbio))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressMessages(library(biovizBase))
suppressMessages(library(org.Hs.eg.db))

######################### FUNCTIONS ###############################

get_TS_from_df = function(DF,start_idx,end_idx) {
  
  transcript_id = c()
  
  for (i in (1:nrow(DF)))
  {
    start_value = as.numeric(DF[i,start_idx])
    end_value = as.numeric(DF[i,end_idx])
    
    if((start_value > (start(current_range) - 1) 
        && start_value < (end(current_range) - 1)) 
       || (end_value > (start(current_range) - 1) 
           && end_value < (end(current_range) - 1))) {
      
      str = strsplit(df13[i,9], split = "; ")[[1]]
      transcript_id_idx = grep(str, pattern = "transcript_id")
      tId = gsub(pattern = "transcript_id ",x = str[transcript_id_idx], replacement = "")
      transcript_id = c(transcript_id, tId)
    }
  }
  
  transcript_id
}

create_gr_from_microTfile = function(DF,start_idx,end_idx,label) {
  label = ""
  left = c()
  right = c()
  base_width = c()
  color = c()
  alpha = c()
  transcript_id = c()
  annotated = c()
  
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
      
      # transcript_id
      str = strsplit(DF[i,9], split = "; ")[[1]]
      transcript_id_idx = grep(str, pattern = "transcript_id")
      tId = gsub(pattern = "transcript_id ",x = str[transcript_id_idx], replacement = "")
      transcript_id = c(transcript_id, tId)
      
      # transcript annotation
      annotation_idx = grep(str, pattern = "tstatus")
      annotation = ! grepl(x = str[annotation_idx], pattern = "unannotated")
      annotated = c(annotated, annotation)
    }
  }
  
  ranges = IRanges(left,right)
  
  if(length(left) > 0) {
    g_ranges <- GRanges(seqnames = current_chr, ranges = ranges, 
                        base_width = base_width, color = color, alpha = rep(0.5,length(left)),
                        annotated = annotated, tsID = transcript_id)
    print(paste0("----Returning a ",length(g_ranges),"-long Gr"))
    invisible(g_ranges)
  }
  else
  {
    invisible(GRanges())
  }
}


######################### PROCESSING ###############################
##################### FILE LOADING ##################################

lncrna_file1 <- paste0("/home/nekomimi/Workspace/COLLAB/mitranscriptome.gtf/mitranscriptome_", current_chr)
df13 <- read.table(lncrna_file1, header=FALSE, stringsAsFactors=FALSE, quote = "\"", sep="\t")

lncrna_file2 <- "/home/nekomimi/Workspace/COLLAB/mitranscriptome.expr.fpkm_select.tsv"
df14 <- read.table(lncrna_file2, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")

snps_df <- read.table(snps_list, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")


##################### SNPS ###############################################

snp_ids <- snps_df$snp_id
snps_ranges <- IRanges(snps_df$location, snps_df$location)
snps <- GRanges(seqnames = current_chr, ranges = snps_ranges, imp = snps_df$metadata)
snps$name <- snp_ids

##################### CURRENT RANGE ##################################

current_range <- IRanges(current_start,current_stop)
current_range <- GRanges(seqnames = current_chr, ranges = current_range)

######################### LNCRNA TRACK ###############################

lnrna <- create_gr_from_microTfile(df13,4,5,"lncrna")

ts = unique(lnrna$tsID)
ts_infos = unique(subset(df14, df14$transcript_id %in% ts))

if(nrow(ts_infos) > 0) {
  #lncrna_df =  data.frame()
  
  
  left = c()
  right = c()
  id = c()
  
  # ranges
  for(i in seq(1:nrow(ts_infos))) {
    id = c(id, ts_infos$transcript_id[i])
    r = reduce(lnrna[lnrna$tsID == ts_infos$transcript_id[i]])
    left = c(left, start(r))
    right = c(right, end(r))
    
  }
  
  # HMEC data
  hmec_data = ts_infos$encode_breast_hmec_cshl_rep1
  
  ranges = IRanges(left,right)
  hmec_lncrna = GRanges(seqnames = current_chr, ranges = ranges, tsID = id, exp = hmec_data)
  hmec_lncrna_track =  ggplot(data = hmec_lncrna) + 
    geom_segment(mapping=aes(x=start, xend=end, color=exp), size = 10) + ylab("") +
    geom_text(aes(x = start, y = 1, label=tsID, vjust=3), size = 2, color = "blue") + 
    theme_bw() + xlim(current_range) + guides(color= FALSE)
  
  # K562 data
  exp = c()
  for(i in seq(1:nrow(ts_infos))) {
    data = mean(c(ts_infos$encode_cml_k562_cshl_rep1[i],
                  ts_infos$encode_cml_k562_cshl_rep2[i], 
                  ts_infos$encode_cml_k562_caltech_rep1[i],
                  ts_infos$encode_cml_k562_caltech_rep2[i]))
    exp = c(exp, data)
  }
  
  ranges = IRanges(left,right)
  k562_lncrna = GRanges(seqnames = current_chr, ranges = ranges, tsID = id, exp = exp)
  k562_lncrna_track =  ggplot(data = k562_lncrna) + 
    geom_segment(mapping=aes(x=start, xend=end, color=exp), size = 10) + ylab("") + 
    geom_text(aes(x = start, y = 1, label=tsID, vjust=3), size = 2, color = "blue") + 
    theme_bw() + xlim(current_range) + guides(color= FALSE)
  
  # MCF7 data
  exp = c()
  for(i in seq(1:nrow(ts_infos))) {
    data = mean(c(ts_infos$encode_breast_mcf7_caltech_rep1[i],
                  ts_infos$encode_breast_mcf7_caltech_rep2[i],
                  ts_infos$encode_breast_mcf7_caltech_rep3[i], 
                  ts_infos$encode_breast_mcf7_cshl_rep1[i],
                  ts_infos$encode_breast_mcf7_cshl_rep2[i]))
    exp = c(exp, data)
  }
  
  ranges = IRanges(left,right)
  mcf7_lncrna = GRanges(seqnames = current_chr, ranges = ranges, tsID = id, exp = exp)
  mcf7_lncrna_track =  ggplot(data = mcf7_lncrna) + 
    geom_segment(mapping=aes(x=start, xend=end, color=exp), size = 10) + ylab("") +
    geom_text(aes(x = start, y = 1, label=tsID, vjust=3), size = 2, color = "blue") + 
    theme_bw() + xlim(current_range) + guides(color= FALSE)
  
  
  #   lnrna_ranges = c(hmec_lncrna,k562_lncrna, mcf7_lncrna)
  #   lncrna_track = ggplot() + 
  #     geom_segment(data=lnrna_ranges[1:5], mapping=aes(x=start, xend=end, color=exp), size=5) + theme_bw() + ylab("") + xlim(current_range)
  #   
  
  print("Track segment_track -> DONE")  
} else { print("No segment_track")}

###################################################################
mcf7_data = c()
k562_data = c()
hmec_data = ts_infos$encode_breast_hmec_cshl_rep1
for(i in seq(1:nrow(ts_infos))) {
  
  
  mcf7_data = c(mcf7_data, mean(c(ts_infos$encode_breast_mcf7_caltech_rep1[i],
                ts_infos$encode_breast_mcf7_caltech_rep2[i],
                ts_infos$encode_breast_mcf7_caltech_rep3[i], 
                ts_infos$encode_breast_mcf7_cshl_rep1[i],
                ts_infos$encode_breast_mcf7_cshl_rep2[i])))
  
  k562_data = c(k562_data, mean(c(ts_infos$encode_cml_k562_cshl_rep1[i],
                ts_infos$encode_cml_k562_cshl_rep2[i], 
                ts_infos$encode_cml_k562_caltech_rep1[i],
                ts_infos$encode_cml_k562_caltech_rep2[i])))
  
}

lnrna_df = data.frame(ts =  rep(times = 3, x = ts_infos$transcript_id),
                      fpkm = c(k562_data, hmec_data, mcf7_data),
                      cell = c(rep(times = 5, x = c("K562")), rep(times = 5, x = c("HMEC")), rep(times = 5, x = c("MCF7"))))


lnrna_track = ggplot(lnrna_df, aes(x = cell, y = fpkm)) + geom_bar(stat = "identity") + facet_grid(. ~ ts)

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

#print(tracks( "K562 ChIAPET" = chiapet_k562_track,"MCF7 ChIAPET" = chiapet_mcf7_track, "SNPS"=snps_track, "Annotation" = genes , label.text.cex = 0.7))


