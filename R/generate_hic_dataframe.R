library(parallel)

chr_list = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
             "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
             "chr18","chr19","chr20","chr21","chr22")

# chr_list = c("chr1","chr2","chr3")

hic_hmec_file_1 <- "/home/nekomimi/Workspace/SNPVIZU/data/HIC_HMEC_Arrowhead_domainlist_"
hic_hmec_file_2 <- "/home/nekomimi/Workspace/SNPVIZU/data/HIC_HMEC_HiCCUPS_looplist_"

hic_k562_file_1 <- "/home/nekomimi/Workspace/SNPVIZU/data/HIC_K562_Arrowhead_domainlist_"
hic_k562_file_2 <- "/home/nekomimi/Workspace/SNPVIZU/data/HIC_K562_HiCCUPS_looplist_"


files_path = c(hic_hmec_file_1, hic_hmec_file_2,hic_k562_file_1, hic_k562_file_2)

names(files_path) = c("hmec_file_1", "hmec_file_2", "k562_file_1", "k562_file_2")

for (current_chr in chr_list) {
  create_df <- function(file_name) {
    my.file <- paste0(files_path[[file_name]], current_chr)
    df <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    saveRDS(df, file=paste0("~/Workspace/SNPVIZU/SNPVizualTools/data/",current_chr,"_",file_name,".Rda"))
  }
  
  mclapply(names(files_path), create_df, mc.cores = 2)
  
}