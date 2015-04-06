library(parallel)

chr_list = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
             "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
             "chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrM")

# chr_list = c("chr1","chr2","chr3")

############################ SINGAPORE DATA ############################  
chiapet_k562_lane13_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_lane13_"
chiapet_k562_lane24_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_lane24_"
chiapet_mcf7_lane11_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_lane11_"
chiapet_mcf7_lane23_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_lane23_"
############################ STANDFORD DATA ############################  
chiapet_k562_rep1_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_rep1_"
chiapet_k562_rep2_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_K562_rep2_"
chiapet_mcf7_rep3_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_rep3_"
chiapet_mcf7_rep4_path <- "/home/nekomimi/Workspace/SNPVIZU/data/chiapet/parsed/CHIAPET_MCF7_rep4_"

sgp_files_path = c(chiapet_k562_lane13_path, chiapet_k562_lane24_path,
                   chiapet_mcf7_lane11_path, chiapet_mcf7_lane23_path)

sdf_file_path = c(chiapet_k562_rep1_path, chiapet_k562_rep2_path,
                  chiapet_mcf7_rep3_path, chiapet_mcf7_rep4_path)

names(sgp_files_path) = c("sgp_k562_lane13", "sgp_k562_lane24", 
                          "sgp_mcf7_lane11", "sgp_mcf7_lane23")

names(sdf_file_path) = c("sfd_k562_rep1", "sfd_k562_rep2",
                         "sfd_mcf7_rep3", "sfd_mcf7_rep4")

for (current_chr in chr_list) {
  create_df_from_sgp <- function(file_name) {
    my.file <- paste0(sgp_files_path[[file_name]], current_chr)
    df <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    saveRDS(df, file=paste0("~/Workspace/SNPVIZU/SNPVizualTools/data/",current_chr,"_",file_name,".Rda"))
  }
  
  create_df_from_sfd <- function(file_name) {
    my.file <- paste0(sdf_file_path[[file_name]], current_chr)
    df <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    data <- strsplit(gsub("\\.\\.", "-", df[,4]), '[:,-]')
    df =  matrix(unlist(data), ncol = 7, byrow = TRUE)
    saveRDS(df, file=paste0("~/Workspace/SNPVIZU/SNPVizualTools/data/",current_chr,"_",file_name,".Rda"))
  }
  
  create_df_from_sgp_inter <- function(file_name) {
    my.file <- paste0(sgp_files_path[[file_name]], current_chr)
    df <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE)
    data <- strsplit(gsub("\\.\\.", "-", df[,4]), '[:,-]')
    df =  matrix(unlist(data), ncol = 7, byrow = TRUE)
    saveRDS(df, file=paste0("~/Workspace/SNPVIZU/SNPVizualTools/data/",current_chr,"_",file_name,"_inter.Rda"))
  }
  
  mclapply(names(sgp_files_path), create_df_from_sgp, mc.cores = 2)
  mclapply(names(sdf_file_path), create_df_from_sfd, mc.cores = 2)
  mclapply(names(sgp_files_path), create_df_from_sgp_inter, mc.cores = 2)
} 



