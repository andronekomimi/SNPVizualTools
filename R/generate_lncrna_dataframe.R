library(parallel)

chr_list = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
             "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
             "chr18","chr19","chr20","chr21","chr22")

# chr_list = c("chr1","chr2","chr3")



files_path <- "/home/nekomimi/Workspace/COLLAB/mitranscriptome.gtf/mitranscriptome_"

for (current_chr in chr_list) {
  my.file <- paste0(files_path, current_chr)
  df <- read.table(my.file, header=FALSE, stringsAsFactors=FALSE, quote = "\"", sep="\t")
  saveRDS(df, file=paste0("~/Workspace/SNPVIZU/SNPVizualTools/data/",current_chr,"_mitrans.Rda")) 
}

lncrna_expr <- "/home/nekomimi/Workspace/COLLAB/mitranscriptome.expr.fpkm_select.tsv"
df <- read.table(lncrna_expr, header=TRUE, stringsAsFactors=FALSE, quote = "\"", sep="\t")
saveRDS(df, file=paste0("~/Workspace/SNPVIZU/SNPVizualTools/data/lncrna_expr.Rda"))



