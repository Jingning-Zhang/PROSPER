
suppressMessages(library("optparse"))

option_list = list(
  make_option("--packagedir", action="store", default=NA, type='character',
              help=" [required]"),
  make_option("--refgeno", action="store", default=NA, type='character',
              help=" [required]"),
  make_option("--hg", action="store", default=NA, type='character',
              help=" [required]"),
  make_option("--chr", action="store", default=NA, type='character',
              help=" [required]"),
  make_option("--workdir", action="store", default=NA, type='character',
              help=" [required]")
)
opt = parse_args(OptionParser(option_list=option_list))

#opt <- list()
#opt$workdir <- "/dcs04/nilanjan/data/23andme/Analysis/JZ/multi-lassosum"
#opt$packagedir <- "/dcs04/nilanjan/data/jzhang2/MEPRS/PRSdir"
#opt$hg <- 19
#chr <- "22"
#opt$refgeno <- "/dcs04/nilanjan/data/23andme/Analysis/JZ/ref_geno/AFR/any_cvd/ref_chr"


####################################################################################################
####################################################################################################

## reformat to starndard data

library(caret)

dir.create(paste0(opt$workdir,"/standard_data"))

block_info <- read.table(paste0(opt$packagedir,"/Berisa.EUR.hg",opt$hg,".bed"),header = T)

chr <- paste0("chr",opt$chr)

block_info_tmp <- block_info[block_info$chr==chr,]

Nsnps <- integer(length = nrow(block_info_tmp)+2)
snps_list <- vector("list", length = nrow(block_info_tmp)+2)
LD_list <- vector("list", length = nrow(block_info_tmp)+2)

#### start

snps <- character()
tmp.snps <- try(read.table(paste0(opt$workdir, "/tmp/byblock/", chr,"/",
                                  chr,"_start_",block_info_tmp$start[1],".bim"), stringsAsFactors = F), silent=TRUE)
if ('try-error' %in% class(tmp.snps)) {
  Nsnps[1] <- 0
}else{
  n.snp.tmp <- nrow(tmp.snps)
  tmp.LD <- readBin(paste0(opt$workdir, "/tmp/LD/", chr,"/",
                            chr,"_start_", block_info_tmp$start[1],".ld.bin"),
                    what="numeric", size=4, n=(n.snp.tmp)^2)
  #tmp.LD[is.nan(tmp.LD)] <- 1; tmp.LD <- matrix(tmp.LD, ncol=n.snp.tmp)
  #if(n.snp.tmp>1){
  #  drop = findCorrelation(tmp.LD,cutoff = 0.999999999)
  #  Nsnps[1] <- n.snp.tmp - length(drop)
  #  snps_list[[1]] <- tmp.snps$V2[-drop]
  #  LD_list[[1]] <- tmp.LD[-drop, -drop, drop=F]
  #  cat(paste0("Total #of SNP is ",n.snp.tmp,", and ",length(drop)," are removed due to homozygosity or perfect correlation.\n"))
  #}else{
  #  Nsnps[1] <- n.snp.tmp
  #  snps_list[[1]] <- tmp.snps$V2
  #  LD_list[[1]] <- tmp.LD
  #  cat(paste0("Total #of SNP is ",n.snp.tmp,", and 0 are removed due to homozygosity or perfect correlation.\n"))
  #}
  tmp.LD <- matrix(tmp.LD, ncol=n.snp.tmp)
  Nsnps[1] <- n.snp.tmp
  snps_list[[1]] <- tmp.snps$V2
  LD_list[[1]] <- tmp.LD
  cat(paste0("Total #of SNP is ",n.snp.tmp,".\n"))
}


#### Middle

for (i in 1:nrow(block_info_tmp)){
  snps <- character()
  tmp.snps <- try(read.table(paste0(opt$workdir, "/tmp/byblock/", chr,"/",
                                    chr,"_",block_info_tmp$start[i],"_",block_info_tmp$stop[i],".bim"), stringsAsFactors = F), silent=TRUE)
  if ('try-error' %in% class(tmp.snps)) {
    Nsnps[i+1] <- 0
  }else{
    n.snp.tmp <- nrow(tmp.snps)
    tmp.LD <- readBin(paste0(opt$workdir, "/tmp/LD/", chr,"/",
                              chr,"_", block_info_tmp$start[i],"_",block_info_tmp$stop[i],".ld.bin"),
                      what="numeric", size=4, n=(n.snp.tmp)^2)
    #tmp.LD[is.nan(tmp.LD)] <- 1; tmp.LD <- matrix(tmp.LD, ncol=n.snp.tmp)
    #if(n.snp.tmp>1){
    #  drop = findCorrelation(tmp.LD,cutoff = 0.999999999)
    #  Nsnps[i+1] <- n.snp.tmp - length(drop)
    #  snps_list[[i+1]] <- tmp.snps$V2[-drop]
    #  LD_list[[i+1]] <- tmp.LD[-drop, -drop, drop=F]
    #  cat(paste0("Total #of SNP is ",n.snp.tmp,", and ",length(drop)," are removed due to homozygosity or perfect correlation.\n"))
    #}else{
    #  Nsnps[i+1] <- n.snp.tmp
    #  snps_list[[i+1]] <- tmp.snps$V2
    #  LD_list[[i+1]] <- tmp.LD
    #  cat(paste0("Total #of SNP is ",n.snp.tmp,", and 0 are removed due to homozygosity or perfect correlation.\n"))
    #}
    tmp.LD <- matrix(tmp.LD, ncol=n.snp.tmp)
    Nsnps[i+1] <- n.snp.tmp
    snps_list[[i+1]] <- tmp.snps$V2
    LD_list[[i+1]] <- tmp.LD
    cat(paste0("Total #of SNP is ",n.snp.tmp,".\n"))

  }
  cat(paste0(chr,": ",i,"/",nrow(block_info_tmp),"\n"))
}
#warnings()

#### end

snps <- character()
tmp.snps <- try(read.table(paste0(opt$workdir, "/tmp/byblock/", chr,"/",
                                  chr,"_",block_info_tmp$stop[i],"_end.bim"), stringsAsFactors = F), silent=TRUE)
if ('try-error' %in% class(tmp.snps)) {
  Nsnps[i+2] <- 0
}else{
  n.snp.tmp <- nrow(tmp.snps)
  tmp.LD <- readBin(paste0(opt$workdir, "/tmp/LD/", chr,"/",
                            chr,"_",block_info_tmp$stop[i],"_end.ld.bin"),
                    what="numeric", size=4, n=(n.snp.tmp)^2)
  #tmp.LD[is.nan(tmp.LD)] <- 1; tmp.LD <- matrix(tmp.LD, ncol=n.snp.tmp)
  #if(n.snp.tmp>1){
  #  drop = findCorrelation(tmp.LD,cutoff = 0.999999999)
  #  Nsnps[i+2] <- n.snp.tmp - length(drop)
  #  snps_list[[i+2]] <- tmp.snps$V2[-drop]
  #  LD_list[[i+2]] <- tmp.LD[-drop, -drop, drop=F]
  #  cat(paste0("Total #of SNP is ",n.snp.tmp,", and ",length(drop)," are removed due to homozygosity or perfect correlation.\n"))
  #}else{
  #  Nsnps[i+2] <- n.snp.tmp
  #  snps_list[[i+2]] <- tmp.snps$V2
  #  LD_list[[i+2]] <- tmp.LD
  #  cat(paste0("Total #of SNP is ",n.snp.tmp,", and 0 are removed due to homozygosity or perfect correlation.\n"))
  #}
  tmp.LD <- matrix(tmp.LD, ncol=n.snp.tmp)
  Nsnps[i+2] <- n.snp.tmp
  snps_list[[i+2]] <- tmp.snps$V2
  LD_list[[i+2]] <- tmp.LD
  cat(paste0("Total #of SNP is ",n.snp.tmp,".\n"))

}

cat(paste0("Saving standard data for ",chr ,"...\n"))
save(Nsnps, snps_list,
     file = paste0(opt$workdir,"/standard_data/",chr,"_snps.RData"))

save(Nsnps, snps_list, LD_list,
     file = paste0(opt$workdir,"/standard_data/",chr,"_LD.RData"))

cat(paste0(chr, " completed.\n"))


