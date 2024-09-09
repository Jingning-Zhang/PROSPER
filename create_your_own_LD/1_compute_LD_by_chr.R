
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
#opt$refgeno <- "/dcs04/nilanjan/data/23andme/Analysis/JZ/ref_geno/AFR/any_cvd/ref_chr"


####################################################################################################
####################################################################################################

## compute LD


dir.create(paste0(opt$workdir, "/tmp"))
dir.create(paste0(opt$workdir, "/tmp/byblock"))
dir.create(paste0(opt$workdir, "/tmp/LD"))

chr <- opt$chr


block_info <- read.table(paste0(opt$packagedir,"/Berisa.EUR.hg",opt$hg,".bed"),header = T)

block_info_tmp <- block_info[block_info$chr==paste0("chr",chr),]

dir.create(paste0(opt$workdir, "/tmp/byblock/chr",chr))
dir.create(paste0(opt$workdir, "/tmp/LD/chr",chr))



a <- paste0(opt$packagedir,"/plink2",
            " --threads 1",
            " --bfile ",opt$refgeno, chr,
            " --chr ",chr," --to-bp ",block_info_tmp$start[1]-1,
            " --make-bed",
            " --out ",opt$workdir, "/tmp/byblock/chr",chr,"/chr",chr,"_start_",block_info_tmp$start[1])

system(a)

a <- paste0(opt$packagedir,"/plink1",
            " --keep-allele-order",
            " --threads 1",
            " --bfile ",opt$workdir, "/tmp/byblock/chr",chr,"/chr",chr,"_start_",block_info_tmp$start[1],
            " --r bin4",
            " --out ",opt$workdir, "/tmp/LD/chr",chr,"/chr",chr,"_start_",block_info_tmp$start[1])
  system(a)


for (i in 1:nrow(block_info_tmp)){

  a <- paste0(opt$packagedir,"/plink2",
            " --threads 1",
            " --bfile ",opt$refgeno, chr,
            " --chr ",chr," --from-bp ",block_info_tmp$start[i]," --to-bp ",block_info_tmp$stop[i]-1,
            " --make-bed",
            " --out ",opt$workdir, "/tmp/byblock/chr",chr,"/chr",chr,"_",block_info_tmp$start[i],"_",block_info_tmp$stop[i])

  system(a)

  a <- paste0(opt$packagedir,"/plink1",
            " --keep-allele-order",
            " --threads 1",
            " --bfile ",opt$workdir, "/tmp/byblock/chr",chr,"/chr",chr,"_",block_info_tmp$start[i],"_",block_info_tmp$stop[i],
            " --r bin4",
            " --out ",opt$workdir, "/tmp/LD/chr",chr,"/chr",chr,"_",block_info_tmp$start[i],"_",block_info_tmp$stop[i])
  system(a)

}


a <- paste0(opt$packagedir,"/plink2",
            " --threads 1",
            " --bfile ",opt$refgeno, chr,
            " --chr ",chr," --from-bp ",block_info_tmp$stop[i],
            " --make-bed",
            " --out ",opt$workdir, "/tmp/byblock/chr",chr,"/chr",chr,"_",block_info_tmp$stop[i],"_end")

system(a)

a <- paste0(opt$packagedir,"/plink1",
            " --keep-allele-order",
            " --threads 1",
            " --bfile ",opt$workdir, "/tmp/byblock/chr",chr,"/chr",chr,"_",block_info_tmp$stop[i],"_end",
            " --r bin4",
            " --out ",opt$workdir, "/tmp/LD/chr",chr,"/chr",chr,"_",block_info_tmp$stop[i],"_end")
  system(a)






