# Create your own LD reference panel



1. Please make sure the bfiles include the same set of variants and the same ref alleles. Even if a variant is not common in a certain population, please also include it in the bfiles to make sure the LD matrix can be matched together. And remove variants with duplicate IDs.

These could be done using plink commands

```
--ref-allele
--rm-dup exclude-mismatch
```

Put those bfiles (by chromosomes 1-22) inside the directory ${dir_for_refgeno_bfile}/${ancestry}/, and named them with ref_chr${chr}.bed, ref_chr${chr}.bim, ref_chr${chr}.fam.

2. Please put plink1 (executable file), plink2 (executable file), Berisa.EUR.hg19.bed and Berisa.EUR.hg38.bed (these two can be downloaded in this repo) inside the directory ${dir_for_plink12_and_ldblock}.

3. Please put the two scripts (1_compute_LD_by_chr.R and 2_reformat_to_RData_by_chr.R) inside the directory ${dir_for_scripts}.

4. Please create two repo: ${dir_for_LD_results_temp} and ${dir_for_LD_results}. The first is a temporary repo to generate temporary LD files, and the second one is the final repo that will be used to store your own LD reference panel.

5. Codes used to generate temporary LD files inside ${dir_for_LD_results_temp}:

```
cd ${dir_for_scripts}

mkdir ${dir_for_LD_results_temp}/
mkdir ${dir_for_LD_results_temp}/${ancestry}/

Rscript ./1_compute_LD_by_chr.R \
--packagedir ${dir_for_plink12_and_ldblock} \
--refgeno ${dir_for_refgeno_bfile}/${ancestry}/ref_chr \
--hg ${genome_buil_of_bfile_19_or_38} \
--chr ${chr} \
--workdir ${dir_for_LD_results_temp}/${ancestry}/

Rscript ./2_reformat_to_RData_by_chr.R \
--packagedir ${dir_for_plink12_and_ldblock} \
--refgeno ${dir_for_refgeno_bfile}/${ancestry}/ref_chr \
--hg ${genome_buil_of_bfile_19_or_38} \
--chr ${chr} \
--workdir ${dir_for_LD_results_temp}/${ancestry}/

```

6. After finishing running it for all chromosomes and ancestries, please copy them to the final repo ${dir_for_LD_results}.

```
mkdir ${dir_for_LD_results}/${ancestry}
cd ${dir_for_LD_results}/${ancestry}
for (( i = 1; i < 23; i++ )); do
	cp ${dir_for_LD_results_temp}/${ancestry}/standard_data/chr${chr}_LD.RData ./
done
```

7. Run the R codes to generate ref_bim.txt（a reference text file for all variants and alleles）.

```
library(bigreadr)

dir_for_LD_results = "PUT YOUR dir_for_LD_results HERE."
dir_for_LD_results_temp = "PUT YOUR dir_for_LD_results_temp HERE."
dir_for_refgeno_bfile = "PUT YOUR dir_for_refgeno_bfile HERE."

for(ancestry in c("AFR","AMR","EAS","EUR","SAS")){

  print(ancestry)

  for(chr in 1:22){
    #print(chr)

    load(paste0(dir_for_LD_results_temp,"/",ancestry,"/standard_data/chr",chr,"_snps.RData"))
    a <- unlist(snps_list)
    b <- fread2(paste0(dir_for_refgeno_bfile,"/",ancestry,"/ref_chr",chr,".bim"))
    b <- b[match(a, b$V2),]
    if(chr==1){
      res <- b
    }else{
      res <- rbind(res, b)
    }
    print(nrow(res))
  }
  fwrite2(res, paste0(dir_for_LD_results,"/",ancestry,"/ref_bim.txt"), col.names = F, sep="\t", nThread=1)

}

res <- NULL
for(ancestry in c("AFR","AMR","EAS","EUR","SAS")){
  dat <- fread2(paste0(dir_for_LD_results,"/",ancestry,"/ref_bim.txt"))
  res <- rbind(res,dat)
  res <- unique(res)
}

## Used to check if there is any duplicate variants.
## If max(a) is bigger than 1, meaning there are duplicate variants, please check your bfile and match the alleles.
a = table(res$V2); max(a)

fwrite2(res, paste0(dir_for_LD_results,"/ref_bim.txt"), col.names = F, sep="\t", nThread=1)

```

8. The LD panel in ${dir_for_LD_results} is now ready to use for PROSPER.



