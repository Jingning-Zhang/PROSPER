# PROSPER

We proposed **PROSPER**, a new multi-ancestry PRS method with penalized regression followed by ensemble learning. This software is a command line tool based on R programming language. Large-scale benchmarking study shows that PROSPER could be the leading method to reduce the disparity of PRS performance across ancestry groups. 


## Anouncement 

Some user have reported failures when using multiple threads in the foreach loop (lines 186 and 361 of PROSPER.R; and lines 206 and 312 of lassosum2.R). If you’re experiencing this issue, switch to a regular loop, or run each chromosome separately. 

## Recent Version History

May 7, 2025: Fix some bugs. Changes include: PROSPER.R line 312-320.     
Sept 9, 2024: Added a repository titled 'create_your_own_LD' that contains code for generating a user-specific LD reference panel.    
Sept 8, 2024: Fix some bugs. Changes include: tuning_testing.R line 119 and line 211. PROSPER.R line 296.    
Jan 14, 2024: Upload AFR, EAS, and SAS LD reference constructed by UKB.    
Sept 18, 2023: Update citation biorxiv version.    
Sept 11, 2023: Upload EUR LD reference constructed by UKB.    
Apr 5, 2023: Update example codes.    
Nov 24, 2022: Repository made public.


<br/>

## Getting Started

- Download or clone this GitHub repository by `git clone https://github.com/Jingning-Zhang/PROSPER.git`. 

- We'll refer to the corresponding directory as `${package}` throughout the following steps. The directory will be used to store PROSPER package and its LD refernce panels. After you clone using the command above, there should be a folder named `scripts` inside the directory of `${package}`. 

- Note that the files required for downloaded below are all stored on Google Drive. You can download them either by clicking on the provided link or by running the following command `gdown <file_id>`, using the file_id provided below. `gdown` is much faster for downloading to high performance clusters. To use `gdown`, please check [here](https://github.com/wkentaro/gdown "here"). You will need to add it to PATH for use.

- Download the `ref_bim.txt` and saved in `${package}` from [ref_bim.txt]( https://drive.google.com/file/d/1PtD4qk7EBPxdhkGrKrG8OktKxaSDF9AT/view?usp=sharing "ref_bim.txt"); Google Drive <file_id>: `1PtD4qk7EBPxdhkGrKrG8OktKxaSDF9AT`.
    
- Download the LD reference panels and decompress files in `${package}`:

     **Note that if you only want to try the toy example below, you will only need the first two: EUR and AFR reference.**
   
     [EUR reference]( https://drive.google.com/file/d/1ger1-jsoD73vCez5vN6h4QD5TMdU_WS6/view?usp=sharing "EUR reference") (~8.6G); Google Drive <file_id>: `1ger1-jsoD73vCez5vN6h4QD5TMdU_WS6`; decompress by `tar -zxvf EUR.tar.gz`

     [AFR reference]( https://drive.google.com/file/d/1aGwAGdIeVmTaTskSuUdEPE6Z9-evSlUG/view?usp=sharing "AFR reference") (~11.3G); Google Drive <file_id>: `1aGwAGdIeVmTaTskSuUdEPE6Z9-evSlUG`; decompress by `tar -zxvf AFR.tar.gz`
     
     [AMR reference]( https://drive.google.com/file/d/1XkVeO6ZqMb7Un_HuTeFrQ908FGKx-i1W/view?usp=sharing "AMR reference") (~9.6G); Google Drive <file_id>: `1XkVeO6ZqMb7Un_HuTeFrQ908FGKx-i1W`; decompress by `tar -zxvf AMR.tar.gz`
        
     [EAS reference]( https://drive.google.com/file/d/1NzltrpebQiaHYRIN67nTqmRChJ_0MEze/view?usp=sharing "EAS reference") (~6.2G); Google Drive <file_id>: `1NzltrpebQiaHYRIN67nTqmRChJ_0MEze`; decompress by `tar -zxvf EAS.tar.gz`
             
     [SAS reference]( https://drive.google.com/file/d/1qFKpSLR3i3YikXGInStSQs6HidMdb9lc/view?usp=sharing "SAS reference") (~7.8G); Google Drive <file_id>: `1qFKpSLR3i3YikXGInStSQs6HidMdb9lc`; decompress by `tar -zxvf SAS.tar.gz`
          
     The above LD reference panels were constructed using the 1000 Genomes Project phase 3 samples for variants in either Hapmap3 or MEGA. To create your own LD reference panels, please contact Jingning Zhang (jingningzhang238@gmail.com) for codes.
  
     **New update: We constructed the EUR(~337K), AFR (~3K), EAS (~0.6K), SAS (~1.8K) reference using unrelated UKB samples from the corresponding genetic ancestry.**
        
     [EUR reference constructed using UKB]( https://drive.google.com/file/d/13CFSAN3unf48Z9JnzdAzyJKoF4csjb7K/view?usp=sharing "EUR reference") (~12.5G); Google Drive <file_id>: `13CFSAN3unf48Z9JnzdAzyJKoF4csjb7K`; decompress by `tar -zxvf EUR_ukb.tar.gz`

     [AFR reference constructed using UKB]( https://drive.google.com/file/d/1C-0gXkDP0VCk3tsroPdnQnic0N7KCXau/view?usp=sharing "AFR reference") (~12.4G); Google Drive <file_id>: `1C-0gXkDP0VCk3tsroPdnQnic0N7KCXau`; decompress by `tar -zxvf AFR_ukb.tar.gz`

     [EAS reference constructed using UKB]( https://drive.google.com/file/d/1kFP6I9Tii7B5N2mFcPOaYfmfoSDs847I/view?usp=sharing "EAS reference") (~7.3G); Google Drive <file_id>: `1kFP6I9Tii7B5N2mFcPOaYfmfoSDs847I`; decompress by `tar -zxvf EAS_ukb.tar.gz`

     [SAS reference constructed using UKB]( https://drive.google.com/file/d/106ctOO-ebpiOA-9W-cclmch6h6FAnpUW/view?usp=sharing "SAS reference") (~10.1G); Google Drive <file_id>: `106ctOO-ebpiOA-9W-cclmch6h6FAnpUW`; decompress by `tar -zxvf SAS_ukb.tar.gz`

  
- Download plink 2.0 from [here](https://www.cog-genomics.org/plink/2.0/ "here").

- Launch R and install required libraries:
```
install.packages(c('optparse','bigreadr','readr','stringr', 'caret', 'SuperLearner', 'glmnet', 'MASS', 'Rcpp', 'RcppArmadillo', 'inline', 'doMC', 'foreach'))
```

- Example codes for all steps above:

```
## download PROSPER package
git clone https://github.com/Jingning-Zhang/PROSPER.git

## go to the directory
cd PROSPER

## download reference SNPs
gdown 1PtD4qk7EBPxdhkGrKrG8OktKxaSDF9AT

## download EUR reference LD
gdown 1ger1-jsoD73vCez5vN6h4QD5TMdU_WS6
tar -zxvf EUR.tar.gz
## download AFR reference LD
gdown 1aGwAGdIeVmTaTskSuUdEPE6Z9-evSlUG
tar -zxvf AFR.tar.gz

## You will still need to install plink2 and the required R packages

```

- **If you only want to try the toy example, you can skip the next sections and directly go down to the last section \<Toy Example\>.**


<br/>

## Using PROSPER

There are two scripts `PROSPER.R` and `tuning_testing.R`. The first one is used to get solutions from all tuning parameters setting, and then the second one is used to perform tuning/testing and get the final ensembled PRS. The parameters are explained in this section. At the end of this README file, there is a toy example, including example codes and results.

`
PROSPER.R --PATH_package --PATH_out --FILE_sst --pop --lassosum_param --chrom --Ll --Lc --verbose --NCORES
`

`
tuning_testing.R --PATH_out --PATH_plink --prefix --SL_library --linear_score --bfile_tuning --pheno_tuning --covar_tuning --testing --bfile_testing --pheno_testing --covar_testing --verbose --cleanup --NCORES
`

 - **PATH_package** (required): Full path to the directory mentioned above that contains: 1. a folder scripts downloaded from github 2. LD reference panels downloaded from Google Drive.

 - **PATH_plink** (required): Full path to plink2 executable. Please do not use plink1.9.

 - **PATH_out** (required): Full path to the directory for output results

 - **FILE_sst** (required): Full path and the file name of the GWAS summary statistics from multiple populations, separated by comma. An example of the format
```
    rsid           chr     a1     a0    beta       beta_se    n_eff
    rs3131972      1       G      A     0.00235    0.01269    15000
    rs3131969      1       G      A     -0.01228   0.01206    15000
    rs1048488      1       T      C     0.00989    0.01162    15000

```
Required columns in the files:
  1. rsid: SNP ID in the format of rsXXXX.
  2. chr: chromosome number in the format of 1,2,...,22.
  3. beta: SNP effect. Note that for binary traits, beta is the coefficient in logistic regression, i.e. log(OR). 
  4. beta_se: Standard error of beta.
  5. a1: effective allele (counted allele in regression).
  6. a0: alternative allele (non-A1 allele).
  7. n_eff: Sample size per variant. Note that for binary traits, it is `effective sample sizes = 4 / (1 / N_control + 1 / N_case)`; and for continuous traits, it is simply the sample size.

Note that the summary statistics files are suggested to be cleaned as follows before using:
  1. Keep variants in reference panels to avoid troubles caused by reading huge files. The rsid of variants in reference panels can be found in `ref_bim.txt`.
  2. Alleles are suggested to match to reference panels to avoid flipping strands. The alleles of variants in reference panels can be found in `ref_bim.txt`.
  3. For each population, remove variants whose allele frequencies are less than 1%. If your summary statistics does not contain information of population-specific allele frequencies, we recommanded to use that in 1000G, which can be found in `ref_af.txt`.
  4. Remove variants with small sample size (90% of median sample size per variant).

 - **pop** (required): Population of the GWAS sample (AFR, AMR, EAS, EUR and SAS), in the same order of the GWAS summary statistics files, separated by comma.

 - **lassosum_param** (required): Full path and the file name of the optimal tuning parameters in lassosum2, separated by comma. The file must have the following format
```
    delta0    lambda0
    49.8      0.002
```
The optimal parameters can be obtained by running the official [lassosum2]( https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html#example-ldpred2-and-lassosum2 "lassosum2"), or used [a Simplified version of lassosum2]( https://github.com/Jingning-Zhang/PROSPER/blob/main/lassosum2.md "lassosum2") written in command line using the same data structure as PROSPER. If you have not run the official lassosum2, we suggests using the **Simplified version** instead. It is faster, requires less computation resource, and has a direct output for input parameter `lassosum_param` of PROSPER. Ideally, individual-level tuning data is needed for all populations; however, in the case where the data is only available for a target population, the tuning parameters can be optimized and tuned towards the target population.

 - **chrom** (optional): The chromosome on which the model is fitted, separated by comma or dash for consecutive chromosomes. Default is 1-22.

 - **Ll** (optional): Length of tuning parameter path for the tuning parameter `lambda`. Default is 5.

 - **Lc** (optional): Length of tuning parameter path for the tuning parameter `c`. Default is 5.

 - **prefix** (required): Output filename prefix. Suggested value is name of target population. 

 - **bfile_tuning** (required): Full path to PLINK binary input file prefix (minus bed/bim/fam) for tuning. We suggest using at least 1000 samples for tuning.

 - **pheno_tuning** (optional): Full path to phenotype file for tuning (PLINK format). This is optional, taken from `bfile_tuning` otherwise.

 - **covar_tuning**  (optional): Full path to quantitative covariates for tuning (PLINK format).

 - **testing** (optional): Whether to perform testing in a seperate dataset. If TRUE, testing datasets should be provided (bfile_testing, pheno_testing, covar_testing). Default is FALSE.

 - **bfile_testing** (optional): Should be provided if `testing=TRUE`. Full path to PLINK binary input file prefix (minus bed/bim/fam) for testing.

 - **pheno_testing** (optional): Should be provided if `testing=TRUE`. Full path to phenotype file for testing (PLINK format). This is optional, taken from `bfile_testing` otherwise.

 - **covar_testing**  (optional): Should be provided if `testing=TRUE`. Full path to quantitative covariates for testing (PLINK format).

 - **SL_library** (optional): Models input to SL.library in SuperLearner function, separated by comma. Default is "SL.glmnet,SL.ridge,SL.lm".

 - **linear_score** (optional): Whether to output PRS text file in the format of `rsid, a1, a0, weight`. Default is TRUE; otherwise, only superlearner model will be saved. Note some models available for SL.library are non-linear models, such as SL.nnet. In this case, PRS text file cannot be output.

 - **verbose** (optional): How much chatter to print: 0=nothing; 1=minimal; 2=all. Default is 1.

 - **NCORES** (optional): Number of cores used for parallel computation. Default is 1. Parallel is used on chromosomes in the script PROSPER.R, and used on plink for PRS calculation in the script tuning_testing.R. So the suggested value is 22 for script PROSPER.R, and as many cores as possible for tuning_testing.R. 

 - **cleanup** (optional): Cleanup all temporary files saved in a `tmp` folder, including PRS scores of tuning and testing samples, and all solutions from PROSPER separated by chromosomes. Default is TRUE.


## Output

The script `PROSPER.R` creates a directory `$PATH_out/before_ensemble/`, and writes two files `score_file.txt` and `score_param.txt` inside it.
1. `score_file.txt` contains all solutions from PROSPER
2. `score_param.txt` contains their ancestry origin and corresponding tuning parameter settings

The script `tuning_testing.R` creates a directory `$PATH_out/after_ensemble_$prefix/`, and writes three files `PROSPER_prs_file.txt`, `R2.txt`, and `superlearner_function.RData` inside it. 
1. **`PROSPER_prs_file.txt` is the final ensembled PRS solution from PROSPER**
2. **`R2.txt` is its R2 on tuning and testing samples (if `testing=TRUE`)**
3. `superlearner_function.RData` contains two variables `sl` (super learner model) and `score_drop` (scores in the `$PATH_out/before_ensemble/score_file.txt` that needs to be dropped before input into `sl`). If a non-linear model is specified in `linear_score`, only `superlearner_function.RData` can be saved.


<br/>

## Toy Example

- Please download [example data]( https://drive.google.com/file/d/1IdgiBveqTXuyuQMfWp7ys8ZAmOeB52wa/view?usp=sharing "example data"); Google Drive <file_id>: `1IdgiBveqTXuyuQMfWp7ys8ZAmOeB52wa` and decompress by using codes:

```
gdown 1IdgiBveqTXuyuQMfWp7ys8ZAmOeB52wa
tar -zxvf example.tar.gz
```

- Run the codes as instructed below (direcotries are needed to be changed based on your local directory), and you will get the example output same as [here]( https://drive.google.com/file/d/16CREbH4emYUHudvy5HpHSbwH5111X-9R/view?usp=sharing "example output"); Google Drive <file_id>: `16CREbH4emYUHudvy5HpHSbwH5111X-9R`; decompress by `tar -zxvf PROSPER_example_results.tar.gz` 
    
- Example codes for running lassosum2 and PROSPER:

```
package='/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/try_from_github/PROSPER'
path_example='/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/try_from_github/PROSPER/example/'
path_result='/dcs04/nilanjan/data/jzhang2/MEPRS/pacakge/try_from_github/PROSPER/PROSPER_example_results/'
path_plink='/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2'

mkdir ${path_result}

Rscript ${package}/scripts/lassosum2.R \
--PATH_package ${package} \
--PATH_out ${path_result}/lassosum2 \
--PATH_plink ${path_plink} \
--FILE_sst ${path_example}/summdata/EUR.txt,${path_example}/summdata/AFR.txt \
--pop EUR,AFR \
--chrom 1-22 \
--bfile_tuning ${path_example}/sample_data/EUR/tuning_geno,${path_example}/sample_data/AFR/tuning_geno \
--pheno_tuning ${path_example}/sample_data/EUR/pheno.fam,${path_example}/sample_data/AFR/pheno.fam \
--bfile_testing ${path_example}/sample_data/EUR/testing_geno,${path_example}/sample_data/AFR/testing_geno \
--pheno_testing ${path_example}/sample_data/EUR/pheno.fam,${path_example}/sample_data/AFR/pheno.fam \
--testing TRUE \
--NCORES 5

Rscript ${package}/scripts/PROSPER.R \
--PATH_package ${package} \
--PATH_out ${path_result}/PROSPER \
--FILE_sst ${path_example}/summdata/EUR.txt,${path_example}/summdata/AFR.txt \
--pop EUR,AFR \
--lassosum_param ${path_result}/lassosum2/EUR/optimal_param.txt,${path_result}/lassosum2/AFR/optimal_param.txt \
--chrom 1-22 \
--NCORES 5

target_pop='AFR'

Rscript ${package}/scripts/tuning_testing.R \
--PATH_plink ${path_plink} \
--PATH_out ${path_result}/PROSPER \
--prefix ${target_pop} \
--testing TRUE \
--bfile_tuning ${path_example}/sample_data/${target_pop}/tuning_geno \
--pheno_tuning ${path_example}/sample_data/${target_pop}/pheno.fam \
--bfile_testing ${path_example}/sample_data/${target_pop}/testing_geno \
--pheno_testing ${path_example}/sample_data/${target_pop}/pheno.fam \
--cleanup F \
--NCORES 5

```

The testing R<sup>2</sup> of PROSPER is in `${path_result}/PROSPER/after_ensemble_${target_pop}/R2.txt`; The linear PRS model of PROSPER is in `${path_result}/PROSPER/after_ensemble_${target_pop}/PROSPER_prs_file.txt`. 

To compare with lassosum2 (the single-ancestry PRS method using penalized regression), you could find the testing R<sup>2</sup> of lassosum in `${path_result}/lassosum2/${target_pop}/R2.txt`.

In this toy example, phenotype data is simulated in a previous paper from [Zhang et al.](https://www.biorxiv.org/content/10.1101/2022.03.24.485519v5.abstract). The sample size is 15K and 100K for AFR and EUR, respectively; and the total simulated common-SNP heritability is assumed to be 0.32 and 0.19 for AFR and EUR, respectively. The resulting testing R<sup>2</sup> for AFR is 0.8% using the single-ancestry lassosum2 model, and 1.8% using the multi-ancestry PROSPER model.

<!-- In this example, the single-ancestry lassosum2 generates a testing R<sup>2</sup> of 0.8% and 11.1% for EUR and AFR, respectively; and the multi-ancestry PROSPER generates a testing R<sup>2</sup> of 11.6% and 1.8% for EUR and AFR, respectively. -->

For your reference, as shown in this example, for a two-ancestry analysis (EUR,AFR) on all autosomal chromosomes (1-22) using 5 cores, it takes ~20 minutes and ~25Gb (~5Gb each core) to run PROSPER. Using the same high performance cluster, for a five-ancestry analysis (EUR,AFR,AMR,EAS,SAS) on all autosomal chromosomes (1-22) using 5 cores, it takes ~43 minutes and ~35Gb (~7Gb each core). 


<br/>

## Support

Please direct any problems or questions to Jingning Zhang (jingningzhang238@gmail.com).

## Citation

Zhang, Jingning, et al. "An ensemble penalized regression method for multi-ancestry polygenic risk prediction." Nature Communications 15.1 (2024): 3238.

## Acknowledgment

I would like to acknowledge Dr. Corinna Thoben, Dr. Jacob Williams, Dr. Haoyu Zhang, Ritwiz Kamal, Dr. Manikandan Narayanan, Dr. Julian McClellan for their valuable assistance in identifying bugs and helping to correct the code.

