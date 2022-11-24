# PRS-epr

**PRS-epr** is a command line tool based on R programming language. It is a powerful multi-ancestry PRS method that jointly models GWAS summary statistics from multiple populations by penalized regression and ensemble approach to improve the performance in minority populations.

Preprint manuscript will be put online very soon. Please contact Jingning Zhang (jzhan218@jhu.edu) for citation.


## Getting Started

- Create a directory for storing PRS-epr package and its LD refernce panels. Denoted by `${package}`.

- Download the scripts from `https://github.com/Jingning-Zhang/PRS-epr`, and saved in `${package}` in a folder named **scripts**

- Download the `ref_bim.txt` from [ref_bim]( https://drive.google.com/file/d/1PtD4qk7EBPxdhkGrKrG8OktKxaSDF9AT/view?usp=sharing "ref_bim"), and saved in `${package}`
    
- Download the LD reference panels and decompress files in `${package}`:
   
     [AFR reference]( https://drive.google.com/file/d/1aGwAGdIeVmTaTskSuUdEPE6Z9-evSlUG/view?usp=sharing "AFR reference") (~12.1G); Google Drive <file_id>: `1aGwAGdIeVmTaTskSuUdEPE6Z9-evSlUG`; Decompress by `tar -zxvf AFR.tar.gz`
     
     [AMR reference]( https://drive.google.com/file/d/1XkVeO6ZqMb7Un_HuTeFrQ908FGKx-i1W/view?usp=sharing "AMR reference") (~10.3G); Google Drive <file_id>: `1XkVeO6ZqMb7Un_HuTeFrQ908FGKx-i1W`; Decompress by `tar -zxvf AMR.tar.gz`
        
     [EAS reference]( https://drive.google.com/file/d/1NzltrpebQiaHYRIN67nTqmRChJ_0MEze/view?usp=sharing "EAS reference") (~6.6G); Google Drive <file_id>: `1NzltrpebQiaHYRIN67nTqmRChJ_0MEze`; Decompress by `tar -zxvf EAS.tar.gz`
        
     [EUR reference]( https://drive.google.com/file/d/1ger1-jsoD73vCez5vN6h4QD5TMdU_WS6/view?usp=sharing "EUR reference") (~8.6G); Google Drive <file_id>: `1ger1-jsoD73vCez5vN6h4QD5TMdU_WS6`; Decompress by `tar -zxvf EUR.tar.gz`
     
     [SAS reference]( https://drive.google.com/file/d/1qFKpSLR3i3YikXGInStSQs6HidMdb9lc/view?usp=sharing "SAS reference") (~8.4G); Google Drive <file_id>: `1qFKpSLR3i3YikXGInStSQs6HidMdb9lc`; Decompress by `tar -zxvf SAS.tar.gz`
     
     The files can be downloaded by the provided link. Besides, **gdown** can fast download files from Google Drive to linux by using the command `gdown <file_id>`. **gdown** can be installed using `pip install gdown`.
     
     The LD reference panels were constructed using the 1000 Genomes Project phase 3 samples for variants in either Hapmap3 or MEGA. To create your own LD reference panels, please contact Jingning Zhang (jzhan218@jhu.edu) for codes. 
     
- Launch R and install required libraries:
```
install.packages(c('optparse','bigreadr','readr','stringr', 'caret', 'SuperLearner', 'glmnet', 'MASS', 'Rcpp', 'RcppArmadillo', 'inline', 'doMC', 'foreach'))
```

## Using PRS-epr

There are two scripts `PRS-epr.R` and `tuning_testing.R`. The first one is used to get solutions from all tuning parameters setting, and then the second one is used to perform tuning/testing and get the final ensembled PRS. The parameters are explained in this section. At the end of this README file, there is a toy example, including example codes and results.

`
PRS-epr.R --PATH_package --PATH_out --FILE_sst --pop --lassosum_param --chrom --Ll --Lc --verbose --NCORES
`

`
tuning_testing.R --PATH_out --PATH_plink --prefix --SL_library --linear_score --bfile_tuning --pheno_tuning --covar_tuning --testing --bfile_testing --pheno_testing --covar_testing --verbose --cleanup --NCORES
`

 - **PATH_package** (required): Full path to the directory mentioned above that contains: 1. a folder scripts downloaded from github 2. LD reference panels downloaded from Google Drive.

 - **PATH_plink** (required): Full path to plink2 executable

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
  7. n_eff: Sample size per variant. Note that for binary traits, it is `effective sample sizes=4 / (1 / N_control + 1 / N_case)`; and for continuous traits, it is simply the sample size.

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
The optimal parameters can be obtained by running the official [lassosum2]( https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html#example-ldpred2-and-lassosum2 "lassosum2"), or used [a simplifier version of lassosum2]( https://github.com/Jingning-Zhang/PRS-epr/blob/main/lassosum2.md "lassosum2") written in command line using the same data structure as PRS-epr. If you have not run the official lassosum2, we suggests using the **simplifier version** instead. It is faster, requires less computation resource, and has a direct output for `lassosum_param`.

 - **chrom** (optional): The chromosome on which the model is fitted, separated by comma or dash for consecutive chromosomes. Default is 1-22.

 - **Ll** (optional): Length of tuning parameter path for the tuning parameter `lambda`. Default is 5.

 - **Lc** (optional): Length of tuning parameter path for the tuning parameter `c`. Default is 5.

 - **prefix** (required): Output filename prefix. Suggested value is name of target population. 

 - **bfile_tuning** (required): Full path to PLINK binary input file prefix (minus bed/bim/fam) for tuning. We suggest using at least 1000 samples for tuning.

 - **pheno_tuning** (optional): Full path to phenotype file for tuning (PLINK format). This is optional, taken from `bfile_tuning` otherwise.

 - **covar_tuning**  (optional): Full path to quantitative covariates for tuning (PLINK format).

 - **testing** (optional): Whether to perform testing in a seperate dataset. If TRUE, tuning datasets should be provided (bfile_testing, pheno_testing, covar_testing). Default is FALSE.

 - **bfile_tuning** (optional): Should be provided if `testing=TRUE`. Full path to PLINK binary input file prefix (minus bed/bim/fam) for testing.

 - **pheno_tuning** (optional): Should be provided if `testing=TRUE`. Full path to phenotype file for testing (PLINK format). This is optional, taken from `bfile_tuning` otherwise.

 - **covar_tuning**  (optional): Should be provided if `testing=TRUE`. Full path to quantitative covariates for testing (PLINK format).

 - **SL_library** (optional): Models input to SL.library in SuperLearner function, separated by comma. Default is "SL.glmnet,SL.ridge,SL.lm".

 - **linear_score** (optional): Whether to output PRS text file in the format of `rsid, a1, a0, weight`. Default is TRUE; otherwise, only superlearner model will be saved. Note some models available for SL.library are non-linear models, such as SL.nnet. In this case, PRS text file cannot be output.

 - **verbose** (optional): How much chatter to print: 0=nothing; 1=minimal; 2=all. Default is 1.

 - **NCORES** (optional): Number of cores used for parallel computation. Default is 1. Parallel is used on chromosomes in the script PRS-epr.R, and used on plink for PRS calculation in the script tuning_testing.R. So the suggested value is 22 for script PRS-epr.R, and as many cores as possible for tuning_testing.R. 

 - **cleanup** (optional): Cleanup all temporary files saved in a `tmp` folder, including PRS scores of tuning and testing samples, and all solutions from PRS-epr separated by chromosomes. Default is TRUE.


## Output

The script `PRS-epr.R` creates a directory `$PATH_out/before_ensemble/`, and writes two files `score_file.txt` and `score_param.txt` inside it.
1. `score_file.txt` contains all solutions from PRS-epr
2. `score_param.txt` contains their ancestry origin and corresponding tuning parameter settings

The script `tuning_testing.R` creates a directory `$PATH_out/after_ensemble_$prefix/`, and writes three files `PRSepr_prs_file.txt`, `R2.txt`, and `superlearner_function.RData` inside it. 
1. **`PRSepr_prs_file.txt` is the final ensembled PRS solution from PRS-epr**
2. **`R2.txt` is its R2 on tuning and testing samples (if `testing=TRUE`)**
3. `superlearner_function.RData` contains two variables `sl` (super learner model) and `score_drop` (scores in the `$PATH_out/before_ensemble/score_file.txt` that needs to be dropped before input into `sl`). If a non-linear model is specified in `linear_score`, only `superlearner_function.RData` can be saved.

## Toy Example

Please download [example data]( https://drive.google.com/file/d/1IdgiBveqTXuyuQMfWp7ys8ZAmOeB52wa/view?usp=sharing "example data"); Google Drive <file_id>: `1IdgiBveqTXuyuQMfWp7ys8ZAmOeB52wa`; Decompress by `tar -zxvf example.tar.gz`. Run the codes as instructed below (changing the direcotries), and you will get the example output same as [here]( https://drive.google.com/file/d/1SWUycZ1oxOV9mGv_NEQaRR1Zl9OGgRi6/view?usp=sharing "example output"); Google Drive <file_id>: `1SWUycZ1oxOV9mGv_NEQaRR1Zl9OGgRi6`; Decompress by `tar -zxvf PRS-epr_example_results.tar.gz` 
    
The files can be downloaded by the provided link, or the command `gdown <file_id>`.

Example codes:

```
package='/dcs04/nilanjan/data/jzhang2/PRS-epr'
path_example='/dcs04/nilanjan/data/jzhang2/example/'
path_plink='/dcs04/nilanjan/data/jzhang2/TOOLS/plink/plink2'
target_pop='AFR'

Rscript ${package}/scripts/lassosum2.R \
--PATH_package ${package} \
--PATH_out ${path_example}/PRS-epr_example_results/lassosum2 \
--PATH_plink ${path_plink} \
--FILE_sst ${path_example}/summdata/EUR.txt,${path_example}/summdata/AFR.txt \
--pop EUR,AFR \
--chrom 1-22 \
--bfile_tuning ${path_example}/sample_data/EUR/tuning_geno,${path_example}/sample_data/AFR/tuning_geno \
--pheno_tuning ${path_example}/sample_data/EUR/pheno.fam,${path_example}/sample_data/AFR/pheno.fam \
--bfile_testing ${path_example}/sample_data/EUR/testing_geno,${path_example}/sample_data/AFR/testing_geno \
--pheno_testing ${path_example}/sample_data/EUR/pheno.fam,${path_example}/sample_data/AFR/pheno.fam \
--testing TRUE \
--NCORES 22

Rscript ${package}/scripts/PRS-epr.R \
--PATH_package ${package} \
--PATH_out ${path_example}/PRS-epr_example_results/PRSepr \
--FILE_sst ${path_example}/summdata/EUR.txt,${path_example}/summdata/AFR.txt \
--pop EUR,AFR \
--lassosum_param ${path_example}/PRS-epr_example_results/lassosum2/EUR/optimal_param.txt,${path_example}/PRS-epr_example_results/lassosum2/AFR/optimal_param.txt \
--chrom 1-22 \
--NCORES 22

Rscript ${package}/scripts/tuning_testing.R \
--PATH_plink ${path_plink} \
--PATH_out ${path_example}/PRS-epr_example_results/PRSepr \
--prefix ${target_pop} \
--testing TRUE \
--bfile_tuning ${path_example}/sample_data/${target_pop}/tuning_geno \
--pheno_tuning ${path_example}/sample_data/${target_pop}/pheno.fam \
--bfile_testing ${path_example}/sample_data/${target_pop}/testing_geno \
--pheno_testing ${path_example}/sample_data/${target_pop}/pheno.fam \
--NCORES 50

```


## Support

Please direct any problems or questions to Jingning Zhang (jzhan218@jhu.edu).


