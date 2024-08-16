## Codes for the GenomicSEM across 14 cognitive traits ##

## The GWAS summary statistics are available at https://yanglab.westlake.edu.cn/data/ukb_fastgwa/imp/
## the index of the data are list in the GWAS_names.txt

library(data.table)
file_name = fread('GWAS_names.txt')

## GenomicSEM are available at https://github.com/GenomicSEM/GenomicSEM

library(GenomicSEM)

setwd("the path of GWAS summary data")

#list the full name of the GWAS summary files (maybe different if the data was updated) 
filenames = paste0(file_name$name,'.v1.1.fastGWA.gz')

# hapmap3 reference: https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v
snplist   = 'L:/GWAS/SpeAcc/Batch/gSEM/w_hm3.snplist'

# munge the GWAS summary data #
for(i in 1:length(filenames)){
  munge_file = filenames[i]
  munge(munge_file, 
        snplist,
        trait.names = file_name$name,
        info.filter = 0.9, 
        maf.filter = 0.01)
}

## run LDSC ##
# list the name of output of munge function
sumfiles = paste0(file_name$name,".sumstats.gz")

# download the reference of ld and wld: eur_w_ld_chr at https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v

LDSCoutp <- ldsc(traits = sumfiles, trait.names = file_name$abbr, 
                 sample.prev = file_name$prev, population.prev = file_name$prev,
                 ld = "eur_w_ld_chr/",
                 wld = "eur_w_ld_chr/")
save(LDSCoutp,file = "LDSCoutp.RData")

############################# Model ################################
# Exploratory Factor Analysis (EFA)
library(Matrix)

load('LDSCoutp.RData')
Ssmooth<-as.matrix((nearPD(LDSCoutp$S, corr = FALSE))$mat)
EFA<-factanal(covmat = Ssmooth, factors = 2, rotation = "promax")
EFA$loadings


# Confirmatory Factor Analysis (CFA) according to the above EFA results （loading > 0.4)
model <- 'F1 =~ Pari_spd + PM_spd + RT_acc + RT_spd + SDS_spd + SDS_acc + TM_err1 + TM_spd1 + TM_spd2
           F2 =~ FI_acc + NM_acc + Pari_acc + PM_acc + TM_err1 + TM_err2
           F1 ~~ F2
           Pari_spd ~~ Pari_acc
           PM_spd ~~ PM_acc
           RT_spd ~~ RT_acc
           SDS_spd ~~ SDS_acc
           TM_spd1 ~~ TM_err1
           TM_spd2 ~~ TM_err2'
CFA <- usermodel(LDSCoutp, estimation = "DWLS", model = model, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
CFA 

## model with SNPs
N = file_name$N
ref   = 'reference.1000G.maf.0.005.txt' #https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v
logit <-c(F,F,F,F,F,F,F,F,F,F,F,F,F,F)
cog_sumstats <- sumstats(files = filenames,ref=ref,trait.names=file_name$abbr,se.logit=logit,OLS=NULL,N=N)

model <- 'F1 =~ Pari_spd + PM_spd + RT_acc + RT_spd + SDS_spd + SDS_acc + TM_err1 + TM_spd1 + TM_spd2
          F2 =~ FI_acc + NM_acc + Pari_acc + PM_acc + TM_err1 + TM_err2
          F1 ~~ F2
          Pari_spd ~~ Pari_acc
          PM_spd ~~ PM_acc
          RT_spd ~~ RT_acc
          SDS_spd ~~ SDS_acc
          TM_spd1 ~~ TM_err1
          TM_spd2 ~~ TM_err2
          F1 ~ SNP
          F2 ~ SNP'
CFAFactors <- userGWAS(covstruc=LDSCoutp,SNPs=cog_sumstats,model=model,sub=c("F1~SNP", "F2~SNP"),Q_SNP=TRUE)

## The results of the two factors (e.g. CPS and CPA) have been deposited in the GWAS Catalog (CPS: GCST90446168 and CPA: GCST90446169)

## Calculate effective Sample Size for the factors (https://github.com/GenomicSEM/GenomicSEM/wiki/5.-Multivariate-GWAS)
# restrict to MAF of 40% and 10%
CFAFactors<-subset(CFAFactors, CFAFactors$MAF <= .4 & CFAFactors$MAF >= .1)
N_hat_F1<-mean(1/((2*CFAFactors[[1]]$MAF*(1-CFAFactors[[1]]$MAF))*CFAFactors[[1]]$SE^2))