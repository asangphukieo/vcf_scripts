#to convert vcf file from HLA imputation server to table style

args<-commandArgs(TRUE)
vcf<-args[1]
out<-args[2]

#library(vcfR)
#install.packages("vcfR")
#vcf="/data/genDOC/work/sangphukieo/HLA_impute_typing/13_MERGE_IMPUTED_VCF/check_imputation_performance.vcf"
#vcf="/data/genDOC/work/sangphukieo/HLA_impute_typing/13_MERGE_IMPUTED_VCF/check_imputation_performance_test.vcf"

library(vcfR)

snps <- vcfR::read.vcfR(vcf,convertNA = F) #TODO
                             
snps_num <- vcfR::extract.gt(snps, # TODO
           element = "GT",
           IDtoRowNames  = T,
           as.numeric = F,
           convertNA = F,
           return.alleles = F)

tab = t(snps_num)
tab = gsub('\\|','\\/',tab)

for(gen in c("A", "B", "C", "DPB1", "DRB1", "DPA1", "DQA1", "DQB1")){
    A_allele_all=data.frame()
    for(i in rownames(tab)){
      A=tab[i, grepl( paste0("HLA_",gen,".*:") , colnames( tab ) )]
      A_allele=""
      if(any(grepl("1", A))) {
        if(any(grepl("1/1", A))) {
          A_allele=rep(names(A[grep("1/1", A)]),2)
        }
        else{
          A_allele=names(A[grep("1", A)])
        }
      }
      if(length(A_allele)==1){
        A_allele=c(A_allele,".")
      }else if(length(A_allele)==0){
        A_allele=c(".", ".")
      }
      A_allele_tab=data.frame(A1=A_allele[1],A2=A_allele[2])
      rownames(A_allele_tab)=i
      colnames(A_allele_tab)=c(paste0(gen,"1"),paste0(gen,"2"))
      A_allele_all=rbind(A_allele_all,A_allele_tab)
    }
    if(gen=="A"){
        sum_tab=A_allele_all
    }else{
        sum_tab=cbind(sum_tab,A_allele_all)
    }   
}

write.table(sum_tab,out,col.name = T)
