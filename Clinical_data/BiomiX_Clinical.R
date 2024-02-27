#ADD THE CORRELATION ANALYSIS WITH CLINICAL DATA
library(dplyr)
library(tidyverse)
library(vroom)

# MANUAL INPUT
# # #
# args = as.list(c("Neutrophils","PAPS"))
# args[1] <-"C4"
# args[2] <-"Control"
# args[3] <-"/home/henry/Desktop/BiomiX2.2"
# #
# directory <- args[3]


COMMAND_ADVANCED <- vroom(paste(directory,"COMMANDS_ADVANCED.tsv",sep="/"), delim = "\t")
COMMAND_MOFA <- vroom(paste(directory,"COMMANDS_MOFA.tsv",sep="/"), delim = "\t")


NUMERICAL_AVAILABLE = as.character(COMMAND_ADVANCED[1,12]) 
BINARY_AVAILABLE = as.character(COMMAND_ADVANCED[2,12]) 

files <- grep(paste("MOFA", "_", args[1] ,"_vs_", args[2],"*", sep=""),list.files(paste(directory,"/MOFA/OUTPUT/",sep="")),value=TRUE)
directory <- args[[3]]

for (fil in files){

directory3 <- paste(directory,"/MOFA/OUTPUT/",fil,sep="") 
setwd(directory3)
print(directory3)

if (file.exists("MOFA_FACTORS_WEIGHTS.tsv")){
FACTORS_WEIGHTS <-read.table("MOFA_FACTORS_WEIGHTS.tsv",sep = "\t", header = TRUE)

#detection significant factors analyzed by BiomiX
files <- grep("*factor_",list.files(directory3),value=TRUE)
files<- files[grep("\\Metabolomics|\\RNAseq|\\Methylomics", files)]
factors<-strsplit(files, "_")
factors<-unlist(factors)
factors<-str_remove(factors,".tsv")
factors<-unique(as.numeric(factors[!is.na(as.numeric(factors))]))


#upload data

directory2 <- paste(directory,"/Clinical_data/INPUT",sep="")
setwd(directory2)


if(NUMERICAL_AVAILABLE == "YES"){
NUMERICAL <-read.table(COMMAND_ADVANCED$ADVANCED_OPTION_CLINIC_DATA_DIRECTORY[1], sep = "\t", header = TRUE)


dir.create(path = paste(directory,"/Clinical_data/OUTPUT/","Clinical", "_", args[1] ,"_vs_", args[2], "_",as.numeric(COMMAND_MOFA[3,2]),"_factors" ,sep="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
directory2 <- paste(directory,"/Clinical_data/OUTPUT/","Clinical", "_", args[1] ,"_vs_", args[2], "_",as.numeric(COMMAND_MOFA[3,2]),"_factors" ,sep="") 
setwd(directory2)

commons<-intersect(NUMERICAL$ID, FACTORS_WEIGHTS$ID)
FACTORS_WEIGHTS<-FACTORS_WEIGHTS[FACTORS_WEIGHTS$ID %in% commons,]
NUMERICAL<-NUMERICAL[NUMERICAL$ID %in% commons,]
FACTORS_WEIGHTS <- FACTORS_WEIGHTS[order(FACTORS_WEIGHTS$ID),]
NUMERICAL <- NUMERICAL[order(NUMERICAL$ID),]


for (fact in factors){
  
        Results = matrix(0, nrow= ncol(NUMERICAL) -1 , ncol= 3)
        for (i in 2:ncol(NUMERICAL)){
          
          try_result <- tryCatch({
            corr <- cor.test(as.numeric(FACTORS_WEIGHTS[,fact + 1]) ,as.numeric(NUMERICAL[,i]))
            print(corr)
            #Results[i-1,]<-c(colnames(NUMERICAL)[i],corr[["p.value"]],corr[["estimate"]][["cor"]])
            tes <- c(colnames(NUMERICAL)[i],corr[["p.value"]],corr[["estimate"]][["cor"]])
            print(tes)
            
          }, error = function(err) {
            message(paste("Not enough observation for", colnames(NUMERICAL)[i], "- Skipping the numeric clinical data, because of", conditionMessage(err)))
            Sys.sleep(2)  # Add a delay before retrying
            tes <- c(colnames(NUMERICAL)[i],"X","X")
            print(tes)

          })
          
          Results[i-1,] <- try_result
        }
        
        colnames(Results) <- c("Numeric_clinical_feature", "p.value", "Correlation(R)")
        Results<-as.data.frame(Results)
        Results<-Results[!Results$Numeric_clinical_feature == 0,]
        Results[,2]<-as.numeric(Results[,2])
        Results[,3]<-as.numeric(Results[,3])
        str(Results)
        real<-length(Results$`Correlation(R)`) - length(Results$`Correlation(R)`[is.na(Results$`Correlation(R)`)])
        Results$p.adj <- p.adjust(Results$p.value,method = "fdr", n = as.numeric(real))
        Results <- Results[order(Results$p.adj),]
        print(fact)
        
        write.table(Results, paste("Correlation_NUMERIC_clinical_data_vs_MOFA_factor_",fact,".tsv",sep=""),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
        dir.create(path = paste(directory,"/MOFA/OUTPUT/",fil,"/Clinical_correlation", sep="")  ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        write.table(Results, paste(directory,"/MOFA/OUTPUT/",fil,"/Clinical_correlation/Correlation_NUMERIC_clinical_data_vs_MOFA_factor_",fact,".tsv",sep=""),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
        
}


directory2 <- paste(directory,"/Clinical_data/INPUT",sep="")
setwd(directory2)

}

if(BINARY_AVAILABLE == "YES"){
BINARY <-vroom(COMMAND_ADVANCED$ADVANCED_OPTION_CLINIC_DATA_DIRECTORY[2], delim = "\t", col_names = TRUE)

dir.create(path = paste(directory,"/Clinical_data/OUTPUT/","Clinical", "_", args[1] ,"_vs_", args[2], "_",as.numeric(COMMAND_MOFA[3,2]),"_factors" ,sep="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
directory2 <- paste(directory,"/Clinical_data/OUTPUT/","Clinical", "_", args[1] ,"_vs_", args[2], "_",as.numeric(COMMAND_MOFA[3,2]),"_factors" ,sep="") 
setwd(directory2)

commons<-intersect(BINARY$ID, FACTORS_WEIGHTS$ID)
FACTORS_WEIGHTS<-FACTORS_WEIGHTS[FACTORS_WEIGHTS$ID %in% commons,]
BINARY<-BINARY[BINARY$ID %in% commons,]


BINARY <- as.data.frame(BINARY)
t<-is.na(BINARY) 
BINARY[t] <- "Unknown"

for (fact in factors){
        
        BINARY$factors<- FACTORS_WEIGHTS[,fact]
        
        pval=NULL
        clinic=NULL
        means=NULL
        for (i2 in 2:length(BINARY)){
                print(i2)
                SLE <-  BINARY[,i2] == "Yes" | BINARY[,i2] == "Present" | BINARY[,i2] == "Past" | BINARY[,i2] == "Male"
                if (sum(SLE) == 0){
                        SLE = NULL
                        }else{
                                SLE <- BINARY[SLE,ncol(BINARY)]
                                }
                
                CTRL <- BINARY[,i2] == "No" | BINARY[,i2] == "Female"
                if (sum(CTRL) == 0){ 
                  CTRL = NULL }else{
                                CTRL <- BINARY[CTRL,ncol(BINARY)]
                                }
                
                if(length(SLE) != 0 && length(CTRL) != 0){
                        
                        res <- wilcox.test(SLE,CTRL, alternative = "two.sided")
                        means <- append(means, (mean(SLE) - mean(CTRL)))
                        pval <-append(pval, res[["p.value"]])
                        clinic<- append(clinic, colnames(BINARY)[i2])
                }else{
                        res <- 1
                        pval <-append(pval, res)
                        means <- append(means, (mean(SLE) - mean(CTRL)))
                        clinic<- append(clinic,colnames(BINARY)[i2])
                }
        }
        
        sum<-as.data.frame(clinic)
        sum$means <- means
        sum$pval <-pval
        real<-length(sum$means) - length(sum$means[is.na(sum$means)])
        outsa<- sum[is.na(sum$means),]
        outsa_2<- sum[!is.na(sum$means),]
        outsa_2$p.adj <-p.adjust(outsa_2$pval,method = "fdr", n = as.numeric(real))
        outsa_2 <- outsa_2[order(outsa_2$p.adj),]
        sum$p.adj = "X"
        sum<-sum[!sum$clinic %in% outsa_2$clinic,] 
        sum<-rbind(outsa_2,sum)

        write.table(sum, paste("Correlation_BINARY_clinical_data_vs_MOFA_factor_",fact,".tsv",sep=""),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
        dir.create(path = paste(directory,"/MOFA/OUTPUT/",fil,"/Clinical_correlation", sep="")  ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        write.table(sum, paste(directory,"/MOFA/OUTPUT/",fil,"/Clinical_correlation/Correlation_BINARY_clinical_data_vs_MOFA_factor_",fact,".tsv",sep=""),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}
}
}else{
        print(paste("No_significant_factor in ", fil, sep = ""))
}
}
