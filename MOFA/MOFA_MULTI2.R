#upload libraries

library(data.table)
library(vroom)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(caret)
library(rlist)
library(MOFA2)
library(reticulate)
#use_python("/usr/bin/python3", required=TRUE)
use_condaenv("BiomiX-env", required = TRUE) 
reticulate::py_config()

#
# # #MANUAL INPUT
# args = as.list(c("BLymphocytes","SJS"))
# args[1] <-"SLE"
# args[2] <-"CTRL"
# args[3] <- "/home/cristia/Scrivania/BiomiX1.7"
# 
# directory <-args[3]
#
# 
Cell_type <- "MOFA"

MART <- vroom(paste(directory,"/MOFA/x_BiomiX_DATABASE/mart_export_37.txt",sep=""), delim = ",")
myList <- list()

COMMAND <- vroom(paste(directory,"COMMANDS.tsv",sep="/"), delim = "\t")
COMMAND_MOFA <- vroom(paste(directory,"COMMANDS_MOFA.tsv",sep="/"), delim = "\t")
COMMAND_ADVANCED <- vroom(paste(directory,"COMMANDS_ADVANCED.tsv",sep="/"), delim = "\t")
Max_features <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY[3])

directory2 <- paste(directory,"/Metabolomics",sep="")

#### METABOLOMICS FUNCTION

Metabolomics_processing <-function(annotation,matrix,mer){
matrix[,-1:-2]<- apply(matrix[,-1:-2],2,log)
annotation <- annotation %>% distinct(NAME.x, .keep_all = TRUE)
annotation$NAMES <- paste(annotation$NAME.x, annotation$Name,sep = "/")
x<- as.data.frame(colnames(matrix))
x$id  <- 1:nrow(x)
colnames(x)[1] <- "Peaks"
xx<-merge(x,annotation,by.x="Peaks", by.y ="NAME.x", all.x = TRUE)
xx <- xx[order(xx$id), ]


for (i in 1:nrow(xx)){
  if (!is.na(xx$NAMES[i])){
    xx$Peaks[i]<- xx$NAMES[i]
  }
}

serum_metabolomics_EASY <- xx$Peaks

list <- list(annotation)


matrix <- matrix %>% filter(CONDITION == args[2] | CONDITION == args[1])
#colnames(matrix) <- gsub("peak","peakS",colnames(matrix))

list[[2]] <- matrix




serum_metabolomics<- reshape(data=list[[2]], idvar="CONDITION",
                             varying = colnames(list[[2]])[3:ncol(list[[2]])],
                             v.name=c("value"),
                             times= colnames(list[[2]])[3:ncol(list[[2]])],
                             new.row.names = 1:1000000,
                             direction="long")
colnames(serum_metabolomics)[3] <- "feature"
serum_metabolomics$view <- paste(mer,"_Metabolomics", sep="")
serum_metabolomics <- serum_metabolomics[,c(1,2,3,5,4)]
colnames(serum_metabolomics)[1] <- "sample"
colnames(serum_metabolomics)[2] <- "group"

#ADD IF STATEMENT IF 2 METABOLOMICS
serum_metabolomics$feature <- gsub("peak",paste("peak_",mer, sep = ""), serum_metabolomics$feature)


list[[1]] <- serum_metabolomics
list[[2]] <- serum_metabolomics_EASY

return(list)

}




#### TRANSCRIPTOMIC FUNCTION


Transcriptomics_processing <-function(annotation,matrix,mer){
        matrix<-merge(matrix, MART, by.x="Gene.name", by.y="Gene.stable.ID")
        matrix<- matrix %>% arrange(desc(variance))
        matrix <-matrix[1:Max_features,]
        
        Bcell_RNAseq_EASY <- paste(matrix$`Gene name`, matrix$`Gene.name`, sep = "/")
        Bcell_RNAseq_VIEW <- as.character(matrix$`Gene name`)
        
        
        matrix <- t(matrix)
        colnames(matrix) <- matrix[1,]
        matrix <- as.data.frame(matrix)
        matrix <- matrix[-1,]
        matrix$ID <- rownames(matrix)
        annotation <- annotation[,c("ID","CONDITION")]
        matrix <- merge(matrix,annotation, by.x = "ID", by.y = "ID")
        
        list <- list(annotation)
        
        
        list[[2]] <- matrix
        
        #Reduction database
        #list[[2]] = list[[2]][1:50,]
        
        Wholeblood_RNAsequi<- reshape(data=list[[2]], idvar="CONDITION",
                                      varying = colnames(list[[2]])[2:(ncol(list[[2]])-1)],
                                      v.name=c("value"),
                                      times= colnames(list[[2]])[2:(ncol(list[[2]])-1)],
                                      new.row.names = 1:100000000,
                                      direction="long")
        colnames(Wholeblood_RNAsequi)[3] <- "feature"
        Wholeblood_RNAsequi$view <-  paste(mer,"_RNAseq", sep="")
        Wholeblood_RNAsequi <- Wholeblood_RNAsequi[,c(1,2,3,5,4)]
        colnames(Wholeblood_RNAsequi)[1] <- "sample"
        colnames(Wholeblood_RNAsequi)[2] <- "group"
        
        
        list[[1]] <- Wholeblood_RNAsequi
        list[[2]] <- Bcell_RNAseq_VIEW
        
        return(list)
}



#### METHYLOMICS FUNCTION


Methylomics_processing <-function(annotation,matrix,metadata,mer){
        annotation <- annotation[,c("gene","CpG_island")]
        matrix<-merge(matrix, annotation, by.x="ID", by.y="CpG_island")
        
        matrix <-matrix[1:Max_features,]
        
        Methylome_WB_VIEW <- as.character(matrix$gene)
        Methylome_WB_EASY <- paste(matrix$ID, matrix$gene, sep = "/")
        
        list <- list(annotation)
        
        matrix <- t(matrix)
        colnames(matrix) <- matrix["ID",]
        matrix <- as.data.frame(matrix)
        y <-rownames(matrix) %in% c("ID","variance")
        matrix <- matrix[!y,]
        matrix$ID <- rownames(matrix)
        
        metadata <- metadata[,c("ID","CONDITION")]
        matrix <- merge(matrix,metadata, by.x = "ID", by.y = "ID")
        
        
        list[[2]] <- matrix
        
        #REDUCTION SAMPLES
        #list[[2]] = list[[2]][1:50,]
        
        Methylome_sequi<- reshape(data=list[[2]], idvar="CONDITION",
                                  varying = colnames(list[[2]])[2:(ncol(list[[2]])-1)],
                                  v.name=c("value"),
                                  times= colnames(list[[2]])[2:(ncol(list[[2]])-1)],
                                  new.row.names = 1:1000000000,
                                  direction="long")
        colnames(Methylome_sequi)[3] <- "feature"
        Methylome_sequi$view <- paste(mer,"_Methylomics", sep="")
        Methylome_sequi <- Methylome_sequi[,c(1,2,3,5,4)]
        colnames(Methylome_sequi)[1] <- "sample"
        colnames(Methylome_sequi)[2] <- "group"
        
        
        list[[1]] <- Methylome_sequi
        list[[2]] <- Methylome_WB_EASY
        
        return(list)
}




##### REARRANGEMENT INPUT1 DATA ----

myList <- list()
for (i in 1:length(COMMAND$INTEGRATION)){


#i <- 3 #N_input
type <- COMMAND$DATA_TYPE[i]
#Prova ad aggiungere questa linea di codice e a sostituire COMMAND con i possibili
#input sulla base del database di comandi di riferimento. 

if(COMMAND$INTEGRATION[i] == "YES"){
if(COMMAND$DATA_TYPE[i] == "Metabolomics"){

#directory2 <- paste(directory,"/Metabolomics",sep="")
directory2 <- paste(directory,"/MOFA/INPUT/", "Metabolomics_", COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], sep ="")
serum_metabolomics <- vroom(paste(directory2,"/Metabolomics_",COMMAND$LABEL[i], "_MOFA.tsv", sep = ""), delim = "\t")
directory2 <- paste(directory,"/Metabolomics/OUTPUT/", COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], sep ="")
serum_annotation <- vroom( paste(directory2,"/",COMMAND$LABEL[i],"_",args[1],"_vs_",args[2],"_results.tsv", sep = ""), delim = "\t")
INPUTX<-Metabolomics_processing(serum_annotation,serum_metabolomics,COMMAND$LABEL[i])
assign(paste("INPUT", i, "_visual", sep=""),INPUTX[[2]])
myList <- list.append(myList,INPUTX[[1]])

}

if(COMMAND$DATA_TYPE[i] == "Transcriptomics"){
        
        print(args[1])
        directory2 <- paste(directory,"/MOFA/INPUT/", COMMAND$LABEL[i],"_",args[1],"_vs_", args[2], sep ="")
        Wholeblood_RNAseq <-  vroom(paste(directory2, "/", COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], "_normalized_vst_variance.tsv",sep = ""), delim = "\t") #read normalization only
        Wholeblood_metadata <-  vroom(paste(directory2, "/","/Metadata_",COMMAND$LABEL[i], "_", args[1],".tsv",sep = ""), delim = "\t")
        INPUTX<-Transcriptomics_processing(Wholeblood_metadata,Wholeblood_RNAseq,COMMAND$LABEL[i])
        assign(paste("INPUT", i, "_visual", sep=""),INPUTX[[2]])
        myList <- list.append(myList,INPUTX[[1]])
        
        
}

if(COMMAND$DATA_TYPE[i] == "Methylomics"){
       

        directory2 <- paste(directory,"/MOFA/INPUT/", "Methylome_",COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], sep ="") 
        Methylome_WB <-  vroom(paste(directory2, "/", COMMAND$LABEL[i], "_matrix_MOFA.tsv",sep = ""), delim = "\t") #read normalization only
        Methylome_metadata <-  vroom(paste(directory2, "/", COMMAND$LABEL[i],"_metadata_MOFA.tsv",sep = "") ,delim = "\t")
        directory2 <- paste(directory,"/Methylomics/OUTPUT/", COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], sep ="")
        Methylome_annotation <- vroom(paste(directory2, "/", "DMP_", COMMAND$LABEL[i], "_Methylome_", args[1] ,"_vs_", args[2],".tsv",sep = ""), delim = "\t", col_names = TRUE)
        INPUTX<-Methylomics_processing(Methylome_annotation,Methylome_WB,Methylome_metadata,COMMAND$LABEL[i])
        assign(paste("INPUT", i, "_visual", sep=""),INPUTX[[2]])
        myList <- list.append(myList,INPUTX[[1]])
        
}
}
}

#### CHOISE MERGING ----

Serum_Urine_Whole_Bcell_Methyl <-rbindlist(myList)
Serum_Urine_Whole_Bcell_Methyl<- as.data.table(Serum_Urine_Whole_Bcell_Methyl)


List <- list()

for (i in 1:length(unique(Serum_Urine_Whole_Bcell_Methyl$view))){
  
  Over<- Serum_Urine_Whole_Bcell_Methyl %>% filter(view == unique(Serum_Urine_Whole_Bcell_Methyl$view)[i])
  Overlap <- as.factor(unique(Over$sample))
  List <- list.append(List,unique(Over$sample))
}

List<- unlist(List)
List <- as.data.frame(List)
#TY<- as.data.frame(table(List) >= 4)
TY<- as.data.frame(table(List) >= as.numeric(COMMAND_MOFA[5,2]))
colnames(TY) <- "Intersection"
saved <- rownames(TY)[TY$Intersection]
out<- Serum_Urine_Whole_Bcell_Methyl$sample %in% saved
Serum_Urine_Whole_Bcell_Methyl <- Serum_Urine_Whole_Bcell_Methyl[out,]



groups<- Serum_Urine_Whole_Bcell_Methyl %>% distinct(sample, .keep_all = TRUE)
groups <- groups[,c(1,2)]
Views<-unique(Serum_Urine_Whole_Bcell_Methyl$view)


Serum_Urine_Whole_Bcell_Methyl

Serum_Urine_Whole_Bcell_Methyl$value<- as.numeric(Serum_Urine_Whole_Bcell_Methyl$value)


#MOFA OBJECT

Serum_Urine_Whole_Bcell_Methyl<- as.data.table(Serum_Urine_Whole_Bcell_Methyl)[,group:=NULL]
unique(Serum_Urine_Whole_Bcell_Methyl$view)[1]

MOFAobject <- create_mofa(Serum_Urine_Whole_Bcell_Methyl)
plot_data_overview(MOFAobject)


#Setting the data options

data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE
head(data_opts)

#Setting the model options

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- as.numeric(COMMAND_MOFA[3,2])
head(model_opts)

#Setting the training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$maxiter <- as.numeric(COMMAND_ADVANCED[1,7])
train_opts$convergence_mode <- as.character(COMMAND_ADVANCED[2,7])
train_opts$freqELBO <- 5
head(train_opts)


#MOFA ANALYSIS
#setwd("/media/cristia/KINGSTON/MOFA/Plot")

MOFAobject@samples_metadata[["condition"]] <- groups$group


MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile = file.path(getwd(),"model.hdf5")
#outfile = "/home/cristia/R/x86_64-pc-linux-gnu-library/4.2/MOFA2/extdata/model.hdf5"

reticulate::py_config()
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = FALSE)



#DOWNSTREAM ANALYSIS

dir.create(path = paste(directory,"/MOFA/OUTPUT/","MOFA", "_", args[1] ,"_vs_", args[2], "_",as.numeric(COMMAND_MOFA[3,2]),"_factors" ,sep="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")

directory2 <- paste(directory,"/MOFA/OUTPUT/","MOFA", "_", args[1] ,"_vs_", args[2], "_",as.numeric(COMMAND_MOFA[3,2]),"_factors" ,sep="") 
setwd(directory2)


pdf(file= paste("MOFA", "_", args[1] ,"_vs_", args[2],".pdf", sep =""), width = 20, height = 9)

#package upload
library(ggplot2)
library(MOFA2)

#upload trained model
#filepath <- system.file("extdata", "model.hdf5", package = "MOFA2")
print(outfile)

model <- load_model(outfile)
p <-plot_data_overview(model)
p <-p + theme(text=element_text(size=18,face = "bold"))
print(p)
model@samples_metadata[["condition"]] <- groups$group

#Add metadata to the model

Nsamples = sum(model@dimensions$N)

sample_metadata <- data.frame(
  sample = samples_names(model)[[1]],
  condition = groups$group)

samples_metadata(model) <- sample_metadata
head(model@samples_metadata, n=3)


#Correlation between factors

#A good sanity check is to verify that the Factors are largely uncorrelated. 
#In MOFA there are no orthogonality constraints such as in Principal Component Analysis, 
#but if there is a lot of correlation between Factors this suggests a poor model fit. Reasons? 
#Perhaps you used too many factors or perhaps the normalisation is not adequate.

p <- plot_factor_cor(MOFAobject.trained)

# Total variance explained per view and group
head(model@cache$variance_explained$r2_total[[1]]) # group 1

# Variance explained for every factor in per view and group
head(model@cache$variance_explained$r2_per_factor[[1]]) # group 1

model@cache$variance_explained$r2_per_factor[[1]]

n_factor<-as.numeric(COMMAND_MOFA[3,2])

setwd(directory2)

line = model@cache$variance_explained$r2_per_factor[[1]]
line<-round(line,2)
total<-apply(line,2,sum)
total<-rbind(line,total)
write.table(total,file=paste("MOFA_MATRIX", "_", args[1] ,"_vs_", args[2],"_",n_factor,"_factors.tsv", sep =""))


# u<-as.data.frame(model@cache$variance_explained$r2_per_factor[[1]])
# u$factor <- rownames(u)
# write.table(u, paste("MOFA_MATRIX", "_", args[2] ,"_vs_", args[1], sep =""),quote = FALSE, row.names = FALSE, sep = "\t")
# 


#Plot variances explained
p<-plot_variance_explained(model, x="view", y="factor", max_r2 = 10,)
p <-p + theme(text=element_text(size=15,face = "bold"))
print(p)

#Plot variances explained total
x<-plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]
x<-x + theme(text=element_text(size=11,face = "bold"))
x<-x + theme(axis.title=element_text(size=15))
x<-x + theme(axis.text.x=element_text(size=12))

print(x)

factors <- get_factors(model, as.data.frame = T)
head(factors, n=3)


#TEST LEVEL OF SEPARATION FOR CANDIDATE FACTOR
factors <- get_factors(model, factors = "all")
meta<-model@samples_metadata
SLE <- meta$condition == args[1]
SLE <- factors$single_group[SLE,]
CTRL <- meta$condition == args[2]
CTRL <- factors$single_group[CTRL,]


u=1
pval =NULL
means = NULL
stdv =NULL

for (u in 1:ncol(factors$single_group)){
        res <- wilcox.test(SLE[,u],CTRL[,u], alternative = "two.sided")
        means <- append(means, (mean(SLE[,u]) - mean(CTRL[,u])))
        stdv <- append(stdv, (sd(SLE[,u]) - sd(CTRL[,u])))
        pval <-append(pval, res[["p.value"]])
}
p <- p.adjust(pval, method ="fdr", length(pval))
D<-as.data.frame(colnames(SLE))
colnames(D) <-"Factors"
D$p.adj <- p
D$means <- means
D$standard_deviation <- stdv

write.table(D,file=paste("MOFA_SEPARATION", "_", args[1] ,"_vs_", args[2],"_",n_factor,"_factors.tsv", sep =""),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


#REMOVAL MOFA TAG FROM TRANSCRIPTOMICS
if(sum(COMMAND$DATA_TYPE == "Transcriptomics") > 1){
position <- which(COMMAND$DATA_TYPE %in% "Transcriptomics")

for (i in position) {       
  l<-list()
  x<-as.data.frame(myList[i])
  l[[1]]<- unique(x$view)
  
  
features_names(model)[[l[[1]]]] <- gsub(paste("_", l[[1]] , sep=""), "" ,features_names(model)[[l[[1]]]])

head(features_names(model)[[l[[1]]]])
}
}

#TEST
n=1

for(i in 1:length(COMMAND$INTEGRATION)){
  print(i)
  
  if(COMMAND$INTEGRATION[i] == "YES"){
    
    if(COMMAND$DATA_TYPE[i] == "Metabolomics"){
      
      o<- list("Serum_metabolomics" = get(paste("INPUT",i,"_visual",sep=""))[-1:-2])
      names(o) <- Views[n]
      print(o)
      features_names(model)[Views[n]] <- o
      n = n+1
      print("OK")
    }else{
      print(n)
      o<- list("Serum_metabolomics" = get(paste("INPUT",i,"_visual",sep="")))
      names(o) <- Views[n]
      print(o)
      features_names(model)[Views[n]] <- o
      n = n+1
      print("OKK")
    }
  }
}

        
# features_names(model)[Views[1]] <- list("Serum_metabolomics" = INPUT1_visual[-1:-2])
# features_names(model)["Urine_metabolomics"] <- list("Urine_metabolomics" = urine_metabolomics_EASY[-1:-2])
# features_names(model)["Methylomics_WB"] <- list("Methylomics_WB" = Methylome_WB_EASY)
#here we can add the gene name associated to each CpG island if required




#FACTORS PLOTTING

x <- plot_factor(model, 
                 factor = 1:n_factor,
                 group_by = "condition",
                 color_by = "condition",
                 shape_by = "condition",
                 legend = FALSE
)

t<-c("CTRL" = "blue", "SLE" = "red")
attr(t, "names")[1]<- args[2]
attr(t, "names")[2]<- args[1]
attr(t, "names")


x <-x + theme(text=element_text(size=14,face = "bold"))


x <- x + 
  scale_color_manual(values=t) +
  scale_fill_manual(values=t)
print(x)

#Adding more options
p <- plot_factor(model, 
                 factors = c(1:n_factor),
                 color_by = "condition",
                 group_by = "condition",
                 dot_size = 2,        # change dot size
                 dodge = T,           # dodge points with different colors
                 legend = F,          # remove legend
                 add_violin = T,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)


# The output of plot_factor is a ggplot2 object that we can edit

t<-c("CTRL" = "blue", "SLE" = "red")
attr(t, "names")[1]<- args[2]
attr(t, "names")[2]<- args[1]
attr(t, "names")


p <-p + theme(text=element_text(size=14,face = "bold"))

p <- p + 
  scale_color_manual(values=t) +
  scale_fill_manual(values=t)

print(p)


#Visualisation of combinations of factors


x<- plot_factors(model,
                factors = 1:n_factor,
                color_by = "condition"
)

print(x)

#Visualisation of feature weights

for(i in Views) {
  
  x<-plot_weights(model,
                  view = i,
                  factor = as.numeric(COMMAND_MOFA[4,2]),
                  nfeatures = 10,     # Number of features to highlight
                  scale = T,          # Scale weights from -1 to 1
                  abs = F, # Take the absolute value?
                  text_size = 3,
  )
  
  print(x)
}

#ADD SAVING OF THE TOP WEIGHTS

nice<-which(D$p.adj < 0.05)
if (length(nice) != 0){
        
        
        pdf(file= paste("Distribution_factors_contributions_", args[1] ,"_vs_", args[2],".pdf", sep =""), width = 10, height = 7)
        
        for ( ir in nice) {
                
                
        y = ir #Factor of interest

        
        print(y)
        
        #USE OF PERCENTILE
        
        weights <- get_weights(model, views = "all", factors = "all")
        fact <- as.data.frame(get_factors(model, factors = "all")[[1]])
        fact$ID <- rownames(fact)
        fact <- fact[, c(ncol(fact), 1:(ncol(fact)-1))]
        write.table(fact, "MOFA_FACTORS_WEIGHTS.tsv",quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
        
        factor = y
        
        for (dimen in 1:model@dimensions[["M"]]){
        
        #=====BCELL=========
        weights_B <- as.data.frame(weights[dimen])
        a1 <- as.data.frame(weights_B[,factor])
        colnames(a1) <- "weights"
        a1$weights<- a1$weights/max(abs(a1$weights))
        a1$features <- rownames(weights_B)
        #print(a1$features)
        
        
        
        a1$abs_weight<-abs(a1$weights)
        nov<-quantile(a1$abs_weight, probs=0.95)
        
        p <-ggplot(a1, aes(x=abs_weight)) +
                geom_density(color="darkblue", fill="lightblue")
        #geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.7, bins=1000) + stat_bin(bins = 1000)
        
        p<- p + labs(title=paste("Factor ", y ," weights in_", names(model@data[dimen]), sep="" ))
        
        p <-p + theme(text=element_text(size=12,face = "bold"), plot.title = element_text(hjust = 0.5, size=25,face = "bold"))
        
        p <- p + geom_vline(aes(xintercept=nov))
        #p <- p + geom_vline(aes(xintercept=cin))
        print(p)
        
        pos<-a1[a1$abs_weight > nov,]
        pos
        pos<-pos[order(pos$abs_weight),]
        
        #setwd("/home/cristia/Scrivania/PhD/Bioinfo/MOFA_integration/Databases_copia/Results_Optimization/MOFA_14_factors/Weights_5_percent") 
        poss<-pos$weights > 0.50
        negg<-pos$weights < -0.50
        possis <- as.data.frame(pos[poss,])
        neggis <- as.data.frame(pos[negg,])
        possi <- as.data.frame(pos[poss,2])
        neggi <- as.data.frame(pos[negg,2])
        
        write.table(possis, paste(names(model@data[dimen]), "_pos_factor_",y,".tsv",sep=""),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
        write.table(neggis, paste(names(model@data[dimen]), "_neg_factor_",y,".tsv",sep=""),quote = FALSE, row.names = FALSE, col.names = TRUE,sep = "\t")
        #write.table(possi, paste(names(model@data[dimen]), "_pos_factor_",y,"_ENRICHR",sep=""),quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
        #write.table(neggi, paste(names(model@data[dimen]), "_neg_factor_",y,"_ENRICHR",sep=""),quote = FALSE, row.names = FALSE, col.names = FALSE,sep = "\t")
        
        
        
        }


        }
        dev.off()
}
        
#Visualisation of feature top weights

for(i in Views) {
  
  x<-plot_top_weights(model,
                      view = i,
                      factor =  as.numeric(COMMAND_MOFA[4,2]),
                      nfeatures = 35
  )
  print(x)
}


#Visualisation of patterns in the input data


dev.off()

pdf(file= paste("MOFA", "_", args[1] ,"_vs_", args[2], "_2",".pdf", sep =""), width = 20, height = 9)


#Heatmaps

for(i in Views) {
  
  x <- plot_data_heatmap(model,
                         view = i,         # view of interest
                         factor =  as.numeric(COMMAND_MOFA[4,2]),             # factor of interest
                         features = 20,          # number of features to plot (they are selected by weight)
                         
                         # extra arguments that are passed to the `pheatmap` function
                         cluster_rows = TRUE, cluster_cols = TRUE,
                         show_rownames = TRUE, show_colnames = TRUE, annotation_samples = "condition")
  
  
  
  x <- x + 
    scale_color_manual(values=t) +
    scale_fill_manual(values=t)
  
  
  print(x)
  
  
}



#Scatter plots

for(i in Views) {
  
  x<- plot_data_scatter(model,
                        view = i,         # view of interest
                        factor =  as.numeric(COMMAND_MOFA[4,2]),             # factor of interest
                        features = 10,           # number of features to plot (they are selected by weight)
                        add_lm = TRUE,          # add linear regression
                        color_by = "condition"
  )
  
  x <- x + 
    scale_color_manual(values=t) +
    scale_fill_manual(values=t)

  print(x)
  
}


dev.off()
#dev.off()



