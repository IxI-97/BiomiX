
library(enrichR)
library(vroom)
library(metpath)
library(tidyverse)
library(dplyr)

# MANUAL INPUT
# args = as.list(c("BLymphocytes","SJS"))
# args[1] <-"C4"
# args[2] <-"Control"
# args[3] <-"/home/henry/Desktop/BiomiX2.2"
# 
# directory <-args[3]
# 
# COMMAND <- vroom(paste(directory,"COMMANDS.tsv",sep="/"), delim = "\t")
# COMMAND_MOFA <- vroom(paste(directory,"COMMANDS_MOFA.tsv",sep="/"), delim = "\t")
# COMMAND_ADVANCED <- vroom(paste(directory,"COMMANDS_ADVANCED.tsv",sep="/"), delim = "\t")


padj_pathways <- as.numeric(COMMAND_ADVANCED[1,11]) 
n_pathways <- as.numeric(COMMAND_ADVANCED[2,11]) 

#BLOCK TRANSCRIPTOMICS
directory <- args[[3]]

files_out <- grep(paste("MOFA", "_", args[1] ,"_vs_", args[2],"*", sep=""),list.files(paste(directory,"/MOFA/OUTPUT/",sep="")),value=TRUE)
dbs <- c("Reactome_2022", "GO_Biological_Process_2023", "CODE_and_ChEA_Consensus_TFs_from_ChIP-X")
websiteLive <- getOption("enrichR.live")


#For loop required
for (fil in files_out){

directory2 <- paste(directory,"/MOFA/OUTPUT/",fil ,sep="") 
setwd(directory2)


#detection total significant factors
files2 <- grep("*factor_",list.files(directory2),value=TRUE)
files2<- files2[grep("\\Metabolomics|\\RNAseq|\\Methylomics", files2)]
factors<-strsplit(files2, "_")
factors<-unlist(factors)
factors<-str_remove(factors,".tsv")
factors<-unique(as.numeric(factors[!is.na(as.numeric(factors))]))

for (numb in factors){

files <- grep(paste("*factor_", numb, sep=""),list.files(directory2),value=TRUE)
files<- files[grep("\\Metabolomics|\\RNAseq|\\Methylomics", files)]
nam<-strsplit(files, "\\Metabolomics|\\RNAseq|\\Methylomics")
nam<-unlist(nam)
nam<-unique(nam[-grep("*factor*", nam)])
nam_adv <- substr(nam, 1, nchar(nam) - 1)
nam_adv<-COMMAND$DATA_TYPE[match(nam_adv, COMMAND$LABEL)]

iter=0

for (omik in (nam)){
        iter <- iter + 1
        print(omik)
        val<-grep(paste(omik,".*_pos_","*", sep=""), files)
        LIST_GENE<- vroom(paste(directory2,"/",files[val], sep=""), delim="\t")
        LIST_GENE<-LIST_GENE %>% filter(abs_weight>0.50) %>% arrange(desc(abs_weight))
        #IF there is only peak the query is removed, peak is removed from the total name
        LIST_GENE$features <-gsub("peak[0-9]*/", "", LIST_GENE$features)
        #Removal peak name
        sav<-grep("peak[0-9]*", LIST_GENE$features)
        if(length(sav) == 0){
                print("selection...")}else{
                        LIST_GENE<- LIST_GENE[-sav,]}
        
        #Retrieve gene name from cpg islands
        LIST_GENE$features <-gsub("cg[0-9]*/", "", LIST_GENE$features)
        sav<-grep("NA", LIST_GENE$features)
        if(length(sav) != 0){
                LIST_GENE<- LIST_GENE[-sav,]}
        
        if(nam_adv[iter] == "Metabolomics"){
        annot<-vroom(paste(directory,"/Metabolomics/OUTPUT/", omik, args[1],"_vs_",args[2], "/", omik, args[1],"_vs_",args[2],"_results.tsv", sep = ""), delim = "\t")
        LIST_GENE<-LIST_GENE[!LIST_GENE$features == "NA",]
        LIST_GENE<-LIST_GENE$features
        if("Name" %in% colnames(annot)){ 
                LIST_GENE<-annot[match(LIST_GENE, annot$Name),]
                LIST_GENE<-LIST_GENE$HMDB}
        
        }        
        
        val<-grep(paste(omik,".*_neg_","*", sep=""), files)
        LIST_GENE2<-vroom(paste(directory2,"/",files[val], sep=""), delim="\t")
        LIST_GENE2<-LIST_GENE2 %>% filter(abs_weight>0.50) %>% arrange(desc(abs_weight))
        LIST_GENE2$features <-gsub("peak[0-9]*/", "", LIST_GENE2$features)
        #Removal peak name
        sav<-grep("peak[0-9]*", LIST_GENE2$features)
        if(length(sav) == 0){
                print("selection...")}else{
                        LIST_GENE2<- LIST_GENE2[-sav,]}
        #Retrieve gene name from cpg islands
        LIST_GENE2$features <-gsub("cg[0-9]*/", "", LIST_GENE2$features)
        sav<-grep("NA", LIST_GENE2$features)
        if(length(sav) != 0){
                LIST_GENE2<- LIST_GENE2[-sav,]}
        
        if(nam_adv[iter] == "Metabolomics"){
                annot<-vroom(paste(directory,"/Metabolomics/OUTPUT/", omik, args[1],"_vs_",args[2], "/", omik, args[1],"_vs_",args[2],"_results.tsv", sep = ""), delim = "\t")
                LIST_GENE2<-LIST_GENE2[!LIST_GENE2$features == "NA",]
                LIST_GENE2<-LIST_GENE2$features
                if("Name" %in% colnames(annot)){ 
                        LIST_GENE2<-annot[match(LIST_GENE2, annot$Name),]
                        LIST_GENE2<-LIST_GENE2$HMDB}
                
        }        
        
        

        

        
        dir.create(path = paste(directory2,"/MOFA_Factor_Pathways",sep="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        pdf(file= paste("MOFA_Factor_Pathways/", omik, "factor_", numb,".pdf", sep =""), width = 20, height = 9)
        directory_path <-  paste(directory2,"/MOFA_Factor_Pathways",sep="")
        
        if(nam_adv[iter] == "Transcriptomics" | nam_adv[iter] == "Methylomics"){
        
                if (websiteLive) {
                enriched <- enrichr(LIST_GENE$features, dbs)}
        
        for(paths in 1:3){
                if (length(enriched[[paths]]) != 0){
                if (nrow(enriched[[paths]]) != 0){
                if (websiteLive) {
                dir.create(paste(directory_path,"/TABLES", sep=""))
                write.table(x=enriched[[paths]], file=paste(directory_path,"/TABLES/",omik,"factor_", numb,"_", names(enriched)[paths],"_positive.tsv", sep="")  ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
                xx<-plotEnrich(enriched[[paths]], showTerms = n_pathways, numChar = 40, y = "Count", orderBy = "P.value", title = paste("Enrichment_analysis_",omik,"positive_contributors_to_factor",numb,"_",dbs[paths], sep=""))
                print(xx)
                }}}}
        
        if (websiteLive) {
                enriched <- enrichr(LIST_GENE2$features, dbs)}

        for(paths in 1:3){
          if (length(enriched[[paths]]) != 0){
                if (nrow(enriched[[paths]]) != 0){
                        if (websiteLive) {
                                dir.create(paste(directory_path,"/TABLES", sep=""))
                                write.table(x=enriched[[paths]], file=paste(directory_path,"/TABLES/",omik,"factor_", numb,"_",names(enriched)[paths],"_negative.tsv", sep="")  ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
                                xx<-plotEnrich(enriched[[paths]], showTerms = n_pathways, numChar = 40, y = "Count", orderBy = "P.value", title = paste("Enrichment_analysis_",omik,"negative_contributors_to_factor",numb, "_",dbs[paths], sep=""))
                                print(xx)
                        }
                }
          }
        }
        }
        
        
        
        
        if(nam_adv[iter] == "Metabolomics"){
        

                
        pathway_class_HMDB = 
                metpath::pathway_class(hmdb_pathway)
        pathway_class_KEGG = 
                metpath::pathway_class(kegg_hsa_pathway)
        
        gc()
        remain_idx = which(unlist(pathway_class_HMDB) == "Metabolic;primary_pathway")
        hmdb_pathway = hmdb_pathway[remain_idx]
        hmdb_pathway
        
        
        result = 
                enrich_hmdb(query_id = unique(LIST_GENE), 
                            query_type = "compound",
                            id_type = "HMDB",
                            pathway_database = hmdb_pathway, 
                            only_primary_pathway = TRUE, 
                            p_cutoff = padj_pathways,  
                            p_adjust_method = "BH",
                            threads = as.numeric(COMMAND_ADVANCED[3,3]))
        
        if (length(result) != 0){
        result
                
        x<-enrich_bar_plot( object = result, x_axis = "p_value_adjust", cutoff = 1.1, top = n_pathways)
        x <- x + ggtitle(paste("Enrichment_analysis_",omik,"positive_contributors_to_factor_",numb, "_HMDB_results", sep=""))
        print(x)
        
        x <-enrich_scatter_plot(object = result)
        x <- x + ggtitle(paste("Enrichment_analysis_",omik,"positive_contributors_to_factor_",numb, "_HMDB_results", sep=""))
        print(x)
        
        write.table(result@result, paste("MOFA_Factor_Pathways/TABLES/", omik, "factor_", numb, "_positive_HMDB_table_results.tsv", sep =""),quote = FALSE, row.names = F, sep = "\t")
        
        }
        
        head(pathway_class_KEGG)
        remain_idx =
                pathway_class_KEGG %>% unlist() %>% stringr::str_detect("Disease") %>% `!`() %>% which()
        
        remain_idx
        
        pathway_database =
                kegg_hsa_pathway[remain_idx]
        
        pathway_database
        
        result = 
                enrich_kegg(query_id = unique(LIST_GENE), 
                            query_type = "compound", 
                            id_type = "KEGG",
                            pathway_database = pathway_database, 
                            p_cutoff = padj_pathways, 
                            p_adjust_method = "BH", 
                            threads = as.numeric(COMMAND_ADVANCED[3,3]))
        
        if (length(result) != 0){
        result
        
        x <-enrich_bar_plot(
                object = result,
                x_axis = "p_value_adjust",
                cutoff = 1.1,
                top = 10
        )
        x <- x + ggtitle(paste("Enrichment_analysis_",omik,"positive_contributors_to_factor_",numb, "_KEGG_results", sep=""))
        print(x)
        
        x<-enrich_scatter_plot(object = result)
        x <- x + ggtitle(paste("Enrichment_analysis_",omik,"positive_contributors_to_factor_",numb, "_KEGG_results", sep=""))
        print(x)
        
        write.table(result@result, paste("MOFA_Factor_Pathways/TABLES/", omik, "factor_", numb, "_positive_KEGG_table_results.tsv", sep =""),quote = FALSE, row.names = F, sep = "\t")
        
        }
        
        
        
        pathway_class_HMDB = 
                metpath::pathway_class(hmdb_pathway)
        pathway_class_KEGG = 
                metpath::pathway_class(kegg_hsa_pathway)
        
        gc()
        remain_idx = which(unlist(pathway_class_HMDB) == "Metabolic;primary_pathway")
        hmdb_pathway = hmdb_pathway[remain_idx]
        hmdb_pathway
        
        
        result = 
                enrich_hmdb(query_id = unique(LIST_GENE2), 
                            query_type = "compound",
                            id_type = "HMDB",
                            pathway_database = hmdb_pathway, 
                            only_primary_pathway = TRUE, 
                            p_cutoff = padj_pathways,  
                            p_adjust_method = "BH",
                            threads = as.numeric(COMMAND_ADVANCED[3,3]))
        
        if (length(result) != 0){
                result
                
                x<-enrich_bar_plot( object = result, x_axis = "p_value_adjust", cutoff = 1.1, top = n_pathways)
                x <- x + ggtitle(paste("Enrichment_analysis_",omik,"negative_contributors_to_factor_",numb, "_HMDB_results", sep=""))
                print(x)
                
                x <-enrich_scatter_plot(object = result)
                x <- x + ggtitle(paste("Enrichment_analysis_",omik,"negative_contributors_to_factor_",numb, "_HMDB_results", sep=""))
                print(x)
                
                write.table(result@result, paste("MOFA_Factor_Pathways/", omik, "factor_", numb, "_negative_HMDB_table_results.tsv", sep =""),quote = FALSE, row.names = F, sep = "\t")
                
        }
        
        head(pathway_class_KEGG)
        remain_idx =
                pathway_class_KEGG %>% unlist() %>% stringr::str_detect("Disease") %>% `!`() %>% which()
        
        remain_idx
        
        pathway_database =
                kegg_hsa_pathway[remain_idx]
        
        pathway_database
        
        result = 
                enrich_kegg(query_id = unique(LIST_GENE2), 
                            query_type = "compound", 
                            id_type = "KEGG",
                            pathway_database = pathway_database, 
                            p_cutoff = padj_pathways, 
                            p_adjust_method = "BH", 
                            threads = as.numeric(COMMAND_ADVANCED[3,3]))
        
        if (length(result) != 0){
                result
                
                x <-enrich_bar_plot(
                        object = result,
                        x_axis = "p_value_adjust",
                        cutoff = 1.1,
                        top = n_pathways
                )
                x <- x + ggtitle(paste("Enrichment_analysis_",omik,"negative_contributors_to_factor_",numb, "_KEGG_results", sep=""))
                print(x)
                
                x<-enrich_scatter_plot(object = result)
                x <- x + ggtitle(paste("Enrichment_analysis_",omik,"negative_contributors_to_factor_",numb, "_KEGG_results", sep=""))
                print(x)
                
                write.table(result@result, paste("MOFA_Factor_Pathways/", omik, "factor_", numb, "_negative_KEGG_table_results.tsv", sep =""),quote = FALSE, row.names = F, sep = "\t")
                
        }
        
        
        }        
        
        dev.off()        
}

}

}

gc()

