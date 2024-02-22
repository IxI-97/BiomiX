# Input list with pathway names
input_list <- read.table("/home/cristia/Scrivania/BiomiX1.6/MOFA/x_BiomiX_DATABASE/Factor_names_selection/Raw_Human_Phenotype_Ontology.txt")
input_list <-input_list[1:6656,1]
# Combine the split lines into single pathway names
combined_list <- character()
i <- 1

while (i <= length(input_list)) {
        if (grepl("^HP_", input_list[i])) {
                # Start of a new pathway name
                pathway_name <- input_list[i]
                i <- i + 1
                
                while (i <= length(input_list) && !grepl("^HP_", input_list[i])) {
                        # Combine the split lines
                        pathway_name <- paste0(pathway_name, input_list[i])
                        i <- i + 1
                }
                
                combined_list <- c(combined_list, pathway_name)
        }
}


write.table(combined_list,file="/home/cristia/Scrivania/BiomiX1.6/MOFA/x_BiomiX_DATABASE/Factor_names_selection/Human_Phenotype_Ontology.txt",append=TRUE, quote = FALSE, row.names = FALSE,col.names = FALSE, sep="\t")     



input_list1 <- read.table("/home/cristia/Scrivania/BiomiX1.6/MOFA/x_BiomiX_DATABASE/Factor_names_selection/Reactome_pathways.txt")
input_list2 <- read.table("/home/cristia/Scrivania/BiomiX1.6/MOFA/x_BiomiX_DATABASE/Factor_names_selection/Human_Phenotype_Ontology.txt")


listing1<-strsplit(input_list1$V1, "_")
listing2<-strsplit(input_list2$V1, "_")

cd1<- unique(unlist(listing1))
cd1 <- cd1[nchar(cd1) > 2]
cd2<- unique(unlist(listing2))
cd2 <- cd2[nchar(cd2) > 2]

cd<- unique(c(cd1,cd2))
write.table(cd,file="/home/cristia/Scrivania/BiomiX1.6/MOFA/x_BiomiX_DATABASE/Factor_names_selection/REA_HPO_vocabolary.txt",append=TRUE, quote = FALSE, row.names = FALSE,col.names = FALSE, sep="\n")     
input_list1 <- read.table("/home/cristia/Scrivania/BiomiX1.6/MOFA/x_BiomiX_DATABASE/Factor_names_selection/REA_HPO_vocabolary.txt")

remove_HP <- function(word) {
        if (substr(word, nchar(word) - 1, nchar(word)) == "HP") {
                return(substr(word, 1, nchar(word) - 2))
        } else {
                return(word)
        }
}

# Apply the function to each element in the vector
processed_vector <- sapply(input_list1$V1, remove_HP)
processed_vector <- unique(processed_vector)
write.table(processed_vector,file="/home/cristia/Scrivania/BiomiX1.6/MOFA/x_BiomiX_DATABASE/Factor_names_selection/REA_HPO_vocabolary_clean.txt",append=TRUE, quote = FALSE, row.names = FALSE,col.names = FALSE, sep="\n")     

processed_vector2 <- sapply(processed_vector, remove_HP)
processed_vector2 <- unique(processed_vector2)
write.table(processed_vector2,file="/home/cristia/Scrivania/BiomiX1.6/MOFA/x_BiomiX_DATABASE/Factor_names_selection/REA_HPO_vocabolary_clean2.txt",append=TRUE, quote = FALSE, row.names = FALSE,col.names = FALSE, sep="\n")     
