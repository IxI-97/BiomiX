
print("looking for the BiomiX_folder")

find_folder <- function(folder_name, depth = 5) {
        # Get the root directory
        root <- normalizePath("/")
        
        # Recursive function to search for the folder within the specified depth
        search_folder <- function(directory, current_depth) {
                if (current_depth <= depth) {
                        files <- list.files(directory, full.names = TRUE)
                        
                        matching_folders <- files[file.info(files)$isdir & basename(files) == folder_name]
                        
                        if (length(matching_folders) > 0) {
                                return(matching_folders[1])
                        } else {
                                subdirectories <- files[file.info(files)$isdir]
                                for (subdir in subdirectories) {
                                        result <- search_folder(subdir, current_depth + 1)
                                        if (!is.null(result)) {
                                                return(result)
                                        }
                                }
                        }
                }
                return(NULL)
        }
        
        # Call the recursive function starting from the root directory
        found_folder <- search_folder(root, 1)
        
        # Return the path to the matching folder (if found)
        if (!is.null(found_folder)) {
                return(found_folder)
        } else {
                return(NULL)
        }
}

setwd(find_folder("BiomiX2.2"))


print("R Package checking..")

chooseCRANmirror(48, ind = TRUE)

library("htmltools")
packageVersion("htmltools")

install.packages("devtools", version="2.4.5")

#NEW BLOCK
#packageurl <- "https://cran.r-project.org/src/contrib/htmltools_0.5.7.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")

packageurl <- "https://cran.r-project.org/src/contrib/Archive/bslib/bslib_0.4.0.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

packageurl <- "https://cran.r-project.org/src/contrib/Archive/shiny/shiny_1.7.4.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

packageurl <- "https://cran.r-project.org/src/contrib/Archive/rmarkdown/rmarkdown_2.17.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

packageurl <- "https://cran.r-project.org/src/contrib/Archive/htmlwidgets/htmlwidgets_1.6.0.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

install.packages("devtools", version="2.4.5")
install.packages("vroom", version="1.6.1")
install.packages("lava", version="1.6.10")
install.packages("recipes", version="1.0.7")
install.packages("dplyr", version="1.1.2")
install.packages("future.apply", version="1.11.0")
install.packages("stringr", version="1.5.0")
install.packages("circlize", version="0.4.15")
install.packages("ggplot2", version="3.4.3")
install.packages("ggrepel", version="0.9.3")
install.packages("enrichR", version="3.2")
install.packages("rlist", version="0.4.6.2")
install.packages("tidyverse", version="2.0.0")
install.packages("data.table", version="1.14.8")
install.packages("caret", version="6.0.94")
install.packages("reticulate", version="1.28")
install.packages("XML", version="3.99.0.14")
install.packages("xml2", version="1.3.3")
install.packages("rentrez", version="1.2.3")
install.packages("remotes", version="2.4.2.1")
install.packages("GGally", version="2.1.2")
install.packages("ggpubr", version="0.6.0")
install.packages("BiocManager", version="1.30.22")


print("Installing Bioconductor packages")

BiocManager::install("DESeq2", version ="3.14")
BiocManager::install("ComplexHeatmap", version ="3.14")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", version ="3.14")
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", version ="3.14")
BiocManager::install("ChAMP", version ="3.14")
BiocManager::install("MOFA2", version ="3.14")
BiocManager::install("MSnbase", version ="3.14")
BiocManager::install("Rdisop", version ="3.14")

library(devtools)
library(remotes)
library(igraph)


print("Installing github packages")
devtools::install_github("lzyacht/cmmr")
remotes::install_gitlab("jaspershen/masstools", dependencies = TRUE)
remotes::install_gitlab("tidymass/metpath", dependencies = TRUE)
remotes::install_github("elizagrames/litsearchr")
devtools::install_github("strengejacke/sjmisc")

#BiocManager::install("ChAMP", version ="3.14")
