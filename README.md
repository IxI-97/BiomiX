## **BiomiX webpage:** https://ixi-97.github.io/
## **BiomiX package/source:** https://github.com/IxI-97/BiomiX
## **News**

### *20 December 2024*

Hi everyone!

Our work in BMC Bioinformatics is coming soon. After the peer review, we implemented many aspects of BiomiX, introducing version 2.4 here.

_**new shiny interfaces**_

_Preview-QC interface:_ The first interface, named preview-QC, opens in a browser and provides an overview of the data through scatterplots (for two selected variables), PCA, correlation heatmaps, and UMAP. You can apply normalization and transformation methods directly in the browser, visualizing changes and proceed once they are satisfied. In the manuscript, we refer to these methods as 'transformation' due to the varying definitions of normalization across different omics. This interface also includes a tab for removing highly variable features that are likely driven by noise for eliminating outliers based on their distances from PCA centroids. Additionally, users also have the option to manually select and remove specific samples.

_Imputation interface:_ The second Shiny interface included in the BiomiX toolkit allows users to select the type of imputation method (NIPALS, random forest, lasso, mean/median and 0 values) along with any potential bias introduced into the data. The BiomiX version 2.4 version can handle both imputed and unimputed data, as MOFA can manage missing data during the integration.

_**Installation process execution**_ 
- Many MacOS installation bugs were fixed and the installation tutorial was updated: [Installation Page](https://ixi-97.github.io/Installation.html#MAC-OS) 
- The BiomiX execution file BiomiX.sh was splitted into Biomix_linux_exe.sh and Biomix_macos_exe.sh.

_**main and advanced option interface**_ 
- Feature renaming in the metabolomics section to clarify the input technology. 
- Each advance option and main interface section has included many "help" sections to guide the user to a wise usage of the parameters. The novel in depth webpage is accessible at this link: [Advanced Parameters Page](https://ixi-97.github.io/Parameters_guidelines.html)  

_**Single omics pipeline**_ 
- Unlock the .mgf file usage for the metabolomics annotation. As they are lighter than .mzML files this optimize the RAM memory usage in the inclusion of high number of files. 
- Although BiomiX does not support instrument normalization for metabolomics, users can still preview QC samples in the Preview-QC interface by labelling the “QC” in the CONDITION column and the loading order by adding order columns in the metadata. QC samples and their loading can in this way, be displayed in PCA and UMAP before being removed automatically afterwards.
- We added a section in the omics area called **“Undefined**”. This new “omics” allows users to perform general QC, transformation, and statistical analysis on any type of data, including omics without a specialized pipeline. Over time, we plan to incorporate specific pipelines based on state-of-the-art methods commonly used for each type of omics.

_**MOFA factor intepretation**_ 
- Bug fixed: the bibliography and correlation analysis are now properly recognizing the undefined omics contributors to correlate and interpret them. 



### *21 August 2024*

Hi everyone!

Based on community feedback, we are excited to announce the release of version 2.3.

This version introduces several **new features**:

- Excel File Format Support: Now handle Excel files seamlessly within the application.

- Undetermined Omics Support: You can now analyze any type of biological data not previously included in the listed omics. This feature performs t-tests and Wilcoxon tests between two groups, and prepares the data for MOFA integration. Integrating any type of data into your multi-omics analysis is now easier than ever!

**Bug Fixes**:

- The BiomiX interface now adapts smoothly to all screen sizes, even after resizing or expanding.

- Package conflicts on Linux and macOS have been resolved.



### *18 June 2024* 
Hi to everyone! 

We have great news.

The BiomiX project has finally started, and our Pre-print is out in Biorxiv (doi: https://doi.org/10.1101/2024.06.14.599059). Please let us know your feedback! and our satisfaction is to know, 
BiomiX can ease your life and your research. Our first Biomix version is the version 2.2.

We are always looking for new members of our community. Feel free to join us by contacting cristian.iperi@univ-brest.fr!

 
 <div align="center">
    <img src="https://github.com/IxI-97/IxI-97.github.io/blob/main/BiomiX_logo3.png?raw=true" width="300" height="300">
</div>



