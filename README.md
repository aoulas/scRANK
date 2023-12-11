# scRANK - Single Cell Ranking Analysis Toolkit

![Alt text](https://github.com/aoulas/scRANK/blob/main/scRANKLogo4.png "scRANK")

## Update
Jun 27, 2023 (Version 1.0.0)

scRANK is now available via devtools installation. 

## Capabilities
scRANK utilizes prior knowledge in combination with expert-user information to guide the choice of cell types from a scRNA-seq analysis that yield the most biologically meaningful results. Prior knowledge is incorporated in a standardized, structured manner, whereby a checklist is attained by querying MalaCards human disease database with a disease of interest. The checklist is comprised of pathways and drugs and optionally drug mode of actions (MOAs), associated with the disease. The user is prompted to “edit” this checklist by removing or adding terms (in the form of keywords) from the list of predefined terms. The user may also define de-novo, a set of keywords that best suit a hypothesis the user is interested in investigating (hypothesis-driven approach). Once the checklist is finalized a “mapping” step is performed. This is done initially against pathway enrichment results attained from analysing the scRNA-seq data. In addition, the user-selected checklist of drug names is mapped against DR result derived from the data. The analysis of scRNA-seq data has the capability to obtain pathways and repurposed drugs by comparing disease-control conditions for every cell type in the analysis. scRANK then uses the user defined information to highlight the cell types that generate the results that best “map” to the predefined prior knowledge provided by the expert user. The methodology is fully automated and a ranking is generated for all cell types in the analysis allowing the user to pinpoint the specific cells that are most prominently affected by the disease under study in accordance to the provided prior knowledge.  The output of our methodology provides an automated validation of the result and further makes it easier for researchers to interpret their findings. In addition, the results provide greater credence to de novo information as it also backed-up by prior knowledge for the disease under study. Novel information is obtained in the form of predicted pathways that are enriched for each cell type, as well as previously unreported repurposed drugs that are obtained from the differentially expressed genes (DEGs) between disease-control conditions.   

## Installation
scRANK R package can be easily installed from Github using devtools:

#install.packages("devtools")\
devtools::install_github("aoulas/scRANK")

Please make sure you have installed all the dependencies. See instructions below.

## Installation of dependencies
### CRAN packages 
Seurat, dplyr, patchwork, metap, ggplot2, cowplot, enrichR, gridExtra, ggpubr, RColorBrewer, crank, riverplot, rvest, stringr.

These will be installed automatically togther with scRANK.\
With the exception of "riverplot" which is no longer available via cran and needs to be downloaded from archives using the link found [here](https://cran.r-project.org/src/contrib/Archive/riverplot/riverplot_0.10.tar.gz)
and then installed using:\
install.packages(path_to_file, repos = NULL, type="source")

### Bioconductor packages
KEGGREST, GO.db, rWikiPathways, ReactomeContentService4R, multtest, msigdb.
### GitHub packages
SeuratDisk:\
remotes::install_github("mojaveazure/seurat-disk")\
BiocFileCache version 2.11.1 (or above):\
devtools::install_github("Bioconductor/BiocFileCache")\
CellChat:\
devtools::install_github("sqjin/CellChat")

## Tutorial
### Downlaod test data 
A test datsets from lymphangioleiomyomatosis (LAM) disease and control (Donor) samples is available for download [here](https://bioinformatics.cing.ac.cy/downloads/scRNA/LAM.tar.gz). Extract data in a local directory.

### Extract prior knowledge from MalaCards database
Paste the following link in your borwser.
```
https://www.malacards.org/card/lymphangioleiomyomatosis?showAll=TRUE
```
For a differnet disease you can change the name of the disease in the url above. Make sure the disease exists with the same name in the database. Ensure to add the 'showAll=TRUE' flag to expand all tables in the web page (this may take some time to load). Once the page has loaded right-click and click save-as to download the html content to the same directory as the test data downloaded above.

### Run scRANK
```
library(scRANK)
#Extract the relevant files from the MalaCards .html file downloaded in the above step
#Define common arguments for scRANK extractMalacards() and runBasicAnalysis() functions
#(note* the path and disease name has
#to be common for all scRANK funnctions to work properly)
path<-"full-path-to-where-data-was-extracted"
disease<-"LAM"
extractMalacards(disease = disease,files = c("name-of-html-file"),path = path)

#Define extra arguments for scRANK functions
annotate<-TRUE
userlabel<-"label"
usercelltype<-"celltype"
checkdrug<-TRUE
scenario<-"Malacards"
scan<-"Cell"


if(scenario =="Hypothesis"){
  keywordsWiki<-c("MTOR","PI3K","MAPK","apoptosis","NF-k","TNF")
  keywordsKEGG<-c("MTOR","PI3K","MAPK","apoptosis","NF-k","TNF")
  keywordsGO<-c("MTOR","PI3K","MAPK","apoptosis","NF-k","TNF")
  keywordsMSIG<-c("MTOR","PI3K","MAPK","apoptosis","NF-k","TNF")
  keywordsReact<-c("MTOR signalling","PI3K","MAPK","apoptosis","NF-k","TNF")
  keywordsMOA<-c("CDK inhibitor","MTOR inhibitor","MEK inhibitor")
}

#Search databases with the terms extracted from Malacards (checks also that files generated from extractMalacards()
#are not empty) or search databases using Hypothesis-driven keywords\
if(scenario =="Hypothesis"){
  listofoutput<-searchDatabases(disease = disease,path=path,scenario=scenario,checkdrug=checkdrug,keywordsWikiUser = keywordsWiki, keywordsKEGGUser =keywordsKEGG,keywordsGOUser =keywordsGO, keywordsMSIGUser = keywordsMSIG,keywordsReactUser = keywordsReact, keywordsMOAUser = keywordsMOA)
}else{
  listofoutput<-searchDatabases(disease = disease,path=path,scenario=scenario,checkdrug=checkdrug)
}

#Run basic analysis 
seuratObject<-runBasicAnalysis(disease = disease,path=path ,annotate = annotate,userlabel = userlabel,
          usercelltype = usercelltype)

#Define extra arguments for scRANK rankCells() function
priorknowledgePathsKEGG<-listofoutput[[1]]
priorknowledgePathsGO<-listofoutput[[2]]
priorknowledgePathsMSIG<-listofoutput[[3]]
priorknowledgePathsWiki<-listofoutput[[4]]
priorknowledgePathsReact<-listofoutput[[5]]
priorknowledgeDRUGSMOA<-listofoutput[[6]]

#Perform mapping and ranking steps - you can use the output from the runBasicAnalysis()
#directly in the rankCells() function.
listofscRANKs<-rankCells(seuratObject,path,scan=scan,priorknowledgePathsKEGG,priorknowledgePathsGO,priorknowledgePathsMSIG,
priorknowledgePathsWiki,priorknowledgePathsReact,priorknowledgeDRUGSMOA,userlabel,usercelltype,checkdrug,scenario=scenario)

#Run CellChat - note the first label is considered as the reference (e.g., control)
foldchangeInterMat<-runCellChat(seuratObject,userlabel,usercelltype)

#Peform basic plots
plotRanks("filename-of-Ranking-results")
plotCellChat("filename-of-CellChat-results")
plotTotalNumberDEGs("filename-of-Total-DEG-results")
plotProportions("filename-of-Cell-Proportion-results")
```




