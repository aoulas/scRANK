# CellRank
## Update
Jun 27, 2023 (Version 1.0.0)

CellRank is now available via devtools installation. 

## Capabilities
CellRank utilizes prior knowledge in combination with expert-user information to guide the choice of cell types from a scRNA-seq analysis that yield the most biologically meaningful results. Prior knowledge is incorporated in a standardized, structured manner, whereby a checklist is attained by querying MalaCards human disease database with a disease of interest. The checklist is comprised of pathways and drugs and optionally drug mode of actions (MOAs), associated with the disease. The user is prompted to “edit” this checklist by removing or adding terms (in the form of keywords) from the list of predefined terms. The user may also define de-novo, a set of keywords that best suit a hypothesis the user is interested in investigating (hypothesis-driven approach). Once the checklist is finalized a “mapping” step is performed. This is done initially against pathway enrichment results attained from analysing the scRNA-seq data. In addition, the user-selected checklist of drug names is mapped against DR result derived from the data. The analysis of scRNA-seq data has the capability to obtain pathways and repurposed drugs by comparing disease-control conditions for every cell type in the analysis. CellRank then uses the user defined information to highlight the cell types that generate the results that best “map” to the predefined prior knowledge provided by the expert user. The methodology is fully automated and a ranking is generated for all cell types in the analysis allowing the user to pinpoint the specific cells that are most prominently affected by the disease under study in accordance to the provided prior knowledge.  The output of our methodology provides an automated validation of the result and further makes it easier for researchers to interpret their findings. In addition, the results provide greater credence to de novo information as it also backed-up by prior knowledge for the disease under study. Novel information is obtained in the form of predicted pathways that are enriched for each cell type, as well as previously unreported repurposed drugs that are obtained from the differentially expressed genes (DEGs) between disease-control conditions.   

## Installation
CellRank R package can be easily installed from Github using devtools:

#install.packages("devtools")\
devtools::install_github("aoulas/CellRank")

Please make sure you have installed all the dependencies. See instructions below.

## Installation of dependencies
### CRAN packages
Seurat, dplyr, patchwork, metap, ggplot2, cowplot, enrichR, gridExtra, ggpubr, RColorBrewer, crank, riverplot, rvest, stringr.
### Bioconductor packages
KEGGREST, GO.db, rWikiPathways, ReactomeContentService4R, multtest, msigdb.
### GitHub packages
SeuratDisk: remotes::install_github("mojaveazure/seurat-disk")\
CellChat: devtools::install_github("sqjin/CellChat")

## Tutorial
### Downlaod test data 
A test datsets from lymphangioleiomyomatosis (LAM) disease and control (Donor) samples is available for download [here](https://bioinformatics.cing.ac.cy/downloads/scRNA/LAM.tar.gz). Extract data in a local directory.

### Extract prior knowledge from MalaCards database
Paste the following link in your borwser.
```
https://www.malacards.org/card/lymphangioleiomyomatosis?showAll=TRUE
```
For a differnet disease you can change the name of the disease in the url above. Make sure the disease exists with the same name in the database. Ensure to add the 'showAll=TRUE' flag to expand all tables in the web page (this may take some time to load). Once the page has loaded right-click and click save-as to download the html content to the same directory as the test data downloaded above.

### Run CellRank
```
library(CellRank)
#Extract the relevant files from the MalaCards .html file downloaded in the above step
#Define common arguments for CellRank extractMalacards() and runBasicAnalysis() functions
#(note* the path and disease name has
#to be common for all CellRank funnctions to work properly)
path<-"full-path-to-where-data-was-extracted"
disease<-"LAM"
extractMalacards(disease = disease,files = c("name-of-html-file"),path = path)

#Define extra arguments for CellRank runBasicAnalysis() function
annotate<-TRUE
userlabel<-"label"
usercelltype<-"celltype"
checkdrug<-TRUE
scenario<-"Malacards"
scan<-"Cell"

#Run basic analysis and search databases
listofoutput<-runBasicAnalysis(disease = disease,path=path ,annotate = annotate,userlabel = userlabel,
          usercelltype = usercelltype,scenario=scenario)

#Define extra arguments for CellRank rankCells() function
seuratObject<-listofoutput[[1]]
priorknowledgePathsKEGG<-listofoutput[[2]]
priorknowledgePathsGO<-listofoutput[[3]]
priorknowledgePathsMSIG<-listofoutput[[4]]
priorknowledgePathsWiki<-listofoutput[[5]]
priorknowledgePathsReact<-listofoutput[[6]]
priorknowledgeDRUGSMOA<-listofoutput[[7]]

#Perform mapping and ranking steps - you can use the output from the runBasicAnalysis()
#directly in the rankCells() function.
listofCellRanks<-rankCells(seuratObject,path,scan=scan,priorknowledgePathsKEGG,priorknowledgePathsGO,priorknowledgePathsMSIG,
priorknowledgePathsWiki,priorknowledgePathsReact,priorknowledgeDRUGSMOA,userlabel,usercelltype,checkdrug,scenario=scenario)

#Run CellChat - note the first label is considered as the reference (control)
foldchangeInterMat<-runCellChat(seuratObject,userlabel,usercelltype)

#Peform basic plots
plotRanks("filename-of-Ranking-results")
plotCellChat("filename-of-CellChat-results")
plotTotalNumberDEGs("filename-of-Total-DEG-results")
plotProportions("filename-of-Cell-Proportion-results")
```




