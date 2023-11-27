searchDatabases<-function(disease,path,scenario="Malacards",checkdrug=TRUE,keywordsWikiUser,keywordsKEGGUser,keywordsGOUser,keywordsMSIGUser,keywordsReactUser,keywordsMOAUser){
  library(enrichR)
  library(ReactomeContentService4R)
  library(GO.db)
  if (missing(disease)) cat("Argument disease is missing") else cat(paste("Argument disease =", disease));
  cat ("\n");
  if (missing(path)) cat("Argument path is missing") else cat(paste("Argument path =", path));
  cat("\n\n");
  
  if(scenario !="Hypothesis" && scenario !="Malacards"){
    cat("Argument scenario can only take values 'Malacards' or 'Hypothesis'")
    cat ("\n");
    stop("Execution terminated")
  }else{
    cat(paste("Argument Scenario =", scenario));
    cat ("\n");
  }
  
  if(scenario=="Hypothesis"){
    if (missing(keywordsWikiUser)) {
      cat("If Hypothesis scenario is selected the argument keywordsWikiUser must also be provided")
      cat ("\n");
      stop("Execution terminated")
    }else{
      cat(paste("Argument keywordsWikiUser =", keywordsWikiUser));
      cat ("\n");
    }
    if (missing(keywordsKEGGUser)){
      cat("If Hypothesis scenario is selected the argument keywordsKEGGUser must also be provided")
      cat ("\n");
      stop("Execution terminated")
    }else {
      cat(paste("Argument keywordsKEGGUser =", keywordsKEGGUser));
      cat("\n\n");
    }
    
    if (missing(keywordsGOUser)){
      cat("If Hypothesis scenario is selected the argument keywordsGOUser must also be provided")
      cat ("\n");
      stop("Execution terminated")
    }else {
      cat(paste("Argument keywordsGOUser =", keywordsGOUser));
      cat("\n\n");
    }
    
    if (missing(keywordsMSIGUser)){
      cat("If Hypothesis scenario is selected the argument keywordsMSIGUser must also be provided")
      cat ("\n");
      stop("Execution terminated")
    }else {
      cat(paste("Argument keywordsMSIGUser =", keywordsMSIGUser));
      cat("\n\n");
    }
    
    if (missing(keywordsReactUser)){
      cat("If Hypothesis scenario is selected the argument keywordsReactUser must also be provided")
      cat ("\n");
      stop("Execution terminated")
    }else {
      cat(paste("Argument keywordsReactUser =", keywordsReactUser));
      cat("\n\n");
    }
    
    if (missing(keywordsMOAUser)){
      cat("If Hypothesis scenario is selected the argument keywordsMOAUser must also be provided")
      cat ("\n");
      stop("Execution terminated")
    }else {
      cat(paste("Argument keywordsMOAUser =", keywordsMOAUser));
      cat("\n\n");
    }
  }
  setwd(path)
  if(checkdrug==FALSE){ 
    if(!file.exists("drug_repurposing_hub.txt")){
      cat("For the option of checkdrug==FALSE, you need to supply a drug repusporing database file, you can download one with the scRANK test data here: https://bioinformatics.cing.ac.cy/downloads/scRNA/LAM.tar.gz")
      cat ("\n");
      stop("Execution terminated")
    }else{
      cat(paste("Argument checkdrug =", checkdrug));
      cat("\n\n");
    }
  }
  xx <- as.list(GOTERM)
  
  if(scenario=="Malacards"){
    #Decided to include all the pathways for KEGG to increase chances of finding KEGG and MSIG pathways
    KEGG<-list.files(pattern = paste(disease,"PathwaysKEGG.txt$",sep=""),)
    GOs<-list.files(pattern = paste(disease,"PathwaysGO.txt$",sep=""))
    React<-list.files(pattern = paste(disease,"PathwaysReactome.txt$",sep=""))
    Wiki<-list.files(pattern = paste(disease,"PathwaysWiki.txt$",sep=""))
    filenamesSD<-c(KEGG,GOs,React,Wiki)
    for(file in filenamesSD){
      if(file.size(file) == 0L){
        print(paste(file," is empty, either insert keywords in the file or search Malacards using additonal terms.",sep=""))
      }
    }
    
    keywordsWiki<-read.table(Wiki, header=F, sep="\t")
    keywordsWiki<-as.array(keywordsWiki[,1])
    keywordsOther<-read.table(KEGG, header=F, sep="\t")
    keywordsOther<-as.array(keywordsOther[,1])
    keywordsGO<-read.table(GOs, header=F, sep="\t")
    keywordsGO<-as.array(keywordsGO[,1])
    keywordsMSIG<-keywordsOther
    keywordsKEGG<-keywordsOther
    termsWiki<-keywordsWiki
    termsWiki<-unique(tolower(termsWiki))
    keywordsReact<-c("")
    tryCatch(
      #try to do this
      {
        keywordsReact<-read.table(React, header=F, sep="\t")
        keywordsReact<-as.array(keywordsReact[,1])
      },
      #if an error occurs, print the error
      error=function(e) {
        message('An Error Occurred... empty file')
        print(e)
      },
      #if a warning occurs, print the warning
      warning=function(w) {
        message('A Warning Occurred')
        print(w)
        return(NA)
      }
    )
    
  }else{
    keywordsWiki<-keywordsWikiUser
    keywordsKEGG<-keywordsKEGGUser
    keywordsGO<-keywordsGOUser
    keywordsMSIG<-keywordsMSIGUser
    keywordsReact<-keywordsReactUser
    
    termsWiki<-c()
    indexWiki<-1
    for(WI in keywordsWiki){
      print(WI)
      termsWikiFirstSearch<-findPathwayNamesByText(WI)
      
      if(indexWiki !=1){
        foundinsearch<-grepl(WI, termsWikiFirstSearch,ignore.case = T)#not sure here double check changed from previous runs
        termsWiki<-c(termsWiki,termsWikiFirstSearch[foundinsearch])
      }else{
        termsWiki<-c(termsWiki,termsWikiFirstSearch)
      }
      indexWiki<-indexWiki+1
    }
    termsWiki<-termsWiki[!is.na(termsWiki)]
    termsWiki<-unique(tolower(termsWiki))
  }
  
  
  if(scenario=="Hypothesis"){
    termsGO<-c()
    for(KI in 1:length(keywordsGO)){
      indexsGO<-grep(keywordsGO[KI], xx,ignore.case = T)
      if(length(indexsGO)!=0){
        for(GOI in 1:length(indexsGO)){
          print(Term(xx[[indexsGO[GOI]]]))
          print(Ontology(xx[[indexsGO[GOI]]]))
          print(GOID(xx[[indexsGO[GOI]]]))
          if(Ontology(xx[[indexsGO[GOI]]])=="BP"){
            termsGO<-c(termsGO,Term(xx[[indexsGO[GOI]]]))
          }
        }
      }
    }
    termsGO<-unique(termsGO)
  }else{
    termsGO<-keywordsGO
  }
  termsKEGG<-c()
  for(KI in 1:length(keywordsKEGG)){
    tryCatch(
      #try to do this
      {
        termsKEGG<-c(termsKEGG,as.character(keggFind("pathway", c(keywordsKEGG[KI])))[1])
        termsKEGG<-c(termsKEGG,as.character(keggFind("disease", c(keywordsKEGG[KI])))[1])
      },
      #if an error occurs, print the error
      error=function(e) {
        message('Warning... term was not found in KEGG')
        print(e)
      },
      #if a warning occurs, print the warning
      warning=function(w) {
        message('A Warning Occurred')
        print(w)
        return(NA)
      }
    )
  }
  
  termsKEGG<-termsKEGG[!is.na(termsKEGG)]
  termsKEGG<-unique(termsKEGG)
  
  misig<-getMsigdb( org = c("hs"), id = c("SYM"), version = getMsigdbVersions() )
  findMSIG<-function(misigdb,collection){
    #retrieeve the hallmarks gene sets
    if(collection==1){
      hallmark<-subsetCollection(misig, 'h')
    }else{
      hallmark<-subsetCollection(misig, 'c2','CP')
    }
    termsMSIG<-c()
    
    for(KI in 1:length(keywordsMSIG)){
      #old approach
      #indexsMSIG<-grep(keywordsMSIG[KI],hallmark,ignore.case = T)
      #print(keywordsMSIG[KI])
      tokens<-strsplit(keywordsMSIG[KI]," ")
      counttokmatch<-0
      for(h in 1:length(hallmark)){
        bool<-c()
        boolDes<-c()
        for(toks in 1:length(tokens[[1]])){
          termfromtoks<-gsub("[\\(\\)]", "", tokens[[1]][toks])
          bool<-c(bool,grepl(termfromtoks,hallmark[[h]]@setName,ignore.case = T))
          boolDes<-c(boolDes,grepl(termfromtoks,hallmark[[h]]@shortDescription,ignore.case = T))
        }
        #print(bool)
        bool<-unique(bool)
        boolDes<-unique(bool)
        if(length(bool)==1){
          if(bool){
            print(keywordsMSIG[KI])
            print(hallmark[[h]]@setName)
            MSIGPath<-setName(hallmark[[h]])
            MSIGPath<-gsub("HALLMARK_","",MSIGPath)
            termsMSIG<-c(termsMSIG,MSIGPath)
          }
        }
        if(length(boolDes)==1){
          if(boolDes){
            print(keywordsMSIG[KI])
            print(hallmark[[h]]@setName)
            MSIGPath<-setName(hallmark[[h]])
            MSIGPath<-gsub("HALLMARK_","",MSIGPath)
            termsMSIG<-c(termsMSIG,MSIGPath)
            
          }
        }
      }
    }
    return(termsMSIG)
  }
  # 
  #old approach
  # if(length(indexsMSIG)!=0){
  #     print(keywordsMSIG[KI])
  #     for(MSIGI in 1:length(indexsMSIG)){
  #     print(setName(hallmark[[indexsMSIG[MSIGI]]]))
  #     print(hallmark[[indexsMSIG[MSIGI]]]@shortDescription)
  #     MSIGPath<-setName(hallmark[[indexsMSIG[MSIGI]]])
  #     MSIGPath<-gsub("HALLMARK_","",MSIGPath)
  #     termsMSIG<-c(termsMSIG,MSIGPath)
  #     }
  #   }
  # }
  
  termsMSIG<-findMSIG(misig,1)
  termsMSIG<-unique(termsMSIG)
  # if(is.null(termsMSIG )){
  #   termsMSIG<-findMSIG(misig,2)
  #   termsMSIG<-unique(termsMSIG)
  # }
  
  termsReact<-c()
  if(keywordsReact[1]!=""){
    for(KI in 1:length(keywordsReact)){
      bdd.search <- searchQuery(query = keywordsReact[KI], species = "human",types = "Pathway")
      bdd.search$results$entries[[1]]$name[1]
      id<-query(id = paste("R-HSA-",bdd.search$results$entries[[1]]$dbId[1],sep=""))
      print(id$displayName)
      termsReact<-c(termsReact,id$displayName)
    }
  }
  
  if(checkdrug==TRUE){
    Drugs<-list.files(pattern = paste(disease,"DrugsSorted.txt$",sep=""))
    if(file.size(Drugs) == 0L){
      print(paste(Drugs," is empty, either insert keywords in the file or search Malacards using additonal terms.",sep=""))
    }
    keywordsMOA<-read.table(Drugs, header=F, sep="\t")
    keywordsMOA<-as.array(keywordsMOA[,1])
    termsMOA<-keywordsMOA
  }else{
    drugInfo<-read.delim("drug_repurposing_hub.txt")
    keywordsMOA<-keywordsMOAUser
    
    termsMOA<-c()
    for(KI in 1:length(keywordsMOA)){
      indexsMOA<-grep(keywordsMOA[KI], drugInfo$moa,ignore.case = T)
      if(length(indexsMOA)!=0){
        for(MMOA in indexsMOA){
          if(grepl(keywordsMOA[KI], drugInfo$moa[MMOA],ignore.case = T)){
            termsMOA<-c(termsMOA,drugInfo$moa[MMOA])
          }
        }
      }
    }
    termsMOA<-unique(termsMOA)
  }
  listofoutput<-list(termsKEGG,termsGO,termsMSIG,termsWiki,termsReact,termsMOA)
  lens<-as.data.frame(c(length(termsKEGG),length(termsGO),length(termsMSIG),length(termsWiki),length(termsReact),length(termsMOA)))
  rownames(lens)<-c("KEGG","GO","MSIG","WIKI","REACT","DRUGS")
  
  lens<-cbind(lens,rownames(lens))
  colnames(lens)<-c("Number","Source")
  source<-rownames(lens)
  p<-ggplot(data = lens,aes(x =  reorder(source,Number,FUN = max,decreasing = FALSE), y=Number)) +
    
    geom_bar(stat="identity",position=position_dodge())+#,
    
    theme(
      plot.title = element_text(size=11,)
    ) +
    ggtitle("Matched prior knowledge") +#cellID_i
    xlab("")
  p<-p+theme(plot.title = element_text(size = 15, face = "bold"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=15),
             axis.title=element_text(size=15,face="bold"),legend.text=element_text(size=15),legend.title=element_text(size=15))
  p<-p + coord_flip()
  print(p)
  return(listofoutput)
}
