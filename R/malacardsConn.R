extractMalacards<-function (disease,files,path){
  
  setwd(path)
  
  
  PathwaysAllNew<-c()
  GOAllNew<-c()
  DrugsAllNew<-c()
  PathwaysAllNewWiki<-c()
  PathwaysAllNewReact<-c()
  PublicationsAllNew<-c()
  MOAsNew<-c()
  
  for(i in 1:length(files)){
    csv<-read.csv2(files[i]) 
    if(grepl("Pathways",files[i])){
      for(j in 1:nrow(csv)){
        superpathway<-str_split_1(csv[j,],",")[3]
        if(!is.na(superpathway) && !superpathway=="SuperPathway"){
          superpathwayclean<-gsub(" [0-9]+$","",superpathway)
          superpathwayclean<-gsub(" ","_",superpathwayclean)
          superpathwayclean<-tolower(superpathwayclean)
          
          xhtml <- tryCatch({
            read_html(paste("https://pathcards.genecards.org/card/",superpathwayclean,sep=""))
          }, error = function(e) {
            message("Error: ", e)
            NULL
          })
          
          if (!is.null(xhtml)) {
            # tbls_ls <- xhtml %>%
            #   html_nodes(".table")%>%
            #   html_table(fill = TRUE)
            
            pathway_source <- xhtml %>%
              html_nodes(".pathway-source")%>%
              html_nodes("img") %>%
              html_attr("title")
            
            subpathways <- xhtml %>%
              html_nodes(".pathway-link")%>%
              html_nodes("a") %>%
              html_text()
            
            
            # subpathways<-as.matrix(tbls_ls[[1]])
            # subpathways<-as.vector(subpathways)
            subpathways<-subpathways[which(!is.na(subpathways))]
            
            
            pathway_source[grep("Wiki",pathway_source)]<-74
            pathway_source[grep("Reactome",pathway_source)]<-66
            pathway_source<-as.numeric(pathway_source)
            pathway_source[which(is.na(pathway_source))]<-1000
            subpathways<-paste(subpathways,pathway_source,sep=" ")
            
            allpaths<-c(superpathway,subpathways)
            
            
            
            
            for(ap in 1:length(allpaths)){
              indexesReactome<-grep("66",allpaths[ap])
              indexesWiki<-grep("74",allpaths[ap])
              indexesother<-grep("74",allpaths[ap],invert = TRUE) && grep("66",allpaths[ap],invert = TRUE)
              
              if(length(indexesReactome) !=0){
                PathwaysAllNewReact<-append(PathwaysAllNewReact, allpaths[ap])
              }
              if(length(indexesWiki) !=0){
                PathwaysAllNewWiki<-append(PathwaysAllNewWiki, allpaths[ap])
              }
              
              if(!is.na(indexesother)){ 
                if(indexesother==TRUE){
                  PathwaysAllNew<-append(PathwaysAllNew, allpaths[ap])
                }
              }
            }
          }
        }
      }
    }else if(grepl("Biological",files[i])){
      for(j in 1:nrow(csv)){
        GO<-str_split_1(csv[j,],",")[2]
        if(!is.na(GO) && !GO=="Name"){
          #print(GO)
          GOAllNew<-append(GOAllNew, GO)
        }
      }
    }else if(grepl("Drugs",files[i])){
      for(j in 1:nrow(csv)){
        DRUGS<-str_split_1(csv[j,],",")[3]
        Status<-str_split_1(csv[j,],",")[4]
        Phase<-str_split_1(csv[j,],",")[5]
        
        if(!is.na(DRUGS) && !DRUGS=="Name"){
          DrugsAllNew<-rbind(DrugsAllNew,cbind(DRUGS,Status,Phase))
        }
      }
    }
    else if(grepl("Publications",files[i])){
      for(j in 1:nrow(csv)){
        Title<-str_split_1(csv[j,],",")[3]
        Year<-str_split_1(csv[j,],",")[6]
        if(!is.na(Title) && !Title=="Title"){
          #print(GO)
          PublicationsAllNew<-rbind(PublicationsAllNew,cbind(Title,Year))
        }
      }
    }
    else if(grepl("Text Mined",files[i])){
      for(j in 1:nrow(csv)){
        Title<-str_split_1(csv[j,],",")[3]
        Year<-str_split_1(csv[j,],",")[6]
        if(!is.na(Title) && !Title=="Title"){
          #print(GO)
          PublicationsAllNew<-rbind(PublicationsAllNew,cbind(Title,Year))
        }
      }
    }
  }
  PublicationsAllNew<-as.data.frame(PublicationsAllNew)
  DrugsAllNew<-as.data.frame(DrugsAllNew)
  #shift approved drugs to the top
  indexApproved<-which(DrugsAllNew$Status=="Approved")
  Drugstoshiftup<-DrugsAllNew$DRUGS[indexApproved]
  RestofDrugs<-DrugsAllNew$DRUGS[-indexApproved]
  DrugsAllNew<-c(Drugstoshiftup,RestofDrugs)
  
  drugscores<-c(1:length(DrugsAllNew))
  for(i in 1:length(DrugsAllNew)){
    nooccur<-length(grep(DrugsAllNew[i],DrugsAllNew,fixed = TRUE))
    if(nooccur !=0){
      drugscores[i]<-drugscores[i]/nooccur
    }
    pubs<-PublicationsAllNew$Year[grep(DrugsAllNew[i],PublicationsAllNew$Title,fixed = TRUE)]
    nupub<-length(pubs)
    if(nupub !=0){
      drugscores[i]<-drugscores[i]/nupub
      latestpub<-pubs[1]
      drugscores[i]<-drugscores[i]/as.numeric(latestpub)
    }
    
  }
  
  drugscores<-as.data.frame(drugscores)
  drugscores$Drugs<-DrugsAllNew
  #rownames(drugscores)<-DrugsAllNew
  drugscores <- drugscores[order(drugscores$drugscores),]
  DrugsAllNew<-drugscores$Drugs
  DrugsAllNew<-unique(DrugsAllNew)
  PathwaysAllNew<-gsub(" [0-9]+$","",PathwaysAllNew)
  PathwaysAllNew<-unique(PathwaysAllNew)
  PathwaysAllNewReact<-gsub(" [0-9]+$","",PathwaysAllNewReact)
  PathwaysAllNewReact<-unique(PathwaysAllNewReact)
  PathwaysAllNewWiki<-gsub(" [0-9]+$","",PathwaysAllNewWiki)
  PathwaysAllNewWiki<-unique(PathwaysAllNewWiki)
  
  write.table(PathwaysAllNew,paste(disease,"PathwaysKEGG.txt",sep=""),quote = F,row.names = F,sep = "\t",col.names = F)
  write.table(PathwaysAllNewReact,paste(disease,"PathwaysReactome.txt",sep=""),quote = F,row.names = F,sep = "\t",col.names = F)
  write.table(PathwaysAllNewWiki,paste(disease,"PathwaysWiki.txt",sep=""),quote = F,row.names = F,sep = "\t",col.names = F)
  write.table(GOAllNew,paste(disease,"PathwaysGO.txt",sep=""),quote = F,row.names = F,sep = "\t",col.names = F)
  write.table(DrugsAllNew,paste(disease,"DrugsSorted.txt",sep=""),quote = F,row.names = F,sep = "\t",col.names = F)
  # write.table(MOAsNew,paste(disease,"MOAs.txt",sep=""),quote = F,row.names = F,sep = "\t",col.names = F)
}