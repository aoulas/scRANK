extractMalacards<-function (disease,files,path){

  #disease<-"Lymphangioleiomyomatosis"
  #files<-c("Lymphangioleiomyomatosis")#"Myeloma_Multiple" "Autism","Autism_Spectrum_Disorder"

  includeRelDis<-"no"
  setwd(path)

  PathwaysAllDirectlyAllKeywords<-c()
  PathwaysAllAllKeywords<-c()
  GOAllAllKeywords<-c()
  DrugsAllAllKeywords<-c()
  PublicationsAllAllKeywords<-c()
  RelDisAllKeywords<-c()

  for(key in 1:length(files)){
    xhtml<-read_html(files[key]) #Autism_disease_ Malacards.html

    tbls_ls <- xhtml %>%
      html_nodes(".section-data")%>%
      html_table(fill = TRUE)

    head_ls<-xhtml %>%
      html_nodes(".section-data")

    indexPathDir<-0
    for(h3 in 1:length(head_ls)){
      if(grepl("Pathways directly related to",head_ls[h3],fixed=TRUE)){
        indexPathDir<-h3
      }
      if(grepl("Pathways related to",head_ls[h3])){
        indexPathRel<-h3
      }
      if(grepl("Biological processes related",head_ls[h3])){
        indexPathGOBP<-h3
      }
      if(grepl("Drugs for.+DrugBank",head_ls[h3])){
        indexDrugs<-h3
      }
      if(grepl("Articles related to",head_ls[h3])){
        indexPub<-h3
      }
      if(grepl("<h3>Diseases related to",head_ls[h3],ignore.case = F)){
        indexRelDis<-h3
      }

    }

    PathwaysAllDirectly<-c()
    if(indexPathDir!=0){
      PathwaysAllDirectly<-as.data.frame(tbls_ls[[indexPathDir]])
      PathwaysAll<-as.data.frame(tbls_ls[[indexPathRel]])
      GOAll<-as.data.frame(tbls_ls[[indexPathGOBP]])
      DrugsAll<-as.data.frame(tbls_ls[[indexDrugs]])
      PublicationsAll<-as.data.frame(tbls_ls[[indexPub]])
      RelDis<-as.data.frame(tbls_ls[[3]])
      PathwaysAll<-PathwaysAll[,-1]
      PathwaysAllDirectly<-PathwaysAllDirectly[,-1]
    }else{
      PathwaysAll<-as.data.frame(tbls_ls[[indexPathRel]])
      GOAll<-as.data.frame(tbls_ls[[indexPathGOBP]])
      DrugsAll<-as.data.frame(tbls_ls[[indexDrugs]])
      PublicationsAll<-as.data.frame(tbls_ls[[indexPub]])
      RelDis<-as.data.frame(tbls_ls[[indexRelDis]])
      PathwaysAll<-PathwaysAll[,-1]
    }
    PathwaysAllDirectlyAllKeywords<-rbind(PathwaysAllDirectlyAllKeywords,PathwaysAllDirectly)
    PathwaysAllAllKeywords<-rbind(PathwaysAllAllKeywords,PathwaysAll)
    GOAllAllKeywords<-rbind(GOAllAllKeywords,GOAll)
    DrugsAllAllKeywords<-rbind(DrugsAllAllKeywords,DrugsAll)
    PublicationsAllAllKeywords<-rbind(PublicationsAllAllKeywords,PublicationsAll)
    indexstartreldis<-grep("Related Disease",RelDis[,2])
    RelDisToAdd<-RelDis[indexstartreldis+1:20,c(2,3)]
    colnames(RelDisToAdd)<-c("Related Disease","Score")
    RelDisAllKeywords<-rbind(RelDisAllKeywords,RelDisToAdd)
  }

  PathwaysAll<-PathwaysAllAllKeywords
  #Adds the directly related pathways to the rest for processing (careful adds to the top so ranking is affected)
  if(length(PathwaysAllDirectlyAllKeywords)!=0){
    PathwaysAllDirectlyAllKeywords<-paste(PathwaysAllDirectlyAllKeywords[,1],str_extract(PathwaysAllDirectlyAllKeywords[,2], "[0-9]+"),sep="\n")
    maxscore<-max(as.numeric(PathwaysAll$Score))
    PathwaysAllDirectlyAllKeywords<-cbind(PathwaysAllDirectlyAllKeywords,maxscore+1,"Genes")
    colnames(PathwaysAllDirectlyAllKeywords)<-c("Super pathways","Score","Top Affiliating Genes")
    PathwaysAll<-rbind(PathwaysAll,PathwaysAllDirectlyAllKeywords)
  }

  PathwaysAll<-PathwaysAll[order(as.numeric(PathwaysAll$Score),decreasing = TRUE),]
  PathwaysAll<-PathwaysAll %>% distinct(`Super pathways`, .keep_all=TRUE)

  GOAll<-GOAllAllKeywords
  GOAll<-GOAll[order(GOAll$Score,decreasing = TRUE),]
  GOAll<-GOAll %>% distinct(Name, .keep_all=TRUE)

  DrugsAll<-DrugsAllAllKeywords
  PublicationsAll<-PublicationsAllAllKeywords

  RelDisAllKeywords<-RelDisAllKeywords[order(as.numeric(RelDisAllKeywords$Score),decreasing = TRUE),]
  RelDisAllKeywords<-RelDisAllKeywords %>% distinct(`Related Disease`, .keep_all=TRUE)

  indexessuperpaths<-grep("Show member pathways",PathwaysAll[,1])

  #If super pathways exist in related pathways expand them and sort them according to source
  if(length(indexessuperpaths)!=0){
    #firstappend<-TRUE
    PathwaysAllNew<-c()
    PathwaysAllNewReact<-c()
    PathwaysAllNewWiki<-c()
    for(i in 1:length(PathwaysAll[,1])){
      if(i %in% indexessuperpaths){
        Pathessplit<-str_split_1(PathwaysAll[i,1],"Show member pathways")
        Pathessplit1<-str_split_1(Pathessplit[1],"\n")
        Pathessplit1<-Pathessplit1[Pathessplit1 != ""]
        Pathessplit1<-Pathessplit1[Pathessplit1 != " "]
        Pathessplit1<-trimws(Pathessplit1)

        Pathessplit2<-str_split_1(Pathessplit[2],"\n")
        Pathessplit2<-Pathessplit2[Pathessplit2 != ""]
        Pathessplit2<-Pathessplit2[Pathessplit2 != " "]
        Pathessplit2<-trimws(Pathessplit2)


        Pathessplit2numeric<-as.numeric(Pathessplit2)

        Pathessplit2<-c(Pathessplit1,Pathessplit2)
        indexesReactome<-which(Pathessplit2=="66")
        indexesWiki<-which(Pathessplit2=="74")
        indexesother<-which(Pathessplit2numeric !=74 & Pathessplit2numeric !=66)
        if(length(indexesReactome) !=0){
          PathwaysAllNewReact<-append(PathwaysAllNewReact, Pathessplit2[indexesReactome-1], after=length(PathwaysAllNewReact))
        }
        if(length(indexesWiki) !=0){
          PathwaysAllNewWiki<-append(PathwaysAllNewWiki, Pathessplit2[indexesWiki-1], after=length(PathwaysAllNewWiki))
        }

        if(length(indexesother) !=0){
          PathwaysAllNew<-append(PathwaysAllNew, Pathessplit2[indexesother-1], after=length(PathwaysAllNew))
        }
      }else{
        PathessplitNoSubPaths<-str_split_1(PathwaysAll[i,1],"\n")
        PathessplitNoSubPaths<-PathessplitNoSubPaths[PathwaysAll[i,1] != ""]
        PathessplitNoSubPaths<-PathessplitNoSubPaths[PathwaysAll[i,1] != " "]
        PathessplitNoSubPaths<-trimws(PathessplitNoSubPaths)
        indexesReactome<-which(PathessplitNoSubPaths=="66")
        indexesWiki<-which(PathessplitNoSubPaths=="74")
        Pathessplit2numeric2<-as.numeric(PathessplitNoSubPaths)
        indexesother<-which(Pathessplit2numeric2 !=74 & Pathessplit2numeric2 !=66)
        if(length(indexesReactome) !=0){
          PathwaysAllNewReact<-append(PathwaysAllNewReact, PathessplitNoSubPaths[indexesReactome-1], after=length(PathwaysAllNewReact))
        }
        if(length(indexesWiki) !=0){
          PathwaysAllNewWiki<-append(PathwaysAllNewWiki, PathessplitNoSubPaths[indexesWiki-1], after=length(PathwaysAllNewWiki))
        }

        if(length(indexesother) !=0){
          PathwaysAllNew<-append(PathwaysAllNew, PathessplitNoSubPaths[indexesother-1], after=length(PathwaysAllNew))
        }
      }
    }
  }else{
    PathwaysAllNew<-c()
    PathwaysAllNewReact<-c()
    PathwaysAllNewWiki<-c()
    for(i in 1:length(PathwaysAll[,1])){
      PathessplitNoSubPaths<-str_split_1(PathwaysAll[i,1],"\n")
      PathessplitNoSubPaths<-PathessplitNoSubPaths[PathwaysAll[i,1] != ""]
      PathessplitNoSubPaths<-PathessplitNoSubPaths[PathwaysAll[i,1] != " "]
      PathessplitNoSubPaths<-trimws(PathessplitNoSubPaths)
      indexesReactome<-which(PathessplitNoSubPaths=="66")
      indexesWiki<-which(PathessplitNoSubPaths=="74")
      Pathessplit2numeric2<-as.numeric(PathessplitNoSubPaths)
      indexesother<-which(Pathessplit2numeric2 !=74 & Pathessplit2numeric2 !=66)
      if(length(indexesReactome) !=0){
        PathwaysAllNewReact<-append(PathwaysAllNewReact, PathessplitNoSubPaths[indexesReactome-1], after=length(PathwaysAllNewReact))
      }
      if(length(indexesWiki) !=0){
        PathwaysAllNewWiki<-append(PathwaysAllNewWiki, PathessplitNoSubPaths[indexesWiki-1], after=length(PathwaysAllNewWiki))
      }

      if(length(indexesother) !=0){
        PathwaysAllNew<-append(PathwaysAllNew, PathessplitNoSubPaths[indexesother-1], after=length(PathwaysAllNew))
      }
    }
  }

  PathwaysAllNew<-gsub("\n.+","",PathwaysAllNew)
  PathwaysAllNew<-unique(PathwaysAllNew)
  PathwaysAllNew<-gsub("\n+","",PathwaysAllNew)

  #Get GO terms
  GOAll<-GOAll[,-1]
  GOAll[,1]<-gsub("\n","",GOAll[,1])
  GOAllNew<-GOAll[,1][GOAll[,1] != "Name"]
  GOAllNew<-gsub("[0-9]+$","",GOAllNew)


  #Get Drugs and Drug MOAs
  indexesSyn<-grep("Synonyms",DrugsAll$X3)
  indexesDrugs<-indexesSyn-1
  Synonyms<-DrugsAll$X3[indexesSyn]
  DrugsAllNames<-unique(DrugsAll$X3[indexesDrugs])
  MOAs<-DrugsAll$X3[-c(indexesDrugs,indexesSyn)]
  MOAs<-MOAs[-(grep("Status",MOAs):length(MOAs))]
  MOAs<-MOAs[-(which((is.na(MOAs))))]
  MOAsNew<-MOAs[MOAs != "Name"]

  #Calculate Drug scores based on occurrence in bibliography and year of latest publication (may take a while)
  drugscores<-c(1:length(DrugsAllNames))
  for(i in 1:length(DrugsAllNames)){
    nooccur<-length(grep(DrugsAllNames[i],DrugsAll$X6))
    if(nooccur !=0){
      drugscores[i]<-drugscores[i]/nooccur
    }
    pubs<-PublicationsAll$Year[grep(DrugsAllNames[i],PublicationsAll$Title)]
    nupub<-length(pubs)
    if(nupub !=0){
      drugscores[i]<-drugscores[i]/nupub
      latestpub<-pubs[1]
      drugscores[i]<-drugscores[i]/latestpub
    }

  }
  drugscores<-as.data.frame(drugscores)
  drugscores$Drugs<-DrugsAllNames
  rownames(drugscores)<-DrugsAllNames
  drugscores <- drugscores[order(drugscores$drugscores),]

  if(includeRelDis=="yes"){
    PathwaysAllNew<-c(PathwaysAllNew,RelDisAllKeywords[,1])
    PathwaysAllNewWiki<-c(PathwaysAllNewWiki,RelDisAllKeywords[,1])
  }
  write.table(PathwaysAllNew,paste(disease,"PathwaysKEGG.txt",sep=""),quote = F,row.names = F,sep = "\t",col.names = F)
  write.table(PathwaysAllNewReact,paste(disease,"PathwaysReactome.txt",sep=""),quote = F,row.names = F,sep = "\t",col.names = F)
  write.table(PathwaysAllNewWiki,paste(disease,"PathwaysWiki.txt",sep=""),quote = F,row.names = F,sep = "\t",col.names = F)
  write.table(GOAllNew,paste(disease,"PathwaysGO.txt",sep=""),quote = F,row.names = F,sep = "\t",col.names = F)
  write.table(drugscores$Drugs,paste(disease,"DrugsSorted.txt",sep=""),quote = F,row.names = F,sep = "\t",col.names = F)
  write.table(MOAsNew,paste(disease,"MOAs.txt",sep=""),quote = F,row.names = F,sep = "\t",col.names = F)
}
