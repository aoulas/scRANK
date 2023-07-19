rankCells<-function (seuratObject,path,scan,priorknowledgePathsKEGG,priorknowledgePathsGO,priorknowledgePathsMSIG,priorknowledgePathsWiki,priorknowledgePathsReact,priorknowledgeMOA,labels,cellIDs,checkdrug,scenario){

  seuratObject$labels.cellIDs <- paste(as.character(seuratObject[[labels]][,1]), as.character(seuratObject[[cellIDs]][,1]), sep = "_")
  Idents(seuratObject) <- "labels.cellIDs"

  LablesUniq<-unique(as.character(seuratObject[[labels]][,1]))
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes
  websiteLive <- TRUE

  #Make sure Assay is RNA
  DefaultAssay(seuratObject) <- "RNA"
  tablecellclounts <-table(seuratObject@meta.data$nCount_RNA, seuratObject@meta.data$labels.cellIDs)
  tablecellclounts<-as.data.frame(colSums(tablecellclounts))
  #seuratObjectcellIDmatchesMOA<-c()
  seuratObjectcellIDmatchesMOAEuc<-c()
  #seuratObjectcellIDmatchesPATHS<-c()
  seuratObjectcellIDmatchesPATHSEuc<-c()
  seuratObjectcellIDmatchesPATHSGOEuc<-c()
  seuratObjectcellIDmatchesPATHSMSIGEuc<-c()
  seuratObjectcellIDmatchesPATHSWikiEuc<-c()
  seuratObjectcellIDmatchesPATHSReactEuc<-c()

  tableMOAsAll<-c()
  mlist<-c()

  #Control normally
  indexcountsCondition1<-0
  #Case normally
  indexcountsCondition2<-0

  AllEnrichedPathsKEGGandGO<-c()
  Allavelog2FC<-c()
  TotalNumberDEGs<-c()

  for(cellID in as.character(unique(seuratObject[[cellIDs]][,1]))){#
    if(scan=="Bulk"){
      print("Scanning Bulk RNA")
      Condition1<-seuratObject$labels.cellIDs[grep(LablesUniq[1],seuratObject$labels.cellIDs)]
      Condition2<-seuratObject$labels.cellIDs[grep(LablesUniq[2],seuratObject$labels.cellIDs)]
      options(future.globals.maxSize = 8000 * 1024^2)
    }else{
      print(cellID)
      Condition1<-seuratObject$labels.cellIDs[grep(paste(LablesUniq[1],"_",cellID,sep=""),seuratObject$labels.cellIDs)]
      Condition2<-seuratObject$labels.cellIDs[grep(paste(LablesUniq[2],"_",cellID,sep=""),seuratObject$labels.cellIDs)]

      indexcountsCondition1<-which(rownames(tablecellclounts)==paste(LablesUniq[1],"_",cellID,sep = ""))
      indexcountsCondition2<-which(rownames(tablecellclounts)==paste(LablesUniq[2],"_",cellID,sep = ""))
    }

    if(length(indexcountsCondition1) !=0 && length(indexcountsCondition2) !=0){
      if((scan=="Cell" && tablecellclounts[indexcountsCondition1,] > 3 && tablecellclounts[indexcountsCondition2,] > 3) || (scan=="Bulk")){
        seuratObject.markersCellIDs <- FindMarkers(seuratObject, ident.1 = Condition1, ident.2 = Condition2, verbose = FALSE)

        indexsigclCellIDs<-which(seuratObject.markersCellIDs$p_val_adj<=0.05)
        seuratObject.markersCellIDs_sig<-seuratObject.markersCellIDs[indexsigclCellIDs,]
        Up<-rownames(seuratObject.markersCellIDs_sig[which(seuratObject.markersCellIDs_sig$avg_log2FC>=0),])
        Down<-rownames(seuratObject.markersCellIDs_sig[which(seuratObject.markersCellIDs_sig$avg_log2FC<0),])


        UpGenes<-seuratObject.markersCellIDs_sig[which(seuratObject.markersCellIDs_sig$avg_log2FC>=0),]
        avelog2FC<-mean(UpGenes$avg_log2FC)
        TotalNumberDEGs<-c(TotalNumberDEGs,(length(Up)+length(Down)))
        Allavelog2FC<-c(Allavelog2FC,avelog2FC)

        dbstouse <- c("Old_CMAP_up","Old_CMAP_down")
        dbstousePaths<-c("KEGG_2021_Human","GO_Biological_Process_2021","Reactome_2022","MSigDB_Hallmark_2020","WikiPathway_2021_Human")

        #Reducing number of sig genes worked better in differentiating between cell-specific and bulk rna-seq
        maxnumofsiggenes<-200

        if(length(Up) > maxnumofsiggenes){
          Up<-Up[1:maxnumofsiggenes]
        }
        if(length(Down) > maxnumofsiggenes){
          Down<-Down[1:maxnumofsiggenes]
        }

        if (websiteLive) {
          enrichedUp <- enrichr(Up, dbstouse)
        }

        if (websiteLive) {
          enrichedDown <- enrichr(Down, dbstouse)
        }

        if(length(enrichedDown) !=0 && length(enrichedUp) !=0){

          if (websiteLive) {
            enrichedPaths <- enrichr(c(Up,Down), dbstousePaths)#
          }
          if(scan=="Cell"){
            if (websiteLive) plot(plotEnrich(enrichedPaths[["KEGG_2021_Human"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = paste("KEGG",cellID,sep=" "))+plotEnrich(enrichedPaths[["GO_Biological_Process_2021"]], showTerms = 40,numChar = 40, y = "Count", orderBy = "P.value",title = paste("GO_Bio_Pro",cellID))+plotEnrich(enrichedPaths[["MSigDB_Hallmark_2020"]], showTerms = 40,numChar = 40, y = "Count", orderBy = "P.value",title = paste("MSigDB",cellID))+plotEnrich(enrichedPaths[["WikiPathway_2021_Human"]], showTerms = 40,numChar = 40, y = "Count", orderBy = "P.value",title = paste("WIKI",cellID))+plotEnrich(enrichedPaths[["Reactome_2022"]], showTerms = 40,numChar = 40, y = "Count", orderBy = "P.value",title = paste("Reactome",cellID)))
          }else{
            if (websiteLive) plot(plotEnrich(enrichedPaths[["KEGG_2021_Human"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = paste("KEGG","Bulk",sep=" "))+plotEnrich(enrichedPaths[["GO_Biological_Process_2021"]], showTerms = 40,numChar = 40, y = "Count", orderBy = "P.value",title = paste("GO_Bio_Pro","Bulk"))+plotEnrich(enrichedPaths[["MSigDB_Hallmark_2020"]], showTerms = 40,numChar = 40, y = "Count", orderBy = "P.value",title = paste("MSigDB","Bulk")))
          }

          paths<-enrichedPaths[["KEGG_2021_Human"]]
          paths<-paths[which(paths$P.value<0.05),]
          pathsGO<-enrichedPaths[["GO_Biological_Process_2021"]]
          pathsGO<-pathsGO[which(pathsGO$P.value<0.05),]
          pathsMSIG<-enrichedPaths[["MSigDB_Hallmark_2020"]]
          pathsWiki<-enrichedPaths[["WikiPathway_2021_Human"]]
          pathsWiki<-pathsWiki[which(pathsWiki$P.value<0.05),]
          pathsReact<-enrichedPaths[["Reactome_2022"]]
          pathsReact<-pathsReact[which(pathsReact$P.value<0.05),]

          if(scan=="Cell"){
            AllEnrichedPathsKEGGandGO<-rbind(AllEnrichedPathsKEGGandGO,cellID,paths,pathsGO,pathsMSIG)
          }else{
            AllEnrichedPathsKEGGandGO<-rbind(AllEnrichedPathsKEGGandGO,"Bulk_RNA",paths,pathsGO)
          }



          matchKEGG<-unique(grep(paste(priorknowledgePathsKEGG,collapse="|"),
                                 paths[,1], value=TRUE,ignore.case = T))


          matchGO<-unique(grep(paste(priorknowledgePathsGO,collapse="|"),
                               pathsGO[,1], value=TRUE,ignore.case = T,fixed = T))

          matchMSIG<-unique(grep(paste(priorknowledgePathsMSIG,collapse="|"),
                                 pathsMSIG[,1], value=TRUE,ignore.case = T))


          #KEGG Euclidean
          matchKEGGIndexesOrdered<-match(tolower(priorknowledgePathsKEGG),tolower(iconv(as.character(paths[,1]),"ISO-8859-1")))
          if(length(unique(matchKEGGIndexesOrdered)) == 1 && is.na(unique(matchKEGGIndexesOrdered))){
            matchKEGGIndexesOrdered<-tidyr::replace_na(matchKEGGIndexesOrdered,1000)
          }else{
            matchKEGGIndexesOrdered<-tidyr::replace_na(matchKEGGIndexesOrdered,(length(priorknowledgePathsKEGG)+1000))
          }

          euclidean <- function(a, b) sqrt(sum((a - b)^2))
          if(is.null(priorknowledgePathsKEGG)){
            celleucpathKEGG<-1000
          }else{
            celleucpathKEGG<-euclidean(c(1:length(priorknowledgePathsKEGG)),matchKEGGIndexesOrdered)
          }
          seuratObjectcellIDmatchesPATHSEuc<-rbind(seuratObjectcellIDmatchesPATHSEuc,cbind(cellID,celleucpathKEGG))


          matchKEGGIndexes<-unique(grep(paste(priorknowledgePathsKEGG,collapse="|"),
                                        paths[,1], ignore.case = T))

          #GO Terms Euclidean
          matchGOIndexesOrdered<-match(tolower(priorknowledgePathsGO),tolower(iconv(as.character(gsub(" +\\(GO.[0-9]+\\)","",pathsGO[,1])),"ISO-8859-1")))
          if(length(unique(matchGOIndexesOrdered)) == 1 && is.na(unique(matchGOIndexesOrdered))){
            matchGOIndexesOrdered<-tidyr::replace_na(matchGOIndexesOrdered,1000)
          }else{
            matchGOIndexesOrdered<-tidyr::replace_na(matchGOIndexesOrdered,(length(priorknowledgePathsGO)+1000))
          }

          if(is.null(priorknowledgePathsGO)){
            celleucpathGO<-1000
          }else{
            celleucpathGO<-euclidean(c(1:length(priorknowledgePathsGO)),matchGOIndexesOrdered)
          }
          seuratObjectcellIDmatchesPATHSGOEuc<-rbind(seuratObjectcellIDmatchesPATHSGOEuc,cbind(cellID,celleucpathGO))

          matchGOIndexes<-unique(grep(paste(priorknowledgePathsGO,collapse="|"),
                                      pathsGO[,1], ignore.case = T,fixed = T))

          #MSIG Euclidean
          formattedMSIG<-gsub(" ","_",tolower(iconv(as.character(pathsMSIG[,1]),"ISO-8859-1")))
          formattedMSIG<-gsub("_+","_",formattedMSIG)
          formattedMSIG<-gsub("/","_",formattedMSIG)
          findfuzz<-agrep(tolower(priorknowledgePathsMSIG)[3], formattedMSIG, max.distance = 4, value = TRUE)
          if(length(findfuzz)!=0){
            indexfindfuzz<-which(formattedMSIG==findfuzz)
          }
          matchMSIGIndexesOrdered<-match(tolower(priorknowledgePathsMSIG),formattedMSIG)
          if(length(unique(matchMSIGIndexesOrdered)) == 1 && is.na(unique(matchMSIGIndexesOrdered))){
            matchMSIGIndexesOrdered<-tidyr::replace_na(matchMSIGIndexesOrdered,1000)
          }else{
            matchMSIGIndexesOrdered<-tidyr::replace_na(matchMSIGIndexesOrdered,(length(priorknowledgePathsMSIG)+1000))
          }

          if(is.null(priorknowledgePathsMSIG)){
            celleucpathMSIG<-1000
          }else{
            celleucpathMSIG<-euclidean(c(1:length(priorknowledgePathsMSIG)),matchMSIGIndexesOrdered)
          }
          seuratObjectcellIDmatchesPATHSMSIGEuc<-rbind(seuratObjectcellIDmatchesPATHSMSIGEuc,cbind(cellID,celleucpathMSIG))

          matchMSIGIndexes<-unique(grep(paste(priorknowledgePathsMSIG,collapse="|"),
                                        pathsMSIG[,1], ignore.case = T))



          #Wiki Euclidean perform fuzzy search on wiki pathways because there are some discrepancies
          indexfindfuzz<-c()
          for(wikii2 in 1:length(priorknowledgePathsWiki)){
            findfuzzWiki<-agrep(tolower(priorknowledgePathsWiki)[wikii2], tolower(iconv(as.character(gsub(" +WP[0-9]+","",pathsWiki[,1])),"ISO-8859-1")), max.distance = 4)
            if(length(findfuzzWiki)!=0){
              indexfindfuzz<-c(indexfindfuzz,findfuzzWiki)
            }else{
              indexfindfuzz<-c(indexfindfuzz,NA)
            }
          }
          matchWikiIndexesOrdered<-match(tolower(priorknowledgePathsWiki),tolower(iconv(as.character(gsub(" +WP[0-9]+","",pathsWiki[,1])),"ISO-8859-1")))
          if(length(unique(matchWikiIndexesOrdered)) == 1 && is.na(unique(matchWikiIndexesOrdered))){
            matchWikiIndexesOrdered<-rep(1000,length(matchWikiIndexesOrdered))
          }else{
            matchWikiIndexesOrdered<-tidyr::replace_na(matchWikiIndexesOrdered,(length(priorknowledgePathsWiki)+1000))
          }

          if(is.null(priorknowledgePathsWiki)){
            celleucpathWiki<-1000
          }else{
            celleucpathWiki<-euclidean(c(1:length(priorknowledgePathsWiki)),matchWikiIndexesOrdered)
          }

          seuratObjectcellIDmatchesPATHSWikiEuc<-rbind(seuratObjectcellIDmatchesPATHSWikiEuc,cbind(cellID,celleucpathWiki))

          indexespriorknowledgePathsWikifound<- which(!is.na(indexfindfuzz))
          indexespathsWikifound<-indexfindfuzz[indexespriorknowledgePathsWikifound]

          #Reactome Euclidean
          matchReactIndexesOrdered<-match(tolower(priorknowledgePathsReact),tolower(iconv(as.character(gsub(" +R-HSA-[0-9]+","",pathsReact[,1])),"ISO-8859-1")))
          if(length(unique(matchReactIndexesOrdered)) == 1 && is.na(unique(matchReactIndexesOrdered))){
            matchReactIndexesOrdered<-tidyr::replace_na(matchReactIndexesOrdered,1000)
          }else{
            matchReactIndexesOrdered<-tidyr::replace_na(matchReactIndexesOrdered,(length(priorknowledgePathsReact)+1000))
          }

          if(is.null(priorknowledgePathsReact)){
            celleucpathReact<-1000
          }else{
            celleucpathReact<-euclidean(c(1:length(priorknowledgePathsReact)),matchReactIndexesOrdered)
          }
          seuratObjectcellIDmatchesPATHSReactEuc<-rbind(seuratObjectcellIDmatchesPATHSReactEuc,cbind(cellID,celleucpathReact))

          #Create river plot of matches##############################################################################
          # testriver<-c()
          # testriver<-data.frame(cbind(priorknowledgePathsWiki[indexespriorknowledgePathsWikifound],pathsWiki[,1][indexespathsWikifound],rep(1,each=length(indexespriorknowledgePathsWikifound))))
          # colnames(testriver)<-c("N1","N2","Value")
          # testriver$Value<-as.numeric(testriver$Value)
          # testriver<-testriver[!duplicated(tolower(testriver$N1)), ]
          # nodes<-data.frame(ID=make.unique(c(priorknowledgePathsWiki,pathsWiki[,1])),x=c(rep(1,each=length(priorknowledgePathsWiki)),rep(2,each=length(pathsWiki[,1]))),y=c(rev(1:length(priorknowledgePathsWiki)),rev(1:length(pathsWiki[,1]))))
          #
          # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
          # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
          #
          #
          # style<-list()
          # for(nodei in 1:length(nodes$ID)){
          #   styletemp<-sapply(nodes$ID[nodei],function(id) list(col=sample(col_vector, 1)),simplify=FALSE)
          #   style<-append(style, styletemp, after = length(style))
          # }
          #
          # #dev.off()
          # len<-1000
          # plot(NULL,xlim=c(1,len),ylim=c(1,len),bty="n",xlab="",ylab="",xaxt='n',yaxt='n')
          # r<-makeRiver(nodes,testriver,styles=style)
          # d<-list(srt=0,textcex=1.2)#,default_style=d)
          # #plot(r,plot_area=1,nodewidth=20,default_style=d,yscale=1)
          # plot(r,plot_area=1,nodewidth=18,default_style=d,yscale=1,add=TRUE,usr=c(1,len,1,len))
          ############################################################################################################

          # if(nposPathKEGG==0){
          #   nposPathKEGG<-nrow(paths)
          # }
          # if(nposPathGO==0){
          #   nposPathGO<-nrow(pathsGO)
          # }
          #
          # if((length(matchKEGG) != 0 && length(which(c(1:nposPathKEGG)%in%matchKEGGIndexes)) !=0) || (length(matchGO) !=0 && length(which(c(1:nposPathGO)%in%matchGOIndexes)))){
          #   seuratObjectcellIDmatchesPATHS<-rbind(seuratObjectcellIDmatchesPATHS,cbind(cellID,"1"))
          # }else{
          #   seuratObjectcellIDmatchesPATHS<-rbind(seuratObjectcellIDmatchesPATHS,cbind(cellID,"0"))
          # }

          drugsUp<-enrichedUp[["Old_CMAP_down"]]$Term[which(enrichedUp[["Old_CMAP_down"]]$P.value <0.05)]
          drugsUp<-gsub("-[0-9]+$","",drugsUp)

          drugsDown<-enrichedDown[["Old_CMAP_up"]]$Term[which(enrichedDown[["Old_CMAP_up"]]$P.value <0.05)]
          drugsDown<-gsub("-[0-9]+$","",drugsDown)

          if(scan=="Cell"){
            if(nrow(enrichedUp[["Old_CMAP_down"]]) != 0 && nrow(enrichedDown[["Old_CMAP_up"]]) != 0){
              if (websiteLive) plot(plotEnrich(enrichedUp[["Old_CMAP_down"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = paste("CMAP", cellID, "UP",sep=" "))+plotEnrich(enrichedDown[["Old_CMAP_up"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = paste("CMAP", cellID, "Down",sep=" ")))
            }else if(nrow(enrichedUp[["Old_CMAP_down"]]) != 0 && nrow(enrichedDown[["Old_CMAP_up"]]) == 0){
              if (websiteLive) plot(plotEnrich(enrichedUp[["Old_CMAP_down"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = paste("CMAP", cellID, "UP",sep=" ")))
            }else if(nrow(enrichedUp[["Old_CMAP_down"]]) == 0 && nrow(enrichedDown[["Old_CMAP_up"]]) != 0){
              if (websiteLive) plotEnrich(enrichedDown[["Old_CMAP_up"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = paste("CMAP", cellID, "Down",sep=" "))
            }
          }else{
            if (websiteLive) plot(plotEnrich(enrichedUp[["Old_CMAP_down"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = paste("CMAP", "Bulk", "UP",sep=" "))+plotEnrich(enrichedDown[["Old_CMAP_up"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = paste("CMAP", "Bulk", "Down",sep=" ")))
          }
          if(!endsWith(path,"/")){
            path<-paste(path,"/",sep="")
          }
          drugInfo<-read.delim(paste(path,"drug_repurposing_hub.txt",sep=""))
          AlldrugsCellIDs<-c()
          if(length(drugsUp) > 50){
            AlldrugsCellIDs<-c(AlldrugsCellIDs,drugsUp[1:50])
          }else{
            AlldrugsCellIDs<-c(AlldrugsCellIDs,drugsUp)
          }

          if(length(drugsDown) > 50){
            AlldrugsCellIDs<-c(AlldrugsCellIDs,drugsDown[1:50])
          }else{
            AlldrugsCellIDs<-c(AlldrugsCellIDs,drugsDown)
          }

          if(length(AlldrugsCellIDs) < 100){
            extradrugs<-100-length(AlldrugsCellIDs)-1
            indexmax<-which(c(length(drugsUp),length(drugsDown))==max(length(drugsUp),length(drugsDown)))
            if(indexmax==2){
              AlldrugsCellIDs<-c(AlldrugsCellIDs,drugsDown[51:(51+extradrugs)])
            }else{
              AlldrugsCellIDs<-c(AlldrugsCellIDs,drugsUp[51:(51+extradrugs)])
            }
          }
          foundds<-0
          MOA<-c()
          if(length(AlldrugsCellIDs) > 100){
            AlldrugsCellIDs<-AlldrugsCellIDs[1:100]
          }
          for(d in 1:length(AlldrugsCellIDs)){
            indexfoundd<-grep(AlldrugsCellIDs[d],drugInfo$pert_iname)
            if(!length(indexfoundd)==0){
              foundds<-foundds+1
              MOA<-rbind(MOA,cbind(drugInfo$pert_iname[indexfoundd],drugInfo$moa[indexfoundd],drugInfo$disease_area[indexfoundd]))
            }
          }
          indexesnull1<-which(MOA[,2] == "")
          if(!length(indexesnull1)==0){
            MOA<-MOA[-indexesnull1,]
          }
          indexesnull2<-which(MOA[,3] == "")
          if(!length(indexesnull2)==0){
            MOA<-MOA[-indexesnull2,]
          }

          # sort by Freq
          tableMOAs<-as.data.frame(table(MOA[,2]))

          #tableMOAs[,1][grep("sera",tableMOAs[,1])]
          matchMOA<-unique(grep(paste(priorknowledgeMOA,collapse="|"),
                                as.character(tableMOAs$Var1), value=TRUE,ignore.case = T))
          tableMOAs <- tableMOAs[order(tableMOAs$Freq,decreasing = T),]
          if(scan=="Cell"){
            tableMOAsAll<-rbind(tableMOAsAll,cellID,tableMOAs)
          }else{
            tableMOAsAll<-rbind(tableMOAsAll,"Bulk_RNA",tableMOAs)
          }

          if(checkdrug==FALSE){
            #Checks MOAs directly
            matchMOAIndexesOrdered<-match(tolower(priorknowledgeMOA),tolower(iconv(as.character(tableMOAs$Var1),"ISO-8859-1")))
          }else{
            #Checks drugs directly
            matchMOAIndexesOrdered<-match(tolower(priorknowledgeMOA),tolower(iconv(as.character(AlldrugsCellIDs),"ISO-8859-1")))
          }
          AlldrugsCellIDs[matchMOAIndexesOrdered[!isNA(matchMOAIndexesOrdered)]]

          sumoffreq<-sum(tableMOAs$Freq[matchMOAIndexesOrdered][!is.na(tableMOAs$Freq[matchMOAIndexesOrdered])])
          if(length(unique(matchMOAIndexesOrdered)) == 1 && is.na(unique(matchMOAIndexesOrdered))){
            matchMOAIndexesOrdered<-tidyr::replace_na(matchMOAIndexesOrdered,1000)
          }else{
            matchMOAIndexesOrdered<-tidyr::replace_na(matchMOAIndexesOrdered,(length(AlldrugsCellIDs)+1))
          }

          matchMOAIndexes<-unique(grep(paste(priorknowledgeMOA,collapse="|"),
                                       as.character(tableMOAs$Var1), ignore.case = T))


          if(!is.na(sumoffreq) && sumoffreq !=0 && checkdrug==FALSE){
            celleuc<-euclidean(c(1:length(priorknowledgeMOA)),matchMOAIndexesOrdered)
            celleuc<-celleuc/sumoffreq
          }else{
            celleuc<-euclidean(c(1:length(priorknowledgeMOA)),matchMOAIndexesOrdered)
          }
          seuratObjectcellIDmatchesMOAEuc<-rbind(seuratObjectcellIDmatchesMOAEuc,cbind(cellID,celleuc))
          # if(nposMOA==0){
          #   nposMOA<-nrow(tableMOAs)
          # }
          # if(length(matchMOA) != 0 && length(which(c(1:nposMOA)%in%matchMOAIndexes)) !=0){
          #   seuratObjectcellIDmatchesMOA<-rbind(seuratObjectcellIDmatchesMOA,cbind(cellID,"1"))
          # }else{
          #   seuratObjectcellIDmatchesMOA<-rbind(seuratObjectcellIDmatchesMOA,cbind(cellID,"0"))
          # }
          mlist<-rbind(mlist,rep(c(cellID),times=ncol(paths)),paths,pathsGO)
        }
      }else{
        Allavelog2FC<-c(Allavelog2FC,0)
        TotalNumberDEGs<-c(TotalNumberDEGs,0)
      }
    }else{
      Allavelog2FC<-c(Allavelog2FC,0)
      TotalNumberDEGs<-c(TotalNumberDEGs,0)
    }
    if(scan=="Bulk"){
      break;
    }
  }
  if(scan=="Bulk"){
    #colnames(seuratObjectcellIDmatchesMOA)<-c("Bulk_RNA","Matched_MOA")
    #seuratObjectcellIDmatchesMOA[1,1]<-"Bulk.RNA"
    seuratObjectcellIDmatchesMOAEuc[1,1]<-"Bulk.RNA"
    seuratObjectcellIDmatchesPATHSEuc[1,1]<-"Bulk.RNA"
    seuratObjectcellIDmatchesPATHSGOEuc[1,1]<-"Bulk.RNA"
    seuratObjectcellIDmatchesPATHSMSIGEuc[1,1]<-"Bulk.RNA"
    seuratObjectcellIDmatchesPATHSWikiEuc[1,1]<-"Bulk.RNA"
    #colnames(seuratObjectcellIDmatchesPATHS)<-c("Bulk)RNA","Matched_Paths")
    #seuratObjectcellIDmatchesPATHS[1,1]<-"Bulk.RNA"
  }else{
    #colnames(seuratObjectcellIDmatchesMOA)<-c("cellIDs","Matched_MOA")
    #colnames(seuratObjectcellIDmatchesPATHS)<-c("cellIDs","Matched_Paths")
    Allavelog2FC<-as.data.frame(Allavelog2FC)
    rownames(Allavelog2FC)<-as.character(unique(seuratObject[[cellIDs]][,1]))
    TotalNumberDEGs<-as.data.frame(TotalNumberDEGs)
    rownames(TotalNumberDEGs)<-as.character(unique(seuratObject[[cellIDs]][,1]))
  }

  listofCellRanks <- list(seuratObjectcellIDmatchesMOAEuc,seuratObjectcellIDmatchesPATHSEuc,seuratObjectcellIDmatchesPATHSGOEuc,seuratObjectcellIDmatchesPATHSMSIGEuc,seuratObjectcellIDmatchesPATHSWikiEuc,seuratObjectcellIDmatchesPATHSReactEuc,tableMOAsAll,AllEnrichedPathsKEGGandGO,Allavelog2FC,TotalNumberDEGs)#mlist,


  Ranks<-c()
  EucValues<-c()
  for(li in c(1,2,3,4,5,6)){
    print(li)
    fromlist<-as.data.frame(listofCellRanks[[li]])
    fromlist[,2]<-as.numeric(fromlist[,2])
    fromlist <- arrange(fromlist, fromlist[,2])
    fromlist2 <- arrange(fromlist, fromlist[,1])
    Ranks<-cbind(Ranks,rank(fromlist2[,2]))
    EucValues<-cbind(EucValues,fromlist2[,2])
  }
  colnames(Ranks)<-c("DRUG","KEGG","GOBP","MSIG","WIKI","REACT")
  if(scan=="Cell"){
    rownames(Ranks)<-fromlist2$cellID
  }else{
    rownames(Ranks)<-"Bulk"
  }
  colnames(EucValues)<-c("DRUG","KEGG","GOBP","MSIG","WIKI","REACT")
  if(scan=="Cell"){
    rownames(EucValues)<-fromlist2$cellID
  }else{
    rownames(EucValues)<-"Bulk"
  }
  write.table(Ranks,paste("Ranks",scenario,"PlusWikiReactFinalPac.txt",sep=""),quote = F,row.names = T,sep = "\t")
  write.table(EucValues,paste("EucValues",scenario,"PlusWikiReactFinalPac.txt",sep=""),quote = F,row.names = T,sep = "\t")

  MeanRanks<-as.data.frame(meanranks(t(Ranks))$mean.ranks)
  colnames(MeanRanks)<-c("Mean Ranking")
  MeanRanks<-arrange(MeanRanks,MeanRanks$`Mean Ranking`)

  write.table(MeanRanks,paste("MeanRanks",scenario,"PlusWikiReactFinalPac.txt",sep=""),quote = F,row.names = T,sep = "\t")


  Allavelog2FC<-as.data.frame(listofCellRanks[[9]])
  write.table(Allavelog2FC,"Allavelog2FCPac.txt",quote = F,row.names = T,sep = "\t")

  TotalNumberDEGs<-as.data.frame(listofCellRanks[[10]])
  TotalNumberDEGs$Var1<-rownames(TotalNumberDEGs)
  TotalNumberDEGs<-as.data.frame(TotalNumberDEGs)
  TotalNumberDEGs <-TotalNumberDEGs[order(TotalNumberDEGs$Var1), ]
  write.table(TotalNumberDEGs,"TotalNumberDEGsPac.txt",quote = F,row.names = T,sep = "\t")

  if(scan=="Cell"){
    # # How many cells are in each cell type or condition?
    Nocellspercelltype<-as.data.frame(table(seuratObject[[usercelltype]][,1]))
    Nocellspercelltype<-arrange(Nocellspercelltype, Nocellspercelltype$Var1)

    TotalNumberDEGsNorm<-as.data.frame(TotalNumberDEGs$TotalNumberDEGs/Nocellspercelltype$Freq)
    rownames(TotalNumberDEGsNorm)<-Nocellspercelltype$Var1
    colnames(TotalNumberDEGsNorm)[1]<-c("TotalNumberDEGsNorm")
    write.table(TotalNumberDEGsNorm,"TotalNumberDEGsNormPac.txt",quote = F,row.names = T,sep = "\t")
  }
  return(listofCellRanks)
}
