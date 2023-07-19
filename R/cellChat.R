runCellChat<-function(seuratObject,labelsC,cellIDs){
  library(CellChat)
  options(stringsAsFactors = FALSE)
  LablesUniq<-unique(as.character(seuratObject[[labelsC]][,1]))
  indexdf.netlist<-1
  df.net<-vector("list", 2)

  for(labelIndex in LablesUniq){
    expr <- FetchData(seuratObject, vars = labelsC)
    seuratObjectSub <- seuratObject[, which(expr == labelIndex)]

    data.input<-GetAssayData(object = seuratObjectSub, slot = "data",assay = "RNA")#Does this need to be RNA check

    cell.use = names(seuratObjectSub$orig.ident)

    meta<-c()
    # Prepare input data for CellChat analysis
    data.input = data.input[, cell.use]
    #meta = meta[cell.use, ]
    #or this
    seuratObjectSub$celltype2<-seuratObjectSub[[cellIDs]][,1]

    meta = data.frame(labels = seuratObjectSub$celltype2[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
    unique(meta$labels) # check the cell labels
    #Drops any unused levels
    meta<-droplevels(meta)
    if(is.null(levels(meta$labels))){
      print("adding levels")
      levels(meta$labels)<-unique(meta$labels)
    }

    #Create a new CellChat object from a data matrix, Seurat or SingleCellExperiment object
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

    #Add cell information into meta slot of the object (Optional)
    #If cell mata information is not added when creating CellChat object, USERS can also add it later using addMeta
    cellchat <- addMeta(cellchat, meta = meta)
    cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
    levels(cellchat@idents) # show factor levels of the cell labels
    groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group


    #Set the ligand-receptor interaction database
    #CellChatDB is a manually curated database of literature-supported ligand-receptor interactions
    #in both human and mouse

    CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
    showDatabaseCategory(CellChatDB)

    # Show the structure of the database
    #dplyr::glimpse(CellChatDB$interaction)


    # use a subset or all (default) of CellChatDB for cell-cell communication analysis
    CellChatDB.use <- CellChatDB
    #subsetDB(CellChatDB, search = "ECM-Receptor") # use Secreted Signaling

    # set the used database in the object
    cellchat@DB <- CellChatDB.use

    unique(cellchat@idents)
    #Preprocessing the expression data for cell-cell communication analysis
    #To infer the cell state-specific communications, we identify over-expressed ligands or receptors in one cell group
    #and then identify over-expressed ligand-receptor interactions if either ligand or receptor is over-expressed

    # subset the expression data of signaling genes for saving computation cost
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    #future::plan("multiprocess", workers = 4) # do parallel

    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)

    # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
    # cellchat <- projectData(cellchat, PPI.human)


    #Part II: Inference of cell-cell communication network
    #CellChat infers the biologically significant cell-cell communication by assigning each interaction with a
    #probability value and peforming a permutation test. CellChat models the probability of cell-cell
    #communication by integrating gene expression with prior known knowledge of the interactions between signaling ligands, receptors and their cofactors using the law of mass action.
    #Compute the communication probability and infer cellular communication network
    #If well-known signaling pathways in the studied biological process are not predicted, USER can try truncatedMean
    #to change the method for calculating the average gene expression per cell group.

    cellchat <- computeCommunProb(cellchat)


    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 10)


    #Extract the inferred cellular communication network as a data frame
    #Function subsetCommunication is used to easily access the inferred cell-cell communications of interest. For example,

    #df.net <- subsetCommunication(cellchat) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
    df.net[[indexdf.netlist]]<-subsetCommunication(cellchat)
    indexdf.netlist<-indexdf.netlist+1

    #Infer the cell-cell communication at a signaling pathway level
    cellchat <- computeCommunProbPathway(cellchat)

    #Calculate the aggregated cell-cell communication network
    #We can calculate the aggregated cell-cell communication network by counting the number of links or
    #summarizing the communication probability.
    #USER can also calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use.

    cellchat <- aggregateNet(cellchat)

    groupSize <- as.numeric(table(cellchat@idents))
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = paste(labelIndex, "Number of interactions",sep=" "))
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = paste(labelIndex, "Interaction weights/strength",sep=" "))
  }

  #comparing number of interactions between disease/control
  foldchangeInterMat<-c()
  mat <- cellchat@net$weight
  for(cellIDindex in rownames(mat)){
    print(cellIDindex)
    numsourceC<-length(grep(cellIDindex,df.net[[1]]$source))
    numtargetC<-length(grep(cellIDindex,df.net[[1]]$target))
    TIC<-nrow(df.net[[1]])
    totalnumINterC<-(numsourceC+numtargetC)/TIC
    numsourceD<-length(grep(cellIDindex,df.net[[2]]$source))
    numtargetD<-length(grep(cellIDindex,df.net[[2]]$target))
    TID<-nrow(df.net[[2]])
    totalnumINterD<-(numsourceD+numtargetD)/TID

    if(totalnumINterC==0){
      totalnumINterC<-1
    }
    foldChangeInter<-log2(totalnumINterD/totalnumINterC)
    print(foldChangeInter)
    foldchangeInterMat<-rbind(foldchangeInterMat,c(cellIDindex,foldChangeInter))
  }

  BulkfoldChangeInter<-log2(TID/TIC)
  foldchangeInterMat<-rbind(foldchangeInterMat,c("Bulk",BulkfoldChangeInter))
  foldchangeInterMat<-as.data.frame(foldchangeInterMat)
  foldchangeInterMat$V2<-as.numeric(foldchangeInterMat$V2)
  foldchangeInterMat <- foldchangeInterMat[order(foldchangeInterMat$V2,decreasing = T),]
  colnames(foldchangeInterMat)<-c("CellID","Fold Diff. Inter.")
  #write.table(df.netDonor,"Donordf.net",quote = F,row.names = F,sep = "\t")
  write.table(foldchangeInterMat,"CellChatPac.txt",quote = F,row.names = F,sep = "\t")
  return(foldchangeInterMat)
}
