runBasicAnalysis<-function(disease,path,annotate=TRUE,userlabel,usercelltype){
  library(enrichR)
  
  if (missing(disease)) cat("Argument disease is missing") else cat(paste("Argument disease =", disease));
  cat ("\n");
  if (missing(path)) cat("Argument path is missing") else cat(paste("Argument path =", path));
  cat("\n\n");
  
  if(annotate==TRUE){
    userlabel<-"label"
    usercelltype<-"celltype"
  }else{
    if (missing(userlabel)) {
      cat("If annotate is FALSE Argument userlabel must also be provided")
      cat ("\n");
      stop("Execution terminated")
    }else{
      cat(paste("Argument userlabel =", userlabel));
      cat ("\n");
    }
    if (missing(usercelltype)) {
      cat("If annotate is FALSE Argument usercelltype must also be provided")
      cat ("\n");
      stop("Execution terminated")
    }else{
      cat(paste("Argument usercelltype =", usercelltype));
      cat ("\n");
    }
  }
  
  
  setwd(path)
  subDir <- "Figures"
  if (!file.exists(subDir)){
    dir.create(file.path(path, subDir))
    print("Directory 'Figures' created")
  }
  
  print(paste("Starting Analysis for",disease,sep=" "))
  
  
  ##########################################Run analysis#############################################################################
  dirs <- list.dirs()
  
  loaded.dataSO.list<-c()
  
  plot.list<-vector("list", length(dirs)-1)
  index<-1
  
  
  #remove directory with html files
  indexhtmlfiles<-grep("*_files",dirs)
  if(length(indexhtmlfiles) !=0){
    dirs<-dirs[-indexhtmlfiles]
  }
  #start from 2 to avoid home directory
  #options(Seurat.object.assay.version = "v3")
  for(i in 2:length(dirs)){
    data_dir<-dirs[i]
    
    filesindir<-list.files(data_dir)
    #check if dir contains seurat raw data using different formats (rds, H5, ect)
    foundfalse<-which(c("barcodes.tsv","genes.tsv","matrix.mtx") %in% filesindir==FALSE)
    
    foundH5<-grep("\\.h5",filesindir)
    foundrds<-grep("\\.rds",filesindir)
    foundmeta<-"meta.txt"%in% filesindir
    if(length(foundfalse)==0){
      # Load the dataset
      loaded.data <- Read10X(data.dir = data_dir)
      project_name<-gsub("[\\.\\/]","",data_dir,)
      print(project_name)
      # Initialize the Seurat object with the raw (non-normalized data).
      if(foundmeta){
        meta <- read.table(paste(project_name,"/meta.txt",sep = ""), header=T, sep="\t", as.is=T, row.names=1)
        loaded.dataSO <- CreateSeuratObject(counts = loaded.data, project = project_name, min.cells = 3, min.features = 200,meta.data=meta)
      }else{
        loaded.dataSO <- CreateSeuratObject(counts = loaded.data, project = project_name, min.cells = 3, min.features = 200)
      }
      #QC and selecting cells for further analysis
      # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
      loaded.dataSO[["percent.mt"]] <- PercentageFeatureSet(loaded.dataSO, pattern = "^MT-")
      # Visualize QC metrics as a violin plot
      # VlnPlot(loaded.dataSO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
      # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
      # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
      
      #plot1 <- FeatureScatter(loaded.dataSO, feature1 = "nCount_RNA", feature2 = "percent.mt")
      plot2 <- FeatureScatter(loaded.dataSO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
      
      plot.list[[index]]<-plot2
      index<-index+1
      
      #Filtering data
      #loaded.dataSO <- subset(loaded.dataSO, subset = nFeature_RNA > 100 & nFeature_RNA < 6500 & percent.mt < 75)
      
      #Normalizing the data. Normalized values are stored in UnWounded1[["RNA"]]@data. *
      loaded.dataSO <- NormalizeData(loaded.dataSO)#, normalization.method = "LogNormalize", scale.factor = 10000)
      #Identification of highly variable features (feature selection) *
      loaded.dataSO <- FindVariableFeatures(loaded.dataSO, selection.method = "vst", nfeatures = 2000) 
      
      loaded.dataSO.list<-c(loaded.dataSO.list,loaded.dataSO)
    }else if(length(foundH5) !=0){
      loaded.data <- Read10X_h5(paste(data_dir,"/",filesindir,sep=""), use.names = TRUE, unique.features = TRUE)
      project_name<-gsub("[\\.\\/]","",data_dir,)
      print(project_name)
      # Initialize the Seurat object with the raw (non-normalized data).
      if(foundmeta){
        meta <- read.table(paste(project_name,"/meta.txt",sep = ""), header=T, sep="\t", as.is=T, row.names=1)
        loaded.dataSO <- CreateSeuratObject(counts = loaded.data, project = project_name, min.cells = 3, min.features = 200,meta.data=meta)
      }else{
        loaded.dataSO <- CreateSeuratObject(counts = loaded.data, project = project_name, min.cells = 3, min.features = 200)
      }
      #QC and selecting cells for further analysis
      # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
      loaded.dataSO[["percent.mt"]] <- PercentageFeatureSet(loaded.dataSO, pattern = "^MT-")
      # Visualize QC metrics as a violin plot
      # VlnPlot(loaded.dataSO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
      # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
      # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
      
      #plot1 <- FeatureScatter(loaded.dataSO, feature1 = "nCount_RNA", feature2 = "percent.mt")
      plot2 <- FeatureScatter(loaded.dataSO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
      
      plot.list[[index]]<-plot2
      index<-index+1
      
      #Filtering data
      #loaded.dataSO <- subset(loaded.dataSO, subset = nFeature_RNA > 100 & nFeature_RNA < 6500 & percent.mt < 75) #changed this
      
      #Normalizing the data. 
      loaded.dataSO <- NormalizeData(loaded.dataSO)#, normalization.method = "LogNormalize", scale.factor = 10000)
      #Identification of highly variable features (feature selection) *
      loaded.dataSO <- FindVariableFeatures(loaded.dataSO, selection.method = "vst", nfeatures = 2000)
      loaded.dataSO.list<-c(loaded.dataSO.list,loaded.dataSO)
    }else if(length(foundrds) !=0){
      
      loaded.dataSO<- readRDS(paste(data_dir,"/",filesindir[foundrds],sep=""))
      project_name<-gsub("[\\.\\/]","",data_dir,)
      print(project_name)
      # Initialize the Seurat object with the raw (non-normalized data).
      # if(foundmeta){
      #   meta <- read.table(paste(project_name,"/meta.txt",sep = ""), header=T, sep="\t", as.is=T, row.names=1)
      #   loaded.dataSO <- CreateSeuratObject(counts = loaded.data, project = project_name, min.cells = 3, min.features = 200,meta.data=meta)
      # }else{
      #   loaded.dataSO <- CreateSeuratObject(counts = loaded.data, project = project_name, min.cells = 3, min.features = 200)
      # }
      #QC and selecting cells for further analysis
      # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
      # loaded.dataSO[["percent.mt"]] <- PercentageFeatureSet(loaded.dataSO, pattern = "^MT-")
      # Visualize QC metrics as a violin plot
      # VlnPlot(loaded.dataSO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
      # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
      # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
      
      #plot1 <- FeatureScatter(loaded.dataSO, feature1 = "nCount_RNA", feature2 = "percent.mt")
      plot2 <- FeatureScatter(loaded.dataSO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
      
      plot.list[[index]]<-plot2
      index<-index+1
      
      #Filtering data
      #loaded.dataSO <- subset(loaded.dataSO, subset = nFeature_RNA > 100 & nFeature_RNA < 6500 & percent.mt < 75) #changed this
      
      DefaultAssay(loaded.dataSO)<-"RNA"
      
      #Normalizing the data. *
      #loaded.dataSO <- NormalizeData(loaded.dataSO)#, normalization.method = "LogNormalize", scale.factor = 10000)
      #Identification of highly variable features (feature selection) *
      loaded.dataSO <- FindVariableFeatures(loaded.dataSO, selection.method = "vst", nfeatures = 2000)
      loaded.dataSO.list<-c(loaded.dataSO.list,loaded.dataSO)
    }
  }
  
  #increase memory size to hold larger objects
  options(future.globals.maxSize = 8000 * 1024^2)
  if(index==2){
    loaded.dataSO.combined <-loaded.dataSO
    jpeg(file=paste(subDir,"/FEATURE_PLOT_",disease,".jpg",sep=""),
         width=1000, height=800,res=100)
    plot(plot2)
    dev.off()
  }else{
    #loaded.dataSO.combined <-loaded.dataSO
    jpeg(file=paste(subDir,"/FEATURE_PLOT_",disease,".jpg",sep=""),
         width=1000, height=800,res=100)
    do.call("grid.arrange", c(plot.list))#, ncol = (index-1))
    dev.off()
    features <- SelectIntegrationFeatures(object.list = loaded.dataSO.list)
    loaded.dataSO.anchors <- FindIntegrationAnchors(object.list = loaded.dataSO.list, dims = 1:20)#,reference = 1,k.filter = NA,anchor.features = features,verbose = TRUE
    #options(Seurat.object.assay.version = "v4")
    loaded.dataSO.combined <- IntegrateData(anchorset = loaded.dataSO.anchors, dims = 1:20)
    
    DefaultAssay(loaded.dataSO.combined) <- "integrated"
  }
  
  # Run the standard workflow for visualization and clustering
  loaded.dataSO.combined <- ScaleData(loaded.dataSO.combined, verbose = FALSE)
  loaded.dataSO.combined <- RunPCA(loaded.dataSO.combined, npcs = 30, verbose = FALSE)
  # t-SNE and Clustering
  loaded.dataSO.combined <- RunUMAP(loaded.dataSO.combined, reduction = "pca", dims = 1:20)
  #loaded.dataSO.combined <- RunTSNE(loaded.dataSO.combined,reduction = "pca",dims = 1:10)
  #Louvein method or improved Leiden
  loaded.dataSO.combined <- FindNeighbors(loaded.dataSO.combined, reduction = "pca", dims = 1:20)
  loaded.dataSO.combined <- FindClusters(loaded.dataSO.combined, resolution  = 0.3)#was at 0.3
  #loaded.dataSO.combined[["RNA"]] <- as(object = loaded.dataSO.combined[["RNA"]], Class = "Assay")
  #Can view object
  #loaded.dataSO.combined[[]]
  if(annotate==FALSE){
    p1combined <- DimPlot(loaded.dataSO.combined, reduction = "umap", group.by = userlabel)
  }
  p2combined <- DimPlot(loaded.dataSO.combined, reduction = "umap", label = TRUE)
  if(annotate==FALSE){
    jpeg(file=paste(subDir,"/UMAP_ANNOT_",disease,".jpg",sep=""),
         width=1000, height=800,res=100)
    plot(p1combined+p2combined)
    dev.off()
  }else{
    jpeg(file=paste(subDir,"/UMAP_",disease,".jpg",sep=""),
         width=1000, height=800,res=100)
    plot(p2combined)
    dev.off()
  }
  
  if(annotate==FALSE){
    jpeg(file=paste(subDir,"/UMAP_ANNOT2_",disease,".jpg",sep=""),
         width=1000, height=800,res=100)
    plot(DimPlot(loaded.dataSO.combined, reduction = "umap", split.by = userlabel,group.by = usercelltype,raster=FALSE))
    dev.off()
  }
  
  if(annotate==TRUE){
    # Identify conserved cell type markers
    # For performing differential expression after integration, we switch back to the original
    # data
    DefaultAssay(loaded.dataSO.combined) <- "RNA"
    loaded.dataSO.combined[["RNA"]]<-JoinLayers(loaded.dataSO.combined[["RNA"]])
    # loaded.dataSO.combined <- IntegrateLayers(object = loaded.dataSO.combined, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
    #                         verbose = FALSE)
    #loaded.dataSO.combined[["RNA"]] <- as(object = loaded.dataSO.combined[["RNA"]], Class = "Assay")
    loaded.dataSO.combined.markers <- FindAllMarkers(loaded.dataSO.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)#,test.use = "wilcox_limma"
    
    loaded.dataSO.combined.markerstop1<-loaded.dataSO.combined.markers %>%
      group_by(cluster) %>%
      slice_max(n = 1, order_by = avg_log2FC) #100 gives best results
    
    loaded.dataSO.combined.markerstop2<-loaded.dataSO.combined.markers %>%
      group_by(cluster) %>%
      slice_max(n = 2, order_by = avg_log2FC) #100 gives best results
    
    loaded.dataSO.combined.markerstop100<-loaded.dataSO.combined.markers %>%
      group_by(cluster) %>%
      slice_max(n = 100, order_by = avg_log2FC) #100  gives best results
    
    loaded.dataSO.combined.markerstop200<-loaded.dataSO.combined.markers %>%
      group_by(cluster) %>%
      slice_max(n = 200, order_by = avg_log2FC) #200 gives second best results
    
    loaded.dataSO.combined.markerstop6<-loaded.dataSO.combined.markers %>%
      group_by(cluster) %>%
      slice_max(n = 6, order_by = avg_log2FC) #100 gives best results
    
    genesAll<-as.data.frame(loaded.dataSO.combined.markerstop6$gene)
    unique(loaded.dataSO.combined$seurat_clusters)
    
    #Visualize top 4 markers
    jpeg(file=paste(subDir,"/TOP4-MARKERS_",disease,".jpg",sep=""),
         width=1000, height=800,res=100)
    plot(FeaturePlot(loaded.dataSO.combined, features = loaded.dataSO.combined.markerstop1$gene[c(1,2,3,4)]))
    dev.off()
    
    #Visualize volcano plot top 1 markers
    jpeg(file=paste(subDir,"/TOP1-MARKERS_VOL_",disease,".jpg",sep=""),
         width=1000, height=800,res=100)
    plot(VlnPlot(loaded.dataSO.combined, features = loaded.dataSO.combined.markerstop1$gene[1]))
    dev.off()
    #FeaturePlot(loaded.dataSO.combined, features = loaded.dataSO.combined.markerstop1$gene[1],label = T)& theme(legend.position = c(0.1,0.2))
    
    #loaded.dataSO.combined <- ScaleData(loaded.dataSO.combined, verbose = FALSE) #this needs to be removed for JoinLayers
    
    jpeg(file=paste(subDir,"/TOP2-MARKERS_HEAT_",disease,".jpg",sep=""),
         width=1000, height=800,res=100)
    plot(DoHeatmap(loaded.dataSO.combined, features = loaded.dataSO.combined.markerstop2$gene) + NoLegend())
    dev.off()
    #Run Enrich R on top 100 Markers
    listEnrichrSites()
    setEnrichrSite("Enrichr") # Human genes
    websiteLive <- TRUE
    #dbs <- listEnrichrDbs()
    #if (is.null(dbs)) websiteLive <- FALSE
    #f (websiteLive) head(dbs)
    annotatedclusters<-c()
    listofresults<-c()
    dbstouse <- c("PanglaoDB_Augmented_2021","CellMarker_Augmented_2021","Tabula_Sapiens")
    for(m in as.numeric(as.character(unique(loaded.dataSO.combined.markerstop100$cluster)))){
      markers_to_search<-loaded.dataSO.combined.markerstop100$gene[loaded.dataSO.combined.markerstop100$cluster==m]
      if (websiteLive) {
        enriched <- enrichr(markers_to_search, dbstouse)
      }
      celltype<-c(enriched[[1]]$Term[1],enriched[[2]]$Term[1],enriched[[3]]$Term[1])
      print(paste("Cluster",m,celltype,sep = " "))
      annotatedclusters<-c(annotatedclusters,toString(celltype[3]))#changed to Tabula Sapiens changed to Panglao for MM
      names(enriched)<-paste(names(enriched),m,sep = "_")
      listofresults<-c(listofresults,enriched)
      
      jpeg(file=paste(subDir,"/CELL_TYPES_CLUSTER_",m,"_",disease,".jpg",sep=""),
           width=1200, height=800)
      if (websiteLive) plot(plotEnrich(enriched[[paste("PanglaoDB_Augmented_2021_",m,sep="")]],  numChar = 40, y = "Count", orderBy = "P.value",title = paste("PanglaoDB",m))+plotEnrich(enriched[[paste("CellMarker_Augmented_2021_",m,sep="")]], numChar = 40, y = "Count", orderBy = "P.value",title = paste("CellMarker",m))+plotEnrich(enriched[[paste("Tabula_Sapiens_",m,sep="")]], numChar = 40, y = "Count", orderBy = "P.value",title = paste("Tabula_Sapiens",m)))
      dev.off()
    }
    
    #Rename indents
    new.cluster.ids <- make.names(annotatedclusters,unique = T)
    names(new.cluster.ids) <- levels(loaded.dataSO.combined)
    loaded.dataSO.combined <- RenameIdents(loaded.dataSO.combined, new.cluster.ids)
    jpeg(file=paste(subDir,"/UMAP_ANNOT3_",disease,".jpg",sep=""),
         width=1000, height=800,res=100)
    plot(DimPlot(loaded.dataSO.combined, reduction = "umap",label = TRUE))
    dev.off()
    #Change the Labels #issue here if no number after label recheck!!!!
    loaded.dataSO.combined$celltype.label <- paste(Idents(loaded.dataSO.combined), gsub("[0-9]+","",loaded.dataSO.combined$orig.ident), sep = "_")
    loaded.dataSO.combined$label <- gsub("[0-9]+","",loaded.dataSO.combined$orig.ident)
    loaded.dataSO.combined$celltype <- Idents(loaded.dataSO.combined)
    
    jpeg(file=paste(subDir,"/UMAP_ANNOT3_SPLIT_",disease,".jpg",sep=""),
         width=1000, height=800,res=100)  
    plot(DimPlot(loaded.dataSO.combined, reduction = "umap", split.by = "label",label = TRUE))
    dev.off()
    
    n <- length(unique(loaded.dataSO.combined$label))
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    cols_to_use<-sample(col_vector, n)
    markers.to.plot <- unique(loaded.dataSO.combined.markerstop2$gene)
    
    jpeg(file=paste(subDir,"/TOP2_MARKERS_DOT_",disease,".jpg",sep=""),
         width=1000, height=800,res=100)
    plot(DotPlot(loaded.dataSO.combined, features = markers.to.plot, cols = cols_to_use, dot.scale = 8, split.by = "label")) +
      RotatedAxis()
    dev.off()
    jpeg(file=paste(subDir,"/TOP1_MARKERS_DOT_",disease,".jpg",sep=""),
         width=1000, height=800,res=100)
    plot(DotPlot(loaded.dataSO.combined, features = c(loaded.dataSO.combined.markerstop1$gene[1]), cols = cols_to_use, dot.scale = 8, split.by = "label") )+
      RotatedAxis()
    dev.off()
    jpeg(file=paste(subDir,"/TOP1_MARKERS_VOL_",disease,".jpg",sep=""),
         width=1000, height=800,res=100)
    plot(VlnPlot(loaded.dataSO.combined, features = c(loaded.dataSO.combined.markerstop1$gene[1]),split.by = "label",split.plot = TRUE))
    dev.off()
  }
  
  
  # What proportion of cells are in each cluster?
  jpeg(file=paste(subDir,"/ALL_CELLS_PROP_",disease,".jpg",sep=""),
       width=1000, height=800,res=100)
  pie(prop.table(table(loaded.dataSO.combined[[usercelltype]][,1])),main = "Proportion of cells per cell-type")
  dev.off()
  #prop.table(table(Idents(loaded.dataSO.combined))) #same result
  
  # How does cluster membership vary by condition?
  tableofcellconts<-table(loaded.dataSO.combined[[usercelltype]][,1], loaded.dataSO.combined[[userlabel]][,1])
  
  # What proportion of cells are in each cluster by condition
  proportions<-as.data.frame(prop.table(table(loaded.dataSO.combined[[usercelltype]][,1], loaded.dataSO.combined[[userlabel]][,1]), margin = 2))
  #proportions1<-as.data.frame(prop.table(table(loaded.dataSO.combined$celltype, loaded.dataSO.combined$label), margin = 2))
  
  write.table(proportions,"proportionsFinalPac.txt",quote = F,row.names = T,sep = "\t")
  
  #listofoutput<-list(loaded.dataSO.combined,termsKEGG,termsGO,termsMSIG,termsWiki,termsReact,termsMOA)
  return(loaded.dataSO.combined)
}