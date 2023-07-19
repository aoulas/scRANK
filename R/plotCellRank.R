# library(ggplot2)
# library(ggpubr)
# library(gridExtra)
#
# library(reshape2)
# library(tidyr)
# library(tidyverse)
# library(dbplyr)
#filename<-"/data/oulas/scRNA-Seq/autism/RanksMalacardsPlusWikiReactFinalPac.txt"
#title<-"Rankings proportions"
#plotRanks(filename)

plotCellChat <-function (filename,title="CellChat Rankings"){

    CellChat <- read.table(filename, sep = "\t",header=T, as.is=T, row.names=1)
    CellChat$CellID<-rownames(CellChat)
    indexBulk<-which(CellChat$CellID=="Bulk")
    colsbulk<-c(rep('black',length(CellChat$CellID)))
    colsbulk[indexBulk]<-"red"
    CellChat <- CellChat[order(CellChat$Fold.Diff..Inter.),]


    #CellChat$CellID<-gsub("Lymph.Node.cd1c.positive.Myeloid.Dendritic.Cell","Myeloid.Dendritic.Cell",CellChat$CellID)

    p<-ggplot(data = CellChat,aes(x =  reorder(CellID,Fold.Diff..Inter.,FUN = max,decreasing = TRUE), y=Fold.Diff..Inter.)) +

      geom_bar(stat="identity",position=position_dodge())+#,

      geom_hline(yintercept = 0,color = "red", size=1)+

      theme(

        plot.title = element_text(size=11,)
      ) +
      ggtitle(title) +
      ylab("LogFDI")+
      xlab("")
    p<-p+theme(plot.title = element_text(size = 15, face = "bold"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=15),
               axis.title=element_text(size=15,face="bold"),legend.text=element_text(size=15),legend.title=element_text(size=15),axis.text.y = element_text(color=colsbulk))
    p<-p + coord_flip()
    print(p)

}

plotRanks <-function (filename,title="Rankings"){

  Ranks <- read.table(filename, sep = "\t",header=T, as.is=T, row.names=1)#Allavelog2FC.txt
  Ranks<-round(Ranks,2)
  #Keep only informative parameters
  Ranks<-Ranks[vapply(Ranks, function(x) length(unique(x)) > 1, logical(1L))]

  CellIDS<-rownames(Ranks)
  Ranks<-as.data.frame(sapply(Ranks, rank))
  MeanRanks<-as.data.frame(apply(Ranks, 1,mean))
  rownames(MeanRanks)<-CellIDS
  MeanRanks<-cbind(CellIDS,MeanRanks)
  MeanRanks<-MeanRanks[order(MeanRanks$`apply(Ranks, 1, mean)`,decreasing = T),]
  #if(Disease=="MM"){
   # indexBulk<-15
  #}else{
  indexBulk<-which(MeanRanks$CellIDS=="Bulk")
  #}
  colsbulk<-c(rep('black',length(CellIDS)))
  colsbulk[indexBulk]<-"red"

  Ranks<-Ranks %>%
    gather(Parameter, Rank)

  Ranks<-cbind(CellIDS,Ranks)

  Ranks$CellIDS<-gsub("Lymph.Node.cd1c.positive.Myeloid.Dendritic.Cell","Myeloid.Dendritic.Cell",Ranks$CellIDS)

  p<-ggplot(data = Ranks,aes(x =  reorder(CellIDS,Rank,FUN = mean,decreasing = TRUE), y = Rank)) +#, color=grey
    geom_boxplot(outlier.colour = NA,color=colsbulk) +

    guides(color = "none")  +

    geom_dotplot(aes(color=Parameter,fill=Parameter),binaxis = "y", stackdir = "center",dotsize=0.3,position = position_dodge(),binwidth=0.5)+ # dotsize=0.3

    stat_summary(fun.y=mean, aes(shape="average Rank"),geom="point",size=2, color="red", fill="red") +
    scale_shape_manual("", values=c("average Rank"=24))+
    scale_y_continuous(breaks=seq(0, max(Ranks$Rank),2))+

  theme_bw()+
    theme(axis.text=element_text(size=15),axis.text.y = element_text(color=colsbulk),
          axis.title=element_text(size=15,face="bold"),legend.text=element_text(size=15),legend.title=element_text(size=15))+


    ggtitle(title) +#cellID_i
    xlab("")+
    ylab("Rank")

  p<-p + coord_flip()
  print(p)

}


plotTotalNumberDEGs <-function (filename,title="Total DEGs Rankings"){
  TotalNumberDEGs <- read.table(filename, sep = "\t",header=T, as.is=T, row.names=1)#Allavelog2FC.txt #TotalNumberDEGs.txt
  TotalNumberDEGs$CellID<-rownames(TotalNumberDEGs)
  p<-ggplot(data = TotalNumberDEGs,aes(x =  reorder(CellID,TotalNumberDEGs,FUN = max,decreasing = FALSE), y=TotalNumberDEGs)) +

    geom_bar(stat="identity",position=position_dodge())+#,

  theme(
    plot.title = element_text(size=11,)
  ) +
    ggtitle(title) +#cellID_i
    xlab("")
  p<-p+theme(plot.title = element_text(size = 15, face = "bold"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.text=element_text(size=15),
             axis.title=element_text(size=15,face="bold"),legend.text=element_text(size=15),legend.title=element_text(size=15))
  p<-p + coord_flip()
  print(p)
}

plotProportions <-function (filename,title="Cell Proportion Rankings"){
  proportions <- read.table(filename, sep = "\t",header=T, as.is=T, row.names=1)#Allavelog2FC.txt
  colnames(proportions)<-c("CellIDs","Condition","Proportion")

  proportions2<-proportions%>%
    arrange(CellIDs) %>%
    group_by(CellIDs) %>%
    mutate(rel_inc= abs(Proportion-lag(Proportion, default=first(Proportion))))

  #proportions$CellID<-rownames(proportions)
  p<-ggplot(data = proportions2,aes(x=reorder(CellIDs,rel_inc,FUN = max,decreasing = TRUE), y=Proportion, color=Condition ,fill=Condition )) +
    #geom_boxplot(position = position_dodge()) +
    #facet_wrap(~Var2)+
    geom_bar(stat="identity",position=position_dodge())+#,
    geom_segment(aes(x = 8, y = max(Proportion), xend = 10, yend = max(Proportion)),
                 arrow = arrow(length = unit(0.5, "cm"),ends = "last"))+
    #geom_text(aes(label=rel_inc), position=position_dodge(width=0.9), vjust=-0.25)+
    geom_label(aes(x=9,y=max(Proportion)-0.04,label="Dec. Diff. in Cell Prop."), colour = "white", fontface = "bold",vjust=-0.25,show.legend = FALSE)+
    #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    #geom_jitter(aes(col=Var2), size=1.5, alpha=0.9) +
    #stat_compare_means(label.y = max(proportions$Freq),label.x = 1.5,size=2,method="wilcox.test")+
    #theme_ipsum() +
    theme(
      #legend.position="none",
      plot.title = element_text(size=11,)
    ) +
    ggtitle(title) +#cellID_i
    xlab("")+ylab("Proportion of cells")
  p<-p+theme(plot.title = element_text(size = 15, face = "bold"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.text=element_text(size=15),
             axis.title=element_text(size=15,face="bold"),legend.text=element_text(size=15),legend.title=element_text(size=15))
  print(p)
  #pie <- p + coord_polar("y", start=0)
  #print(pie)
}

