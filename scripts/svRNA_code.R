rm(list = ls())
#GSE166635_RAW

suppressMessages(library(clusterProfiler))
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(magrittr))
suppressMessages(library(gtools))
suppressMessages(library(stringr))
suppressMessages(library(Matrix))
suppressMessages(library(tidyverse))
suppressMessages(library(patchwork))
suppressMessages(library(data.table))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggpubr))
suppressMessages(library(scales))
suppressMessages(library(ggsci))
suppressMessages(library(sctransform))
suppressMessages(library(harmony))
suppressMessages(library(tidydr))
suppressMessages(library(celldex))
suppressMessages(library(pheatmap))
suppressMessages(library(clustree))
suppressMessages(library(xlsx))
suppressMessages(library(gridExtra))
suppressMessages(library(ggthemes))
suppressMessages(library(ggnewscale))
suppressMessages(library(CellChat))
suppressMessages(library(ggpubr))
suppressMessages(library(patchwork))
suppressMessages(library(monocle))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(ggrepel))
options(stringsAsFactors = F)


##########
dir_name=list.files('00_origin_datas/GSE166635_RAW/')
datalist=list()
for (i in 1:length(dir_name)){
  dir.10x = paste0("00_origin_datas/GSE166635_RAW/",dir_name[i])
  list.files(dir.10x)
  my.data <- Read10X(data.dir = dir.10x)
  datalist[[i]] = CreateSeuratObject(counts = my.data, project = dir_name[i], min.cells = 3, min.features = 200)
  Samples1=dir_name[i]
  datalist[[i]] = AddMetaData(datalist[[i]] , Samples1,col.name = "Samples")
}
names(datalist)=dir_name
rm(my.data)


####
for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  #sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
  datalist[[i]] <- sce
  rm(sce)
}
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
rm(datalist)
unique(sce@meta.data$Samples)
head(sce@meta.data)
#sce  <- subset(sce,Samples %in% s_id1) 
############
saveRDS(sce,"00_origin_datas/sce.rds")
#sce <- readRDS("00_origin_datas/sce.rds")
scRNA = subset(sce, subset=nFeature_RNA>200&nFeature_RNA<8000&nCount_RNA <= 100000&percent.mt<15)
sum(table(scRNA@meta.data$Samples))
p2= VlnPlot(scRNA, features=c("nFeature_RNA", "nCount_RNA","percent.mt"),  pt.size=0)+NoLegend()+theme(text=element_text(family="Times"))+theme(axis.text.x = element_text(angle = 45, hjust = 1))
sc <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 1000)
#top10 <- head(VariableFeatures(sc), 10)
#plot1 <- VariableFeaturePlot(sc)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
sc <- ScaleData(sc, features = rownames(sc))
#scRNA = SCTransform(scRNA, vars.to.regress="percent.mt", verbose=FALSE)
###############################################################################
scRNA <- RunPCA(sc, features = VariableFeatures(object = sc))
sum(table(scRNA@meta.data$Samples))

DimPlot(scRNA , reduction = "pca" )
VizDimLoadings(scRNA , dims = 1:6, reduction = "pca")
DimHeatmap(scRNA , dims = 1:20, cells = 500, balanced = TRUE)
scRNA = RunHarmony(scRNA, group.by.vars="Samples")

scRNA_umap <- RunUMAP(scRNA, dims=1:20, reduction="harmony")
scRNA_umap <- FindNeighbors(scRNA_umap, dims=1:20, reduction="harmony")
scRNA_umap <- FindClusters(scRNA_umap  , resolution=0.1)
mydata<-scRNA_umap

#scRNA_tsne <- RunTSNE(scRNA, dims=1:15, reduction="harmony")
#scRNA_tsne <- FindNeighbors(scRNA_tsne, dims=1:15, reduction="harmony")
#scRNA_tsne <- FindClusters(scRNA_tsne  , resolution=0.1)
#mydata<-scRNA_tsne

######################################
p1 = DimPlot(mydata , reduction="umap", group.by="Samples", pt.size=1)+theme(legend.position="right", plot.title=element_blank())
#DimPlot(mydata , reduction="umap", label=T,group.by="seurat_clusters", pt.size=1)+theme(legend.position="right", plot.title=element_blank())
#DimPlot(scRNA_tsne , reduction="tsne", label=T,group.by="seurat_clusters", pt.size=1)+theme(legend.position="right", plot.title=element_blank())
FigureS1=mg_merge_plot(p2,p1,ncol=2,nrow=1,labels = c('A','B'),legend ="bottom",common.legend =T,widths = c(2.5,1))
ggsave('PDFs/Figure S1.pdf',FigureS1,height = 5,width = 12)
ggsave('PDFs/Figure S1.jpg',FigureS1,height = 5,width = 12)



colors=c("#E89DA0","#B2D3A4","#88CEE6","#F5D2A8","#80C1C4")
colors=c("#B383B9","#FCED82","#3C77AF","#D1352B","#8FA4AE","#AECDE1","#E89DA0","#F5D2A8","#BBDD78")

u_plot_clusters=DimPlot(mydata, reduction="umap",label = T,label.size = 8,group.by="seurat_clusters", pt.size=0.1)+
  theme(legend.position="right", plot.title=element_blank())+
  theme(text=element_text(family="Times"))+
  #scale_color_manual(values=colors)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed"))
u_plot_clusters

markers <- FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(markers, "scRNA/All_celltype_DEG.txt", col.names=T, row.names=T, quote=F, sep="\t")
#markers <- read.table("01_scRNA/All_celltype_DEG.txt",header = T,check.names = F,fill=T,sep = "\t")
markers["pct.diff"]=markers$pct.1-markers$pct.2

#################################################################################
###
Top5 <- markers %>% group_by(cluster) %>% slice_max(n =5, order_by =avg_logFC )
Top51 <- markers %>% group_by(cluster) %>% slice_max(n =5, order_by =pct.diff )
Cellmarker <- read.table
Cellmarker <- Cellmarker[,c(9,7)]
colnames(Top51)[7]="marker"
Cellmarker2 <- merge(Cellmarker,Top51,by="marker")
write.table(Cellmarker2, "scRNA/Cellmarker2.txt", col.names=T, row.names=F, quote=F, sep="\t")



diff.marker.dotplot= 
  DotPlot(object = mydata, features = as.vector(unique(Top51$marker)),
          # dot.scale =8,dot.min = 0,
          scale =T)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('')+ylab('')+coord_flip()+scale_color_gradientn(colors=c( "grey","white", "hotpink"))
diff.marker.dotplot
VlnPlot(mydata, features=c("CPA3","GATA2","TPSAB1","TPSB2"), pt.size=0)+NoLegend()+theme(axis.title.x=element_blank())
u_plot_clusters

#0 T cell:"CD2","CD3E","IL7R"
#1 Macrophage："AIF1","FCER1G"
#2 HPCs："EPCAM","KRT19","CD24"
#3 Macrophage："AIF1","FCER1G"
#4 CAFs:"COL1A2","BGN","ACTA2"
#5 B cell:"CD79A","MS4A1"
#6 HPCs："EPCAM","KRT19","CD24"
#7 Plasma cell:"DERL3","JCHAIN","MZB1"
#8 TECs:"PECAM1","ENG"
#9 Mast cell:"CPA3","GATA2","TPSAB1","TPSB2"



genes = unique(c(
  "CD2","CD3E","IL7R"
 ,"AIF1","FCER1G"
 ,"EPCAM","KRT19","CD24"
 ,"AIF1","FCER1G"
,"COL1A2","BGN","ACTA2"
 ,"CD79A","MS4A1"
 ,"EPCAM","KRT19","CD24"
 ,"DERL3","JCHAIN","MZB1"
,"PECAM1","ENG"
 ,"CPA3","GATA2","TPSAB1","TPSB2"
  
))


cell_label = c(
 "T cell",
  "Macrophage",
  "HPCs",
 "Macrophage",
  "CAFs",
  "B cell",
 "HPCs",
  "Plasma cell",
 "TECs",
 "Mast cell"
)
#################################################################################
##
names(cell_label) <- levels(mydata)
mydata1 <- RenameIdents(mydata, cell_label)
mydata1[["cell_type"]] = Idents(mydata1)

u_plot=DimPlot(mydata1, pt.size=0.2, label=T, label.size=3)+
  theme( plot.title=element_blank())+
  scale_color_manual(values=colors)+
  theme(text=element_text(family="Times"))+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed"))


diff.marker.dotplot1= DotPlot(object = mydata1, features = unique(genes),
                              dot.scale =6,
                              dot.min = 0,
                              scale =T)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size = 8)) +xlab('')+ylab('')+coord_flip()+scale_color_gradientn(colors=c( "snow", "#9998FF"))
diff.marker.dotplot1
cell.prop<-as.data.frame(prop.table(table(Idents(mydata1), as.vector(mydata1$Samples))))
colnames(cell.prop)<-c("Cell_type","Samples","proportion")
unique(cell.prop$Samples)

bili=ggplot(cell.prop,aes(Samples,proportion,fill=Cell_type))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=colors)+
  ggtitle("")+
  theme_bw()+
  theme(text=element_text(family="Times"), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(face="bold",angle=45, size=12, hjust=1), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank(), legend.position="right")


hub=c( "G6PD","AKR1B15", "S100A9",  "ADH4")

umap_sig2= DotPlot(object = mydata1, features =hub ,
                   dot.scale =6,
                   dot.min = 0,
                   scale =T)+
  RotatedAxis()+ ggtitle("Hub Genes")+
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('')+ylab('')+coord_flip()+scale_color_gradientn(colors=c( "snow", "blue"))

VlnPlot(mydata1, features=hub, pt.size=0)+NoLegend()+theme(axis.title.x=element_blank())

Figure1a=mg_merge_plot(u_plot_clusters,u_plot,ncol=2,nrow=1,labels = c('A','B'))
Figure1d=mg_merge_plot(diff.marker.dotplot1,bili,umap_sig2,ncol=3,nrow=1, legend = "right",labels = c("C",'D',"E"))

Figure1=mg_merge_plot(Figure1a,Figure1d,ncol=1,nrow=2, legend = "right")
ggsave('PDFs/Fig6.pdf',Figure1,height = 15,width = 16)
ggsave('PDFs/Fig6.jpg',Figure1,height = 15,width = 16)


