library(Seurat)
library(SoupX)

#define SoupX function
#tod: raw matrix without filter
#toc:filtered matrix
#rho:cotamination fraction
run_soupx<-function(toc,tod,rho=NULL){
  #according to SoupX protocol, Simple clustering is necessary for providing cluster information to SoupX
  all<-toc
  all<-CreateSeuratObject(all)
  all<-NormalizeData(all,normalization.method="LogNormalize",scale.factor=10000)
  all<-FindVariableFeatures(all,selection.method="vst",nfeatures=3000)
  all.genes<-rownames(all)
  all<-ScaleData(all,features=all.genes)
  all<-RunPCA(all,features=VariableFeatures(all),npcs=40,verbose=F)
  all<-FindNeighbors(all,dims=1:30)
  all<-FindClusters(all,resolution=1.0)
  all<-RunUMAP(all,dims=1:30)
  #get cluster information after clustering
  matx<-all@meta.data
  sc=SoupChannel(tod,toc)
  sc=setClusters(sc,setNames(matx$seurat_clusters,rownames(matx)))
  #calculate contamination fraction automaticaly
  if (is.null(rho)) {
    tryCatch(
      {sc=autoEstCont(sc)},
      error=function(e){
        print("autoEstCont Error !")
        sc=setContaminationFraction(sc,0.2)}
    )
  }
  else{ #if there is rho value, use set rho
    sc=setContaminationFraction(sc,rho)
  }
  #adjust scRNA expression matrix
  out=adjustCounts(sc)
  #save two matrix 
  saveRDS(sc,"sc.rds")
  #save adjusted matrix and output it as 10X format
  DropletUtils::write10xCounts("./soupX_matrix",out,version="3")
}

##S14
toc_S14 <-Read10X("F:/神经退行斑马鱼脑基因代谢分析/2024_10_18结果/4.SoupX去除环境污染/data/S14/filtered_feature_bc_matrix")
tod_S14 <- Read10X("F:/神经退行斑马鱼脑基因代谢分析/2024_10_18结果/4.SoupX去除环境污染/data/S14/raw_feature_bc_matrix")
tod_S14<-tod_S14[rownames(toc_S14),]
setwd("F:/神经退行斑马鱼脑基因代谢分析/2024_10_18结果/4.SoupX去除环境污染/data/S14")
run_soupx(toc=toc_S14, tod=tod_S14)
#S18
toc_S18 <-Read10X("F:/神经退行斑马鱼脑基因代谢分析/2024_10_18结果/4.SoupX去除环境污染/data/S18/filtered_feature_bc_matrix")
tod_S18 <- Read10X("F:/神经退行斑马鱼脑基因代谢分析/2024_10_18结果/4.SoupX去除环境污染/data/S18/raw_feature_bc_matrix")
tod_S18<-tod_S18[rownames(toc_S18),]
setwd("F:/神经退行斑马鱼脑基因代谢分析/2024_10_18结果/4.SoupX去除环境污染/data/S18")
run_soupx(toc=toc_S18, tod=tod_S18)
#S30
toc_S30 <-Read10X("F:/神经退行斑马鱼脑基因代谢分析/2024_10_18结果/4.SoupX去除环境污染/data/S30/filtered_feature_bc_matrix")
tod_S30 <- Read10X("F:/神经退行斑马鱼脑基因代谢分析/2024_10_18结果/4.SoupX去除环境污染/data/S30/raw_feature_bc_matrix")
tod_S30<-tod_S30[rownames(toc_S30),]
setwd("F:/神经退行斑马鱼脑基因代谢分析/2024_10_18结果/4.SoupX去除环境污染/data/S30")
run_soupx(toc=toc_S30, tod=tod_S30)
