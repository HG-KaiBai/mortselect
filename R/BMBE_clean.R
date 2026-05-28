#' Aggregate plot
#'
#' @param AggregateLogRate.reshape Sum of log rate for certain series.
#' @param Seires Female/Male/Total
#'
#' @return The plot of aggregatelogRate
#' @import ggplot2
#'
Plot_aggregatelog = function(AggregateLogRate.reshape, Series){
  Gplot1=ggplot(AggregateLogRate.reshape,aes(x=t,y=AggregateLogRates,group=population))+
    geom_line(lwd=1.2)+ggtitle(paste0("Series: ",Series))
  pdf(file=paste0(Series,"_AggregateLogRate.pdf"))
  Gplot1
  dev.off()
}


#' cluster aggregate plot
#'
#' @param AggregateLogRate.reshape Sum of log rate for certain series.
#' @param Seires Female/Male/Total
#'
#' @return The plot of aggregatelogRate
#' @import ggplot2
#'
Plot_cluster_aggregatelog = function(AggregateLogRate.reshape2, Series){
  Gplot1=ggplot(AggregateLogRate.reshape2,aes(x=t,y=AggregateLogRates,group=population,color=as.factor(cluster)))+
    geom_line(lwd=1.2)+ggtitle(paste0("Series: ",Series))
  pdf(file=paste0(Series,"_cluster_AggregateLogRate.pdf"))
  Gplot1
  dev.off()
}


#' MDS cluster plot
#'
#' @param MDSresult MDS fit on data
#' @param Series Female/Male/Total
#'
#' @return MDS plot file
#' @import ggplot2
#' @import ggrepel
#'
Plot_MDS = function(MDSresult, Series){
  Gplot2=ggplot(MDSresult,aes(x=x,y=y,color=cluster))+geom_point()+
    ggrepel::geom_label_repel(aes(label=population,color=cluster),check_overlap = FALSE)+ggtitle(paste0("Series: ",Series))
  pdf(file=paste0(Series,"_MDSresult.pdf"))
  Gplot2
  dev.off()
}


#' Group the data with MDS method
#'
#' @param Series male/female
#' @param Dataset_training traingin set
#' @import reshape2
#' @import MASS
#'
#'
BMBE_clean = function(Series, Dataset_training){
  AggregateLogRate=matrix(nrow=length(Dataset_training),
                          ncol=length(Dataset_training[[1]]$year))
  for(i in 1:length(Dataset_training)){
    temp=Dataset_training[[i]]$rate[[which(names(Dataset_training[[i]]$rate)==Series)]]
    temp1=colSums(log(temp))
    AggregateLogRate[i,]=temp1
  }
  rownames(AggregateLogRate)=names(Dataset_training)
  colnames(AggregateLogRate)=Dataset_training[[1]]$year
  ## Plot for AggregateLogRate
  AggregateLogRate.reshape=(reshape2::melt(AggregateLogRate, value.name = "AggregateLogRates",
                                           varnames = c("population", "t")))
  # Save plot
  #Plot_aggregatelog(AggregateLogRate.reshape, Series)



  ## MDS and Clustering
  DisimilarityMatrix=matrix(0,nrow=length(Dataset_training),ncol=length(Dataset_training))
  rownames(DisimilarityMatrix)=colnames(DisimilarityMatrix)=names(Dataset_training)
  for (i in 1:(length(Dataset_training)-1))
  {
    for(j in (i+1):length(Dataset_training))
      DisimilarityMatrix[i,j]=DisimilarityMatrix[j,i]=var(AggregateLogRate[i,]-AggregateLogRate[j,])
  }


  # MDS
  fit <- isoMDS(DisimilarityMatrix, k=2) # k is the number of dim
  MDSresult=data.frame(population=names(Dataset_training),x=fit$points[,1],y=fit$points[,2])
  #Cluster
  Cluster=kmeans(fit$points,2,nstart=1000)$cluster
  MDSresult$cluster=as.factor(Cluster)
  #Plot_MDS(MDSresult, Series)


  AggregateLogRate.reshape2=AggregateLogRate.reshape
  AggregateLogRate.reshape2$cluster=rep(Cluster,ncol(AggregateLogRate))
  #Plot_cluster_aggregatelog(AggregateLogRate.reshape2, Series)
}
