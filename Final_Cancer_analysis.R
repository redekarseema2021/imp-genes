install.packages(c("readr","tidyr","RColorBrewer","ggpubr","dplyr"))
install.packages("msu","infotheo","data.table","igraph")



library('msu')
library("dplyr")
library('infotheo')
library('data.table')
library("igraph")
library("readr")
library("tidyr")
library("RColorBrewer")

#set the path here according to your directory
setwd("C:\\Cancer")
cdf=fread('Significant genes for network construction.csv')


# Function to build the data frame
# parameters involve
# cdf=The genome Dataset
# option=Choice of relation.Possible choices are('info_gain','spearman','pearson')
# returns a reduced dataset with the format
#            from    to   choice

df_builder<-function(cdf,option){
  if (option=='info_gain'){
    info_gain<-c()
    for (x in 1:ncol(cdf)){
      ans<-c()
      for (y in 1:ncol(cdf)){
        ans<-c(ans,mutinformation(discretize(cdf[,x,with=FALSE],),discretize(cdf[,y,with=FALSE])))
      }
      info_gain<-c(info_gain,ans)
    }
    
    info_gain<-array(info_gain,dim=c(ncol(cdf),ncol(cdf)))
    
    
    
    info_gain<-matrix(info_gain,nrow=ncol(cdf),ncol=ncol(cdf))
    rownames(info_gain)=colnames(cdf)
    colnames(info_gain)=colnames(cdf)
    info_gain
    
    for (i in 1:ncol(info_gain)){
      info_gain[i,i]<-NA
    }
    
    emp<-data.frame(matrix(vector(), 0, 3,
                           dimnames=list(c(), c("from", "to", "ig"))),
                    stringsAsFactors=F)
    ind_ig<-data.frame(matrix(vector(), 0, 3,
                              dimnames=list(c(), c("from", "to", "ig"))),
                       stringsAsFactors=F)
    for (x in 1:155){
      f<-as.vector(info_gain[x,])
      ind_ig<-rbind(ind_ig,transpose(rbind(emp,from=rep(rownames(info_gain)[x],(156-x)),to=colnames(info_gain)[(x+1):156],ig=f[(x+1):156])))
      
    }
    colnames(ind_ig)<-c('from','to','ig')
    ind_ig$ig<-as.numeric(ind_ig$ig)
    thresh_ig<-arrange(ind_ig,desc(ig))[1:100,]
    return (thresh_ig)
  }
  else if(option=='spearman')
  {
    spearman_corr<-cor(cdf,cdf,method=c("spearman"))
    
    
    for (i in 1:ncol(spearman_corr)){
      spearman_corr[i,i]<-NA
    }
    
    emp<-data.frame(matrix(vector(), 0, 3,
                           dimnames=list(c(), c("from", "to", "s_corr"))),
                    stringsAsFactors=F)
    ind_sc<-data.frame(matrix(vector(), 0, 3,
                              dimnames=list(c(), c("from", "to", "s_corr"))),
                       stringsAsFactors=F)
    for (x in 1:155){
      f<-as.vector(spearman_corr[x,])
      ind_sc<-rbind(ind_sc,transpose(rbind(emp,from=rep(rownames(spearman_corr)[x],(156-x)),to=colnames(spearman_corr)[(x+1):156],s_corr=f[(x+1):156])))
      
    }
    colnames(ind_sc)<-c('from','to','s_corr')
    ind_sc$s_corr<-as.numeric(ind_sc$s_corr)
    thresh_sc<-arrange(ind_sc,desc(s_corr))[1:100,]
    return (thresh_sc)
  }
  else if(option=='pearson')
  {
    pearson_corr<-cor(cdf,cdf,method=c("pearson"))
    
    
    for (i in 1:ncol(pearson_corr)){
      pearson_corr[i,i]<-NA
    }
    
    emp<-data.frame(matrix(vector(), 0, 3,
                           dimnames=list(c(), c("from", "to", "p_corr"))),
                    stringsAsFactors=F)
    ind_sc<-data.frame(matrix(vector(), 0, 3,
                              dimnames=list(c(), c("from", "to", "p_corr"))),
                       stringsAsFactors=F)
    for (x in 1:155){
      f<-as.vector(pearson_corr[x,])
      ind_sc<-rbind(ind_sc,transpose(rbind(emp,from=rep(rownames(pearson_corr)[x],(156-x)),to=colnames(pearson_corr)[(x+1):156],p_corr=f[(x+1):156])))
      
    }
    colnames(ind_sc)<-c('from','to','p_corr')
    ind_sc$p_corr<-as.numeric(ind_sc$p_corr)
    thresh_sc<-arrange(ind_sc,desc(p_corr))[1:100,]
    return (thresh_sc)
  }
}


#Function for plotting the graph
#Arguements include
# df=the dataframe obtained from df_builder function
# metric=the relationship between graphs;i.e('ig','s_corr','p_corr') 
total_plot<-function(df,metric){
  sc_graph<- graph_from_data_frame(df, directed = FALSE)
  sc_deg<-degree(sc_graph,mode=c("All"))
  V(sc_graph)$degree<-sc_deg
  V(sc_graph)$evc<-eigen_centrality(sc_graph)$vector
  V(sc_graph)$bt<-betweenness(sc_graph)
  V_df<-data.frame(list(Vertex=as_ids(V(sc_graph)),Degree=V(sc_graph)$degree),stringsAsFactors = FALSE)
  V_df$EVC<-V(sc_graph)$evc
  V_df$bt<-V(sc_graph)$bt
  V_df$color<-rep('red',dim(V_df)[1])
  
  #This code here selects only the top n(in this case 20) vertices for each of the following 
  #1.Degree Centrality
  #2.EigenVector Centrality
  #3.Betweeness
  deg<-order(V_df$Degree,decreasing = TRUE)[1:20]
  eig<-order(V_df$EVC,decreasing = TRUE)[1:20]
  btw<-order(V_df$bt,decreasing = TRUE)[1:20]
  
  #Use the below code for color coding
  tryCatch({V_df[intersect(btw,intersect(deg,eig)),]$color<-'yellow'},error=function(cond){})
  tryCatch({V_df[setdiff(intersect(deg,eig),btw),]$color<-'orange'},error=function(cond){})
  tryCatch({V_df[setdiff(intersect(deg,btw),eig),]$color<-'pink'},error=function(cond){})
  tryCatch({V_df[setdiff(intersect(eig,btw),deg),]$color<-'gray'},error=function(cond){})
  tryCatch({V_df[setdiff(deg,union(eig,btw)),]$color<-'green'},error=function(cond){})
  tryCatch({V_df[setdiff(eig,union(deg,btw)),]$color<-'purple'},error=function(cond){})
  tryCatch({V_df[setdiff(btw,union(deg,eig)),]$color<-'brown'},error=function(cond){})
  
  V(sc_graph)$color<-V_df$color
  
  
  if (metric=='ig')
  {
    E(sc_graph)$weight<-E(sc_graph)$ig
  }
  else if(metric=='s_corr'){
    E(sc_graph)$weight<-E(sc_graph)$s_corr
  }
  else if(metric=='p_corr'){
    E(sc_graph)$weight<-E(sc_graph)$p_corr
  }
  
  #plotting the  i graph
  set.seed(1001)
  library(RColorBrewer) # This is the color library
  pal<-brewer.pal(length(unique(V(sc_graph))), "Set3")
  par(mar=c(1,1,1,1))
  plot(sc_graph,edge.color = 'black',vertex.color=(sc_graph)$color,vertex.label.cex =0.8,vertex.size=8
       ,edge.width=sqrt(E(sc_graph)$weight/800),
       layout = layout.sphere,asp=0.4)
}


#############################################################################################









df<-df_builder(cdf=cdf,'spearman')
total_plot(df,'s_corr')

