# Some helper functions



procrustes_r <- function(A, B, normalize=F ){
  # center and normalize A 
  # A.centered <- t(scale(t(A), center = TRUE, scale = FALSE))
  A.size <- norm(A, type = "F") 
  A.normalized <- A
  # 
  # # center and normalize B
  # B.centered <- t(scale(t(B), center = TRUE, scale = FALSE))
  B.size <- norm(B, type = "F")
  if (normalize == T){
    B.normalized <- B /B.size
  } else {B.normalized <- B }
  
  # Rotation matrix T 
  svd.results <- svd(B.normalized %*% t(A.normalized))
  U <- svd.results$u
  V <- svd.results$v
  T <- V %*% t(U)
  
  # B transformed
  B.transformed <- T %*% B.normalized
  
  # Error after superimposition
  RSS <- norm(A.normalized - B.normalized,  type = "F")
  
  # Return
  return(list(rotation.mtx = T, B.transformed = B.transformed, RSS = RSS))
}


distance_squared_inner_prod = function(mu1,mu2,Sigma1,Sigma2){
  value = t(mu1)%*% mu2 %*% t(mu2) %*%mu1 +sum(diag(Sigma1 %*% Sigma2)) +  t(mu1)%*% Sigma2 %*%mu1 + t(mu2)%*% Sigma1 %*%mu2
  
  return(value)
}

positions_to_edges_Gaussian = function(Xt,beta,sigma){
  n = length(Xt[,1])
  edge_mat = matrix(rep(0,n*n),nrow=n)
  t=1
  for (i in 1:n){
    for(j in 1:n){
      if(j > i)
      {edge_mat[j,i] = rnorm(n=1,mean=beta+t(Xt[i,])%*%Xt[j,],sd=sigma)
      t= t+1}
    }
  }
  edge_mat = edge_mat+ t(edge_mat)
  return(edge_mat)
}

positions_to_edges = function(Xt,beta){
  n = length(Xt[,1])  
  edge_mat = matrix(rep(0,n*n),nrow=n)
  t=1
  for (i in 1:n){
    for(j in 1:n){
      if(j > i)
      {
        binary_mean = 1/(1+exp(-beta-t(Xt[i,])%*%Xt[j,]))
        edge_mat[j,i] = rbinom(n=1,size=1, prob = binary_mean)
        t= t+1}
    }
  }
  edge_mat = edge_mat+ t(edge_mat)
  return(edge_mat)
}


ggplot_dis= function(Xm,Ym){
  n0= length(Xm)
  data_fr = as.data.frame(cbind(c(1:n0),Xm,Ym))
  names(data_fr) = c("index","neg_inner_prod","connected")
  data_fr$connected = as.factor(data_fr$connected)
  levels(data_fr$connected) = c("unconnected","connected")
  gp = ggplot(data_fr,aes(x=index,y=neg_inner_prod))+geom_point(size=4,aes(shape=connected,color=connected))+ 
    theme(legend.position = "none")+theme(axis.text=element_text(size=16),
                                          axis.title=element_text(size=18,face="bold"))+theme(legend.title = element_blank())+theme(legend.text = element_text( size=16))
  return(gp)
}

radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

plot_clus = function(Xm1,Y1,clus,t){
  n = length(Xm1[,1])
  diag(Y1) = 0
  Y1[lower.tri(Y1, diag = FALSE)] = 0
  NodeList <- data.frame((1:n), x=Xm1[,1] ,y=Xm1[,2])
#    Edge_list = data.frame(matrix(ncol = 2, nrow = 0))
 Edge_list = igraph::as_data_frame	(igraph::graph_from_adjacency_matrix(Y1,mode="undirected"),'edge')
  pal <- RColorBrewer::brewer.pal(4,"Accent")
  vertex.col = pal[factor(clus)]
  v_shape= setdiff(igraph::shapes(), "")[c(1,1)]
  v_size = c(5,5)
  a<- igraph::graph_from_data_frame(vertices = NodeList, d= Edge_list, directed = FALSE)
  #lab.locs <- radian.rescale(x=1:n, direction=-1, start=0)
#  igraph::plot.igraph(a,vertex.color=vertex.col,vertex.shape=v_shape[clus], vertex.label=1:n,vertex.size=v_size[clus],
 #             edge.width=0.1,vertex.label.dist=2,vertex.label.color=pal[3])
  Lay<- ggraph::create_layout(a, layout = "nicely")
  Lay$x<- NodeList$x
  Lay$y<- NodeList$y
 g_plot = ggraph::ggraph(Lay)+
    geom_edge_link(colour = "black", alpha = 0.8, show.legend = F) + 
    geom_node_point()

#  title(main=paste("t =",t,sep = ' '),cex.main=1.5)
#  legend('topleft',legend=levels(factor(clus)),pch=c(1,0),bty = "n",cex=1.6)
#  legend('bottom',inset = c(0, -.1),legend=paste("t =",t,sep = ' '),cex = 1.5,bty = "n")
  return( g_plot)
}

plot_clus_igraph = function(Xm1,Y1,clus,t,label,NA_ind){
  n = length(Xm1[,1])
  diag(Y1) = 0
  Y1[lower.tri(Y1, diag = FALSE)] = 0
  NodeList <- data.frame((1:n), x=Xm1[,1] ,y=Xm1[,2])
  #    Edge_list = data.frame(matrix(ncol = 2, nrow = 0))
  Edge_list = igraph::as_data_frame	(igraph::graph_from_adjacency_matrix(Y1,mode="undirected"),'edge')
 # pal <- RColorBrewer::brewer.pal(12,"Paired")
  pal <- RColorBrewer::brewer.pal(8,"Accent")
  vertex.col = pal[factor(clus)]
  vertex.col[is.na(label)]=adjustcolor("white", alpha.f = 0)
  v_shape= setdiff(igraph::shapes(), "")[c(1,1)]
  v_size = c(20)
  a<- igraph::graph_from_data_frame(vertices = NodeList, d= Edge_list, directed = FALSE)
  igraph::plot.igraph(a,vertex.color=vertex.col,vertex.shape=v_shape[1], vertex.label=label,vertex.size=v_size,
                      edge.arrow.size=0.001,vertex.label.cex =1,vertex.label.color = "black"  ,vertex.frame.color = adjustcolor("black", alpha.f = 0),
            vertex.color = adjustcolor("white", alpha.f = 0),edge.color=adjustcolor("pink", alpha.f = 1),display.isolates=FALSE,xlim=c(-1,1),ylim=c(-1,1), asp = 0)
  #  title(main=paste("t =",t,sep = ' '),cex.main=1.5)
  #  legend('topleft',legend=levels(factor(clus)),pch=c(1,0),bty = "n",cex=1.6)
 legend('bottom',inset = c(0, -.1),legend=paste("t =",t,sep = ' '),cex = 1.5,bty = "n")
  return( NULL)
}


plot_clus_igraph_sub = function(Xm1,Y1,clus,t,label,axis_lim){
  NA_ind = which(is.na(label))
  Xm1 = Xm1[-NA_ind,]
  Y1= Y1[-NA_ind,-NA_ind]
  
  pal <- RColorBrewer::brewer.pal(8,"Accent")
  vertex.col_full = pal[factor(clus)]

  clus=clus[-NA_ind]
  label = label[-NA_ind]
  
  vertex.col=vertex.col_full[-NA_ind]
  
  n = length(Xm1[,1])
  diag(Y1) = 0
  Y1[lower.tri(Y1, diag = FALSE)] = 0
  NodeList <- data.frame((1:n), x=Xm1[,1] ,y=Xm1[,2])
  #    Edge_list = data.frame(matrix(ncol = 2, nrow = 0))
  Edge_list = igraph::as_data_frame	(igraph::graph_from_adjacency_matrix(Y1,mode="undirected"),'edge')
  

  v_shape= setdiff(igraph::shapes(), "")[c(1,1)]
  v_size = c(40)
  a<- igraph::graph_from_data_frame(vertices = NodeList, d= Edge_list, directed = FALSE)
  igraph::plot.igraph(a,vertex.shape=v_shape[1],vertex.color=vertex.col, vertex.label=label,vertex.size=v_size,
                      edge.arrow.size=0.001,vertex.label.cex = 1.5,vertex.label.color = adjustcolor("black", alpha.f = 1) ,vertex.frame.color = adjustcolor("white", alpha.f = 0),
                      edge.color=adjustcolor("pink", alpha.f = 1),rescale=F,axes = T,
                      xlim=c(axis_lim[1],axis_lim[2]),ylim=c(axis_lim[3],axis_lim[4]), asp = 0)
  #  title(main=paste("t =",t,sep = ' '),cex.main=1.5)
  #  legend('topleft',legend=levels(factor(clus)),pch=c(1,0),bty = "n",cex=1.6)
  legend('top',inset = c(0, -0.1),legend=paste("year =",t,sep = ' '),cex = 2.5,bty = "n")
  return( NULL)
}

plot_sub_preprocess = function(Xm,Y,t1,t2,data_label){
  non_NA_ind = NULL
  x_min = 0
  x_max = 0
  y_min = 0
  y_max = 0
  for (t in t1:t2){
    non_NA_ind = union(non_NA_ind, which(apply(Y[[t]],1,sum)!=0))
    x_min =min(x_min,min(Xm[[t]][,1]-0.5))
    x_max =max(x_max,max(Xm[[t]][,1]+0.5))
    y_min =min(y_min,min(Xm[[t]][,2]-0.5))
    y_max =max(y_max,max(Xm[[t]][,2]+0.5))
  }
  axis_lim = c(x_min,x_max,y_min,y_max)
  label.plot_non_NA = rep(NA,n)
  label.plot_non_NA[non_NA_ind] = as.character(data_label[non_NA_ind])
  return(list(axis_lim =axis_lim, label.plot_non_NA = label.plot_non_NA))
}


label_check = function(Xm,data_label,t,node,k,NA_list){
  distance_inner = Xm[[t]]%*%Xm[[t]][node,]
  distance_inner[NA_list] = NA
  top_k = as.character(data_label[sort(distance_inner,index.return=TRUE,decreasing = TRUE,na.last = TRUE)$ix[1:k]])
  bottom_k =   as.character(data_label[sort(-distance_inner,index.return=TRUE,decreasing = TRUE,na.last = TRUE)$ix[1:k]])
  return(list(country= as.character(data_label[node]),top_close=top_k,top_far=bottom_k))
}
