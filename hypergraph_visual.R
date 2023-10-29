library(ggplot2)
library(ggalt) 
library(Thresher)



Xm = Xm2[[10]]
  km = kmeans(Xm,centers=4,nstart=10)
  label_base = factor(km$cluster)
  df1 = NULL
  df1$x = Xm[,1]
  df1$y = Xm[,2]
  df1$hyper_cluster = factor(km$cluster)
  df1$index = 1:20
  df1$status = factor(c('active','active',rep('non-active',n-2)))
  df1= as.data.frame(df1)

  


  Xm = Xm2[[25]]
  

  km = kmeans(Xm,centers=4,nstart=10)
  df2 = NULL
  df2$x = Xm[,1]
  df2$y = Xm[,2]
  label_new = factor(km$cluster)
  R <- factor(remap(label_base, label_new))
  df2$hyper_cluster = R
  df2$index = 1:20
  df2$status = factor(c('active','active',rep('non-active',n-2)))
  df2= as.data.frame(df2)

  
  Xm = Xm2[[50]]
  
  km = kmeans(Xm,centers=4,nstart=10)

  df3 = NULL
  df3$x = Xm[,1]
  df3$y = Xm[,2]
  label_new = factor(km$cluster)
  R <- factor(remap(label_base, label_new))
  df3$hyper_cluster = R
  df3$index = 1:20
  df3$status = factor(c('active','active',rep('non-active',n-2)))
  df3= as.data.frame(df3)

  
  
  
  df = rbind(df1,df2,df3)
  
  df$subfigure = factor(c(rep('t1',20),rep('t2',20),rep('t3',20)))
  
  
  
  ggplot(df, aes(x, y, colour = status)) + geom_point() +
    geom_encircle(aes(group=hyper_cluster,fill=hyper_cluster),alpha=0.3)+ 
    facet_wrap(~subfigure)+theme(text = element_text(size=15)) 
