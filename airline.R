rm(list=ls())
set.seed(1234)
file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste(file_location,'/function',sep=""))
source('MP_gibbs_multi_Sigma_adaptive.R')
source('MP_gibbs_tensor_adaptive_stochastic.R')
library(MASS)
library(mvtnorm)
library(igraph)
library(DiagrammeR)
library(pROC)
library(rTensor)
library(foreach)
library(doParallel)
library(animation)
library(stringr)
library(gridExtra)


#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]) 
registerDoParallel(cl)


setwd(paste(file_location,'/data',sep=""))


airline_2018 = read.csv('2018.csv',na.strings = "") 
airline_2019 = read.csv('2019.csv',na.strings = "") 
airline_2020 = read.csv('2020.csv',na.strings = "") 
airline_2021 = read.csv('2021.csv',na.strings = "")  



airline = as.data.frame(rbind(airline_2018,airline_2019,airline_2020,airline_2021))

airline = na.omit(airline)



carrier = names(sort(table(airline$CARRIER_NAME),decreasing = TRUE)[1:12])



city = names(sort(table(c(airline$ORIGIN_CITY_NAME,airline$DEST_CITY_NAME)),decreasing = TRUE)[1:100])

airline = airline[airline$CARRIER_NAME %in% carrier,]

airline = airline[airline$ORIGIN_CITY_NAME %in% city,]

airline = airline[airline$DEST_CITY_NAME %in% city,]


n1 = 100

n2 = n1

n3 = length(unique(airline$CARRIER_NAME))

T = 4*12

Y = array(rep(0,),dim=c(n1,n2,n3,T))

for (i in 1:length(airline$CARRIER_NAME)){
  row_value = airline[i,]
  row_time = (row_value[4]-2018)*12+row_value[5]
  Y[which(city==as.character(row_value[2])),which(city==as.character(row_value[3])), 
    which(carrier==as.character(row_value[1])),as.numeric(row_time)] =  1
    #Y[which(city==as.character(row_value[2])),which(city==as.character(row_value[3])), 
     #  which(carrier==as.character(row_value[1])),as.numeric(row_time)]+1
 # Y[which(city==as.character(row_value[3])),which(city==as.character(row_value[2])), 
#    which(carrier==as.character(row_value[1])),as.numeric(row_time)]  = 1
    #Y[which(city==as.character(row_value[3])),which(city==as.character(row_value[2])), 
    #   which(carrier==as.character(row_value[1])),as.numeric(row_time)] +1
}



d = 2

MF_list = MP_Gaussain_tensor_adaptive (Y, gap =1e-2, gap_per_iter=1e-4,max_iter=50,d=d,global_prior='Cauthy',loss = 'Gaussian')
stopCluster(cl)

data_cp = array(rep(0,),dim=c(n1,n2,n3,T))
factor_cp = array(rep(0,),dim=c(n3,d,T))
for (t in 1:T){
  CPD <- cp(as.tensor(Y[,,,t]),num_components =d,max_iter = 500, tol = 1e-3)
  for(k in 1:n3){
    data_cp[,,k,t] = CPD$est@data[,,k]
  }
  factor_cp[,,t] = CPD$U[[3]]
}

CPD_all <- cp(as.tensor(Y),num_components =d,max_iter = 500, tol = 1e-3)

data_cp_all = CPD_all$est@data
#MF_diff = rep(0,T)


# for (t in 1:T){
#   MF_diff[t] = norm(as.matrix(Y[,,,t], nrow=1),'F')/sqrt(n1*n2*n3)
#   }


myDate = seq(as.Date('2018-1-1'),to=as.Date('2021-12-1'),by='1 month')
#myDate = format(myDate,"%m-%Y")

#df_plot = data.frame(list(time=as.Date(myDate),diff=c(MF_diff),method =c(rep('FFS',T)) ))
#,rep('raw data',T-1),rep

#,data_diff,data_f_diff

carrier = word(carrier,1)


library(ggplot2)
# ggplot(data=df_plot, aes(x=time, y=diff)) +
#   geom_point(size=1.5,color=1) + ylab("Mean Connection Intensities")+scale_linetype_manual(values=c("twodash", "solid"))+
#   theme(plot.title = element_text(hjust = 0.5))+
#   scale_x_date(date_labels = "%m-%Y")+
#   theme(axis.text=element_text(size=14),
#         axis.title=element_text(size=16),
#         legend.title = element_text( size = 14),
#         legend.text = element_text( size = 14))







pred_X = MF_list$Mean_X




#base_month = as.Date("2018-01")

#index = c(28,29,30,31,32,48)

index = c(16,17,18,19,20,36)+12

locations = rbind(pred_X[1:n3,,3,index[1]],pred_X[1:n3,,3,index[2]],
                  pred_X[1:n3,,3,index[3]],pred_X[1:n3,,3,index[4]],pred_X[1:n3,,3,index[5]],
                  pred_X[1:n3,,3,index[6]])

#locations_norm = pred_X

# locations = rbind(factor_cp[1:n3,,index[1]],factor_cp[1:n3,,index[2]],
#                   factor_cp[1:n3,,index[3]],factor_cp[1:n3,,index[4]],factor_cp[1:n3,,index[5]],
#                   factor_cp[1:n3,,index[6]])

myDate2=as.character(format(myDate,"%m-%Y"))

carrier2 = carrier

carrier2[carrier2!='Southwest']=NA

df_plot2 = data.frame(list(time=c(rep(myDate2[index[1]],n3) ,rep(myDate2[index[2]],n3),
                                  rep(myDate2[index[3]],n3),rep(myDate2[index[4]],n3),rep(myDate2[index[5]],n3),
                                  rep(myDate2[index[6]],n3)),
           x1=locations[,1],x2=locations[,2],carrier=rep(carrier[1:n3],6),carrier2=rep(carrier2[1:n3],6)))



p1 = ggplot(data=df_plot2,aes(x=x1, y=x2, group=time))+
  geom_point(size=1.5,aes(color = time, shape = time))+
  ggtitle("Latent Factors for Effects of Carriers") + theme(plot.title = element_text(hjust = 0.5))+
  scale_shape_manual(values=seq(0,7))+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.position = "none")


p2 = ggplot(data=df_plot2,aes(x=x1, y=x2, group=time))+
  geom_point(size=1.5,aes(color = time, shape = time))+geom_text(aes(label = carrier),size = 2)+
ggtitle("Latent Factors for Effects of Carriers") + theme(plot.title = element_text(hjust = 0.5))+
  scale_shape_manual(values=seq(0,7))+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.title = element_text( size = 14),
        legend.text = element_text( size = 14))

#grid.arrange(p1, p2, nrow = 1)

library(cowplot)

plot_grid(p1, p2, nrow = 1, rel_widths = c(10.5/24, 13.5/24))


# ggplot(data=df_plot2,aes(x=x1, y=x2))+
#   geom_point(size=1.5)+geom_text(aes(label = carrier),size = 2)+
#   facet_wrap(~time,ncol=3, scales = "free")+
#   ggtitle("Latent Space for Effects of Carriers at different time")+ theme(plot.title = element_text(hjust = 0.5))
# 
library(vegan)
library(gridGraphics)
library(grid)
# 
# 
# # load package
# library(pheatmap)
# grab_grob = function(){
#   grid.echo()
#   grid.grab()
# }
# gl = lapply(1:4, function(i){
#   op <- par(cex.main = 0.7, font.main = 1)
#  heatmap(as.matrix(vegdist(t(pred_X[i,,3,]),method='euclidean')),Rowv=NA,
#           Colv=NA,symm = T,main=carrier[i], 
#           key=FALSE, keysize=1.0, symkey=FALSE, density.info='none'
#   )
#   par(op)
#   grab_grob()
# })
# dev.new()
# grid.arrange(grobs=gl,ncol=1,clip=TRUE)
# 
# 

carrier0 = names(sort(table(airline$CARRIER_NAME),decreasing = TRUE)[1:n3])


carrier0[carrier0=='Federal Express Corporation']='Federal Express Co.'

library(magick)

n_length = 12
data_dist = NULL
data_dist_prob = NULL
time_index=1:48
T_len = length(time_index)
for (i in 1:n_length){
  data_dist_each = NULL
  data_dist_mean = NULL
  #lf_matrix = t(matrix(Y[1:n_length,1:n_length,i,],n_length*n_length,T))
  #lf_matrix =  t(pred_X[i,,3,])
  #lf_matrix = t(matrix(data_cp_all[1:n_length,1:n_length,i,],n_length*n_length,T))
  lf_matrix = t(matrix(MF_list$pred_mean[1:n_length,1:n_length,i,time_index],n_length*n_length,T_len))
  lf_matrix_mean = rowMeans(lf_matrix)
  for (t1 in 1:T_len){
    for (t2 in 1:T_len){
      data_new= sqrt(sum((lf_matrix[t1,]-lf_matrix[t2,])^2))/sqrt(sum((lf_matrix[t1,])^2))/sqrt(sum((lf_matrix[t2,])^2))
      if (sum((lf_matrix[t1,])^2)!= 0 & sum((lf_matrix[t2,])^2)!= 0){
      data_dist_each = c(data_dist_each,data_new)
      } else{
        data_dist_each = c(data_dist_each,0)
      }
      
     data_dist_mean= c(data_dist_mean,sqrt(sum((lf_matrix_mean[t1]-lf_matrix_mean[t2])^2))/
                         abs(lf_matrix_mean[t1])/abs(lf_matrix_mean[t2]))
    }
  }
#  data_dist_each = data_dist_each/ max(data_dist_each)
  data_dist =c(data_dist, data_dist_each)
  data_dist_prob = c(data_dist_prob,data_dist_mean)
}
data_dist = data_dist/max(data_dist)
data_dist_prob = data_dist_prob/max(data_dist_prob)
data_date_1 = rep(as.Date(myDate[time_index]),each=T_len)
data_date_2 = as.Date(rep(myDate[time_index],T_len))
       
data_plot_heat = as.data.frame(list(date_1=data_date_1,date_2= data_date_2,dist=data_dist, prob_diff = data_dist_prob,
                                    carrier=rep(carrier0[1:n_length],each=T_len*T_len)))

ggplot(data_plot_heat, aes(date_1, date_2, fill= dist)) + 
  facet_wrap(~carrier,scales='free')+scale_x_date(date_labels = "%Y")+
  scale_fill_gradient(low="blue", high="red") +
  geom_tile()



ggplot(data_plot_heat, aes(date_1, date_2, fill= prob_diff)) + 
  facet_wrap(~carrier,scales='free')+scale_x_date(date_labels = "%Y")+
  scale_fill_gradient(low="blue", high="red") +
  geom_tile()


#matrix(CPD_all$est@data[,,i,],n1*n2,T)
#pred_X[i,,3,],2,function(x)(x/sqrt(sum(x^2)))))
#/sqrt(apply(Y[,,i,],3,sum)
#apply(data_cp[i,,3,],2,function(x)(x/sqrt(sum(x^2))))


#par(mfrow=c(4,4))
#MF_list$pred_mean[,,1,]

# for (i in 1:n3) {
#   ifile <- paste0(i,'_heatmap.pdf')
#   pdf(ifile)
#  # cc <- rainbow(T,start = 0,end = 0.18)
#  # lf_matrix = t(matrix(MF_list$pred_mean[,,i,],n1*n2,T))
#   lf_matrix = t(pred_X[i,,3,])
#   par(cex.main=1.4)
#   heatmap(as.matrix(vegdist(lf_matrix,
#                              method='euclidean')),Rowv=NA,
#           Colv=NA,symm = T,main=carrier0[i],
#           #RowSideColors=cc,
#           key=FALSE, keysize=1.0, symkey=FALSE, density.info='none'
#   )
#   dev.off()
# }
# system('montage -geometry 100% -tile 4x4 ./*_heatmap.pdf outfile.pdf')


# df_plot_all = data.frame(list(time=rep(myDate2,each=n3),
#                            x1=locations[,1],x2=locations[,2],carrier=rep(carrier,6),carrier2=rep(carrier2,6)))


