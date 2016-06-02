X_Values <- c(6, 19, 15, 11, 18, 9, 19, 18, 5, 4, 7, 21, 1, 1, 0, 5)
Y_Values <- c(12, 7, 4, 0, 12, 20, 22, 17, 11, 18, 15, 18, 19, 4, 9, 11)
dataPoints = data.frame(X_Values, Y_Values)

distance1 <- function(dp,centers,m)
{  

  f_distance1=c(0,10000,10000,10000)
  for(i in 1:m)
  {
    f_distance1[i]= sqrt((dp[1,1]-centers[i,1])^2+(dp[1,2]-centers[i,2])^2)
  }
 # distance1[i]=a[1]
 
  #min distance calc
  min_distance = min(f_distance1,na.rm = TRUE);
  #center calac
  cent=which(f_distance1==min_distance)
  C_Center=centers[cent,]
  x= c(min_distance,as.matrix(C_Center))
   
  return(x) 
}




BSAS <- function(dataPoints,theta,max_cluster,centers)
{
  clustersAndCenter=matrix(nrow=16, ncol=5)
  clustersAndCenter[1,]= c(dataPoints[1,1],dataPoints[1,2],dataPoints[1,1],dataPoints[1,2],1) 
  #centers= dataPoints[1,]
  
  m=1
  
  for(i in  2 : nrow(dataPoints)) 
  {  
    dp= dataPoints[i,]
    
    # x=as.data.frame(x)
    x=distance1(dp,centers,m)
    if(x[1]>theta && m < max_cluster)
    {
      m= m+1
      centers[m,] = dp 
      
     # cluster_cen[m]= 1
      clustersAndCenter[i,1]=dp[1,1]
      clustersAndCenter[i,2]=dp[1,2]
      clustersAndCenter[i,3]=centers[m,1]
      clustersAndCenter[i,4]=centers[m,2]
      clustersAndCenter[i,5]=m
    } 
    else
    { 
      
      clustersAndCenter[i,1]=dp[1,1]
      clustersAndCenter[i,2]=dp[1,2]
      clustersAndCenter[i,3]=x[2]
      clustersAndCenter[i,4]=x[3]
     
      for(j in 1:m)
      { 
        if(centers[j,1]==x[2] && centers[j,2]== x[3])
        { 
          clustersAndCenter[i,5]=j
        #update cluster center
          a1=which(clustersAndCenter[,5]==j)
        temp_x=0
        temp_y=0
        for(temp_ce in a1)
        {
        temp_x =  temp_x + clustersAndCenter[temp_ce,1]
        temp_y =temp_y + clustersAndCenter[temp_ce,2]
        } 
        centers[j,1]=temp_x/length(a1)
        centers[j,2]=temp_y/length(a1)
        
       }
      }
    }
    
  }
  return(clustersAndCenter)
}



max_cluster=4
theta = 12 
centers=dataPoints[1,]
clustersAndCenter = BSAS(dataPoints,theta,max_cluster,centers)
library(cluster)
scaled =scale(dataPoints)
clusplot(scaled, clustersAndCenter[,5], color=TRUE, shade=TRUE, labels=2,lines=0)

#1.b
X_Values=rev(X_Values)
Y_Values=rev(Y_Values)

dataPoints_rev = data.frame(X_Values, Y_Values)
centers=dataPoints_rev[1,]
clustersAndCenter_rev = BSAS(dataPoints_rev,theta,max_cluster,centers)

library(cluster)
scaled =scale(dataPoints_rev)

clusplot(scaled, clustersAndCenter_rev[,5], color=TRUE, shade=TRUE, labels=2,lines=0)



#1.c

library(fossil)
rand.index(clustersAndCenter_rev[,5],clustersAndCenter[,5])
#BSAS algorithm depends on the order of the data it process. First point is treated as first cluster center
# then next point is distance <theta then it is added to clusters otherwise it is treated as  custer center.
# this is the reason we have different cluster centers and different clusters if same data is processed in different order






#2


scaled =scale(dataPoints)

distance_sing <- dist(as.matrix(dataPoints))
hcluster_sing <- hclust(distance_sing, method='single')

plot(hcluster_sing)
rect.hclust(hcluster_sing, k = 3, border = 3)
cut_tree=cutree(hcluster_sing, k = 3)
library(cluster)
clusplot(scaled, cut_tree, color=TRUE, shade=TRUE, labels=2,lines=0)


#b
distance <- dist(as.matrix(dataPoints))

hcluster_compl <- hclust(distance)
plot(hcluster_compl)
rect.hclust(hcluster_compl, k = 3, border = 3)
cut_tree=cutree(hcluster_compl, k = 3)
library(cluster)
clusplot(scaled, cut_tree, color=TRUE, shade=TRUE, labels=2, lines=0)



#C
library(GMD)
SSE=css(distance,cut_tree)
WSE=css(distance,cut_tree)$wss

#D
proximity_m <- dist(as.matrix(dataPoints), upper=TRUE, diag=TRUE)
inc_m <- as.matrix(dist(as.matrix(cut_tree_sing), upper=TRUE, diag=TRUE, method="manhattan"))
inc_m[inc_m==1]<-2
inc_m[inc_m==0]<-1
inc_m[inc_m==2]<-0
cor(c(inc_m), c(as.matrix(proximity_m)))

comp_im <- as.matrix(dist(as.matrix(cut_tree_comp), upper=TRUE, diag=TRUE, method="manhattan"))
comp_im[comp_im==1]<-2
comp_im[comp_im==0]<-1
comp_im[comp_im==2]<-0

cor(c(comp_im), c(as.matrix(proximity_m)))


#3
xValues_DB= c(1, 3, 5, 6, 8, 11, 12, 13, 14, 15, 16, 22, 28, 32, 33, 34, 35, 36, 37, 42, 58)
library('fpc')
DB_scan_ep4=dbscan(xValues_DB, eps=4, MinPts=3, seeds = TRUE, showplot = 1, countmode=NULL)
#clusters obtained
DB_scan_ep4$cluster
length(which(DB_scan_ep4$cluster==0)) 
length(which(DB_scan_ep4$cluster==1))
length(which(DB_scan_ep4$cluster==2))


DB_scan_ep6=dbscan(xValues_DB, eps=6, MinPts=3, seeds = TRUE, showplot = 1, countmode=NULL)
#clusters obtained
DB_scan_ep6$cluster
length(which(DB_scan_ep6$cluster==0)) 
length(which(DB_scan_ep6$cluster==1))
length(which(DB_scan_ep6$cluster==2))

library(fossil)
rand.index(DB_scan_ep6$cluster,DB_scan_ep4$cluster)
