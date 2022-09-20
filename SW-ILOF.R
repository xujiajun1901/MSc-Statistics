
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(matrixStats)
set.seed(12345)



#Simulate a univariate data using normal distribution
N=200
mean=0
sd=1
#dat2=my.rnorm(200,10,0,1)
set.seed(123)
data1<- rnorm(N, mean = 0, sd = 1)
ind <- sample(1:N, 10, replace=FALSE )
data1[ind] <- (abs(data1[ind]) + 3*sd) * sign(data1[ind])

#local outlier factor
#install.packages("DescTools")
library(DescTools)
library(Rlof)
library(FNN)
library(dbscan)


data11<-as.data.frame(data1)
data1_lof<-LOF(data1,10)
data1_Rlof<-lof(data11,10)
which(data1_lof > 2)
length(which(data1_lof > 2))
outliers2<-order(data1_lof, decreasing=T)[1:length(which(data1_lof > 2))]
sort(outliers2)
sort(ind)
data1_lof[outliers2]
length(which(outliers2%in% ind))


#lof decomposition
data=data1
data=data1[151:200]
data=as.matrix(data)


knneigh.vect <- function(x, data, k) {
  temp = as.matrix(data)
  numrow = dim(data)[1]
  dimnames(temp) = NULL
  difference <- scale(temp, x, FALSE)
  dtemp <- drop(difference^2 %*% rep(1, ncol(data)))
  dtemp = sqrt(dtemp)
  order.dist <- order(dtemp)
  nndist = dtemp[order.dist]
  knndist = nndist[k + 1]
  neighborhood = drop(nndist[nndist <= knndist])
  neighborhood = neighborhood[-1]
  numneigh = length(neighborhood)
  index.neigh = order.dist[1:numneigh + 1]
  num1 = length(index.neigh) + 3
  num2 = length(index.neigh) + numneigh + 2
  neigh.dist = c(num1, num2, index.neigh, neighborhood)
  return(neigh.dist)
}

#dataset=data
#neighbors=10

dist.to.knn <- function(dataset, neighbors) {
  numrow = dim(dataset)[1]
  knndist = rep(0, 0)
  for (i in 1:numrow) {
    neighdist = knneigh.vect(dataset[i, ], dataset, neighbors)
    if (i == 2) {
      if (length(knndist) < length(neighdist)) {
        z = length(neighdist) - length(knndist)
        zeros = rep(0, z)
        knndist = c(knndist, zeros)
      }
      else if (length(knndist) > length(neighdist)) {
        z = length(knndist) - length(neighdist)
        zeros = rep(0, z)
        neighdist = c(neighdist, zeros)
      }
    }
    else {
      if (i != 1) {
        if (dim(knndist)[1] < length(neighdist)) {
          z = (length(neighdist) - dim(knndist)[1])
          zeros = rep(0, z * dim(knndist)[2])
          zeros = matrix(zeros, z, dim(knndist)[2])
          knndist = rbind(knndist, zeros)
        }
        else if (dim(knndist)[1] > length(neighdist)) {
          z = (dim(knndist)[1] - length(neighdist))
          zeros = rep(0, z)
          neighdist = c(neighdist, zeros)
        }
      }
    }
    knndist = cbind(knndist, neighdist)
  }
  return(knndist)
}


reachability_dist <- function(distdata, k) {
  p = dim(distdata)[2]
  reach_dist<-distdata
  #lrd = rep(0, p)
  #reach_dist = rep(0, p)
  for (i in 1:p) {
    j = seq(3, 3 + (distdata[2, i] - distdata[1, i]))
    numneigh = distdata[2, i] - distdata[1, i] + 1
    temp = rbind(diag(distdata[distdata[2, distdata[j, 
                                                    i]], distdata[j, i]]), distdata[j + numneigh, 
                                                                                    i])
    #reach = 1/(sum(apply(temp, 2, max))/numneigh)
    reach_dist[j + numneigh,i] = apply(temp, 2, max)
  }
  reach_dist
}



reachability <- function(distdata, k) {
  p = dim(distdata)[2]
  #reach_dist<-distdata
  lrd = rep(0, p)
  #reach_dist = rep(0, p)
  for (i in 1:p) {
    j = seq(3, 3 + (distdata[2, i] - distdata[1, i]))
    numneigh = distdata[2, i] - distdata[1, i] + 1
    temp = rbind(diag(distdata[distdata[2, distdata[j, 
                                                    i]], distdata[j, i]]), distdata[j + numneigh, 
                                                                                    i])
    reach = 1/(sum(apply(temp, 2, max))/numneigh)
    lrd[i] = reach
  }
  lrd
}






k=10
data = as.matrix(data1[151:200])
distdata = dist.to.knn(data, k)
p = dim(distdata)[2]
reachdist=reachability_dist(distdata, k)
lrddata = reachability(distdata, k)
lof = rep(0, p)
for (i in 1:p) {
  nneigh = distdata[2, i] - distdata[1, i] + 1
  j = seq(0, (nneigh - 1))
  local.factor = sum(lrddata[distdata[3 + j, i]]/lrddata[i])/nneigh
  lof[i] = local.factor
}
lof


rNN <- sapply(1:n, function(i) {
  as.vector(which(dist.obj$id[, 1:max_nb] == i, arr.ind = TRUE)[, 
                                                                1])
})


data2<-data
distdata1<-distdata
#Deletion
#reverse neighbours of deletion point
dist.obj1 <- dbscan::kNN(data2, k)
#n <- nrow(data2)

rNN_x_1 <- sapply(1, function(i) {
  as.vector(which(dist.obj1$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                           1])
})

data2<-as.matrix(data2[-c(1),])

#Update k-distance

j1=seq(3, 3 + (distdata1[2, 1] - distdata1[1, 1]))
j2=seq(distdata1[1, 1],distdata1[2, 1])
distdata1[j1,]=distdata1[j1,]-1
dist.obj2 <- dbscan::kNN(data2, k)

for (i in 1:length(rNN_x_1)) {
  if (length(rNN_x_1)==1){
    if (length(rNN_x_1[[1]])!=0){
      x1<-rNN_x_1[i]
      numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
      if (numneigh == k) {
        j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
        j2=seq(distdata1[1, x1],distdata1[2, x1])
        distdata1[j1,x1]=dist.obj2$id[x1-1,]
        distdata1[j2,x1]=dist.obj2$dist[x1-1,]
      } else if (numneigh!=k) {
        print("warning numneigh!=k")
      }
    }
  } else {
    x1<-rNN_x_1[i,1]
    numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
    if (numneigh == k) {
      j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
      j2=seq(distdata1[1, x1],distdata1[2, x1])
      distdata1[j1,x1]=dist.obj2$id[x1-1,]
      distdata1[j2,x1]=dist.obj2$dist[x1-1,]
    } else if (numneigh!=k) {
      print("warning numneigh!=k")
    }
  }
}





distdata1<-distdata1[,-c(1)]
s_update_lrd=rNN_x_1
j32=seq(3, 3 + (distdata1[2, 1] - distdata1[1, 1]))
j4=seq(distdata1[1, 1],distdata1[2, 1])
reachdist[j32,]=reachdist[j32,]-1
#update reach_dist
for (i in 1:length(rNN_x_1)) {
  if (length(rNN_x_1)==1) {
    if (length(rNN_x_1[[1]])!=0) {
      x2<-rNN_x_1[i]-1
      knn_x2<-as.matrix(dist.obj2$id[x2,-c(length(dist.obj2$id[x2,]))])
      numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
      if (numneigh1 == k) {
        j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
        j4=seq(distdata1[1, x2],distdata1[2, x2])
        inds<-reachdist[j3,x2+1]
        inds1<-reachdist[j4,x2+1]
        inds2<-match(dist.obj2$id[x2,],reachdist[j3,x2+1])
        reachdist[j3,x2+1]=dist.obj2$id[x2,]
        reachdist[j4,x2+1]=inds1[inds2]
        x3<-length(which(is.na(reachdist[,x2+1])))
        x4<-which(is.na(reachdist[,x2+1]))
        if (x3!=0){
          temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                              x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                   x2])
          #reach = 1/(sum(apply(temp, 2, max))/numneigh)
          
          reach_dist_x2 = apply(temp1, 2, max)
          reachdist[x4,x2+1]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
        }
      } else if (numneigh1!=k) {
        print("warning numneigh!=k")
      }
      for (i in 1:dim(knn_x2)[1]) {
        x5<-knn_x2[i,1]
        x6<-which(reachdist[,x5+1]== x2)+numneigh1
        reachdist[x6,x5+1]=distdata1[length(distdata1[,x2]),x2]
        #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
        #x8<-which(reachdist[,x2]== dim(data2)[1])
        if (x2 %in% dist.obj2$id[x5,]) {
          s_update_lrd=as.matrix(union(s_update_lrd,(x5+1)))
        }
        
      }
    }
    
  } else {
    x2<-rNN_x_1[i,1]-1
    knn_x2<-as.matrix(dist.obj2$id[x2,-c(length(dist.obj2$id[x2,]))])
    numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
    if (numneigh1 == k) {
      j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
      j4=seq(distdata1[1, x2],distdata1[2, x2])
      inds<-reachdist[j3,x2+1]
      inds1<-reachdist[j4,x2+1]
      inds2<-match(dist.obj2$id[x2,],reachdist[j3,x2+1])
      reachdist[j3,x2+1]=dist.obj2$id[x2,]
      reachdist[j4,x2+1]=inds1[inds2]
      x3<-length(which(is.na(reachdist[,x2+1])))
      x4<-which(is.na(reachdist[,x2+1]))
      if (x3!=0){
        temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                            x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                 x2])
        #reach = 1/(sum(apply(temp, 2, max))/numneigh)
        
        reach_dist_x2 = apply(temp1, 2, max)
        reachdist[x4,x2+1]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
      }
    } else if (numneigh1!=k) {
      print("warning numneigh!=k")
    }
    for (i in 1:dim(knn_x2)[1]) {
      x5<-knn_x2[i,1]
      x6<-which(reachdist[,x5+1]== x2)+numneigh1
      reachdist[x6,x5+1]=distdata1[length(distdata1[,x2]),x2]
      #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
      #x8<-which(reachdist[,x2]== dim(data2)[1])
      if (x2 %in% dist.obj2$id[x5,]) {
        s_update_lrd=as.matrix(union(s_update_lrd,(x5+1)))
      }
      
    }
  }
}



#distdata1<-distdata1[,-c(1)]
reachdist<-reachdist[,-c(1)]

#update lrd
s_update_lof_d<-s_update_lrd
for (i in 1:length(s_update_lrd)) {
  if (length(s_update_lrd)==1){
    if (length(s_update_lrd[[1]])!=0){
      x7<-s_update_lrd[i]-1
      numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
      if (numneigh2 == k) {
          j5=seq(reachdist[1, x7],reachdist[2, x7])
          reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
          lrddata[x7+1] = reach2
          #dist.obj <- dbscan::kNN(data2, k)
          #n <- nrow(data2)
          rNN_x7 <- sapply(x7, function(i) {
            as.vector(which(dist.obj2$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                      1])
          })
          if (length(rNN_x7)==1){
            if (length(rNN_x7[[1]])!=0){
              s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
            }
          } else {
            s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
          }
          
        } else if (numneigh2!=k) {
          print("warning numneigh!=k")
        }
    }
  } else {
    x7<-s_update_lrd[i,1]-1
    numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
    if (numneigh2 == k) {
      j5=seq(reachdist[1, x7],reachdist[2, x7])
      reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
      lrddata[x7+1] = reach2
      #dist.obj <- dbscan::kNN(data2, k)
      #n <- nrow(data2)
      rNN_x7 <- sapply(x7, function(i) {
        as.vector(which(dist.obj2$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                 1])
      })
      if (length(rNN_x7)==1){
        if (length(rNN_x7[[1]])!=0){
          s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
        }
      } else {
        s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
      }
    } else if (numneigh2!=k) {
      print("warning numneigh!=k")
    }
  }
  
  
}

lrddata<-lrddata[-c(1)]


#update lof
if (length(s_update_lof_d)==1){
  if (length(s_update_lof_d[[1]])!=0){
    #s_update_lof_d<-as.matrix(s_update_lof[-c(which(s_update_lof[]== dim(data2)[1])),])
    for (i in 1:length(s_update_lof_d)) {
      x8<-s_update_lof_d[i]-1
      nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
      if (nneigh1 == k) {
        j8 = seq(0, (nneigh1 - 1))
        local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
        lof[x8+1] = local.factor1
      } else if (nneigh1!=k) {
        print("warning numneigh!=k")
      }
    }
    
  }
} else {
  #s_update_lof<-as.matrix(s_update_lof[-c(which(s_update_lof[,1]== dim(data2)[1])),])
  for (i in 1:length(s_update_lof_d)) {
    x8<-s_update_lof_d[i,1]-1
    nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
    if (nneigh1 == k) {
      j8 = seq(0, (nneigh1 - 1))
      local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
      lof[x8+1] = local.factor1
    } else if (nneigh1!=k) {
      print("warning numneigh!=k")
    }
  }
  
}


#lrddata<-lrddata[-c(1)]
lof<-lof[-c(1)]









#Normal distribution 
k<-10
data2<-data
distdata = dist.to.knn(data, k)
distdata1<-distdata
reachdist=reachability_dist(distdata, k)
lrddata = reachability(distdata, k)
p = dim(distdata)[2]
lof = rep(0, p)
for (i in 1:p) {
  nneigh = distdata[2, i] - distdata[1, i] + 1
  j = seq(0, (nneigh - 1))
  local.factor = sum(lrddata[distdata[3 + j, i]]/lrddata[i])/nneigh
  lof[i] = local.factor
}
lof


outlier_detection9 = matrix(nrow=200, ncol=4)
for (i in 1:200){
  m=i%%11
  m1=i
  if (m == 0){
    x1<- rnorm(1, mean = 0, sd = 1)
    x_1 <- (abs(x1) + 3*1) * sign(x1)
    
  } else {
    x_1=rnorm(1,mean=0, sd = 1)
  }
  outlier_detection9[i,1]<-x_1
  
  #Deletion
  #reverse neighbours of deletion point
  dist.obj1 <- dbscan::kNN(data2, k)
  #n <- nrow(data2)
  
  rNN_x_1 <- sapply(1, function(i) {
    as.vector(which(dist.obj1$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                              1])
  })
  
  data2<-as.matrix(data2[-c(1),])
  
  #Update k-distance
  
  j1=seq(3, 3 + (distdata1[2, 1] - distdata1[1, 1]))
  j2=seq(distdata1[1, 1],distdata1[2, 1])
  distdata1[j1,]=distdata1[j1,]-1
  dist.obj2 <- dbscan::kNN(data2, k)
  
  for (i in 1:length(rNN_x_1)) {
    if (length(rNN_x_1)==1){
      if (length(rNN_x_1[[1]])!=0){
        x1<-rNN_x_1[i]
        numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
        if (numneigh == k) {
          j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
          j2=seq(distdata1[1, x1],distdata1[2, x1])
          distdata1[j1,x1]=dist.obj2$id[x1-1,]
          distdata1[j2,x1]=dist.obj2$dist[x1-1,]
        } else if (numneigh!=k) {
          print("warning numneigh!=k")
        }
      }
    } else {
      x1<-rNN_x_1[i,1]
      numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
      if (numneigh == k) {
        j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
        j2=seq(distdata1[1, x1],distdata1[2, x1])
        distdata1[j1,x1]=dist.obj2$id[x1-1,]
        distdata1[j2,x1]=dist.obj2$dist[x1-1,]
      } else if (numneigh!=k) {
        print("warning numneigh!=k")
      }
    }
  }
  
  
  distdata1<-distdata1[,-c(1)]
  s_update_lrd=rNN_x_1
  j32=seq(3, 3 + (distdata1[2, 1] - distdata1[1, 1]))
  j4=seq(distdata1[1, 1],distdata1[2, 1])
  reachdist[j32,]=reachdist[j32,]-1
  #update reach_dist
  for (i in 1:length(rNN_x_1)) {
    if (length(rNN_x_1)==1) {
      if (length(rNN_x_1[[1]])!=0) {
        x2<-rNN_x_1[i]-1
        knn_x2<-as.matrix(dist.obj2$id[x2,-c(length(dist.obj2$id[x2,]))])
        numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
        if (numneigh1 == k) {
          j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
          j4=seq(distdata1[1, x2],distdata1[2, x2])
          inds<-reachdist[j3,x2+1]
          inds1<-reachdist[j4,x2+1]
          inds2<-match(dist.obj2$id[x2,],reachdist[j3,x2+1])
          reachdist[j3,x2+1]=dist.obj2$id[x2,]
          reachdist[j4,x2+1]=inds1[inds2]
          x3<-length(which(is.na(reachdist[,x2+1])))
          x4<-which(is.na(reachdist[,x2+1]))
          if (x3!=0){
            temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                                x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                     x2])
            #reach = 1/(sum(apply(temp, 2, max))/numneigh)
            
            reach_dist_x2 = apply(temp1, 2, max)
            reachdist[x4,x2+1]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
          }
        } else if (numneigh1!=k) {
          print("warning numneigh!=k")
        }
        for (i in 1:dim(knn_x2)[1]) {
          x5<-knn_x2[i,1]
          x6<-which(reachdist[,x5+1]== x2)+numneigh1
          reachdist[x6,x5+1]=distdata1[length(distdata1[,x2]),x2]
          #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
          #x8<-which(reachdist[,x2]== dim(data2)[1])
          if (x2 %in% dist.obj2$id[x5,]) {
            s_update_lrd=as.matrix(union(s_update_lrd,(x5+1)))
          }
          
        }
      }
      
    } else {
      x2<-rNN_x_1[i,1]-1
      knn_x2<-as.matrix(dist.obj2$id[x2,-c(length(dist.obj2$id[x2,]))])
      numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
      if (numneigh1 == k) {
        j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
        j4=seq(distdata1[1, x2],distdata1[2, x2])
        inds<-reachdist[j3,x2+1]
        inds1<-reachdist[j4,x2+1]
        inds2<-match(dist.obj2$id[x2,],reachdist[j3,x2+1])
        reachdist[j3,x2+1]=dist.obj2$id[x2,]
        reachdist[j4,x2+1]=inds1[inds2]
        x3<-length(which(is.na(reachdist[,x2+1])))
        x4<-which(is.na(reachdist[,x2+1]))
        if (x3!=0){
          temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                              x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                   x2])
          #reach = 1/(sum(apply(temp, 2, max))/numneigh)
          
          reach_dist_x2 = apply(temp1, 2, max)
          reachdist[x4,x2+1]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
        }
      } else if (numneigh1!=k) {
        print("warning numneigh!=k")
      }
      for (i in 1:dim(knn_x2)[1]) {
        x5<-knn_x2[i,1]
        x6<-which(reachdist[,x5+1]== x2)+numneigh1
        reachdist[x6,x5+1]=distdata1[length(distdata1[,x2]),x2]
        #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
        #x8<-which(reachdist[,x2]== dim(data2)[1])
        if (x2 %in% dist.obj2$id[x5,]) {
          s_update_lrd=as.matrix(union(s_update_lrd,(x5+1)))
        }
        
      }
    }
  }
  
  
  
  #distdata1<-distdata1[,-c(1)]
  reachdist<-reachdist[,-c(1)]
  
  #update lrd
  s_update_lof_d<-s_update_lrd
  for (i in 1:length(s_update_lrd)) {
    if (length(s_update_lrd)==1){
      if (length(s_update_lrd[[1]])!=0){
        x7<-s_update_lrd[i]-1
        numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
        if (numneigh2 == k) {
          j5=seq(reachdist[1, x7],reachdist[2, x7])
          reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
          lrddata[x7+1] = reach2
          #dist.obj <- dbscan::kNN(data2, k)
          #n <- nrow(data2)
          rNN_x7 <- sapply(x7, function(i) {
            as.vector(which(dist.obj2$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                      1])
          })
          if (length(rNN_x7)==1){
            if (length(rNN_x7[[1]])!=0){
              s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
            }
          } else {
            s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
          }
          
        } else if (numneigh2!=k) {
          print("warning numneigh!=k")
        }
      }
    } else {
      x7<-s_update_lrd[i,1]-1
      numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
      if (numneigh2 == k) {
        j5=seq(reachdist[1, x7],reachdist[2, x7])
        reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
        lrddata[x7+1] = reach2
        #dist.obj <- dbscan::kNN(data2, k)
        #n <- nrow(data2)
        rNN_x7 <- sapply(x7, function(i) {
          as.vector(which(dist.obj2$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                    1])
        })
        if (length(rNN_x7)==1){
          if (length(rNN_x7[[1]])!=0){
            s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
          }
        } else {
          s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
        }
      } else if (numneigh2!=k) {
        print("warning numneigh!=k")
      }
    }
    
    
  }
  
  lrddata<-lrddata[-c(1)]
  
  
  #update lof
  if (length(s_update_lof_d)==1){
    if (length(s_update_lof_d[[1]])!=0){
      #s_update_lof_d<-as.matrix(s_update_lof[-c(which(s_update_lof[]== dim(data2)[1])),])
      for (i in 1:length(s_update_lof_d)) {
        x8<-s_update_lof_d[i]-1
        nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
        if (nneigh1 == k) {
          j8 = seq(0, (nneigh1 - 1))
          local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
          lof[x8+1] = local.factor1
        } else if (nneigh1!=k) {
          print("warning numneigh!=k")
        }
      }
      
    }
  } else {
    #s_update_lof<-as.matrix(s_update_lof[-c(which(s_update_lof[,1]== dim(data2)[1])),])
    for (i in 1:length(s_update_lof_d)) {
      x8<-s_update_lof_d[i,1]-1
      nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
      if (nneigh1 == k) {
        j8 = seq(0, (nneigh1 - 1))
        local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
        lof[x8+1] = local.factor1
      } else if (nneigh1!=k) {
        print("warning numneigh!=k")
      }
    }
    
  }
  
  
  #lrddata<-lrddata[-c(1)]
  lof<-lof[-c(1)]
  
  
  #knn of the inserted observation
  data2=as.matrix(c(data2,x_1))
  neighdist = knneigh.vect(data2[dim(data2)[1], ], data2, k)
  distdata1= cbind(distdata1, neighdist)
  
  #reach_dist of x_1
  i=dim(data2)[1]
  j = seq(3, 3 + (distdata1[2, i] - distdata1[1, i]))
  numneigh = distdata1[2, i] - distdata1[1, i] + 1
  temp = rbind(diag(distdata1[distdata1[2, distdata1[j, 
                                                     i]], distdata1[j, i]]), distdata1[j + numneigh, 
                                                                                       i])
  #reach = 1/(sum(apply(temp, 2, max))/numneigh)
  
  reach_dist_x1 = apply(temp, 2, max) #reach_dist(x_1, knn)
  #1/(sum(reach_dist_x1)/numneigh)
  
  
  
  #update neighbors of x_1
  dist.obj <- dbscan::kNN(data2, k)
  n <- nrow(data2)
  
  rNN_x1 <- sapply(n, function(i) {
    as.vector(which(dist.obj$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                             1])
  })
  
  
  #update k_distance
  for (i in 1:length(rNN_x1)) {
    if (length(rNN_x1)==1){
      if (length(rNN_x1[[1]])!=0){
        x1<-rNN_x1[i]
        numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
        if (numneigh == k) {
          j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
          j2=seq(distdata1[1, x1],distdata1[2, x1])
          distdata1[j1,x1]=dist.obj$id[x1,]
          distdata1[j2,x1]=dist.obj$dist[x1,]
        } else if (numneigh!=k) {
          print("warning numneigh!=k")
        }
      }
    } else {
      x1<-rNN_x1[i,1]
      numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
      if (numneigh == k) {
        j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
        j2=seq(distdata1[1, x1],distdata1[2, x1])
        distdata1[j1,x1]=dist.obj$id[x1,]
        distdata1[j2,x1]=dist.obj$dist[x1,]
      } else if (numneigh!=k) {
        print("warning numneigh!=k")
      }
    }
  }
  
  #update reach_dist
  s_update=rNN_x1
  for (i in 1:length(rNN_x1)) {
    if (length(rNN_x1)==1) {
      if (length(rNN_x1[[1]])!=0) {
        x2<-rNN_x1[i]
        knn_x2<-as.matrix(dist.obj$id[x2,])
        knn_x2<-as.matrix(knn_x2[-c(which(knn_x2[,1]== dim(data2)[1])),])
        numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
        if (numneigh1 == k) {
          j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
          j4=seq(distdata1[1, x2],distdata1[2, x2])
          inds<-reachdist[j3,x2]
          inds1<-reachdist[j4,x2]
          inds2<-match(dist.obj$id[x2,],reachdist[j3,x2])
          reachdist[j3,x2]=dist.obj$id[x2,]
          reachdist[j4,x2]=inds1[inds2]
          x3<-length(which(is.na(reachdist[,x2])))
          x4<-which(is.na(reachdist[,x2]))
          if (x3!=0){
            temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                                x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                     x2])
            #reach = 1/(sum(apply(temp, 2, max))/numneigh)
            
            reach_dist_x2 = apply(temp1, 2, max)
            reachdist[x4,x2]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
          }
        } else if (numneigh1!=k) {
          print("warning numneigh!=k")
        }
        for (i in 1:dim(knn_x2)[1]) {
          x5<-knn_x2[i,1]
          x6<-which(reachdist[,x5]== x2)+numneigh1
          reachdist[x6,x5]=distdata1[length(distdata1[,x2]),x2]
          #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
          #x8<-which(reachdist[,x2]== dim(data2)[1])
          if (x2 %in% dist.obj$id[x5,]) {
            s_update=as.matrix(union(s_update,x5))
          }
          
        }
      }
      
    } else {
      x2<-rNN_x1[i,1]
      knn_x2<-as.matrix(dist.obj$id[x2,])
      knn_x2<-as.matrix(knn_x2[-c(which(knn_x2[,1]== dim(data2)[1])),])
      numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
      if (numneigh1 == k) {
        j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
        j4=seq(distdata1[1, x2],distdata1[2, x2])
        inds<-reachdist[j3,x2]
        inds1<-reachdist[j4,x2]
        inds2<-match(dist.obj$id[x2,],reachdist[j3,x2])
        reachdist[j3,x2]=dist.obj$id[x2,]
        reachdist[j4,x2]=inds1[inds2]
        x3<-length(which(is.na(reachdist[,x2])))
        x4<-which(is.na(reachdist[,x2]))
        if (x3!=0){
          temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                              x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                   x2])
          #reach = 1/(sum(apply(temp, 2, max))/numneigh)
          
          reach_dist_x2 = apply(temp1, 2, max)
          reachdist[x4,x2]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
        }
      } else if (numneigh1!=k) {
        print("warning numneigh!=k")
      }
      for (i in 1:dim(knn_x2)[1]) {
        x5<-knn_x2[i,1]
        x6<-which(reachdist[,x5]== x2)+numneigh1
        reachdist[x6,x5]=distdata1[length(distdata1[,x2]),x2]
        #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
        #x8<-which(reachdist[,x2]== dim(data2)[1])
        if (x2 %in% dist.obj$id[x5,]) {
          s_update=as.matrix(union(s_update,x5))
        }
        
      }
    }
  }
  
  if (numneigh == k) {
    ind<-dim(data2)[1]
    #j = seq(3, 3 + (distdata1[2, i] - distdata1[1, i]))
    #numneigh = distdata1[2, i] - distdata1[1, i] + 1
    #temp = rbind(diag(distdata1[distdata1[2, distdata1[j, 
    #                                                   i]], distdata1[j, i]]), distdata1[j + numneigh, 
    #                                                                                     i])
    #reach = 1/(sum(apply(temp, 2, max))/numneigh)
    reachdist=cbind(reachdist,neighdist)
    ind1<-seq(reachdist[1, ind],reachdist[2, ind])
    #reach_dist_x1 = apply(temp, 2, max) #reach_dist(x_1, knn)
    reachdist[ind1,ind]=reach_dist_x1
  } else if (numneigh1!=k) {
    print("warning numneigh!=k")
  }
  
  
  
  
  #update lrd
  s_update_lof<-s_update
  for (i in 1:length(s_update)) {
    if (length(s_update)==1){
      if (length(s_update[[1]])!=0){
        x7<-s_update[i]
        numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
        if (numneigh2 == k) {
          j5=seq(reachdist[1, x7],reachdist[2, x7])
          reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
          lrddata[x7] = reach2
          #dist.obj <- dbscan::kNN(data2, k)
          #n <- nrow(data2)
          rNN_x7 <- sapply(x7, function(i) {
            as.vector(which(dist.obj$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                     1])
          })
          if (length(rNN_x7)==1){
            if (length(rNN_x7[[1]])!=0){
              s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
            }
          } else {
            s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
          }
          
        } else if (numneigh2!=k) {
          print("warning numneigh!=k")
        }
      }
    } else {
      x7<-s_update[i,1]
      numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
      if (numneigh2 == k) {
        j5=seq(reachdist[1, x7],reachdist[2, x7])
        reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
        lrddata[x7] = reach2
        #dist.obj <- dbscan::kNN(data2, k)
        #n <- nrow(data2)
        rNN_x7 <- sapply(x7, function(i) {
          as.vector(which(dist.obj$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                   1])
        })
        if (length(rNN_x7)==1){
          if (length(rNN_x7[[1]])!=0){
            s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
          }
        } else {
          s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
        }
      } else if (numneigh2!=k) {
        print("warning numneigh!=k")
      }
    }
    
    
  }
  
  #calculate lrd for x_1
  x9<-dim(data2)[1]
  numneigh3 = reachdist[2, x9] - reachdist[1, x9] + 1
  if (numneigh3 == k) {
    j6=seq(reachdist[1, x9],reachdist[2, x9])
    reach3 = 1/(sum(reachdist[j6,x9])/numneigh3)
    lrddata[x9] = reach3
    #dist.obj <- dbscan::kNN(data2, k)
    #n <- nrow(data2)
  } else if (numneigh3!=k) {
    print("warning numneigh!=k")
  }
  
  
  
  #update lof
  if (length(s_update_lof)==1){
    if (length(s_update_lof[[1]])!=0){
      if (length(which(s_update_lof[]== dim(data2)[1]))!=0){
        s_update_lof<-as.matrix(s_update_lof[-c(which(s_update_lof[]== dim(data2)[1])),])
      }
      
      for (i in 1:length(s_update_lof)) {
        x8<-s_update_lof[i]
        nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
        if (nneigh1 == k) {
          j8 = seq(0, (nneigh1 - 1))
          local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
          lof[x8] = local.factor1
        } else if (nneigh1!=k) {
          print("warning numneigh!=k")
        }
      }
      
    }
  } else if (length(s_update_lof)>1){
    if (length(which(s_update_lof[,1]== dim(data2)[1]))!=0){
      s_update_lof<-as.matrix(s_update_lof[-c(which(s_update_lof[,1]== dim(data2)[1])),])
    }
    
    for (i in 1:length(s_update_lof)) {
      x8<-s_update_lof[i,1]
      nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
      if (nneigh1 == k) {
        j8 = seq(0, (nneigh1 - 1))
        local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
        lof[x8] = local.factor1
      } else if (nneigh1!=k) {
        print("warning numneigh!=k")
      }
    }
    
  }
  
  outlier_detection9[m1,4]<-length(union(s_update_lof, s_update_lof_d))
  #calculate lof for x_1
  x9<-dim(data2)[1]
  nneigh2 = distdata1[2, x9] - distdata1[1, x9] + 1
  if (nneigh2 == k) {
    j7 = seq(0, (nneigh2 - 1))
    local.factor2 = sum(lrddata[distdata1[3 + j7, x9]]/lrddata[x9])/nneigh2
    lof[x9] = local.factor2
    if (lof[x9]>1.7){
      outlier_detection9[m1,2]<-1
      outlier_detection9[m1,3]<-lof[x9]
      
    } else {
      outlier_detection9[m1,2]<-0
      outlier_detection9[m1,3]<-lof[x9]
    }
  } else if (nneigh2!=k) {
    print("warning numneigh!=k")
  }
  
  
}

mean(outlier_detection9[,4])
length(which(outlier_detection9[,2]==1))
which(outlier_detection9[,2]==1)
#which(lof[201:400] > 1.35)
#outlier3<-union(which(outlier_detection2[,2]==1),which(lof[201:400] > 1.35) )
outlier3<-which(outlier_detection9[,2]==1)
which(c(1:200)%%11==0)
#length(which(lof[201:400] > 1.2))

TP_V<-length(which(outlier3 %in% which(c(1:200)%%11==0)))
FN_V<-18-TP_V
FP_V<-length(which(outlier_detection9[,2]==1))-TP_V
TN_V<-200-length(which(c(1:200)%%11==0))-FP_V
SE_V<-TP_V/(TP_V+FN_V) #true positive rate
SP_V<-TN_V/(FP_V+TN_V) #true negative rate
G_mean_lof_SWLIOF_S<-sqrt(SE_V*SP_V)
H_mean_lof_SWILOF_S<-2*SE_V*SP_V/(SE_V+SP_V)


outliers2<-order(data1_lof, decreasing=T)[1:length(which(data1_lof > 2))]
sort(outliers2)
sort(ind)

outlier_detection12<-as.data.frame(outlier_detection9)
outlier_detection12[,5]=(1:200)
outlier_detection12$V2<-as.factor(outlier_detection12$V2)
colnames(outlier_detection12)[2] <- "Outliers"
# plot
p5 <- ggplot(outlier_detection12, aes(x=V5,y=V1, label=V5)) +
  geom_point(aes(colour=Outliers)) + 
  #geom_point(aes(y=deaths), colour='red') +
  #geom_errorbar(aes(ymin=V4, ymax=V3)) +
  geom_text(aes(label=ifelse(V3>1.7,as.character(V5),'')),hjust=0,vjust=0)+
  theme_bw() +            
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))+
  #scale_y_continuous(expand=c(0,0), lim=c(0,1e3)) +
  labs(x='timepoint of streaming data', y='value of observation', title='SW-ILOF, k=10, alpha=1.7, W=50')
#coord_cartesian(ylim=c(0,300)) +
#facet_grid(~age_group)
p5

distdata1<-distdata
reachdist=reachability_dist(distdata, k)
lrddata = reachability(distdata, k)
p = dim(distdata)[2]
lof = rep(0, p)
for (i in 1:p) {
  nneigh = distdata[2, i] - distdata[1, i] + 1
  j = seq(0, (nneigh - 1))
  local.factor = sum(lrddata[distdata[3 + j, i]]/lrddata[i])/nneigh
  lof[i] = local.factor
}
lof

#performance measures
my_rnorm_SW_ILOF_S <- function(data=data1, k=10,w=151, a1=2.9 ){
  data2<-as.matrix(data[w:200])
  distdata = dist.to.knn(data2, k)
  distdata1<-distdata
  reachdist=reachability_dist(distdata, k)
  lrddata = reachability(distdata, k)
  p = dim(distdata)[2]
  lof = rep(0, p)
  for (i in 1:p) {
    nneigh = distdata[2, i] - distdata[1, i] + 1
    j = seq(0, (nneigh - 1))
    local.factor = sum(lrddata[distdata[3 + j, i]]/lrddata[i])/nneigh
    lof[i] = local.factor
  }
  
  
  outlier_detection9 = matrix(nrow=200, ncol=4)
  for (i in 1:200){
    m=i%%11
    m1=i
    if (m == 0){
      x1<- rnorm(1, mean = 0, sd = 1)
      x_1 <- (abs(x1) + 3*1) * sign(x1)
      
    } else {
      x_1=rnorm(1,mean=0, sd = 1)
    }
    outlier_detection9[i,1]<-x_1
    
    #Deletion
    #reverse neighbours of deletion point
    dist.obj1 <- dbscan::kNN(data2, k)
    #n <- nrow(data2)
    
    rNN_x_1 <- sapply(1, function(i) {
      as.vector(which(dist.obj1$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                1])
    })
    
    data2<-as.matrix(data2[-c(1),])
    
    #Update k-distance
    
    j1=seq(3, 3 + (distdata1[2, 1] - distdata1[1, 1]))
    j2=seq(distdata1[1, 1],distdata1[2, 1])
    distdata1[j1,]=distdata1[j1,]-1
    dist.obj2 <- dbscan::kNN(data2, k)
    
    for (i in 1:length(rNN_x_1)) {
      if (length(rNN_x_1)==1){
        if (length(rNN_x_1[[1]])!=0){
          x1<-rNN_x_1[i]
          numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
          if (numneigh == k) {
            j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
            j2=seq(distdata1[1, x1],distdata1[2, x1])
            distdata1[j1,x1]=dist.obj2$id[x1-1,]
            distdata1[j2,x1]=dist.obj2$dist[x1-1,]
          } else if (numneigh!=k) {
            print("warning numneigh!=k")
          }
        }
      } else {
        x1<-rNN_x_1[i,1]
        numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
        if (numneigh == k) {
          j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
          j2=seq(distdata1[1, x1],distdata1[2, x1])
          distdata1[j1,x1]=dist.obj2$id[x1-1,]
          distdata1[j2,x1]=dist.obj2$dist[x1-1,]
        } else if (numneigh!=k) {
          print("warning numneigh!=k")
        }
      }
    }
    
    
    distdata1<-distdata1[,-c(1)]
    s_update_lrd=rNN_x_1
    j32=seq(3, 3 + (distdata1[2, 1] - distdata1[1, 1]))
    j4=seq(distdata1[1, 1],distdata1[2, 1])
    reachdist[j32,]=reachdist[j32,]-1
    #update reach_dist
    for (i in 1:length(rNN_x_1)) {
      if (length(rNN_x_1)==1) {
        if (length(rNN_x_1[[1]])!=0) {
          x2<-rNN_x_1[i]-1
          knn_x2<-as.matrix(dist.obj2$id[x2,-c(length(dist.obj2$id[x2,]))])
          numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
          if (numneigh1 == k) {
            j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
            j4=seq(distdata1[1, x2],distdata1[2, x2])
            inds<-reachdist[j3,x2+1]
            inds1<-reachdist[j4,x2+1]
            inds2<-match(dist.obj2$id[x2,],reachdist[j3,x2+1])
            reachdist[j3,x2+1]=dist.obj2$id[x2,]
            reachdist[j4,x2+1]=inds1[inds2]
            x3<-length(which(is.na(reachdist[,x2+1])))
            x4<-which(is.na(reachdist[,x2+1]))
            if (x3!=0){
              temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                                  x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                       x2])
              #reach = 1/(sum(apply(temp, 2, max))/numneigh)
              
              reach_dist_x2 = apply(temp1, 2, max)
              reachdist[x4,x2+1]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
            }
          } else if (numneigh1!=k) {
            print("warning numneigh!=k")
          }
          for (i in 1:dim(knn_x2)[1]) {
            x5<-knn_x2[i,1]
            x6<-which(reachdist[,x5+1]== x2)+numneigh1
            reachdist[x6,x5+1]=distdata1[length(distdata1[,x2]),x2]
            #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
            #x8<-which(reachdist[,x2]== dim(data2)[1])
            if (x2 %in% dist.obj2$id[x5,]) {
              s_update_lrd=as.matrix(union(s_update_lrd,(x5+1)))
            }
            
          }
        }
        
      } else {
        x2<-rNN_x_1[i,1]-1
        knn_x2<-as.matrix(dist.obj2$id[x2,-c(length(dist.obj2$id[x2,]))])
        numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
        if (numneigh1 == k) {
          j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
          j4=seq(distdata1[1, x2],distdata1[2, x2])
          inds<-reachdist[j3,x2+1]
          inds1<-reachdist[j4,x2+1]
          inds2<-match(dist.obj2$id[x2,],reachdist[j3,x2+1])
          reachdist[j3,x2+1]=dist.obj2$id[x2,]
          reachdist[j4,x2+1]=inds1[inds2]
          x3<-length(which(is.na(reachdist[,x2+1])))
          x4<-which(is.na(reachdist[,x2+1]))
          if (x3!=0){
            temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                                x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                     x2])
            #reach = 1/(sum(apply(temp, 2, max))/numneigh)
            
            reach_dist_x2 = apply(temp1, 2, max)
            reachdist[x4,x2+1]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
          }
        } else if (numneigh1!=k) {
          print("warning numneigh!=k")
        }
        for (i in 1:dim(knn_x2)[1]) {
          x5<-knn_x2[i,1]
          x6<-which(reachdist[,x5+1]== x2)+numneigh1
          reachdist[x6,x5+1]=distdata1[length(distdata1[,x2]),x2]
          #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
          #x8<-which(reachdist[,x2]== dim(data2)[1])
          if (x2 %in% dist.obj2$id[x5,]) {
            s_update_lrd=as.matrix(union(s_update_lrd,(x5+1)))
          }
          
        }
      }
    }
    
    
    
    #distdata1<-distdata1[,-c(1)]
    reachdist<-reachdist[,-c(1)]
    
    #update lrd
    s_update_lof_d<-s_update_lrd
    for (i in 1:length(s_update_lrd)) {
      if (length(s_update_lrd)==1){
        if (length(s_update_lrd[[1]])!=0){
          x7<-s_update_lrd[i]-1
          numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
          if (numneigh2 == k) {
            j5=seq(reachdist[1, x7],reachdist[2, x7])
            reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
            lrddata[x7+1] = reach2
            #dist.obj <- dbscan::kNN(data2, k)
            #n <- nrow(data2)
            rNN_x7 <- sapply(x7, function(i) {
              as.vector(which(dist.obj2$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                        1])
            })
            if (length(rNN_x7)==1){
              if (length(rNN_x7[[1]])!=0){
                s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
              }
            } else {
              s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
            }
            
          } else if (numneigh2!=k) {
            print("warning numneigh!=k")
          }
        }
      } else {
        x7<-s_update_lrd[i,1]-1
        numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
        if (numneigh2 == k) {
          j5=seq(reachdist[1, x7],reachdist[2, x7])
          reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
          lrddata[x7+1] = reach2
          #dist.obj <- dbscan::kNN(data2, k)
          #n <- nrow(data2)
          rNN_x7 <- sapply(x7, function(i) {
            as.vector(which(dist.obj2$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                      1])
          })
          if (length(rNN_x7)==1){
            if (length(rNN_x7[[1]])!=0){
              s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
            }
          } else {
            s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
          }
        } else if (numneigh2!=k) {
          print("warning numneigh!=k")
        }
      }
      
      
    }
    
    lrddata<-lrddata[-c(1)]
    
    
    #update lof
    if (length(s_update_lof_d)==1){
      if (length(s_update_lof_d[[1]])!=0){
        #s_update_lof_d<-as.matrix(s_update_lof[-c(which(s_update_lof[]== dim(data2)[1])),])
        for (i in 1:length(s_update_lof_d)) {
          x8<-s_update_lof_d[i]-1
          nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
          if (nneigh1 == k) {
            j8 = seq(0, (nneigh1 - 1))
            local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
            lof[x8+1] = local.factor1
          } else if (nneigh1!=k) {
            print("warning numneigh!=k")
          }
        }
        
      }
    } else {
      #s_update_lof<-as.matrix(s_update_lof[-c(which(s_update_lof[,1]== dim(data2)[1])),])
      for (i in 1:length(s_update_lof_d)) {
        x8<-s_update_lof_d[i,1]-1
        nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
        if (nneigh1 == k) {
          j8 = seq(0, (nneigh1 - 1))
          local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
          lof[x8+1] = local.factor1
        } else if (nneigh1!=k) {
          print("warning numneigh!=k")
        }
      }
      
    }
    
    
    #lrddata<-lrddata[-c(1)]
    lof<-lof[-c(1)]
    
    
    #knn of the inserted observation
    data2=as.matrix(c(data2,x_1))
    neighdist = knneigh.vect(data2[dim(data2)[1], ], data2, k)
    distdata1= cbind(distdata1, neighdist)
    
    #reach_dist of x_1
    i=dim(data2)[1]
    j = seq(3, 3 + (distdata1[2, i] - distdata1[1, i]))
    numneigh = distdata1[2, i] - distdata1[1, i] + 1
    temp = rbind(diag(distdata1[distdata1[2, distdata1[j, 
                                                       i]], distdata1[j, i]]), distdata1[j + numneigh, 
                                                                                         i])
    #reach = 1/(sum(apply(temp, 2, max))/numneigh)
    
    reach_dist_x1 = apply(temp, 2, max) #reach_dist(x_1, knn)
    #1/(sum(reach_dist_x1)/numneigh)
    
    
    
    #update neighbors of x_1
    dist.obj <- dbscan::kNN(data2, k)
    n <- nrow(data2)
    
    rNN_x1 <- sapply(n, function(i) {
      as.vector(which(dist.obj$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                               1])
    })
    
    
    #update k_distance
    for (i in 1:length(rNN_x1)) {
      if (length(rNN_x1)==1){
        if (length(rNN_x1[[1]])!=0){
          x1<-rNN_x1[i]
          numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
          if (numneigh == k) {
            j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
            j2=seq(distdata1[1, x1],distdata1[2, x1])
            distdata1[j1,x1]=dist.obj$id[x1,]
            distdata1[j2,x1]=dist.obj$dist[x1,]
          } else if (numneigh!=k) {
            print("warning numneigh!=k")
          }
        }
      } else {
        x1<-rNN_x1[i,1]
        numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
        if (numneigh == k) {
          j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
          j2=seq(distdata1[1, x1],distdata1[2, x1])
          distdata1[j1,x1]=dist.obj$id[x1,]
          distdata1[j2,x1]=dist.obj$dist[x1,]
        } else if (numneigh!=k) {
          print("warning numneigh!=k")
        }
      }
    }
    
    #update reach_dist
    s_update=rNN_x1
    for (i in 1:length(rNN_x1)) {
      if (length(rNN_x1)==1) {
        if (length(rNN_x1[[1]])!=0) {
          x2<-rNN_x1[i]
          knn_x2<-as.matrix(dist.obj$id[x2,])
          knn_x2<-as.matrix(knn_x2[-c(which(knn_x2[,1]== dim(data2)[1])),])
          numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
          if (numneigh1 == k) {
            j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
            j4=seq(distdata1[1, x2],distdata1[2, x2])
            inds<-reachdist[j3,x2]
            inds1<-reachdist[j4,x2]
            inds2<-match(dist.obj$id[x2,],reachdist[j3,x2])
            reachdist[j3,x2]=dist.obj$id[x2,]
            reachdist[j4,x2]=inds1[inds2]
            x3<-length(which(is.na(reachdist[,x2])))
            x4<-which(is.na(reachdist[,x2]))
            if (x3!=0){
              temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                                  x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                       x2])
              #reach = 1/(sum(apply(temp, 2, max))/numneigh)
              
              reach_dist_x2 = apply(temp1, 2, max)
              reachdist[x4,x2]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
            }
          } else if (numneigh1!=k) {
            print("warning numneigh!=k")
          }
          for (i in 1:dim(knn_x2)[1]) {
            x5<-knn_x2[i,1]
            x6<-which(reachdist[,x5]== x2)+numneigh1
            reachdist[x6,x5]=distdata1[length(distdata1[,x2]),x2]
            #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
            #x8<-which(reachdist[,x2]== dim(data2)[1])
            if (x2 %in% dist.obj$id[x5,]) {
              s_update=as.matrix(union(s_update,x5))
            }
            
          }
        }
        
      } else {
        x2<-rNN_x1[i,1]
        knn_x2<-as.matrix(dist.obj$id[x2,])
        knn_x2<-as.matrix(knn_x2[-c(which(knn_x2[,1]== dim(data2)[1])),])
        numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
        if (numneigh1 == k) {
          j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
          j4=seq(distdata1[1, x2],distdata1[2, x2])
          inds<-reachdist[j3,x2]
          inds1<-reachdist[j4,x2]
          inds2<-match(dist.obj$id[x2,],reachdist[j3,x2])
          reachdist[j3,x2]=dist.obj$id[x2,]
          reachdist[j4,x2]=inds1[inds2]
          x3<-length(which(is.na(reachdist[,x2])))
          x4<-which(is.na(reachdist[,x2]))
          if (x3!=0){
            temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                                x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                     x2])
            #reach = 1/(sum(apply(temp, 2, max))/numneigh)
            
            reach_dist_x2 = apply(temp1, 2, max)
            reachdist[x4,x2]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
          }
        } else if (numneigh1!=k) {
          print("warning numneigh!=k")
        }
        for (i in 1:dim(knn_x2)[1]) {
          x5<-knn_x2[i,1]
          x6<-which(reachdist[,x5]== x2)+numneigh1
          reachdist[x6,x5]=distdata1[length(distdata1[,x2]),x2]
          #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
          #x8<-which(reachdist[,x2]== dim(data2)[1])
          if (x2 %in% dist.obj$id[x5,]) {
            s_update=as.matrix(union(s_update,x5))
          }
          
        }
      }
    }
    
    if (numneigh == k) {
      ind<-dim(data2)[1]
      #j = seq(3, 3 + (distdata1[2, i] - distdata1[1, i]))
      #numneigh = distdata1[2, i] - distdata1[1, i] + 1
      #temp = rbind(diag(distdata1[distdata1[2, distdata1[j, 
      #                                                   i]], distdata1[j, i]]), distdata1[j + numneigh, 
      #                                                                                     i])
      #reach = 1/(sum(apply(temp, 2, max))/numneigh)
      reachdist=cbind(reachdist,neighdist)
      ind1<-seq(reachdist[1, ind],reachdist[2, ind])
      #reach_dist_x1 = apply(temp, 2, max) #reach_dist(x_1, knn)
      reachdist[ind1,ind]=reach_dist_x1
    } else if (numneigh1!=k) {
      print("warning numneigh!=k")
    }
    
    
    
    
    #update lrd
    s_update_lof<-s_update
    for (i in 1:length(s_update)) {
      if (length(s_update)==1){
        if (length(s_update[[1]])!=0){
          x7<-s_update[i]
          numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
          if (numneigh2 == k) {
            j5=seq(reachdist[1, x7],reachdist[2, x7])
            reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
            lrddata[x7] = reach2
            #dist.obj <- dbscan::kNN(data2, k)
            #n <- nrow(data2)
            rNN_x7 <- sapply(x7, function(i) {
              as.vector(which(dist.obj$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                       1])
            })
            if (length(rNN_x7)==1){
              if (length(rNN_x7[[1]])!=0){
                s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
              }
            } else {
              s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
            }
            
          } else if (numneigh2!=k) {
            print("warning numneigh!=k")
          }
        }
      } else {
        x7<-s_update[i,1]
        numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
        if (numneigh2 == k) {
          j5=seq(reachdist[1, x7],reachdist[2, x7])
          reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
          lrddata[x7] = reach2
          #dist.obj <- dbscan::kNN(data2, k)
          #n <- nrow(data2)
          rNN_x7 <- sapply(x7, function(i) {
            as.vector(which(dist.obj$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                     1])
          })
          if (length(rNN_x7)==1){
            if (length(rNN_x7[[1]])!=0){
              s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
            }
          } else {
            s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
          }
        } else if (numneigh2!=k) {
          print("warning numneigh!=k")
        }
      }
      
      
    }
    
    #calculate lrd for x_1
    x9<-dim(data2)[1]
    numneigh3 = reachdist[2, x9] - reachdist[1, x9] + 1
    if (numneigh3 == k) {
      j6=seq(reachdist[1, x9],reachdist[2, x9])
      reach3 = 1/(sum(reachdist[j6,x9])/numneigh3)
      lrddata[x9] = reach3
      #dist.obj <- dbscan::kNN(data2, k)
      #n <- nrow(data2)
    } else if (numneigh3!=k) {
      print("warning numneigh!=k")
    }
    
    
    
    #update lof
    if (length(s_update_lof)==1){
      if (length(s_update_lof[[1]])!=0){
        if (length(which(s_update_lof[]== dim(data2)[1]))!=0){
          s_update_lof<-as.matrix(s_update_lof[-c(which(s_update_lof[]== dim(data2)[1])),])
        }
        
        for (i in 1:length(s_update_lof)) {
          x8<-s_update_lof[i]
          nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
          if (nneigh1 == k) {
            j8 = seq(0, (nneigh1 - 1))
            local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
            lof[x8] = local.factor1
          } else if (nneigh1!=k) {
            print("warning numneigh!=k")
          }
        }
        
      }
    } else if (length(s_update_lof)>1){
      if (length(which(s_update_lof[,1]== dim(data2)[1]))!=0){
        s_update_lof<-as.matrix(s_update_lof[-c(which(s_update_lof[,1]== dim(data2)[1])),])
      }
      
      for (i in 1:length(s_update_lof)) {
        x8<-s_update_lof[i,1]
        nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
        if (nneigh1 == k) {
          j8 = seq(0, (nneigh1 - 1))
          local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
          lof[x8] = local.factor1
        } else if (nneigh1!=k) {
          print("warning numneigh!=k")
        }
      }
      
    }
    
    outlier_detection9[m1,4]<-length(union(s_update_lof, s_update_lof_d))
    #calculate lof for x_1
    x9<-dim(data2)[1]
    nneigh2 = distdata1[2, x9] - distdata1[1, x9] + 1
    if (nneigh2 == k) {
      j7 = seq(0, (nneigh2 - 1))
      local.factor2 = sum(lrddata[distdata1[3 + j7, x9]]/lrddata[x9])/nneigh2
      lof[x9] = local.factor2
      if (lof[x9]>a1){
        outlier_detection9[m1,2]<-1
        outlier_detection9[m1,3]<-lof[x9]
        
      } else {
        outlier_detection9[m1,2]<-0
        outlier_detection9[m1,3]<-lof[x9]
      }
    } else if (nneigh2!=k) {
      print("warning numneigh!=k")
    }
    
    
  }
  
  mean(outlier_detection9[,4])
  length(which(outlier_detection9[,2]==1))
  which(outlier_detection9[,2]==1)
  #which(lof[201:400] > 1.35)
  #outlier3<-union(which(outlier_detection2[,2]==1),which(lof[201:400] > 1.35) )
  outlier3<-which(outlier_detection9[,2]==1)
  which(c(1:200)%%11==0)
  #length(which(lof[201:400] > 1.2))
  
  TP_V<-length(which(outlier3 %in% which(c(1:200)%%11==0)))
  FN_V<-18-TP_V
  FP_V<-length(which(outlier_detection9[,2]==1))-TP_V
  TN_V<-200-length(which(c(1:200)%%11==0))-FP_V
  SE_V<-TP_V/(TP_V+FN_V) #true positive rate
  SP_V<-TN_V/(FP_V+TN_V) #true negative rate
  G_mean_lof_SWLIOF_S<-sqrt(SE_V*SP_V)
  H_mean_lof_SWILOF_S<-2*SE_V*SP_V/(SE_V+SP_V)
  
  
  
  
  c(G_mean_lof_SWLIOF_S, H_mean_lof_SWILOF_S)
}


# grid over which we will perform the hyperparameter search:
hparam_grid1 <- as.data.frame(expand.grid(k=seq(10, 20, by=2), a1=seq(1.5, 4, by=0.4),w=seq(141, 181, by=10)))

# to store the OOB estimates of the MSE
oob_mses1 <- rep(0.0, nrow(hparam_grid1))

# perform the gridsearch
for(hparam_idx in 1:nrow(hparam_grid1)) {
  # train candidate model
  this_k <- hparam_grid1[hparam_idx, 1]
  this_a1 <- hparam_grid1[hparam_idx, 2]
  this_w <- hparam_grid1[hparam_idx, 3]
  test2<-my_rnorm_SW_ILOF_S(data=data1, k=this_k,w=this_w, a1=this_a1)[2]
  #rf <- randomForest(x_train, y_train, mtry=this_mtry, maxnodes=this_maxnodes)
  
  # calculate H-mean
  oob_mses1[hparam_idx] <- test2
}

# select the best model (that which has the minimum OOB MSE)
best_hparam_set2 <- hparam_grid1[which.max(oob_mses1),]



#my_rnorm_streaming_boxplot_S_performance <- function(data1=data1, q_lower=q_lower, n_p_lower=n_p_lower, n_d_lower=n_d_lower, dn_d_lower=dn_d_lower,q_upper=q_upper, n_p_upper=n_p_upper, n_d_upper=n_d_upper, dn_d_upper=dn_d_upper){
outlier_detection91 = matrix(nrow=100, ncol=2)
for (i in 1:100) {
  test1<-my_rnorm_SW_ILOF_S(data=data1, k=20,w=151, a1=1.9)
  outlier_detection91[i,1]<-test1[1]
  outlier_detection91[i,2]<-test1[2]
}

mean(outlier_detection91[,1])
mean(outlier_detection91[,2])
mean(outlier_detection9[,1])
mean(outlier_detection9[,2])





#Normal distribution with variations
k<-10
#data=data1
data=data1[151:200]
data=as.matrix(data)
data2<-data
distdata = dist.to.knn(data, k)
distdata1<-distdata
reachdist=reachability_dist(distdata, k)
lrddata = reachability(distdata, k)
p = dim(distdata)[2]
lof = rep(0, p)
for (i in 1:p) {
  nneigh = distdata[2, i] - distdata[1, i] + 1
  j = seq(0, (nneigh - 1))
  local.factor = sum(lrddata[distdata[3 + j, i]]/lrddata[i])/nneigh
  lof[i] = local.factor
}
lof


outlier_detection13 = matrix(nrow=300, ncol=4)
for (i in 1:300){
  m=i%%15
  m1=i
  
  if (m == 0){
    if (m1%%2==0 & m1>=100){
      x1<- rnorm(1, mean = 0.02*i, sd = 1)
      x_1 <- x1- 5*1
    }else{
      x1<- rnorm(1, mean = 0.02*i, sd = 1)
      x_1 <- (abs(x1) + 5*1) * sign(x1)
    }
  } else {
    x_1=rnorm(1,mean=0.02*i, sd = 1)
  }
  outlier_detection13[i,1]<-x_1
  
  #Deletion
  #reverse neighbours of deletion point
  dist.obj1 <- dbscan::kNN(data2, k)
  #n <- nrow(data2)
  
  rNN_x_1 <- sapply(1, function(i) {
    as.vector(which(dist.obj1$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                              1])
  })
  
  data2<-as.matrix(data2[-c(1),])
  
  #Update k-distance
  
  j1=seq(3, 3 + (distdata1[2, 1] - distdata1[1, 1]))
  j2=seq(distdata1[1, 1],distdata1[2, 1])
  distdata1[j1,]=distdata1[j1,]-1
  dist.obj2 <- dbscan::kNN(data2, k)
  
  for (i in 1:length(rNN_x_1)) {
    if (length(rNN_x_1)==1){
      if (length(rNN_x_1[[1]])!=0){
        x1<-rNN_x_1[i]
        numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
        if (numneigh == k) {
          j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
          j2=seq(distdata1[1, x1],distdata1[2, x1])
          distdata1[j1,x1]=dist.obj2$id[x1-1,]
          distdata1[j2,x1]=dist.obj2$dist[x1-1,]
        } else if (numneigh!=k) {
          print("warning numneigh!=k")
        }
      }
    } else {
      x1<-rNN_x_1[i,1]
      numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
      if (numneigh == k) {
        j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
        j2=seq(distdata1[1, x1],distdata1[2, x1])
        distdata1[j1,x1]=dist.obj2$id[x1-1,]
        distdata1[j2,x1]=dist.obj2$dist[x1-1,]
      } else if (numneigh!=k) {
        print("warning numneigh!=k")
      }
    }
  }
  
  
  distdata1<-distdata1[,-c(1)]
  s_update_lrd=rNN_x_1
  j32=seq(3, 3 + (distdata1[2, 1] - distdata1[1, 1]))
  j4=seq(distdata1[1, 1],distdata1[2, 1])
  reachdist[j32,]=reachdist[j32,]-1
  #update reach_dist
  for (i in 1:length(rNN_x_1)) {
    if (length(rNN_x_1)==1) {
      if (length(rNN_x_1[[1]])!=0) {
        x2<-rNN_x_1[i]-1
        knn_x2<-as.matrix(dist.obj2$id[x2,-c(length(dist.obj2$id[x2,]))])
        numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
        if (numneigh1 == k) {
          j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
          j4=seq(distdata1[1, x2],distdata1[2, x2])
          inds<-reachdist[j3,x2+1]
          inds1<-reachdist[j4,x2+1]
          inds2<-match(dist.obj2$id[x2,],reachdist[j3,x2+1])
          reachdist[j3,x2+1]=dist.obj2$id[x2,]
          reachdist[j4,x2+1]=inds1[inds2]
          x3<-length(which(is.na(reachdist[,x2+1])))
          x4<-which(is.na(reachdist[,x2+1]))
          if (x3!=0){
            temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                                x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                     x2])
            #reach = 1/(sum(apply(temp, 2, max))/numneigh)
            
            reach_dist_x2 = apply(temp1, 2, max)
            reachdist[x4,x2+1]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
          }
        } else if (numneigh1!=k) {
          print("warning numneigh!=k")
        }
        for (i in 1:dim(knn_x2)[1]) {
          x5<-knn_x2[i,1]
          x6<-which(reachdist[,x5+1]== x2)+numneigh1
          reachdist[x6,x5+1]=distdata1[length(distdata1[,x2]),x2]
          #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
          #x8<-which(reachdist[,x2]== dim(data2)[1])
          if (x2 %in% dist.obj2$id[x5,]) {
            s_update_lrd=as.matrix(union(s_update_lrd,(x5+1)))
          }
          
        }
      }
      
    } else {
      x2<-rNN_x_1[i,1]-1
      knn_x2<-as.matrix(dist.obj2$id[x2,-c(length(dist.obj2$id[x2,]))])
      numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
      if (numneigh1 == k) {
        j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
        j4=seq(distdata1[1, x2],distdata1[2, x2])
        inds<-reachdist[j3,x2+1]
        inds1<-reachdist[j4,x2+1]
        inds2<-match(dist.obj2$id[x2,],reachdist[j3,x2+1])
        reachdist[j3,x2+1]=dist.obj2$id[x2,]
        reachdist[j4,x2+1]=inds1[inds2]
        x3<-length(which(is.na(reachdist[,x2+1])))
        x4<-which(is.na(reachdist[,x2+1]))
        if (x3!=0){
          temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                              x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                   x2])
          #reach = 1/(sum(apply(temp, 2, max))/numneigh)
          
          reach_dist_x2 = apply(temp1, 2, max)
          reachdist[x4,x2+1]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
        }
      } else if (numneigh1!=k) {
        print("warning numneigh!=k")
      }
      for (i in 1:dim(knn_x2)[1]) {
        x5<-knn_x2[i,1]
        x6<-which(reachdist[,x5+1]== x2)+numneigh1
        reachdist[x6,x5+1]=distdata1[length(distdata1[,x2]),x2]
        #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
        #x8<-which(reachdist[,x2]== dim(data2)[1])
        if (x2 %in% dist.obj2$id[x5,]) {
          s_update_lrd=as.matrix(union(s_update_lrd,(x5+1)))
        }
        
      }
    }
  }
  
  
  
  #distdata1<-distdata1[,-c(1)]
  reachdist<-reachdist[,-c(1)]
  
  #update lrd
  s_update_lof_d<-s_update_lrd
  for (i in 1:length(s_update_lrd)) {
    if (length(s_update_lrd)==1){
      if (length(s_update_lrd[[1]])!=0){
        x7<-s_update_lrd[i]-1
        numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
        if (numneigh2 == k) {
          j5=seq(reachdist[1, x7],reachdist[2, x7])
          reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
          lrddata[x7+1] = reach2
          #dist.obj <- dbscan::kNN(data2, k)
          #n <- nrow(data2)
          rNN_x7 <- sapply(x7, function(i) {
            as.vector(which(dist.obj2$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                      1])
          })
          if (length(rNN_x7)==1){
            if (length(rNN_x7[[1]])!=0){
              s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
            }
          } else {
            s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
          }
          
        } else if (numneigh2!=k) {
          print("warning numneigh!=k")
        }
      }
    } else {
      x7<-s_update_lrd[i,1]-1
      numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
      if (numneigh2 == k) {
        j5=seq(reachdist[1, x7],reachdist[2, x7])
        reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
        lrddata[x7+1] = reach2
        #dist.obj <- dbscan::kNN(data2, k)
        #n <- nrow(data2)
        rNN_x7 <- sapply(x7, function(i) {
          as.vector(which(dist.obj2$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                    1])
        })
        if (length(rNN_x7)==1){
          if (length(rNN_x7[[1]])!=0){
            s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
          }
        } else {
          s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
        }
      } else if (numneigh2!=k) {
        print("warning numneigh!=k")
      }
    }
    
    
  }
  
  lrddata<-lrddata[-c(1)]
  
  
  #update lof
  if (length(s_update_lof_d)==1){
    if (length(s_update_lof_d[[1]])!=0){
      #s_update_lof_d<-as.matrix(s_update_lof[-c(which(s_update_lof[]== dim(data2)[1])),])
      for (i in 1:length(s_update_lof_d)) {
        x8<-s_update_lof_d[i]-1
        nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
        if (nneigh1 == k) {
          j8 = seq(0, (nneigh1 - 1))
          local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
          lof[x8+1] = local.factor1
        } else if (nneigh1!=k) {
          print("warning numneigh!=k")
        }
      }
      
    }
  } else {
    #s_update_lof<-as.matrix(s_update_lof[-c(which(s_update_lof[,1]== dim(data2)[1])),])
    for (i in 1:length(s_update_lof_d)) {
      x8<-s_update_lof_d[i,1]-1
      nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
      if (nneigh1 == k) {
        j8 = seq(0, (nneigh1 - 1))
        local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
        lof[x8+1] = local.factor1
      } else if (nneigh1!=k) {
        print("warning numneigh!=k")
      }
    }
    
  }
  
  
  #lrddata<-lrddata[-c(1)]
  lof<-lof[-c(1)]
  
  
  #knn of the inserted observation
  data2=as.matrix(c(data2,x_1))
  neighdist = knneigh.vect(data2[dim(data2)[1], ], data2, k)
  distdata1= cbind(distdata1, neighdist)
  
  #reach_dist of x_1
  i=dim(data2)[1]
  j = seq(3, 3 + (distdata1[2, i] - distdata1[1, i]))
  numneigh = distdata1[2, i] - distdata1[1, i] + 1
  temp = rbind(diag(distdata1[distdata1[2, distdata1[j, 
                                                     i]], distdata1[j, i]]), distdata1[j + numneigh, 
                                                                                       i])
  #reach = 1/(sum(apply(temp, 2, max))/numneigh)
  
  reach_dist_x1 = apply(temp, 2, max) #reach_dist(x_1, knn)
  #1/(sum(reach_dist_x1)/numneigh)
  
  
  
  #update neighbors of x_1
  dist.obj <- dbscan::kNN(data2, k)
  n <- nrow(data2)
  
  rNN_x1 <- sapply(n, function(i) {
    as.vector(which(dist.obj$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                             1])
  })
  
  
  #update k_distance
  for (i in 1:length(rNN_x1)) {
    if (length(rNN_x1)==1){
      if (length(rNN_x1[[1]])!=0){
        x1<-rNN_x1[i]
        numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
        if (numneigh == k) {
          j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
          j2=seq(distdata1[1, x1],distdata1[2, x1])
          distdata1[j1,x1]=dist.obj$id[x1,]
          distdata1[j2,x1]=dist.obj$dist[x1,]
        } else if (numneigh!=k) {
          print("warning numneigh!=k")
        }
      }
    } else {
      x1<-rNN_x1[i,1]
      numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
      if (numneigh == k) {
        j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
        j2=seq(distdata1[1, x1],distdata1[2, x1])
        distdata1[j1,x1]=dist.obj$id[x1,]
        distdata1[j2,x1]=dist.obj$dist[x1,]
      } else if (numneigh!=k) {
        print("warning numneigh!=k")
      }
    }
  }
  
  #update reach_dist
  s_update=rNN_x1
  for (i in 1:length(rNN_x1)) {
    if (length(rNN_x1)==1) {
      if (length(rNN_x1[[1]])!=0) {
        x2<-rNN_x1[i]
        knn_x2<-as.matrix(dist.obj$id[x2,])
        knn_x2<-as.matrix(knn_x2[-c(which(knn_x2[,1]== dim(data2)[1])),])
        numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
        if (numneigh1 == k) {
          j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
          j4=seq(distdata1[1, x2],distdata1[2, x2])
          inds<-reachdist[j3,x2]
          inds1<-reachdist[j4,x2]
          inds2<-match(dist.obj$id[x2,],reachdist[j3,x2])
          reachdist[j3,x2]=dist.obj$id[x2,]
          reachdist[j4,x2]=inds1[inds2]
          x3<-length(which(is.na(reachdist[,x2])))
          x4<-which(is.na(reachdist[,x2]))
          if (x3!=0){
            temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                                x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                     x2])
            #reach = 1/(sum(apply(temp, 2, max))/numneigh)
            
            reach_dist_x2 = apply(temp1, 2, max)
            reachdist[x4,x2]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
          }
        } else if (numneigh1!=k) {
          print("warning numneigh!=k")
        }
        for (i in 1:dim(knn_x2)[1]) {
          x5<-knn_x2[i,1]
          x6<-which(reachdist[,x5]== x2)+numneigh1
          reachdist[x6,x5]=distdata1[length(distdata1[,x2]),x2]
          #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
          #x8<-which(reachdist[,x2]== dim(data2)[1])
          if (x2 %in% dist.obj$id[x5,]) {
            s_update=as.matrix(union(s_update,x5))
          }
          
        }
      }
      
    } else {
      x2<-rNN_x1[i,1]
      knn_x2<-as.matrix(dist.obj$id[x2,])
      knn_x2<-as.matrix(knn_x2[-c(which(knn_x2[,1]== dim(data2)[1])),])
      numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
      if (numneigh1 == k) {
        j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
        j4=seq(distdata1[1, x2],distdata1[2, x2])
        inds<-reachdist[j3,x2]
        inds1<-reachdist[j4,x2]
        inds2<-match(dist.obj$id[x2,],reachdist[j3,x2])
        reachdist[j3,x2]=dist.obj$id[x2,]
        reachdist[j4,x2]=inds1[inds2]
        x3<-length(which(is.na(reachdist[,x2])))
        x4<-which(is.na(reachdist[,x2]))
        if (x3!=0){
          temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                              x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                   x2])
          #reach = 1/(sum(apply(temp, 2, max))/numneigh)
          
          reach_dist_x2 = apply(temp1, 2, max)
          reachdist[x4,x2]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
        }
      } else if (numneigh1!=k) {
        print("warning numneigh!=k")
      }
      for (i in 1:dim(knn_x2)[1]) {
        x5<-knn_x2[i,1]
        x6<-which(reachdist[,x5]== x2)+numneigh1
        reachdist[x6,x5]=distdata1[length(distdata1[,x2]),x2]
        #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
        #x8<-which(reachdist[,x2]== dim(data2)[1])
        if (x2 %in% dist.obj$id[x5,]) {
          s_update=as.matrix(union(s_update,x5))
        }
        
      }
    }
  }
  
  if (numneigh == k) {
    ind<-dim(data2)[1]
    #j = seq(3, 3 + (distdata1[2, i] - distdata1[1, i]))
    #numneigh = distdata1[2, i] - distdata1[1, i] + 1
    #temp = rbind(diag(distdata1[distdata1[2, distdata1[j, 
    #                                                   i]], distdata1[j, i]]), distdata1[j + numneigh, 
    #                                                                                     i])
    #reach = 1/(sum(apply(temp, 2, max))/numneigh)
    reachdist=cbind(reachdist,neighdist)
    ind1<-seq(reachdist[1, ind],reachdist[2, ind])
    #reach_dist_x1 = apply(temp, 2, max) #reach_dist(x_1, knn)
    reachdist[ind1,ind]=reach_dist_x1
  } else if (numneigh1!=k) {
    print("warning numneigh!=k")
  }
  
  
  
  
  #update lrd
  s_update_lof<-s_update
  for (i in 1:length(s_update)) {
    if (length(s_update)==1){
      if (length(s_update[[1]])!=0){
        x7<-s_update[i]
        numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
        if (numneigh2 == k) {
          j5=seq(reachdist[1, x7],reachdist[2, x7])
          reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
          lrddata[x7] = reach2
          #dist.obj <- dbscan::kNN(data2, k)
          #n <- nrow(data2)
          rNN_x7 <- sapply(x7, function(i) {
            as.vector(which(dist.obj$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                     1])
          })
          if (length(rNN_x7)==1){
            if (length(rNN_x7[[1]])!=0){
              s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
            }
          } else {
            s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
          }
          
        } else if (numneigh2!=k) {
          print("warning numneigh!=k")
        }
      }
    } else {
      x7<-s_update[i,1]
      numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
      if (numneigh2 == k) {
        j5=seq(reachdist[1, x7],reachdist[2, x7])
        reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
        lrddata[x7] = reach2
        #dist.obj <- dbscan::kNN(data2, k)
        #n <- nrow(data2)
        rNN_x7 <- sapply(x7, function(i) {
          as.vector(which(dist.obj$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                   1])
        })
        if (length(rNN_x7)==1){
          if (length(rNN_x7[[1]])!=0){
            s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
          }
        } else {
          s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
        }
      } else if (numneigh2!=k) {
        print("warning numneigh!=k")
      }
    }
    
    
  }
  
  #calculate lrd for x_1
  x9<-dim(data2)[1]
  numneigh3 = reachdist[2, x9] - reachdist[1, x9] + 1
  if (numneigh3 == k) {
    j6=seq(reachdist[1, x9],reachdist[2, x9])
    reach3 = 1/(sum(reachdist[j6,x9])/numneigh3)
    lrddata[x9] = reach3
    #dist.obj <- dbscan::kNN(data2, k)
    #n <- nrow(data2)
  } else if (numneigh3!=k) {
    print("warning numneigh!=k")
  }
  
  
  
  #update lof
  if (length(s_update_lof)==1){
    if (length(s_update_lof[[1]])!=0){
      s_update_lof<-as.matrix(s_update_lof[-c(which(s_update_lof[]== dim(data2)[1])),])
      for (i in 1:length(s_update_lof)) {
        x8<-s_update_lof[i]
        nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
        if (nneigh1 == k) {
          j8 = seq(0, (nneigh1 - 1))
          local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
          lof[x8] = local.factor1
        } else if (nneigh1!=k) {
          print("warning numneigh!=k")
        }
      }
      
    }
  } else if(length(s_update_lof)>1){
    s_update_lof<-as.matrix(s_update_lof[-c(which(s_update_lof[,1]== dim(data2)[1])),])
    for (i in 1:length(s_update_lof)) {
      x8<-s_update_lof[i,1]
      nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
      if (nneigh1 == k) {
        j8 = seq(0, (nneigh1 - 1))
        local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
        lof[x8] = local.factor1
      } else if (nneigh1!=k) {
        print("warning numneigh!=k")
      }
    }
    
  }
  
  outlier_detection13[m1,4]<-length(union(s_update_lof, s_update_lof_d))
  #calculate lof for x_1
  x9<-dim(data2)[1]
  nneigh2 = distdata1[2, x9] - distdata1[1, x9] + 1
  if (nneigh2 == k) {
    j7 = seq(0, (nneigh2 - 1))
    local.factor2 = sum(lrddata[distdata1[3 + j7, x9]]/lrddata[x9])/nneigh2
    lof[x9] = local.factor2
    if (lof[x9]>2.9){
      outlier_detection13[m1,2]<-1
      outlier_detection13[m1,3]<-lof[x9]
      
    } else {
      outlier_detection13[m1,2]<-0
      outlier_detection13[m1,3]<-lof[x9]
    }
  } else if (nneigh2!=k) {
    print("warning numneigh!=k")
  }
  
  
}

mean(outlier_detection13[,4])
length(which(outlier_detection13[,2]==1))
which(outlier_detection13[,2]==1)
#outlier3<-union(which(outlier_detection2[,2]==1),which(lof[201:400] > 1.35) )
outlier3<-which(outlier_detection13[,2]==1)
which(c(1:300)%%15==0)
#length(which(lof[201:400] > 1.2))

TP_V<-length(which(outlier3 %in% which(c(1:300)%%15==0)))
FN_V<-20-TP_V
FP_V<-length(which(outlier_detection13[,2]==1))-TP_V
TN_V<-300-length(which(c(1:300)%%15==0))-FP_V
SE_V<-TP_V/(TP_V+FN_V) #true positive rate
SP_V<-TN_V/(FP_V+TN_V) #true negative rate
G_mean_lof_SWILOF_V<-sqrt(SE_V*SP_V)
H_mean_lof_SWILOF_V<-2*SE_V*SP_V/(SE_V+SP_V)

outlier_detection14<-as.data.frame(outlier_detection13)
outlier_detection14[,5]=(1:300)
outlier_detection14$V2<-as.factor(outlier_detection14$V2)
colnames(outlier_detection14)[2] <- "Outliers"

# plot
p13 <- ggplot(outlier_detection14, aes(x=V5, y=V1,label=V5)) +
  geom_point(aes(colour= Outliers)) + 
  #geom_point(aes(y=deaths), colour='red') +
  #geom_errorbar(aes(ymin=V4, ymax=V3)) +
  geom_text(aes(label=ifelse(Outliers==1,as.character(V5),'')),hjust=0,vjust=0)+
  theme_bw() +            
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))+
  #scale_y_continuous(expand=c(0,0), lim=c(0,1e3)) +
  labs(x='timepoint of streaming data', y='value of observation', title='SW-ILOF in non-stationary data stream \n k=10, alpha=2.9, W=50') 
#coord_cartesian(ylim=c(0,300)) +
#facet_grid(~age_group)
p13


#performance measures
my_rnorm_SW_ILOF_V <- function(data=data1, k=10,w=151, a1=2.9 ){
  data2<-as.matrix(data[w:200])
  distdata = dist.to.knn(data2, k)
  distdata1<-distdata
  reachdist=reachability_dist(distdata, k)
  lrddata = reachability(distdata, k)
  p = dim(distdata)[2]
  lof = rep(0, p)
  for (i in 1:p) {
    nneigh = distdata[2, i] - distdata[1, i] + 1
    j = seq(0, (nneigh - 1))
    local.factor = sum(lrddata[distdata[3 + j, i]]/lrddata[i])/nneigh
    lof[i] = local.factor
  }
  
  outlier_detection13 = matrix(nrow=300, ncol=3)
  for (i in 1:300){
    m=i%%15
    m1=i
    
    if (m == 0){
      if (m1%%2==0 & m1>=100){
        x1<- rnorm(1, mean = 0.02*i, sd = 1)
        x_1 <- x1- 5*1
      }else{
        x1<- rnorm(1, mean = 0.02*i, sd = 1)
        x_1 <- (abs(x1) + 5*1) * sign(x1)
      }
    } else {
      x_1=rnorm(1,mean=0.02*i, sd = 1)
    }
    outlier_detection13[i,1]<-x_1
    
    #Deletion
    #reverse neighbours of deletion point
    dist.obj1 <- dbscan::kNN(data2, k)
    #n <- nrow(data2)
    
    rNN_x_1 <- sapply(1, function(i) {
      as.vector(which(dist.obj1$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                1])
    })
    
    data2<-as.matrix(data2[-c(1),])
    
    #Update k-distance
    
    j1=seq(3, 3 + (distdata1[2, 1] - distdata1[1, 1]))
    j2=seq(distdata1[1, 1],distdata1[2, 1])
    distdata1[j1,]=distdata1[j1,]-1
    dist.obj2 <- dbscan::kNN(data2, k)
    
    for (i in 1:length(rNN_x_1)) {
      if (length(rNN_x_1)==1){
        if (length(rNN_x_1[[1]])!=0){
          x1<-rNN_x_1[i]
          numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
          if (numneigh == k) {
            j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
            j2=seq(distdata1[1, x1],distdata1[2, x1])
            distdata1[j1,x1]=dist.obj2$id[x1-1,]
            distdata1[j2,x1]=dist.obj2$dist[x1-1,]
          } else if (numneigh!=k) {
            print("warning numneigh!=k")
          }
        }
      } else {
        x1<-rNN_x_1[i,1]
        numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
        if (numneigh == k) {
          j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
          j2=seq(distdata1[1, x1],distdata1[2, x1])
          distdata1[j1,x1]=dist.obj2$id[x1-1,]
          distdata1[j2,x1]=dist.obj2$dist[x1-1,]
        } else if (numneigh!=k) {
          print("warning numneigh!=k")
        }
      }
    }
    
    
    distdata1<-distdata1[,-c(1)]
    s_update_lrd=rNN_x_1
    j32=seq(3, 3 + (distdata1[2, 1] - distdata1[1, 1]))
    j4=seq(distdata1[1, 1],distdata1[2, 1])
    reachdist[j32,]=reachdist[j32,]-1
    #update reach_dist
    for (i in 1:length(rNN_x_1)) {
      if (length(rNN_x_1)==1) {
        if (length(rNN_x_1[[1]])!=0) {
          x2<-rNN_x_1[i]-1
          knn_x2<-as.matrix(dist.obj2$id[x2,-c(length(dist.obj2$id[x2,]))])
          numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
          if (numneigh1 == k) {
            j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
            j4=seq(distdata1[1, x2],distdata1[2, x2])
            inds<-reachdist[j3,x2+1]
            inds1<-reachdist[j4,x2+1]
            inds2<-match(dist.obj2$id[x2,],reachdist[j3,x2+1])
            reachdist[j3,x2+1]=dist.obj2$id[x2,]
            reachdist[j4,x2+1]=inds1[inds2]
            x3<-length(which(is.na(reachdist[,x2+1])))
            x4<-which(is.na(reachdist[,x2+1]))
            if (x3!=0){
              temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                                  x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                       x2])
              #reach = 1/(sum(apply(temp, 2, max))/numneigh)
              
              reach_dist_x2 = apply(temp1, 2, max)
              reachdist[x4,x2+1]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
            }
          } else if (numneigh1!=k) {
            print("warning numneigh!=k")
          }
          for (i in 1:dim(knn_x2)[1]) {
            x5<-knn_x2[i,1]
            x6<-which(reachdist[,x5+1]== x2)+numneigh1
            reachdist[x6,x5+1]=distdata1[length(distdata1[,x2]),x2]
            #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
            #x8<-which(reachdist[,x2]== dim(data2)[1])
            if (x2 %in% dist.obj2$id[x5,]) {
              s_update_lrd=as.matrix(union(s_update_lrd,(x5+1)))
            }
            
          }
        }
        
      } else {
        x2<-rNN_x_1[i,1]-1
        knn_x2<-as.matrix(dist.obj2$id[x2,-c(length(dist.obj2$id[x2,]))])
        numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
        if (numneigh1 == k) {
          j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
          j4=seq(distdata1[1, x2],distdata1[2, x2])
          inds<-reachdist[j3,x2+1]
          inds1<-reachdist[j4,x2+1]
          inds2<-match(dist.obj2$id[x2,],reachdist[j3,x2+1])
          reachdist[j3,x2+1]=dist.obj2$id[x2,]
          reachdist[j4,x2+1]=inds1[inds2]
          x3<-length(which(is.na(reachdist[,x2+1])))
          x4<-which(is.na(reachdist[,x2+1]))
          if (x3!=0){
            temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                                x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                     x2])
            #reach = 1/(sum(apply(temp, 2, max))/numneigh)
            
            reach_dist_x2 = apply(temp1, 2, max)
            reachdist[x4,x2+1]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
          }
        } else if (numneigh1!=k) {
          print("warning numneigh!=k")
        }
        for (i in 1:dim(knn_x2)[1]) {
          x5<-knn_x2[i,1]
          x6<-which(reachdist[,x5+1]== x2)+numneigh1
          reachdist[x6,x5+1]=distdata1[length(distdata1[,x2]),x2]
          #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
          #x8<-which(reachdist[,x2]== dim(data2)[1])
          if (x2 %in% dist.obj2$id[x5,]) {
            s_update_lrd=as.matrix(union(s_update_lrd,(x5+1)))
          }
          
        }
      }
    }
    
    
    
    #distdata1<-distdata1[,-c(1)]
    reachdist<-reachdist[,-c(1)]
    
    #update lrd
    s_update_lof_d<-s_update_lrd
    for (i in 1:length(s_update_lrd)) {
      if (length(s_update_lrd)==1){
        if (length(s_update_lrd[[1]])!=0){
          x7<-s_update_lrd[i]-1
          numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
          if (numneigh2 == k) {
            j5=seq(reachdist[1, x7],reachdist[2, x7])
            reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
            lrddata[x7+1] = reach2
            #dist.obj <- dbscan::kNN(data2, k)
            #n <- nrow(data2)
            rNN_x7 <- sapply(x7, function(i) {
              as.vector(which(dist.obj2$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                        1])
            })
            if (length(rNN_x7)==1){
              if (length(rNN_x7[[1]])!=0){
                s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
              }
            } else {
              s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
            }
            
          } else if (numneigh2!=k) {
            print("warning numneigh!=k")
          }
        }
      } else {
        x7<-s_update_lrd[i,1]-1
        numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
        if (numneigh2 == k) {
          j5=seq(reachdist[1, x7],reachdist[2, x7])
          reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
          lrddata[x7+1] = reach2
          #dist.obj <- dbscan::kNN(data2, k)
          #n <- nrow(data2)
          rNN_x7 <- sapply(x7, function(i) {
            as.vector(which(dist.obj2$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                      1])
          })
          if (length(rNN_x7)==1){
            if (length(rNN_x7[[1]])!=0){
              s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
            }
          } else {
            s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
          }
        } else if (numneigh2!=k) {
          print("warning numneigh!=k")
        }
      }
      
      
    }
    
    lrddata<-lrddata[-c(1)]
    
    
    #update lof
    if (length(s_update_lof_d)==1){
      if (length(s_update_lof_d[[1]])!=0){
        #s_update_lof_d<-as.matrix(s_update_lof[-c(which(s_update_lof[]== dim(data2)[1])),])
        for (i in 1:length(s_update_lof_d)) {
          x8<-s_update_lof_d[i]-1
          nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
          if (nneigh1 == k) {
            j8 = seq(0, (nneigh1 - 1))
            local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
            lof[x8+1] = local.factor1
          } else if (nneigh1!=k) {
            print("warning numneigh!=k")
          }
        }
        
      }
    } else {
      #s_update_lof<-as.matrix(s_update_lof[-c(which(s_update_lof[,1]== dim(data2)[1])),])
      for (i in 1:length(s_update_lof_d)) {
        x8<-s_update_lof_d[i,1]-1
        nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
        if (nneigh1 == k) {
          j8 = seq(0, (nneigh1 - 1))
          local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
          lof[x8+1] = local.factor1
        } else if (nneigh1!=k) {
          print("warning numneigh!=k")
        }
      }
      
    }
    
    
    #lrddata<-lrddata[-c(1)]
    lof<-lof[-c(1)]
    
    
    #knn of the inserted observation
    data2=as.matrix(c(data2,x_1))
    neighdist = knneigh.vect(data2[dim(data2)[1], ], data2, k)
    distdata1= cbind(distdata1, neighdist)
    
    #reach_dist of x_1
    i=dim(data2)[1]
    j = seq(3, 3 + (distdata1[2, i] - distdata1[1, i]))
    numneigh = distdata1[2, i] - distdata1[1, i] + 1
    temp = rbind(diag(distdata1[distdata1[2, distdata1[j, 
                                                       i]], distdata1[j, i]]), distdata1[j + numneigh, 
                                                                                         i])
    #reach = 1/(sum(apply(temp, 2, max))/numneigh)
    
    reach_dist_x1 = apply(temp, 2, max) #reach_dist(x_1, knn)
    #1/(sum(reach_dist_x1)/numneigh)
    
    
    
    #update neighbors of x_1
    dist.obj <- dbscan::kNN(data2, k)
    n <- nrow(data2)
    
    rNN_x1 <- sapply(n, function(i) {
      as.vector(which(dist.obj$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                               1])
    })
    
    
    #update k_distance
    for (i in 1:length(rNN_x1)) {
      if (length(rNN_x1)==1){
        if (length(rNN_x1[[1]])!=0){
          x1<-rNN_x1[i]
          numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
          if (numneigh == k) {
            j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
            j2=seq(distdata1[1, x1],distdata1[2, x1])
            distdata1[j1,x1]=dist.obj$id[x1,]
            distdata1[j2,x1]=dist.obj$dist[x1,]
          } else if (numneigh!=k) {
            print("warning numneigh!=k")
          }
        }
      } else {
        x1<-rNN_x1[i,1]
        numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
        if (numneigh == k) {
          j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
          j2=seq(distdata1[1, x1],distdata1[2, x1])
          distdata1[j1,x1]=dist.obj$id[x1,]
          distdata1[j2,x1]=dist.obj$dist[x1,]
        } else if (numneigh!=k) {
          print("warning numneigh!=k")
        }
      }
    }
    
    #update reach_dist
    s_update=rNN_x1
    for (i in 1:length(rNN_x1)) {
      if (length(rNN_x1)==1) {
        if (length(rNN_x1[[1]])!=0) {
          x2<-rNN_x1[i]
          knn_x2<-as.matrix(dist.obj$id[x2,])
          knn_x2<-as.matrix(knn_x2[-c(which(knn_x2[,1]== dim(data2)[1])),])
          numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
          if (numneigh1 == k) {
            j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
            j4=seq(distdata1[1, x2],distdata1[2, x2])
            inds<-reachdist[j3,x2]
            inds1<-reachdist[j4,x2]
            inds2<-match(dist.obj$id[x2,],reachdist[j3,x2])
            reachdist[j3,x2]=dist.obj$id[x2,]
            reachdist[j4,x2]=inds1[inds2]
            x3<-length(which(is.na(reachdist[,x2])))
            x4<-which(is.na(reachdist[,x2]))
            if (x3!=0){
              temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                                  x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                       x2])
              #reach = 1/(sum(apply(temp, 2, max))/numneigh)
              
              reach_dist_x2 = apply(temp1, 2, max)
              reachdist[x4,x2]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
            }
          } else if (numneigh1!=k) {
            print("warning numneigh!=k")
          }
          for (i in 1:dim(knn_x2)[1]) {
            x5<-knn_x2[i,1]
            x6<-which(reachdist[,x5]== x2)+numneigh1
            reachdist[x6,x5]=distdata1[length(distdata1[,x2]),x2]
            #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
            #x8<-which(reachdist[,x2]== dim(data2)[1])
            if (x2 %in% dist.obj$id[x5,]) {
              s_update=as.matrix(union(s_update,x5))
            }
            
          }
        }
        
      } else {
        x2<-rNN_x1[i,1]
        knn_x2<-as.matrix(dist.obj$id[x2,])
        knn_x2<-as.matrix(knn_x2[-c(which(knn_x2[,1]== dim(data2)[1])),])
        numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
        if (numneigh1 == k) {
          j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
          j4=seq(distdata1[1, x2],distdata1[2, x2])
          inds<-reachdist[j3,x2]
          inds1<-reachdist[j4,x2]
          inds2<-match(dist.obj$id[x2,],reachdist[j3,x2])
          reachdist[j3,x2]=dist.obj$id[x2,]
          reachdist[j4,x2]=inds1[inds2]
          x3<-length(which(is.na(reachdist[,x2])))
          x4<-which(is.na(reachdist[,x2]))
          if (x3!=0){
            temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                                x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                     x2])
            #reach = 1/(sum(apply(temp, 2, max))/numneigh)
            
            reach_dist_x2 = apply(temp1, 2, max)
            reachdist[x4,x2]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
          }
        } else if (numneigh1!=k) {
          print("warning numneigh!=k")
        }
        for (i in 1:dim(knn_x2)[1]) {
          x5<-knn_x2[i,1]
          x6<-which(reachdist[,x5]== x2)+numneigh1
          reachdist[x6,x5]=distdata1[length(distdata1[,x2]),x2]
          #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
          #x8<-which(reachdist[,x2]== dim(data2)[1])
          if (x2 %in% dist.obj$id[x5,]) {
            s_update=as.matrix(union(s_update,x5))
          }
          
        }
      }
    }
    
    if (numneigh == k) {
      ind<-dim(data2)[1]
      #j = seq(3, 3 + (distdata1[2, i] - distdata1[1, i]))
      #numneigh = distdata1[2, i] - distdata1[1, i] + 1
      #temp = rbind(diag(distdata1[distdata1[2, distdata1[j, 
      #                                                   i]], distdata1[j, i]]), distdata1[j + numneigh, 
      #                                                                                     i])
      #reach = 1/(sum(apply(temp, 2, max))/numneigh)
      reachdist=cbind(reachdist,neighdist)
      ind1<-seq(reachdist[1, ind],reachdist[2, ind])
      #reach_dist_x1 = apply(temp, 2, max) #reach_dist(x_1, knn)
      reachdist[ind1,ind]=reach_dist_x1
    } else if (numneigh1!=k) {
      print("warning numneigh!=k")
    }
    
    
    
    
    #update lrd
    s_update_lof<-s_update
    for (i in 1:length(s_update)) {
      if (length(s_update)==1){
        if (length(s_update[[1]])!=0){
          x7<-s_update[i]
          numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
          if (numneigh2 == k) {
            j5=seq(reachdist[1, x7],reachdist[2, x7])
            reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
            lrddata[x7] = reach2
            #dist.obj <- dbscan::kNN(data2, k)
            #n <- nrow(data2)
            rNN_x7 <- sapply(x7, function(i) {
              as.vector(which(dist.obj$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                       1])
            })
            if (length(rNN_x7)==1){
              if (length(rNN_x7[[1]])!=0){
                s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
              }
            } else {
              s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
            }
            
          } else if (numneigh2!=k) {
            print("warning numneigh!=k")
          }
        }
      } else {
        x7<-s_update[i,1]
        numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
        if (numneigh2 == k) {
          j5=seq(reachdist[1, x7],reachdist[2, x7])
          reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
          lrddata[x7] = reach2
          #dist.obj <- dbscan::kNN(data2, k)
          #n <- nrow(data2)
          rNN_x7 <- sapply(x7, function(i) {
            as.vector(which(dist.obj$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                     1])
          })
          if (length(rNN_x7)==1){
            if (length(rNN_x7[[1]])!=0){
              s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
            }
          } else {
            s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
          }
        } else if (numneigh2!=k) {
          print("warning numneigh!=k")
        }
      }
      
      
    }
    
    #calculate lrd for x_1
    x9<-dim(data2)[1]
    numneigh3 = reachdist[2, x9] - reachdist[1, x9] + 1
    if (numneigh3 == k) {
      j6=seq(reachdist[1, x9],reachdist[2, x9])
      reach3 = 1/(sum(reachdist[j6,x9])/numneigh3)
      lrddata[x9] = reach3
      #dist.obj <- dbscan::kNN(data2, k)
      #n <- nrow(data2)
    } else if (numneigh3!=k) {
      print("warning numneigh!=k")
    }
    
    
    
    #update lof
    if (length(s_update_lof)==1){
      if (length(s_update_lof[[1]])!=0){
        if (length(which(s_update_lof[]== dim(data2)[1]))!=0){
          s_update_lof<-as.matrix(s_update_lof[-c(which(s_update_lof[]== dim(data2)[1])),])
        }
        
        for (i in 1:length(s_update_lof)) {
          x8<-s_update_lof[i]
          nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
          if (nneigh1 == k) {
            j8 = seq(0, (nneigh1 - 1))
            local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
            lof[x8] = local.factor1
          } else if (nneigh1!=k) {
            print("warning numneigh!=k")
          }
        }
        
      }
    } else if(length(s_update_lof)>1){
      if (length(which(s_update_lof[,1]== dim(data2)[1]))!=0){
        s_update_lof<-as.matrix(s_update_lof[-c(which(s_update_lof[,1]== dim(data2)[1])),])
      }
      
      for (i in 1:length(s_update_lof)) {
        x8<-s_update_lof[i,1]
        nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
        if (nneigh1 == k) {
          j8 = seq(0, (nneigh1 - 1))
          local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
          lof[x8] = local.factor1
        } else if (nneigh1!=k) {
          print("warning numneigh!=k")
        }
      }
      
    }
    
    
    #calculate lof for x_1
    x9<-dim(data2)[1]
    nneigh2 = distdata1[2, x9] - distdata1[1, x9] + 1
    if (nneigh2 == k) {
      j7 = seq(0, (nneigh2 - 1))
      local.factor2 = sum(lrddata[distdata1[3 + j7, x9]]/lrddata[x9])/nneigh2
      lof[x9] = local.factor2
      if (lof[x9]>a1){
        outlier_detection13[m1,2]<-1
        outlier_detection13[m1,3]<-lof[x9]
        
      } else {
        outlier_detection13[m1,2]<-0
        outlier_detection13[m1,3]<-lof[x9]
      }
    } else if (nneigh2!=k) {
      print("warning numneigh!=k")
    }
    
    
  }
  
  
  length(which(outlier_detection13[,2]==1))
  which(outlier_detection13[,2]==1)
  #outlier3<-union(which(outlier_detection2[,2]==1),which(lof[201:400] > 1.35) )
  outlier3<-which(outlier_detection13[,2]==1)
  which(c(1:300)%%15==0)
  #length(which(lof[201:400] > 1.2))
  
  TP_V<-length(which(outlier3 %in% which(c(1:300)%%15==0)))
  FN_V<-20-TP_V
  FP_V<-length(which(outlier_detection13[,2]==1))-TP_V
  TN_V<-300-length(which(c(1:300)%%15==0))-FP_V
  SE_V<-TP_V/(TP_V+FN_V) #true positive rate
  SP_V<-TN_V/(FP_V+TN_V) #true negative rate
  G_mean_lof_SWILOF_V<-sqrt(SE_V*SP_V)
  H_mean_lof_SWILOF_V<-2*SE_V*SP_V/(SE_V+SP_V)
  
  
  
  
  c(G_mean_lof_SWILOF_V, H_mean_lof_SWILOF_V)
}


# grid over which we will perform the hyperparameter search:
hparam_grid <- as.data.frame(expand.grid(k=seq(5, 15, by=2), a1=seq(1.5, 4, by=0.4),w=seq(141, 181, by=10)))

# to store the OOB estimates of the MSE
oob_mses <- rep(0.0, nrow(hparam_grid))

# perform the gridsearch
for(hparam_idx in 1:nrow(hparam_grid)) {
  # train candidate model
  this_k <- hparam_grid[hparam_idx, 1]
  this_a1 <- hparam_grid[hparam_idx, 2]
  this_w <- hparam_grid[hparam_idx, 3]
  test2<-my_rnorm_SW_ILOF_V(data=data1, k=this_k,w=this_w, a1=this_a1)[2]
  #rf <- randomForest(x_train, y_train, mtry=this_mtry, maxnodes=this_maxnodes)
  
  # calculate H-mean
  oob_mses[hparam_idx] <- test2
}

# select the best model (that which has the minimum OOB MSE)
best_hparam_set1 <- hparam_grid[which.max(oob_mses),]



#my_rnorm_streaming_boxplot_S_performance <- function(data1=data1, q_lower=q_lower, n_p_lower=n_p_lower, n_d_lower=n_d_lower, dn_d_lower=dn_d_lower,q_upper=q_upper, n_p_upper=n_p_upper, n_d_upper=n_d_upper, dn_d_upper=dn_d_upper){
outlier_detection10 = matrix(nrow=100, ncol=2)
for (i in 1:100) {
  test1<-my_rnorm_SW_ILOF_V(data=data1, k=15,w=171, a1=2.3)
  outlier_detection10[i,1]<-test1[1]
  outlier_detection10[i,2]<-test1[2]
}

mean(outlier_detection10[,1])
mean(outlier_detection10[,2])






#Normal distribution set in purpose
k<-15
#data=data1
data=data1[151:200]
data=as.matrix(data)
data2<-data
distdata = dist.to.knn(data, k)
distdata1<-distdata
reachdist=reachability_dist(distdata, k)
lrddata = reachability(distdata, k)
p = dim(distdata)[2]
lof = rep(0, p)
for (i in 1:p) {
  nneigh = distdata[2, i] - distdata[1, i] + 1
  j = seq(0, (nneigh - 1))
  local.factor = sum(lrddata[distdata[3 + j, i]]/lrddata[i])/nneigh
  lof[i] = local.factor
}
lof

m=c(102, 350, 450, 750, 1050, 1450)
outlier_detection15 = matrix(nrow=1500, ncol=4)
for (i in 1:1500){
  #m=i%%15
  m1=i
  #m=i%%11
  i_1=i
  if (i<= 200 ){
    if (i_1 %in% m ){
      x1<- rnorm(1, mean = 0, sd = 1)
      x_1 <- (abs(x1) + 5*1) * sign(x1)
      
    } else {
      x_1=rnorm(1, mean = 0, sd = 1)
    }
  } else if (i>=201 & i<= 250){
    if (i_1 %in% m) {
      x1<- rnorm(1, mean = 0.1*(i-200), sd = 1)
      x_1 <- (abs(x1) + 4*1) * sign(x1)
      
    } else {
      x_1=rnorm(1, mean = 0.1*(i-200), sd = 1)
      
    }
  } else if (i>=251 & i<= 550){
    if (i_1 == m[2]) {
      x1<- rnorm(1, mean = 5, sd = 1)
      x_1 <- (abs(x1) + 5*1) * sign(x1)
      
    } else if (i_1 == m[3]) {
      x1<- rnorm(1, mean = 5, sd = 1)
      x_1 <- (abs(x1) - 5*1) * sign(x1)
    } else {
      x_1=rnorm(1, mean = 5, sd = 1)
      
    }
  } else if (i>=551 & i<= 650){
    if (i_1 %in% m) {
      x1<- rnorm(1, mean = 5+0.05*(i-550), sd = 1)
      x_1 <- (abs(x1) + 5*1) * sign(x1)
      
    } else {
      x_1=rnorm(1, mean = 5+0.05*(i-550), sd = 1)
      
    }
  } else if (i>=651 & i<= 850){
    if (i_1 %in% m) {
      x1<- rnorm(1, mean = 10, sd = 1)
      x_1 <- (abs(x1) - 5*1) * sign(x1)
      
    } else {
      x_1=rnorm(1, mean = 10, sd = 1)
      
    }
  }  else if (i>=851 & i<= 900){
    if (i_1 %in% m) {
      x1<- rnorm(1, mean = 10-0.2*(i-850), sd = 1)
      x_1 <- (abs(x1) + 5*1) * sign(x1)
      
    } else {
      x_1=rnorm(1, mean = 10-0.2*(i-850), sd = 1)
      
    }
  } else if (i>=901 & i<= 1200){
    if (i_1 %in% m) {
      x1<- rnorm(1, mean = -0.015*(i-900), sd = 1)
      x_1 <- (x1 + 5*1) 
      
    } else {
      x_1=rnorm(1, mean = -0.015*(i-900), sd = 1)
      
    }
  } else if (i>=1201 & i<= 1300){
    if (i_1 %in% m) {
      x1<- rnorm(1, mean = -4.5+0.045*(i-1200), sd = 1)
      x_1 <- (abs(x1) + 5*1) * sign(x1)
      
    } else {
      x_1=rnorm(1, mean = -4.5+0.045*(i-1200), sd = 1)
      
    }
  } else if (i>=1301 & i<= 1500){
    if (i_1 %in% m) {
      x1<- rnorm(1, mean = 0, sd = 1)
      x_1 <- (abs(x1) + 5*1) * sign(x1)
      
    } else {
      x_1=rnorm(1, mean = 0, sd = 1)
      
    }
  }    
  
  outlier_detection15[i,1]<-x_1
  
  #Deletion
  #reverse neighbours of deletion point
  dist.obj1 <- dbscan::kNN(data2, k)
  #n <- nrow(data2)
  
  rNN_x_1 <- sapply(1, function(i) {
    as.vector(which(dist.obj1$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                              1])
  })
  
  data2<-as.matrix(data2[-c(1),])
  
  #Update k-distance
  
  j1=seq(3, 3 + (distdata1[2, 1] - distdata1[1, 1]))
  j2=seq(distdata1[1, 1],distdata1[2, 1])
  distdata1[j1,]=distdata1[j1,]-1
  dist.obj2 <- dbscan::kNN(data2, k)
  
  for (i in 1:length(rNN_x_1)) {
    if (length(rNN_x_1)==1){
      if (length(rNN_x_1[[1]])!=0){
        x1<-rNN_x_1[i]
        numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
        if (numneigh == k) {
          j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
          j2=seq(distdata1[1, x1],distdata1[2, x1])
          distdata1[j1,x1]=dist.obj2$id[x1-1,]
          distdata1[j2,x1]=dist.obj2$dist[x1-1,]
        } else if (numneigh!=k) {
          print("warning numneigh!=k")
        }
      }
    } else {
      x1<-rNN_x_1[i,1]
      numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
      if (numneigh == k) {
        j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
        j2=seq(distdata1[1, x1],distdata1[2, x1])
        distdata1[j1,x1]=dist.obj2$id[x1-1,]
        distdata1[j2,x1]=dist.obj2$dist[x1-1,]
      } else if (numneigh!=k) {
        print("warning numneigh!=k")
      }
    }
  }
  
  
  distdata1<-distdata1[,-c(1)]
  s_update_lrd=rNN_x_1
  j32=seq(3, 3 + (distdata1[2, 1] - distdata1[1, 1]))
  j4=seq(distdata1[1, 1],distdata1[2, 1])
  reachdist[j32,]=reachdist[j32,]-1
  #update reach_dist
  for (i in 1:length(rNN_x_1)) {
    if (length(rNN_x_1)==1) {
      if (length(rNN_x_1[[1]])!=0) {
        x2<-rNN_x_1[i]-1
        knn_x2<-as.matrix(dist.obj2$id[x2,-c(length(dist.obj2$id[x2,]))])
        numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
        if (numneigh1 == k) {
          j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
          j4=seq(distdata1[1, x2],distdata1[2, x2])
          inds<-reachdist[j3,x2+1]
          inds1<-reachdist[j4,x2+1]
          inds2<-match(dist.obj2$id[x2,],reachdist[j3,x2+1])
          reachdist[j3,x2+1]=dist.obj2$id[x2,]
          reachdist[j4,x2+1]=inds1[inds2]
          x3<-length(which(is.na(reachdist[,x2+1])))
          x4<-which(is.na(reachdist[,x2+1]))
          if (x3!=0){
            temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                                x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                     x2])
            #reach = 1/(sum(apply(temp, 2, max))/numneigh)
            
            reach_dist_x2 = apply(temp1, 2, max)
            reachdist[x4,x2+1]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
          }
        } else if (numneigh1!=k) {
          print("warning numneigh!=k")
        }
        for (i in 1:dim(knn_x2)[1]) {
          x5<-knn_x2[i,1]
          x6<-which(reachdist[,x5+1]== x2)+numneigh1
          reachdist[x6,x5+1]=distdata1[length(distdata1[,x2]),x2]
          #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
          #x8<-which(reachdist[,x2]== dim(data2)[1])
          if (x2 %in% dist.obj2$id[x5,]) {
            s_update_lrd=as.matrix(union(s_update_lrd,(x5+1)))
          }
          
        }
      }
      
    } else {
      x2<-rNN_x_1[i,1]-1
      knn_x2<-as.matrix(dist.obj2$id[x2,-c(length(dist.obj2$id[x2,]))])
      numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
      if (numneigh1 == k) {
        j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
        j4=seq(distdata1[1, x2],distdata1[2, x2])
        inds<-reachdist[j3,x2+1]
        inds1<-reachdist[j4,x2+1]
        inds2<-match(dist.obj2$id[x2,],reachdist[j3,x2+1])
        reachdist[j3,x2+1]=dist.obj2$id[x2,]
        reachdist[j4,x2+1]=inds1[inds2]
        x3<-length(which(is.na(reachdist[,x2+1])))
        x4<-which(is.na(reachdist[,x2+1]))
        if (x3!=0){
          temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                              x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                   x2])
          #reach = 1/(sum(apply(temp, 2, max))/numneigh)
          
          reach_dist_x2 = apply(temp1, 2, max)
          reachdist[x4,x2+1]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
        }
      } else if (numneigh1!=k) {
        print("warning numneigh!=k")
      }
      for (i in 1:dim(knn_x2)[1]) {
        x5<-knn_x2[i,1]
        x6<-which(reachdist[,x5+1]== x2)+numneigh1
        reachdist[x6,x5+1]=distdata1[length(distdata1[,x2]),x2]
        #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
        #x8<-which(reachdist[,x2]== dim(data2)[1])
        if (x2 %in% dist.obj2$id[x5,]) {
          s_update_lrd=as.matrix(union(s_update_lrd,(x5+1)))
        }
        
      }
    }
  }
  
  
  
  #distdata1<-distdata1[,-c(1)]
  reachdist<-reachdist[,-c(1)]
  
  #update lrd
  s_update_lof_d<-s_update_lrd
  for (i in 1:length(s_update_lrd)) {
    if (length(s_update_lrd)==1){
      if (length(s_update_lrd[[1]])!=0){
        x7<-s_update_lrd[i]-1
        numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
        if (numneigh2 == k) {
          j5=seq(reachdist[1, x7],reachdist[2, x7])
          reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
          lrddata[x7+1] = reach2
          #dist.obj <- dbscan::kNN(data2, k)
          #n <- nrow(data2)
          rNN_x7 <- sapply(x7, function(i) {
            as.vector(which(dist.obj2$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                      1])
          })
          if (length(rNN_x7)==1){
            if (length(rNN_x7[[1]])!=0){
              s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
            }
          } else {
            s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
          }
          
        } else if (numneigh2!=k) {
          print("warning numneigh!=k")
        }
      }
    } else {
      x7<-s_update_lrd[i,1]-1
      numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
      if (numneigh2 == k) {
        j5=seq(reachdist[1, x7],reachdist[2, x7])
        reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
        lrddata[x7+1] = reach2
        #dist.obj <- dbscan::kNN(data2, k)
        #n <- nrow(data2)
        rNN_x7 <- sapply(x7, function(i) {
          as.vector(which(dist.obj2$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                    1])
        })
        if (length(rNN_x7)==1){
          if (length(rNN_x7[[1]])!=0){
            s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
          }
        } else {
          s_update_lof_d=as.matrix(union(s_update_lof_d,(rNN_x7+1)))
        }
      } else if (numneigh2!=k) {
        print("warning numneigh!=k")
      }
    }
    
    
  }
  
  lrddata<-lrddata[-c(1)]
  
  
  #update lof
  if (length(s_update_lof_d)==1){
    if (length(s_update_lof_d[[1]])!=0){
      #s_update_lof_d<-as.matrix(s_update_lof[-c(which(s_update_lof[]== dim(data2)[1])),])
      for (i in 1:length(s_update_lof_d)) {
        x8<-s_update_lof_d[i]-1
        nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
        if (nneigh1 == k) {
          j8 = seq(0, (nneigh1 - 1))
          local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
          lof[x8+1] = local.factor1
        } else if (nneigh1!=k) {
          print("warning numneigh!=k")
        }
      }
      
    }
  } else {
    #s_update_lof<-as.matrix(s_update_lof[-c(which(s_update_lof[,1]== dim(data2)[1])),])
    for (i in 1:length(s_update_lof_d)) {
      x8<-s_update_lof_d[i,1]-1
      nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
      if (nneigh1 == k) {
        j8 = seq(0, (nneigh1 - 1))
        local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
        lof[x8+1] = local.factor1
      } else if (nneigh1!=k) {
        print("warning numneigh!=k")
      }
    }
    
  }
  
  
  #lrddata<-lrddata[-c(1)]
  lof<-lof[-c(1)]
  
  
  #knn of the inserted observation
  data2=as.matrix(c(data2,x_1))
  neighdist = knneigh.vect(data2[dim(data2)[1], ], data2, k)
  distdata1= cbind(distdata1, neighdist)
  
  #reach_dist of x_1
  i=dim(data2)[1]
  j = seq(3, 3 + (distdata1[2, i] - distdata1[1, i]))
  numneigh = distdata1[2, i] - distdata1[1, i] + 1
  temp = rbind(diag(distdata1[distdata1[2, distdata1[j, 
                                                     i]], distdata1[j, i]]), distdata1[j + numneigh, 
                                                                                       i])
  #reach = 1/(sum(apply(temp, 2, max))/numneigh)
  
  reach_dist_x1 = apply(temp, 2, max) #reach_dist(x_1, knn)
  #1/(sum(reach_dist_x1)/numneigh)
  
  
  
  #update neighbors of x_1
  dist.obj <- dbscan::kNN(data2, k)
  n <- nrow(data2)
  
  rNN_x1 <- sapply(n, function(i) {
    as.vector(which(dist.obj$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                             1])
  })
  
  
  #update k_distance
  for (i in 1:length(rNN_x1)) {
    if (length(rNN_x1)==1){
      if (length(rNN_x1[[1]])!=0){
        x1<-rNN_x1[i]
        numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
        if (numneigh == k) {
          j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
          j2=seq(distdata1[1, x1],distdata1[2, x1])
          distdata1[j1,x1]=dist.obj$id[x1,]
          distdata1[j2,x1]=dist.obj$dist[x1,]
        } else if (numneigh!=k) {
          print("warning numneigh!=k")
        }
      }
    } else {
      x1<-rNN_x1[i,1]
      numneigh = distdata1[2, x1] - distdata1[1, x1] + 1
      if (numneigh == k) {
        j1=seq(3, 3 + (distdata1[2, x1] - distdata1[1, x1]))
        j2=seq(distdata1[1, x1],distdata1[2, x1])
        distdata1[j1,x1]=dist.obj$id[x1,]
        distdata1[j2,x1]=dist.obj$dist[x1,]
      } else if (numneigh!=k) {
        print("warning numneigh!=k")
      }
    }
  }
  
  #update reach_dist
  s_update=rNN_x1
  for (i in 1:length(rNN_x1)) {
    if (length(rNN_x1)==1) {
      if (length(rNN_x1[[1]])!=0) {
        x2<-rNN_x1[i]
        knn_x2<-as.matrix(dist.obj$id[x2,])
        knn_x2<-as.matrix(knn_x2[-c(which(knn_x2[,1]== dim(data2)[1])),])
        numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
        if (numneigh1 == k) {
          j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
          j4=seq(distdata1[1, x2],distdata1[2, x2])
          inds<-reachdist[j3,x2]
          inds1<-reachdist[j4,x2]
          inds2<-match(dist.obj$id[x2,],reachdist[j3,x2])
          reachdist[j3,x2]=dist.obj$id[x2,]
          reachdist[j4,x2]=inds1[inds2]
          x3<-length(which(is.na(reachdist[,x2])))
          x4<-which(is.na(reachdist[,x2]))
          if (x3!=0){
            temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                                x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                     x2])
            #reach = 1/(sum(apply(temp, 2, max))/numneigh)
            
            reach_dist_x2 = apply(temp1, 2, max)
            reachdist[x4,x2]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
          }
        } else if (numneigh1!=k) {
          print("warning numneigh!=k")
        }
        for (i in 1:dim(knn_x2)[1]) {
          x5<-knn_x2[i,1]
          x6<-which(reachdist[,x5]== x2)+numneigh1
          reachdist[x6,x5]=distdata1[length(distdata1[,x2]),x2]
          #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
          #x8<-which(reachdist[,x2]== dim(data2)[1])
          if (x2 %in% dist.obj$id[x5,]) {
            s_update=as.matrix(union(s_update,x5))
          }
          
        }
      }
      
    } else {
      x2<-rNN_x1[i,1]
      knn_x2<-as.matrix(dist.obj$id[x2,])
      knn_x2<-as.matrix(knn_x2[-c(which(knn_x2[,1]== dim(data2)[1])),])
      numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
      if (numneigh1 == k) {
        j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
        j4=seq(distdata1[1, x2],distdata1[2, x2])
        inds<-reachdist[j3,x2]
        inds1<-reachdist[j4,x2]
        inds2<-match(dist.obj$id[x2,],reachdist[j3,x2])
        reachdist[j3,x2]=dist.obj$id[x2,]
        reachdist[j4,x2]=inds1[inds2]
        x3<-length(which(is.na(reachdist[,x2])))
        x4<-which(is.na(reachdist[,x2]))
        if (x3!=0){
          temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                              x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                   x2])
          #reach = 1/(sum(apply(temp, 2, max))/numneigh)
          
          reach_dist_x2 = apply(temp1, 2, max)
          reachdist[x4,x2]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
        }
      } else if (numneigh1!=k) {
        print("warning numneigh!=k")
      }
      for (i in 1:dim(knn_x2)[1]) {
        x5<-knn_x2[i,1]
        x6<-which(reachdist[,x5]== x2)+numneigh1
        reachdist[x6,x5]=distdata1[length(distdata1[,x2]),x2]
        #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
        #x8<-which(reachdist[,x2]== dim(data2)[1])
        if (x2 %in% dist.obj$id[x5,]) {
          s_update=as.matrix(union(s_update,x5))
        }
        
      }
    }
  }
  
  if (numneigh == k) {
    ind<-dim(data2)[1]
    #j = seq(3, 3 + (distdata1[2, i] - distdata1[1, i]))
    #numneigh = distdata1[2, i] - distdata1[1, i] + 1
    #temp = rbind(diag(distdata1[distdata1[2, distdata1[j, 
    #                                                   i]], distdata1[j, i]]), distdata1[j + numneigh, 
    #                                                                                     i])
    #reach = 1/(sum(apply(temp, 2, max))/numneigh)
    reachdist=cbind(reachdist,neighdist)
    ind1<-seq(reachdist[1, ind],reachdist[2, ind])
    #reach_dist_x1 = apply(temp, 2, max) #reach_dist(x_1, knn)
    reachdist[ind1,ind]=reach_dist_x1
  } else if (numneigh1!=k) {
    print("warning numneigh!=k")
  }
  
  
  
  
  #update lrd
  s_update_lof<-s_update
  for (i in 1:length(s_update)) {
    if (length(s_update)==1){
      if (length(s_update[[1]])!=0){
        x7<-s_update[i]
        numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
        if (numneigh2 == k) {
          j5=seq(reachdist[1, x7],reachdist[2, x7])
          reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
          lrddata[x7] = reach2
          #dist.obj <- dbscan::kNN(data2, k)
          #n <- nrow(data2)
          rNN_x7 <- sapply(x7, function(i) {
            as.vector(which(dist.obj$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                     1])
          })
          if (length(rNN_x7)==1){
            if (length(rNN_x7[[1]])!=0){
              s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
            }
          } else {
            s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
          }
          
        } else if (numneigh2!=k) {
          print("warning numneigh!=k")
        }
      }
    } else {
      x7<-s_update[i,1]
      numneigh2 = reachdist[2, x7] - reachdist[1, x7] + 1
      if (numneigh2 == k) {
        j5=seq(reachdist[1, x7],reachdist[2, x7])
        reach2 = 1/(sum(reachdist[j5,x7])/numneigh2)
        lrddata[x7] = reach2
        #dist.obj <- dbscan::kNN(data2, k)
        #n <- nrow(data2)
        rNN_x7 <- sapply(x7, function(i) {
          as.vector(which(dist.obj$id[, 1:k] == i, arr.ind = TRUE)[, 
                                                                   1])
        })
        if (length(rNN_x7)==1){
          if (length(rNN_x7[[1]])!=0){
            s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
          }
        } else {
          s_update_lof=as.matrix(union(s_update_lof,rNN_x7))
        }
      } else if (numneigh2!=k) {
        print("warning numneigh!=k")
      }
    }
    
    
  }
  
  #calculate lrd for x_1
  x9<-dim(data2)[1]
  numneigh3 = reachdist[2, x9] - reachdist[1, x9] + 1
  if (numneigh3 == k) {
    j6=seq(reachdist[1, x9],reachdist[2, x9])
    reach3 = 1/(sum(reachdist[j6,x9])/numneigh3)
    lrddata[x9] = reach3
    #dist.obj <- dbscan::kNN(data2, k)
    #n <- nrow(data2)
  } else if (numneigh3!=k) {
    print("warning numneigh!=k")
  }
  
  
  
  #update lof
  if (length(s_update_lof)==1){
    if (length(s_update_lof[[1]])!=0){
      s_update_lof<-as.matrix(s_update_lof[-c(which(s_update_lof[]== dim(data2)[1])),])
      for (i in 1:length(s_update_lof)) {
        x8<-s_update_lof[i]
        nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
        if (nneigh1 == k) {
          j8 = seq(0, (nneigh1 - 1))
          local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
          lof[x8] = local.factor1
        } else if (nneigh1!=k) {
          print("warning numneigh!=k")
        }
      }
      
    }
  } else if(length(s_update_lof)>1){
    s_update_lof<-as.matrix(s_update_lof[-c(which(s_update_lof[,1]== dim(data2)[1])),])
    for (i in 1:length(s_update_lof)) {
      x8<-s_update_lof[i,1]
      nneigh1 = distdata1[2, x8] - distdata1[1, x8] + 1
      if (nneigh1 == k) {
        j8 = seq(0, (nneigh1 - 1))
        local.factor1 = sum(lrddata[distdata1[3 + j8, x8]]/lrddata[x8])/nneigh1
        lof[x8] = local.factor1
      } else if (nneigh1!=k) {
        print("warning numneigh!=k")
      }
    }
    
  }
  
  outlier_detection15[m1,4]<-length(union(s_update_lof, s_update_lof_d))
  #calculate lof for x_1
  x9<-dim(data2)[1]
  nneigh2 = distdata1[2, x9] - distdata1[1, x9] + 1
  if (nneigh2 == k) {
    j7 = seq(0, (nneigh2 - 1))
    local.factor2 = sum(lrddata[distdata1[3 + j7, x9]]/lrddata[x9])/nneigh2
    lof[x9] = local.factor2
    if (lof[x9]>4.5){
      outlier_detection15[m1,2]<-1
      outlier_detection15[m1,3]<-lof[x9]
      
    } else {
      outlier_detection15[m1,2]<-0
      outlier_detection15[m1,3]<-lof[x9]
    }
  } else if (nneigh2!=k) {
    print("warning numneigh!=k")
  }
  
  
}

mean(outlier_detection15[,4])
length(which(outlier_detection15[,2]==1))
which(outlier_detection15[,2]==1)
#outlier3<-union(which(outlier_detection2[,2]==1),which(lof[201:400] > 1.35) )
outlier3<-which(outlier_detection15[,2]==1)
#which(c(1:300)%%15==0)
#length(which(lof[201:400] > 1.2))
m

TP_V<-length(which(outlier3 %in% m))
FN_V<-length(m)-TP_V
FP_V<-length(which(outlier_detection15[,2]==1))-TP_V
TN_V<-1500-length(m)-FP_V
SE_V<-TP_V/(TP_V+FN_V) #true positive rate
SP_V<-TN_V/(FP_V+TN_V) #true negative rate
G_mean_lof_SWILOF_SV<-sqrt(SE_V*SP_V)
H_mean_lof_SWILOF_SV<-2*SE_V*SP_V/(SE_V+SP_V)

outlier_detection16<-as.data.frame(outlier_detection15)
outlier_detection16[,5]=(1:1500)
outlier_detection16$V2<-as.factor(outlier_detection16$V2)
colnames(outlier_detection16)[2] <- "Outliers"

# plot
p15 <- ggplot(outlier_detection16, aes(x=V5, y=V1,label=V5)) +
  geom_point(aes(colour= Outliers)) + 
  #geom_point(aes(y=deaths), colour='red') +
  #geom_errorbar(aes(ymin=V4, ymax=V3)) +
  geom_text(aes(label=ifelse(Outliers==1,as.character(V5),'')),hjust=0,vjust=0)+
  theme_bw() +            
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))+
  #scale_y_continuous(expand=c(0,0), lim=c(0,1e3)) +
  labs(x='timepoint of streaming data', y='value of observation', title='SW-ILOF in non-stationary data stream (Extreme case) \n k=15, alpha=4.5, W=50') 
#coord_cartesian(ylim=c(0,300)) +
#facet_grid(~age_group)
p15







s_update_lrd=rNN_x_1
j3=seq(3, 3 + (distdata1[2, 1] - distdata1[1, 1]))
j4=seq(distdata1[1, 1],distdata1[2, 1])
reachdist[j3,]=reachdist[j3,]-1
#update reach_dist
for (i in 1:length(rNN_x_1)) {
  if (length(rNN_x_1)==1) {
    if (length(rNN_x_1[[1]])!=0) {
      x2<-rNN_x_1[i]
      knn_x2<-as.matrix(dist.obj2$id[x2-1,-c(length(dist.obj2$id[x2-1,]))])
      numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
      if (numneigh1 == k) {
        j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
        j4=seq(distdata1[1, x2],distdata1[2, x2])
        inds<-reachdist[j3,x2]
        inds1<-reachdist[j4,x2]
        inds2<-match(dist.obj2$id[x2-1,],reachdist[j3,x2])
        reachdist[j3,x2]=dist.obj2$id[x2-1,]
        reachdist[j4,x2]=inds1[inds2]
        x3<-length(which(is.na(reachdist[,x2])))
        x4<-which(is.na(reachdist[,x2]))
        if (x3!=0){
          temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                              x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                   x2])
          #reach = 1/(sum(apply(temp, 2, max))/numneigh)
          
          reach_dist_x2 = apply(temp1, 2, max)
          reachdist[x4,x2]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
        }
      } else if (numneigh1!=k) {
        print("warning numneigh!=k")
      }
      for (i in 1:dim(knn_x2)[1]) {
        x5<-knn_x2[i,1]
        x6<-which(reachdist[,x5+1]== x2-1)+numneigh1
        reachdist[x6,x5+1]=distdata1[length(distdata1[,x2]),x2]
        #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
        #x8<-which(reachdist[,x2]== dim(data2)[1])
        if ((x2-1) %in% dist.obj2$id[x5,]) {
          s_update_lrd=as.matrix(union(s_update_lrd,x5+1))
        }
        
      }
    }
    
  } else {
    x2<-rNN_x_1[i,1]
    knn_x2<-as.matrix(dist.obj2$id[x2-1,-c(length(dist.obj2$id[x2-1,]))])
    numneigh1 = distdata1[2, x2] - distdata1[1, x2] + 1
    if (numneigh1 == k) {
      j3=seq(3, 3 + (distdata1[2, x2] - distdata1[1, x2]))
      j4=seq(distdata1[1, x2],distdata1[2, x2])
      inds<-reachdist[j3,x2]
      inds1<-reachdist[j4,x2]
      inds2<-match(dist.obj2$id[x2-1,],reachdist[j3,x2])
      reachdist[j3,x2]=dist.obj2$id[x2-1,]
      reachdist[j4,x2]=inds1[inds2]
      x3<-length(which(is.na(reachdist[,x2])))
      x4<-which(is.na(reachdist[,x2]))
      if (x3!=0){
        temp1 = rbind(diag(distdata1[distdata1[2, distdata1[j3, 
                                                            x2]], distdata1[j3, x2]]), distdata1[j3 + numneigh1, 
                                                                                                 x2])
        #reach = 1/(sum(apply(temp, 2, max))/numneigh)
        
        reach_dist_x2 = apply(temp1, 2, max)
        reachdist[x4,x2]=reach_dist_x2[x4-(distdata1[1, x2]-1)]
      }
    } else if (numneigh1!=k) {
      print("warning numneigh!=k")
    }
    for (i in 1:dim(knn_x2)[1]) {
      x5<-knn_x2[i,1]
      x6<-which(reachdist[,x5+1]== x2-1)+numneigh1
      reachdist[x6,x5+1]=distdata1[length(distdata1[,x2]),x2]
      #x7<-length(which(reachdist[,x2]== dim(data2)[1]))
      #x8<-which(reachdist[,x2]== dim(data2)[1])
      if ((x2-1) %in% dist.obj2$id[x5,]) {
        s_update_lrd=as.matrix(union(s_update_lrd,x5+1))
      }
      
    }
  }
}











