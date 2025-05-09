# Creating synthetic network based complex data input

library(spatstat)
library(tidyverse)
library(extraDistr)
#library(ergm)
library(igraph)
library(netrankr)
library(DiscreteLaplace)
library(lpSolve)
library(ggraph)



# --------------------------- start private network generation --------------------------- # 

# --------------------------- parameter input --------------------------- # 
eps = 2 #privacy parameter
bs=1 #box size
m = 10 #size of partition (we consider the uniform partition of the d-dim cube of size l)
partl = bs/m

# --------------------------- data input --------------------------- #
#Example 1 
#data input
data(florentine_m)
wealth = V(florentine_m)$wealth

#normalize (to be in [0,1])
wealth = wealth/max(wealth)

x = wealth #true data input
n = length(x)
a = 16 #expected number points in network based on observed data
b = 16 #expected number points in network based on private data
c = 2  #constant factor in edge creation (k(x,y) = min(1,x*y*c))


#Example 2 
#n = 100
#x = c(runif(n/2, 0, 0),runif(n/2, 1,1))
#a = 100 #expected number points in network based on observed data
#b = 100 #expected number points in network based on private data
#c = 1   #constant factor in edge creation (k(x,y) = min(1,x*y*c))
#m = round(sqrt(eps)*sqrt(n))
#partl = bs/m

# --------------------------- apply privacy mechanism on node level --------------------------- # 

#partiton

setofboxes = data.frame(c(0,bs))
names(setofboxes) <- "Whole Space"
row.names(setofboxes) <- c("xmin","xmax")#, "ymin","ymax")

for (k in 1:m) {
  #for (j in 1:sqrt(n)){
  partbox <- c((k-1)*partl ,k*partl)#cbind(c((k-1)*partl ,k*partl), c((j-1)*partl , j*partl))
  setofboxes <- setofboxes %>% add_column(c(partbox)) 
  #}
}

#true counts
tcounts = rep(0,m)

for (i in 1:n) {
  for (k in 1:m) {
    #print(k)
    if(x[i]>= setofboxes[[k+1]][1] && x[i]<= setofboxes[[k+1]][2]){
      tcounts[k] = tcounts[k] +1
    }
  }
}

#new dataset
y = rep(0,n)
y_rep = rep(0,n)

i = 0
for(k in 1:m){
  y_rep[k] = runif(1,min = setofboxes[[k+1]][1], setofboxes[[k+1]][2])
  if(tcounts[k]>0){
    y[(i+1):(i+tcounts[k])] = rep(y_rep[k], tcounts[k])
    i = i+tcounts[k] 
  }
}

#add noise
p = exp(-eps) #dLap(1/eps)
noise = rdlaplace(m,p,p) #p=q

v = (tcounts+noise)/n

#linear programming
## objective 
obj = append(rep(1,m),rep(0,m))

##constraints matrix
A1 = cbind(diag(nrow = m),diag(nrow = m))
A2 = cbind(diag(nrow = m),diag(x = -1,nrow = m))
A3 = cbind(diag(x = 0, nrow = m),diag(nrow = m))
A4 = append(rep(0,m),rep(1,m))
A = rbind(A1,A2,A3,A4)

##constr
B1 = v
B2 = -v
B3 = rep(0,m)
B4 = 1
B = append(append(append(B1,B2),B3),B4)

constraints_dir = append(rep(">=",3*m),"=")

##solve
v_hat = lp(direction = "min", objective.in = obj, const.mat = A, const.dir = constraints_dir, const.rhs = B)$solution[(m+1):(2*m)]
v_hat
# --------------------------- jointly generate network data based on Chung-Lu model --------------------------- # 

#sample graph size
Kn = rpois(1,(a-min(a,b))*bs)
Km = rpois(1,(b-min(a,b))*bs)
L = rpois(1,min(a,b)*bs)
N = Kn + L
M = Km + L

#create common vertex counts 
p_vec = pmin(v_hat,tcounts/n)
p = append(p_vec, 1-sum(p_vec)) 
Z_vec = rmultinom(L,1,p)[1:m,] #last entry corresponds to sum lk = 0 
Z = sum(Z_vec)

#create non-common vertex counts
M_vec = rmultinom(M-Z,1, v_hat)
N_vec = rmultinom(N-Z,1,tcounts/n)

#create vertices
xi = data.frame("id" = 1:M," y" = rep(0,M), "k" = rep(0,M), "M_k" = rep(0,M))
eta = data.frame("id" = 1:N, "y" = rep(0,N), "k" = rep(0,N), "N_k" = rep(0,N))

k_xi = 0
k_eta = 0

for(k in 1:m){
  Mk = sum(Z_vec[k,]) + sum(M_vec[k,])
  Nk = sum(Z_vec[k,]) + sum(N_vec[k,])
  #print(paste("k = ",k, "Mk = ",Mk, "Nk = ", Nk))
  for(i in 1:Mk){
    if(Mk >0){
      xi$y[k_xi + i] = y_rep[k]
      xi$M_k[k_xi + i] = i
      xi$k[k_xi + i] = k
    }
  }
  
  for(i in 1:Nk){
    if (Nk > 0) {
      ptsinbox = x[x <= setofboxes[[k+1]][2] & x >= setofboxes[[k+1]][1]]
      eta$y[k_eta + i] = ptsinbox[rdunif(1 ,min = 1,length(ptsinbox))]
      eta$N_k[k_eta + i] = i
      eta$k[k_eta + i] = k
    }
  }
  k_xi = k_xi + Mk
  k_eta = k_eta + Nk
}

#create edges
## define function kappa
kappa <- function(x,y){
  return (min(1,x*y*c))
}

## create edges for common vertex counts
sigma = data.frame("id" = 1:(M*(M-1)/2), from = rep(0, M*(M-1)/2), to = rep(0, M*(M-1)/2))
tau = data.frame("id" = 1:(N*(N-1)/2), from = rep(0, N*(N-1)/2), to = rep(0, N*(N-1)/2))

id_s = 0
id_t = 0

for (k in 1:m) {
  for (l in k:m){
    for(i in 1:sum(Z_vec[k,])){
      if(sum(Z_vec[k,]) == 0){
        next
      }
      for(j in 1:sum(Z_vec[l,])){
        #if(max(tau$from == eta$id[eta$k == k & eta$N_k == i + sum(Z_vec[k,])] & tau$to == eta$id[eta$k == l & eta$N_k == j + sum(Z_vec[l,])])==1){
          #if edge already exists
        #  print("double edge in Z_vec for tau")
        #  next
        #}
        #if(max(sigma$from == xi$id[xi$k == k & xi$M_k == i + sum(Z_vec[k,])] & sigma$to == xi$id[xi$k == l & xi$M_k == j + sum(Z_vec[l,])])==1){
          #if edge already exists
        #  print("double edge in Z_vec for tau")
        #  next
        #}
        if(sum(Z_vec[l,]) == 0){
          next
        }
        if(k == l && i == j){
          next
        }
        # states (tau,sigma) are (1,1), (1,0), (0,1), (0,0)
        kappa_t = max(0, kappa(max(0,eta$y[eta$k == k & eta$N_k == i]),max(0,eta$y[eta$k == l & eta$N_k == j])))
        kappa_s = max(0,kappa(max(0,xi$y[xi$k == k & xi$M_k == i]),max(0,xi$y[xi$k == l & xi$M_k == j])))
        #print(paste("kappa_t ", kappa_t, "k = ", k, "i = ",i,"l = ", l, "j = ", j))
        p_1 = min(kappa_t,kappa_s)
        p_2 = kappa_t - p_1
        p_3 = kappa_s - p_1
        edges = rmultinom(1,1,c(p_1,p_2,p_3,1-p_1-p_2-p_3))
        
        if(edges[1] == 1 || edges[2] ==1 ){
          #add edge to tau
          tau$from[id_t] = eta$id[eta$k == k & eta$N_k == i]
          tau$to[id_t] = eta$id[eta$k == l & eta$N_k == j]
          id_t = id_t + 1
        }
        if(edges[1] == 1 || edges[3] == 1){
          #add edge to sigma
          sigma$from[id_s] = xi$id[xi$k == k & xi$M_k == i]
          sigma$to[id_s] = xi$id[xi$k == l & xi$M_k == j]
          id_s = id_s + 1
        }
      }
    }
      for(i in 1:sum(N_vec[k,])){
        for(j in 1:sum(N_vec[l,])){
          if(sum(N_vec[k,]) == 0){
            next
          }
          if(sum(N_vec[l,]) == 0){
            next
          }
          if(max(tau$from == eta$id[eta$k == k & eta$N_k == i + sum(Z_vec[k,])] & tau$to == eta$id[eta$k == l & eta$N_k == j + sum(Z_vec[l,])])==1){
            #if edge already exists
            next
          }
          if(k == l && j ==i){
            #not self loops
            next
          }
          
          #create edges for non-common vertex counts in tau
          kappa_t = max(0, kappa(max(0,eta$y[eta$k == k & eta$N_k == i + sum(Z_vec[k,])]),max(0,eta$y[eta$k == l & eta$N_k == j + sum(Z_vec[l,])])))
          edge = rbern(1, prob = kappa_t)
          #print(paste(kappa_t, edge))
          if(edge == 1){
            #add edge to tau
            #print("edge to tau")
            tau$from[id_t] = eta$id[eta$k == k & eta$N_k == i + sum(Z_vec[k,])]
            tau$to[id_t] = eta$id[eta$k == l & eta$N_k == j + sum(Z_vec[l,])]
            id_t = id_t + 1
          }
        }
      }
      for(i in 1:sum(M_vec[k,])){
        for(j in 1:sum(M_vec[l,])){
          if(sum(M_vec[k,]) == 0){
            next
          }
          if(sum(M_vec[l,]) == 0){
            next
          }
          if(max(sigma$from == xi$id[xi$k == k & xi$M_k == i + sum(Z_vec[k,])] & sigma$from ==xi$id[xi$k == l & xi$M_k == j + sum(Z_vec[l,])])==1){
            #if edge already exists
            print("double edge detected")
            next
          }
          if(k == l && j ==i){
            #no self loops
            next
          }
          
          #create edges for non-common vertex counts in sigma
          kappa_s = max(0,kappa(max(0,xi$y[xi$k == k & xi$M_k == i + sum(Z_vec[k,])]),max(0,xi$y[xi$k == l & xi$M_k == j + sum(Z_vec[l,])])))
          edge = rbern(1, prob = kappa_s)
          if(edge == 1){
            #add edge to sigma
            #print(paste("Length of from",length(xi$id[xi$k == k & xi$M_k == i + sum(Z_vec[k,])]),"with i,j",i,j))
            sigma$from[id_s] = xi$id[xi$k == k & xi$M_k == i + sum(Z_vec[k,])]
            sigma$to[id_s] = xi$id[xi$k == l & xi$M_k == j + sum(Z_vec[l,])]
            #print(paste(sigma$from[id_s],sigma$to[id_s]))
            id_s = id_s + 1
          }
        
      }
    }
  }
  
}

# --------------------------- plot graphs --------------------------- # 

edges_sigma = cbind("from" =sigma$from[sigma$from>0],"to" = sigma$to[sigma$to>0])
edges_tau = cbind("from" =tau$from[tau$from>0],"to" = tau$to[tau$to>0])

g_s <- igraph::graph_from_data_frame(d = edges_sigma, directed = FALSE, vertices = cbind(xi$id,xi$y))
g_t <- igraph::graph_from_data_frame(d = edges_tau, directed = FALSE, vertices = cbind(eta$id,eta$y))

#plot original graph
plot_gTau <- ggraph(g_t, layout = "graphopt")+ #layout = cbind(eta$coord_x1,eta$coord_x2))+ #, label = pointset$name)+
  geom_edge_link(edge_width = 1, color = "grey")+
  geom_node_point(aes(color = eta$y), size=3, show.legend = FALSE)+
  #scale_color_gradient(low="lightblue3", high="darkblue")+
  #geom_node_point(color = "darkred", size =3)+
  #geom_node_text(aes(label = verticesG1), color = 'white', size =0)+
  #coord_fixed(xlim=c(0.5,2.35), ylim = c(0.4,2.3))+
  #coord_fixed(xlim=c(min(eta$coord_x1),max(eta$coord_x1)), ylim = c(min(eta$coord_x2),max(eta$coord_x2)))+
  theme_void()+
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank())
plot_gTau

#match locations to other graph  
##recover vertex locations in realization based on observed data
eta$coord_x1 = plot_gTau$data$x
eta$coord_x2 = plot_gTau$data$y

xi$coord_x1 = rep(0,length(xi$id))
xi$coord_x2 = rep(0,length(xi$id))
##assign locations to private graph using the coupling
for(k in 1:m){
  for (i in 1:sum(Z_vec[k,])) {
    if(sum(Z_vec[k,]) >0){
      xi$coord_x1[xi$k == k & xi$M_k == i] = eta$coord_x1[xi$k == k & xi$M_k == i]
      xi$coord_x2[xi$k == k & xi$M_k == i] = eta$coord_x2[xi$k == k & xi$M_k == i]
    }
    
  }
  for(i in 1:sum(N_vec[k,])){
    if(sum(N_vec[k,]) >0){
      xi$coord_x1[xi$k == k & xi$M_k == i + sum(Z_vec[k,])] = runif(1,min = min(eta$coord_x1), max = max(eta$coord_x1))
      xi$coord_x2[xi$k == k & xi$M_k == i+ sum(Z_vec[k,])] = runif(1,min = min(eta$coord_x2), max = max(eta$coord_x2))
    }
  }
}
plot_gSigma <- ggraph(g_s, layout = cbind(xi$coord_x1,xi$coord_x2))+#, label = pointset$name)+
  geom_edge_link(edge_width = 1, color = "grey")+
  geom_node_point(aes(color = xi$y), size=3, show.legend = FALSE)+
  #scale_color_gradient(low="lightblue3", high="darkblue")+
  #geom_node_point(color = "darkred", size =3)+
  #geom_node_text(aes(label = verticesG1), color = 'white', size =0)+
  #coord_fixed(xlim=c(min(eta$coord_x1),max(eta$coord_x1)), ylim = c(min(eta$coord_x2),max(eta$coord_x2)))+
  theme_void()+
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.minor = element_blank())
#axis.line = element_line(colour = "black"),
#panel.border = element_rect(colour = "black",fill=NA, linewidth = 1))
plot_gSigma


# --------------------------- computation of FGW distance in python file --------------------------- # 

# export graphs for computation in python

path = getwd()
## create folder
#if (dir.exists(file.path(path, "GraphData"))){
#  setwd(file.path(path, "GraphData"))
#} else {
#  dir.create(file.path(path, "GraphData"))
#  setwd(file.path(path, "GraphData"))
#}

#file_name1 <- paste("NetworkObservedDataEdgesBoxSize=",bs, "intensity=", a," , ", b ,"partition=", m, "eps=",eps, ".csv",col="", sep="")
#file_name2 <- paste("NetworkObservedDataWeightsBoxSize=",bs, "intensity=", a," , ", b ,"partition=", m, "eps=",eps, ".csv",col="", sep="")
#file_name3 <- paste("NetworkPrivacyMechEdgesBoxSize=",bs, "intensity=", a," , ", b ,"partition=", m, "eps=",eps, ".csv",col="", sep="")
#file_name4 <- paste("NetworkPrivacyMechWeightsBoxSize=",bs, "intensity=", a," , ", b ,"partition=", m, "eps=",eps, ".csv",col="", sep="")

if(FALSE){
  write.csv(cbind( tau$from[tau$from>0], tau$to[tau$from>0]), file_name1)
  write.csv(eta$y, file_name2)
  
  write.csv(cbind( sigma$from[sigma$from>0], sigma$to[sigma$from>0]), file_name3)
  write.csv(eta$y, file_name4)
}

setwd(path)

## run python code to calculate FGW
library(reticulate)

ot <- import("ot")
np <- import("numpy")
alpha = 0.5

#pairwise distances of node features
Ma = ot$dist(np$transpose(np$asmatrix(eta$y)),np$transpose(np$asmatrix(xi$y)))

#adjacency matrices
adj_sigma = as.matrix(as_adjacency_matrix(g_s))
adj_tau = as.matrix(as_adjacency_matrix(g_t))

FGW = ot$fused_gromov_wasserstein(M = Ma , C1 = np$asmatrix(adj_tau), C2 = np$asmatrix(adj_sigma), alpha=alpha,features_metric='sqeuclidean')
FGW[1]
N
