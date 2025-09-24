# Creating synthetic network based complex data input

library(spatstat)
#library(tidyverse)
library(extraDistr)
#library(ergm)
library(igraph)
library(netrankr)
library(DiscreteLaplace)
library(lpSolve)
library(ggraph)



# --------------------------- start private network generation --------------------------- # 

# --------------------------- parameter input --------------------------- # 
eps = 0.1 #privacy parameter
bs=1 #box size
maxiter = 1000 #number of iterations
example = 3 #one of 1,2,3
#if example == 3: make n dependent on eps
n = 50000
seed = runif(1,min = 1, max = 100000)

# --------------------------- data input --------------------------- #
# --------------------------- Example 1 --------------------------- #
if(example == 1){
  #data input
  data(florentine_m)
  wealth = V(florentine_m)$wealth
  
  #normalize (to be in [0,1])
  wealth = wealth/max(wealth)
  
  x = matrix(wealth,ncol = d) #true data input
  n = length(x)
  d = 1
  m = 10
  a = 16 #expected number points in network based on observed data
  b = 16 #expected number points in network based on private data
  c = 2  #constant factor in edge creation (k(x,y) = min(1,x*y*c))
  
  ## define function kappa
  kappa <- function(x,y){
    return (min(1,x*y*c))
  }
}

# --------------------------- Example 2 --------------------------- #
if(example == 2){
  n = 1000
  d = 1
  x = c(runif(n/2, 0, 0),runif(n/2, 1,1))
  c = 1   #constant factor in edge creation (k(x,y) = min(1,x*y*c))
  m = round((eps*(n))^(d/(d+1))) #optimal partition size
  a = m^(2/d) #expected number points in network based on observed data
  b = m^(2/d) #expected number points in network based on private data
  
  ## define function kappa
  kappa <- function(x,y){
    return (min(1,x*y*c))
  }
}

# --------------------------- Example 3 --------------------------- #
if(example == 3){
  # data input: census 2011 teaching data 
  # https://www.ons.gov.uk/census/2011census/2011censusdata/censusmicrodata/microdatateachingfile)
  data = read.csv("2011 Census Microdata Teaching File.csv")
  colnames(data) = data[1,]
  data = data[-1,]
  data_wo_students = data[data$Student != 1,] #exclude students as they have many "-9"
  data_wo_students = data_wo_students[data_wo_students$`Economic Activity` != -9,]
  
  x <- data_wo_students[c("Age", "Health", "Economic Activity")]
  x <- x[1:n,]
  
  #normalize data
  x$Age <- as.numeric(as.character(x$Age))
  x$Age <- x$Age/max(x$Age)
  x$Health <- as.numeric(as.character(x$Health))
  x$Health <- x$Health/max(x$Health)
  x$`Economic Activity` <- as.numeric(as.character(x$`Economic Activity`))
  x$`Economic Activity` <- x$`Economic Activity`/max(x$`Economic Activity`)
  
  #parameters
  d=3
  n = length(x$Age)
  m = round((eps*n)^(1/(d+1))) #boxes in 1 dim such that m_overall = m^d is optimal
  m_overall = m^d #overall number of boxes
  a = m_overall^(2/d) #optimal graph size
  b = a
  
  ## define function kappa
  kappa <- function(x,y){
    # x,y 3-dim vectors
    f_age = abs(x[1]-y[1]) #similar age is favorable for edges
    f_health = abs(1-x[2]*y[2]) #the smaller, the "healthier" and thus the higher prob of edges
    f_ea = abs(x[3]-y[3]) #similar economic activity is favourable for edges
    return (1/3*(f_age +f_health + f_ea))
  }
}
# --------------------------- run several iterations --------------------------- # 
set.seed(seed)
timing = rep(0,maxiter)
dist_FGW= rep(0,maxiter)
timing_FGW = rep(0,maxiter)
#print parameters
print(paste0("graph size a = ",a))
print(paste0("partition size m = ",m_overall))


for (iter in 1:maxiter){
  print(iter)
# --------------------------- apply privacy mechanism on node level --------------------------- # 
  #partiton
  m_overall = m^d #overall number of boxes
  partl = (bs/m)
  
  #find beginning of each box in 1 dim
  low_bounds_1d <- seq(0, bs - partl, length.out = m)
  
  # get beginnings and endings of all boxes in d-dim
  lower_bounds <- expand.grid(replicate(d, low_bounds_1d, simplify = FALSE))
  higher_bounds <- lower_bounds + partl
  
  #round up to avoid problems with different rounding in x and lower_bounds
  #lower_bounds <- round(lower_bounds,10) #round to avoid floating problems
  #higher_bounds <- round(higher_bounds,10) #round to avoid floating problems
  
  #add a little bit to boundary bs to use "<" in the following
  higher_bounds[higher_bounds == max(higher_bounds)] <- max(higher_bounds) +0.000001
  
  #-----------------algorithm starts here -------------------#
  timing[iter]<- system.time({
  #true counts
  tcounts = rep(0,m_overall)

  assigned <- rep(FALSE, n)  # track which points have already been counted
  
  # #for each box
  #for (k in 1:m_overall) {
  #  print(k)
  #  # Check which rows of x fall into the k-th box
  #  in_box <- apply(x, 1, function(row) all(row >= lower_bounds[k,] & row < higher_bounds[k,]))
  #  tcounts[k] <- sum(in_box)
  #}
  
  # Vectorized over all boxes
  x_mat <- as.matrix(x)#round_up(x,10))
  for (k in 1:m_overall) {
    lb <- lower_bounds[k, ]
    hb <- higher_bounds[k, ]
    
    # Only check unassigned points
    unassigned_x <- x_mat[!assigned, , drop = FALSE]
    
    # Vectorized logical condition
    eps_box <- 1e-12
    in_box <- rowSums(sweep(unassigned_x, 2, lb-eps_box, `>=`) & sweep(unassigned_x, 2, hb+eps_box, `<`)) == ncol(unassigned_x)
    
    tcounts[k] <- sum(in_box)
    
    # Update assigned mask
    # Need to find row indices in full x_mat
    if (tcounts[k] > 0) {
      matched_rows <- which(!assigned)[in_box]
      assigned[matched_rows] <- TRUE
    }
  }
  sum(tcounts)-n #should be 0
  

  #new dataset
  y = matrix(rep(0,n*d),ncol = d)
  y_rep = matrix(rep(0,m_overall*d),ncol = d)

  i = 0
  for(k in 1:m_overall){
    for (j in 1:d){
      y_rep[k,j] = runif(1,min = lower_bounds[k,j], higher_bounds[k,j])
    }
    if(tcounts[k]>0){
      y[(i+1):(i+tcounts[k]),] = matrix(rep(y_rep[k,], each = tcounts[k]), ncol = d)
      i = i+tcounts[k] 
    }
  }

  #add noise
  p = exp(-eps) #dLap(1/eps)
  noise = rdlaplace(m_overall,p,p) #p=q

  v = (tcounts+noise)/n

  #linear programming
  ## objective 
  obj = append(rep(1,m_overall),rep(0,m_overall))
  
  ##constraints matrix
  A1 = cbind(diag(nrow = m_overall),diag(nrow = m_overall))
  A2 = cbind(diag(nrow = m_overall),diag(x = -1,nrow = m_overall))
  A3 = cbind(diag(x = 0, nrow = m_overall),diag(nrow = m_overall))
  A4 = append(rep(0,m_overall),rep(1,m_overall))
  A = rbind(A1,A2,A3,A4)
  
  ##constr
  B1 = v
  B2 = -v
  B3 = rep(0,m_overall)
  B4 = 1
  B = append(append(append(B1,B2),B3),B4)
  
  constraints_dir = append(rep(">=",3*m_overall),"=")
  
  ##solve
  v_hat = lp(direction = "min", objective.in = obj, const.mat = A, const.dir = constraints_dir, const.rhs = B)$solution[(m_overall+1):(2*m_overall)]
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
  Z_vec = rmultinom(L,1,p)[1:m_overall,] #last entry corresponds to sum lk = 0 
  Z = sum(Z_vec)

  #create non-common vertex counts
  M_vec = rmultinom(M-Z,1, v_hat)
  N_vec = rmultinom(N-Z,1,tcounts/n)

  #create vertices
  xi = data.frame("id" = 1:M, "k" = rep(0,M), "M_k" = rep(0,M))
  eta = data.frame("id" = 1:N, "k" = rep(0,N), "N_k" = rep(0,N))
  
  ## add locations
  for (j in 1:d) {
    xi[[paste0("y", j)]] <- numeric(M)
    eta[[paste0("y", j)]] <- numeric(N)
  }
  
  k_xi = 0
  k_eta = 0

  for(k in 1:m_overall){
    Mk = sum(Z_vec[k,]) + sum(M_vec[k,])
    Nk = sum(Z_vec[k,]) + sum(N_vec[k,])
    #print(paste("k = ",k, "Mk = ",Mk, "Nk = ", Nk))
    for(i in 1:Mk){
      if(Mk >0){
        for(j in 1:d){
          xi[[paste0("y",j)]][k_xi + i] = y_rep[k,j]
        }
        xi$M_k[k_xi + i] = i
        xi$k[k_xi + i] = k
      }
    }
  
    for(i in 1:Nk){
      if (Nk > 0) {
        #get all points in box k (d-dim)
        rows_in_box <- which(apply(x, 1, function(pt) {
          all(pt >= lower_bounds[k, ] & pt <= higher_bounds[k, ])
        }))
        ptsinbox <- as.matrix(x[rows_in_box, , drop = FALSE])
        pt_id <- rdunif(1 ,min = 1,length(ptsinbox[,1]))
        for(j in 1:d){
          eta[[paste0("y",j)]][k_eta + i] = ptsinbox[pt_id,j]
        }
        eta$N_k[k_eta + i] = i
        eta$k[k_eta + i] = k
      }
    }
    k_xi = k_xi + Mk
    k_eta = k_eta + Nk
  }

  #create edges


  ## create edges for common vertex counts
  sigma = data.frame("id" = 1:(M*(M-1)/2), from = rep(0, M*(M-1)/2), to = rep(0, M*(M-1)/2))
  tau = data.frame("id" = 1:(N*(N-1)/2), from = rep(0, N*(N-1)/2), to = rep(0, N*(N-1)/2))

  id_s = 0
  id_t = 0

  for (k in 1:m_overall) {
    for (l in k:m_overall){
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
          kappa_t = max(0, kappa(eta[eta$k == k & eta$N_k == i,paste0("y",1:d)],eta[eta$k == l & eta$N_k == j,paste0("y",1:d)]))
          kappa_s = max(0,kappa(xi[xi$k == k & xi$M_k == i,paste0("y",1:d)],xi[xi$k == l & xi$M_k == j,paste0("y",1:d)]))
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
            kappa_t = max(0, kappa(eta[eta$k == k & eta$N_k == i + sum(Z_vec[k,]), paste0("y",1:d)],eta[eta$k == l & eta$N_k == j + sum(Z_vec[l,]), paste0("y",1:d)]))
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
            kappa_s = max(0,kappa(max(0,xi$y1[xi$k == k & xi$M_k == i + sum(Z_vec[k,])]),max(0,xi$y1[xi$k == l & xi$M_k == j + sum(Z_vec[l,])])))
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

  # --------------------------- create graphs --------------------------- #
  edges_sigma = cbind("from" =sigma$from[sigma$from>0],"to" = sigma$to[sigma$to>0])
  edges_tau = cbind("from" =tau$from[tau$from>0],"to" = tau$to[tau$to>0])
  
  g_s <- igraph::graph_from_data_frame(d = edges_sigma, directed = FALSE, vertices = cbind(xi$id,xi[paste0("y",1:d)]))
  g_t <- igraph::graph_from_data_frame(d = edges_tau, directed = FALSE, vertices = cbind(eta$id,eta[paste0("y",1:d)]))
  })["elapsed"] #end of algorithm
  
  C_sigma <- as_adjacency_matrix(g_s)
  C_tau <- as_adjacency_matrix(g_t)
  A_sigma <- vertex_attr(g_s)[2:(d+1)]
  A_tau <- vertex_attr(g_t)[2:(d+1)]
  
  # --------------------------- plot graphs --------------------------- # 

  #edges_sigma = cbind("from" =sigma$from[sigma$from>0],"to" = sigma$to[sigma$to>0])
  #edges_tau = cbind("from" =tau$from[tau$from>0],"to" = tau$to[tau$to>0])

  #plot original graph
  #plot_gTau <- ggraph(g_t, layout = "graphopt")+ #layout = cbind(eta$coord_x1,eta$coord_x2))+ #, label = pointset$name)+
  #  geom_edge_link(edge_width = 1, color = "grey")+
  #  geom_node_point(aes(color = eta$y1), size=3, show.legend = FALSE)+
  #  #scale_color_gradient(low="lightblue3", high="darkblue")+
  #  #geom_node_point(color = "darkred", size =3)+
  #  #geom_node_text(aes(label = verticesG1), color = 'white', size =0)+
  #  #coord_fixed(xlim=c(0.5,2.35), ylim = c(0.4,2.3))+
  #  #coord_fixed(xlim=c(min(eta$coord_x1),max(eta$coord_x1)), ylim = c(min(eta$coord_x2),max(eta$coord_x2)))+
  #  theme_void()+
  #  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
  #        panel.background = element_blank(),
  #        panel.grid.major = element_blank(),
  #        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
  #        panel.grid.minor = element_blank())
  #plot_gTau
  #
  ##match locations to other graph  
  ###recover vertex locations in realization based on observed data
  #eta$coord_x1 = plot_gTau$data$x
  #eta$coord_x2 = plot_gTau$data$y
  #
  #xi$coord_x1 = rep(0,length(xi$id))
  #xi$coord_x2 = rep(0,length(xi$id))
  ###assign locations to private graph using the coupling
  #for(k in 1:m){
  #  for (i in 1:sum(Z_vec[k,])) {
  #    if(sum(Z_vec[k,]) >0){
  #      xi$coord_x1[xi$k == k & xi$M_k == i] = eta$coord_x1[xi$k == k & xi$M_k == i]
  #      xi$coord_x2[xi$k == k & xi$M_k == i] = eta$coord_x2[xi$k == k & xi$M_k == i]
  #    }
  #    
  #  }
  #  for(i in 1:sum(N_vec[k,])){
  #    if(sum(N_vec[k,]) >0){
  #      xi$coord_x1[xi$k == k & xi$M_k == i + sum(Z_vec[k,])] = runif(1,min = min(eta$coord_x1), max = max(eta$coord_x1))
  #      xi$coord_x2[xi$k == k & xi$M_k == i+ sum(Z_vec[k,])] = runif(1,min = min(eta$coord_x2), max = max(eta$coord_x2))
  #    }
  #  }
  #}
  #plot_gSigma <- ggraph(g_s, layout = cbind(xi$coord_x1,xi$coord_x2))+#, label = pointset$name)+
  #  geom_edge_link(edge_width = 1, color = "grey")+
  #  geom_node_point(aes(color = xi$y1), size=3, show.legend = FALSE)+
  #  #scale_color_gradient(low="lightblue3", high="darkblue")+
  #  #geom_node_point(color = "darkred", size =3)+
  #  #geom_node_text(aes(label = verticesG1), color = 'white', size =0)+
  #  #coord_fixed(xlim=c(min(eta$coord_x1),max(eta$coord_x1)), ylim = c(min(eta$coord_x2),max(eta$coord_x2)))+
  #  theme_void()+
  #  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
  #        panel.background = element_blank(),
  #        panel.grid.major = element_blank(),
  #        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
  #        panel.grid.minor = element_blank())
  ##axis.line = element_line(colour = "black"),
  ##panel.border = element_rect(colour = "black",fill=NA, linewidth = 1))
  #plot_gSigma


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
    write.csv(eta$y1, file_name2)
    
    write.csv(cbind( sigma$from[sigma$from>0], sigma$to[sigma$from>0]), file_name3)
    write.csv(eta$y1, file_name4)
  }
  
  setwd(path)
  
  ## run python code to calculate FGW
  library(reticulate)
  
  ot <- import("ot")
  np <- import("numpy")
  alpha = 0.5
  
  #pairwise distances of node features
  Ma = ot$dist((np$asarray(eta[paste0("y",1:d)])),(np$asarray(xi[paste0("y",1:d)])))
  
  #adjacency matrices
  adj_sigma = as.matrix(as_adjacency_matrix(g_s))
  adj_tau = as.matrix(as_adjacency_matrix(g_t))
  
  timing_FGW[iter] <- system.time({
    FGW = ot$fused_gromov_wasserstein2(M = Ma , C1 = np$asmatrix(adj_tau), C2 = np$asmatrix(adj_sigma), alpha=alpha,features_metric='sqeuclidean')
  })["elapsed"]
  dist_FGW[iter] <- FGW
}

#save data
results <- data.frame("dist_FGW" = dist_FGW, "runtime GG" = timing, "runtime FGW" = timing_FGW )
# ----------------------- Save data ----------------------- #
file_name <- paste(path,"/Results/PSGGExample",example,"withEps=",eps, "seed=",seed,"n=",n,"a=",round(a), "m = ",m ,".csv",col="", sep="")
write.csv(results, file_name)
