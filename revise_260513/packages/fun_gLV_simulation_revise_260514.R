# %% loading library and functions
library(igraph)
library(Matrix)
library(codyn)

# Functions ----------------------------------------------------------------

N0 = \(sp=4)( runif(sp, 0.1, 0.8) )
N0_absence = \(sp=4,min=-0.5,max=0.5){ru <- runif(sp,min,max)
                                      return(ifelse(ru < 0, 0, ru))}

makeA = function(N=100, connectance=0.3, power=1, mu=0, sigma=0.1, output_g=FALSE){
  
  g = igraph::sample_gnp(N, p=connectance, directed = TRUE)
  # g = igraph::sample_pa(N, power = power,
  #                       start.graph = g,  m=1)
  
  A = as_adjacency_matrix(g)
  link_posi = (A!=0)
  A[which(link_posi)] = rnorm(sum(link_posi), mu, sigma)*rlnorm(sum(link_posi), 0, 1)
  
  if(output_g) return(g)
  return(as.matrix(A))
}

makeB = function(N=100, connectance=0.5, power=1, mu=0, sigma=0.01, output_g=FALSE){
  
  g = igraph::sample_gnp(N, p=connectance, directed = TRUE)
  # g = igraph::sample_pa(N, power = power,
  #                       start.graph = g,  m=1)
  
  A = as_adjacency_matrix(g)
  link_posi = (A>0)
  A[which(link_posi)] = rnorm(sum(link_posi), mu, sigma)
  
  if(output_g) return(g)
  return(as.matrix(A))
}


library(deSolve)
N0 = \(sp=4)( runif(sp, 0.1, 0.8) )

lv_fun = function(time, N, parms) {
  with(parms, {
    
    N <- pmax(N, 0)
    dN <- N * (r + A %*% N)
    
    return(list(dN))
  })
}

lv_fun = \(A, r, N) {
  N * (r + A %*% N)+rnorm(length(N), 0, 0.01)*rlnorm(length(N), 0, 2)
}

lv_fun_nonoise = \(A, r, N) {
  N * (r + A %*% N)
}

# lv_fun = \(A, r, N) {
#   N * (r + A %*% N)+ifelse(N>0,rnorm(length(N), 0, 0.01)*rlnorm(length(N), 0, 2),0)
#   }
# ms_ricker = \(N, r, A) exp( r*(1+A%*%N) ) * N

gen_dyn = function(time, dt=1000, A, r, N0, type="nonoise"){
  
  mat = matrix(0, ncol=length(N0), nrow=time)
  mat[1,] = N = N0
  fn = ifelse(type=="nonoise", lv_fun_nonoise, lv_fun)
  
  for(t in 2:time){
    
    for(ddt in 1:dt){
      dn = fn(A, r, N)      
      N = N+dn/dt
      N = pmax(N, 0)
    }      
    mat[t,] = N
  }
  
  #print( tail(mat) )
  return(mat)
}

calc_aburptness = function(before, after,ab_th=0){
  
  b_com = colMeans(before[-5:0+nrow(before),])
  a_com = colMeans(after[-5:0+nrow(after),])
  
  com = rbind(b_com, a_com)
  com[com<=0] = 0
  com_bi = ifelse(com>ab_th, 1, 0)
  
  b = vegan::vegdist(com, method="bray")
  j = vegan::vegdist(com_bi, method="jaccard", binary=TRUE)
  #turn_over = rank(b_com)-rank(a_com)
  c(BC=b, Jaccard=j,#rank_turn_over=turn_over,
    ab_th=ab_th)
}

conv_to_graph = function(A, weight=TRUE){
  diag(A) = 0
  Aw = A
  A[A!=0] = 1
  g = graph_from_adjacency_matrix(A)
  V(g)$name = V(g)
  if(weight){
    E(g)$weight = t(Aw)[which(t(Aw)!=0)]
    return(g)  
  }else{
    return(g)  
  }
  
}

calc_centrality = function(A){
  
  diag(A) = 0
  g = conv_to_graph(A)
  df = data.frame(
    id = as.integer(V(g)),
    degree_out = degree(g, mode="out", normalized = TRUE),
    degree_out_w = colSums(A)/(nrow(A)-1),
    degree_out_abs = colSums(abs(A))/(nrow(A)-1),
    degree_out_p = colSums(ifelse(A>0,A,0))/(nrow(A)-1),
    degree_out_n = colSums(ifelse(A<0,A,0))/(nrow(A)-1),
    degree_in = degree(g, mode="in", normalized = TRUE),
    betweenness = betweenness(g, weights = NA, normalized = TRUE),
    closenes = closeness(g, weights = NA, normalized = TRUE),
    pagerank = page_rank(g, weights = NA)$vector
  )
  return(df)
}

# calc_network_topology = function(A){
#   
#   g = conv_to_graph(A)
#   
#   df = data.frame(
#     edge_density = edge_density(g),
#     orderability =O_FUNC(g),
#     tree = T_FUNC(g)
#   )
#   return(df)
# }

plt_dynamics = function(mat,color=NULL, title=""){
  plot(0, type="n", xlim=c(1,nrow(mat)), ylim=range(mat), main=title)
  if(!is.null(color[1])){
    for(i in 1:ncol(mat)) lines(mat[,i], col=color[i])  
  }else{
    for(i in 1:ncol(mat)) lines(mat[,i], col=i)  
  }
}




silhouette_k <- function(D, 
                         k_max = 10, 
                         method = "ward.D2") {

  require(cluster)  
  ## dist オブジェクトに統一
  if (!inherits(D, "dist")) {
    D <- as.dist(D)
  }
  
  ## 階層クラスタリング
  hc <- hclust(D, method = method)
  
  ks <- 2:k_max
  sil_mean <- numeric(length(ks))
  
  for (i in seq_along(ks)) {
    k <- ks[i]
    cl <- cutree(hc, k = k)
    
    sil <- silhouette(cl, D)
    sil_mean[i] <- mean(sil[, "sil_width"])
  }
  
  ## 最適 k
  k_opt <- ks[which.max(sil_mean)]
  
  return(list(
    k = ks,
    silhouette = sil_mean,
    k_opt = k_opt,
    hclust = hc
  ))
}

edge_asymmetry <- function(A){
  S <- nrow(A)
  I <- numeric(S)
  
  for(i in 1:S){
    for(j in 1:S){
      if(i != j){
        I[i] <- I[i] + abs(A[i, j] + A[j, i])
      }
    }
  }
  I
}

delta_lambda_max <- function(A){
  eig_full <- eigen(A, only.values = TRUE)$values
  lambda_full <- max(Re(eig_full))
  
  S <- nrow(A)
  delta_lambda <- numeric(S)
  
  for(i in 1:S){
    A_i <- A[-i, -i]
    eig_i <- eigen(A_i, only.values = TRUE)$values
    lambda_i <- max(Re(eig_i))
    delta_lambda[i] <- lambda_i - lambda_full
  }
  
  delta_lambda
}

delta_modularity <- function(A){
  S <- nrow(A)
  
  # full network
  A0 <- abs(A)
  diag(A0) <- 0
  
  g0 <- graph_from_adjacency_matrix(
    A0, mode = "directed", weighted = TRUE
  )
  g0u <- as.undirected(g0, mode = "collapse", edge.attr.comb = "sum")
  
  cl0 <- cluster_louvain(g0u, weights = E(g0u)$weight)
  Q0 <- modularity(cl0)
  
  deltaQ <- numeric(S)
  
  for(i in 1:S){
    A_i <- abs(A[-i, -i])
    diag(A_i) <- 0
    
    g_i <- graph_from_adjacency_matrix(
      A_i, mode = "directed", weighted = TRUE
    )
    g_iu <- as.undirected(g_i, mode = "collapse", edge.attr.comb = "sum")
    
    cl_i <- cluster_louvain(g_iu, weights = E(g_iu)$weight)
    Q_i <- modularity(cl_i)
    
    deltaQ[i] <- Q_i - Q0
  }
  
  deltaQ
}

compute_multistability_metrics <- function(A){
  data.frame(
    species = 1:nrow(A),
    delta_lambda_max = delta_lambda_max(A),
    delta_modularity = delta_modularity(A),
    edge_asymmetry = edge_asymmetry(A)
  )
}


