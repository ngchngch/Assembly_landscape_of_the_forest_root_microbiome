
source("packages/fun_gLV_simulation_revise_260514.R")
source("packages/01_1_function.R")
# %% loading library and functions
library(tidyverse)
library(ggplot2)
library(doParallel)
library(igraph)
library(deSolve)
library(vegan)
library(cluster)
library(patchwork)

save.dir <- "05_01_gLV_simulation_Nsp80_nonoise_260514"
dir.create(save.dir)
# Parameters
total_sample <- 200
rep =200
time=100
Nsp=c(80)
connectance_list=c(0.2,0.3,0.4)
A_self = -1
#r0 = 1




#-- parallel
#nthread = detectCores()
cluster = makeCluster(n.core)
registerDoParallel(cluster)  
#on.exit(stopCluster(cluster)) 


# Main --------------------------------------------------------------------

dir.create(sprintf("%s/seed%s",save.dir,seed), showWarnings = FALSE)
set.seed(seed)

connectance=connectance_list[connect]
# %% Make interaction network
set = expand.grid(Nsp, connectance)


rlist = foreach(n=1:nrow(set))%do%{
  N=set[n,1]
  r <- runif(N)
  return(r)
}|> set_names(apply(set,1,paste, collapse="_"))
saveRDS(rlist, sprintf("%s/seed%s/rparameter_seed%s_Nsp%s.rds",save.dir,seed,connectance,Nsp))


Alist = foreach(n=1:nrow(set), .packages = c("igraph"))%do%{
  inf_ab <- TRUE
  count <- 1
  while(inf_ab){
    A = makeA(set[n,1], connectance = (set[n,2]), power=1, mu=-0.1,
              sigma = 0.2#ceiling(abs(A_self)/sqrt(set[n,1]*set[n,2])*100)/100
    )
    diag(A) = A_self
    eigA <- eigen(A, only.values = TRUE)$values
    # 最大実部固有値
    cat("Max Re eigenvalue:", max(Re(eigA)), "\n")
    test_mat = gen_dyn(30, dt=1000, A, rlist[[1]], runif(Nsp), type="nonoise")
    
    if(!any(is.na(test_mat))&max(test_mat)<20){
      inf_ab <- FALSE
    }
    
    if(count>10000){
      inf_ab <- FALSE
    }
    count <- count +1
  }
  return(A)
}|> set_names(apply(set,1,paste, collapse="_"))

saveRDS(Alist, sprintf("%s/seed%s/Amatrix_seed%s.rds",save.dir,seed,connectance))


# Plot
# par(mfrow=c(2,3))
# map(Alist, ~plot( conv_to_graph(.x, weight=FALSE), edge.arrow.size=0.2 ))
# map(Alist, hist)

# %% Run simulation

set = expand.grid(Nsp, connectance, 1:rep)|>
  set_names(c("N", "C", "rep"))|>
  as.matrix()


dynamics = foreach(n=1:nrow(set),
                   .packages = c("igraph", "tidyverse", "deSolve", "doParallel"))%dopar%{
                     
                     
                     s = paste(set[n,1:2], collapse="_")
                     A=Alist[[s]]
                     N=set[n,1]
                     r_set <- rlist[[s]]
                     
                     mat = gen_dyn(time=time, dt=1000, A=A, r=r_set, N0=runif(N), type="nonoise") #ode(y = N0(N), times = seq(1,500), func = lv_fun, parms = list(A=A, r=runif(N)))[,-1]
                     #plt_dynamics(mat)
                     
                     #N_perturbate = mat[time*0.9,]
                     #mat_nonperturbate = gen_dyn(time=time, dt=1000, A=A, r=r_set,
                     #                            N0=N_perturbate, type="nonoise")
                     
                     pert_list2 = pert_list3 = diff = c()
                     for(i in 1:N){
                       init <- c(runif(N-1),0.5)[order(c(c(1:N)[-i],i))]
                       mat_woi = gen_dyn(time=time, dt=1000, A=A[-i,-i], r=r_set[-i],
                                         N0=init[-i], type="nonoise")
                       mat_wi = gen_dyn(time, dt=1000, A, r_set,
                                        init, type="nonoise")
                       #plt_dynamics(mat_wi)
                       #mat_perturbate = gen_dyn(time, dt=1000, A[-i,-i], r_set[-i],
                       #                         N_perturbate[-i], type="gLV")
                       #plt_dynamics(rbind( mat[1:time*0.9,-i], mat_perturbate[,-i]))
                       #abruptness = calc_aburptness(before=mat[,-i], after=mat_perturbate[,-i])
                       
                       # diff = rbind(diff, c(set[n,], id=i, 
                       #                      abruptness))
                       
                       s2 = paste(c(set[n,1:3],"sp",i) ,collapse="_")
                       pert_list3[[s2]] = mat_wi
                       pert_list2[[s2]] = cbind(mat_woi,0)[,order(c(c(1:N)[-i],i))]
                       
                     }
                     
                     return( list(dynamics=mat, #pertabation=pert_list,
                                  dyn_woi=pert_list2,
                                  #non_pertabation=mat_nonperturbate,
                                  dyn_wi=pert_list3) )
                     
                   }|> set_names(apply(set,1,paste, collapse="_"))

picks <- sample(length(dynamics),10)
pdf(sprintf("%s/seed%s/figs_dynamics_conn%s.pdf",save.dir,seed,connectance), onefile = TRUE, w=6, h=6)
map2(dynamics[picks], names(dynamics)[picks], ~.x[[1]]|>plt_dynamics(title=.y) )
dev.off()

# pdf(sprintf("%s/seed%s/figs_dynamics_w_removing_conn%s.pdf",save.dir,seed,connectance),
#     onefile = TRUE, w=6, h=6)
# map2(dynamics[picks], names(dynamics)[picks], ~.x[[2]][[9]]|>plt_dynamics(title=.y) )
# dev.off()


##create community matrix
non_pertabate = map(dynamics, ~.x[[1]])
#pertabate = map(dynamics, ~.x[[2]])
woi = map(dynamics, ~.x[[2]])
wi = map(dynamics, ~.x[[3]])

sim_set <- matrix(unlist(strsplit(names(non_pertabate),"_")),ncol=3,byrow=TRUE)

mat0 <- NULL
for(i in 1:nrow(sim_set[sim_set[,2]==as.character(connectance),])){
  mat0 <- rbind(mat0,non_pertabate[[which(sim_set[,2]==as.character(connectance))[i]]][sample(1:time,total_sample/rep),])
}

d_bray <- vegdist(mat0,method="bray")

qt <- quantile(mat0, probs = c(0.1,0.9))
bin_th <- seq(qt[1],qt[2], length.out = 30)

cr <- c()
for(i in 1:length(bin_th)){
  cat("\r", i, "/", length(bin_th))
  bth <- bin_th[i]
  mat2 <- mat0
  mat2[mat2<=bth] <- 0
  mat2[mat2>0] <- 1
  d_jac <- vegdist(mat2,method="jaccard")
  cr[i] <- cor(d_bray, d_jac,method = "kendall")
}

max_cr <- which.max(cr)
best_bin_th <- bin_th[max_cr[order(max_cr)[1]]]

mat2 <- mat0
mat2[mat2<=best_bin_th] <- 0
mat2[mat2>0] <- 1

colMeans(mat2)
saveRDS(non_pertabate, sprintf("%s/seed%s/without_perturbation_dyn_conn%s.rds",save.dir,seed,connectance))
#saveRDS(pertabate, sprintf("%s/seed%s/with_perturbation_dyn_conn%s.rds",save.dir,seed,connectance))
saveRDS(woi, sprintf("%s/seed%s/dyn_woi_conn%s.rds",save.dir,seed,connectance))
saveRDS(wi, sprintf("%s/seed%s/dyn_wi_conn%s.rds",save.dir,seed,connectance))

####
kness = foreach(n=1:nrow(set),
                .packages = c("igraph", "tidyverse", "deSolve", "doParallel"))%do%{
                  
                  s = paste(set[n,1:2], collapse="_")
                  A=Alist[[s]]
                  cur_N=set[n,1]
                  cur_connectance=set[n,2]
                  cur_rep=set[n,3]
                  r_set <- rlist[[s]]
                  
                  dyn_n <- dynamics[[ paste(set[n,], collapse="_") ]]
                  mat = dyn_n$dynamics
                  #plt_dynamics(mat)
                  
                  diff <- NULL
                  com_wis <- NULL; com_wois <- NULL; b_com <- NULL
                  for(i in 1:N){
                    
                    #mat_perturbate = dyn_n$pertabation[[paste(cur_N,cur_connectance,cur_rep,"sp",i,sep="_")]]
                    #mat_non_perturbate = dyn_n$non_pertabation#[[1]]
                    cwoi = dyn_n$dyn_woi[[paste(cur_N,cur_connectance,cur_rep,"sp",i,sep="_")]]
                    cwi = dyn_n$dyn_wi[[paste(cur_N,cur_connectance,cur_rep,"sp",i,sep="_")]]
                    #plt_dynamics(rbind( mat[1:time*0.9,-i], mat_perturbate[,-i]))
                    #b_com_woi <- mat[,-i]
                    #a_com_woi <- mat_perturbate[((time*0.9)+1):nrow(mat_perturbate),-i]
                    #a_com_wi <- mat_non_perturbate[((time*0.9)+1):nrow(mat_non_perturbate),-i]
                    com_woi <- cwoi[,-i]
                    com_wi <- cwi[,-i]
                    com_wi_woi <- cbind(cwi[,-i],matrix(0,nrow=nrow(cwi),ncol=1))[,order(c(c(1:ncol(cwi))[-i],i))]
                    # abruptness = calc_aburptness(before=b_com_woi,
                    #                              after=a_com_woi,
                    #                              ab_th = best_bin_th)
                    
                    # abruptness_nonpert = calc_aburptness(before=b_com_woi,
                    #                                      after=a_com_wi,
                    #                                      ab_th = best_bin_th)
                    abruptness_diffstart = calc_aburptness(before=com_woi,
                                                           after=com_wi,
                                                           ab_th = best_bin_th)
                    
                    #names(abruptness_short) <- paste0(names(abruptness_short),"_short")
                    #names(abruptness_nonpert) <- paste0(names(abruptness_nonpert),"_nonpert")
                    names(abruptness_diffstart) <- paste0(names(abruptness_diffstart),"_diffstart")
                    
                    diff = rbind(diff, c(set[n,], id=i, 
                                         abruptness_diffstart[1:2]))
                    com_wois = rbind(com_wois,
                                     as.data.frame(matrix(c(set[n,],id=i,colMeans(com_woi[-5:0+nrow(com_woi),])),nrow=1)))
                    com_wis = rbind(com_wis,
                                    as.data.frame(matrix(c(set[n,],id=i,colMeans(com_wi[-5:0+nrow(com_woi),])),nrow=1)))
                    # b_com = rbind(b_com,
                    #               as.data.frame(matrix(c(set[n,],id=i,colMeans(b_com_woi[-5:0+nrow(b_com_woi),])),nrow=1)))
                  }
                  
                  keystoneness = 
                    data.frame(diff,r=r_set, abundance=mat[time*0.9,], ab_rank=rank(mat[time*0.9,])) |>
                    left_join(
                      calc_centrality(A),
                      by="id"
                    )
                  # 
                  return( list(keystoneness=keystoneness,
                               comm_mean_wi=com_wis,
                               comm_mean_woi=com_wois) )
                  
                }|> set_names(apply(set,1,paste, collapse="_"))

############# Analyze keystoneness
keystonness = map_df(kness, ~.x[[1]])

saveRDS(keystonness, sprintf("%s/seed%s/keystonness_vs_index_conn%s.rds",save.dir,seed,connectance))


comm_wi = map_df(kness, ~.x[[2]])
comm_woi = map_df(kness, ~.x[[3]])

saveRDS(comm_wi, sprintf("%s/seed%s/comm_mean_with_Spi_conn%s.rds",
                         save.dir,seed,connectance))
saveRDS(comm_woi, sprintf("%s/seed%s/comm_mean_without_Spi_conn%s.rds",
                          save.dir,seed,connectance))

md_wi <- c(); md_wi_woi <- c();mdb_l2r <- c();mj_wi <- c(); mj_wi_woi <- c(); mj_l2r <- c()
k_bray <- c()
k_jac <- c()
sil_val_bray <- c()
sil_val_jac <- c()
shannon_bray <- c()
shannon_jac <- c()
kn_wi <- c();kn_woi <- c();kn_wi_jac <- c();kn_woi_jac <- c()
for(id in 1:Nsp){
  b = vegan::vegdist(comm_wi[comm_wi[,4]==id,-c(1:4)], method="bray")
  b[is.na(b)] <- 0
  abcomm <- rbind(comm_wi[comm_wi[,4]==id,-c(1:4)],
                  comm_woi[comm_woi[,4]==id,-c(1:4)])
  binmat <- ifelse(abcomm>best_bin_th, 1, 0)
  b = vegan::vegdist(abcomm, method="bray")
  b[is.na(b)] <- 0
  b_wi <- as.matrix(b)[1:sum(comm_wi[,4]==id),
                       1:sum(comm_wi[,4]==id)]
  md_wi[id] <- mean(b_wi[upper.tri(b_wi)])
  md_wi_woi[id] <- mean(as.matrix(b)[1:sum(comm_wi[,4]==id),
                                     sum(comm_wi[,4]==id)+1:sum(comm_woi[,4]==id)])
  
  mdb_l2r[id] <- log2(md_wi_woi[id]/md_wi[id])
  
  sil_bray <- silhouette_k(D=b,k_max=rep-1)
  k_bray[id] <- sil_bray$k_opt
  cluster_wi <- cutree(sil_bray$hc, k=k_bray[id])[1:sum(comm_wi[,4]==id)]
  pcwi <- table(cluster_wi)/length(cluster_wi)
  kn_wi[id] <- length(pcwi)
  cluster_woi <- cutree(sil_bray$hc, k=k_bray[id])[1:sum(comm_wi[,4]==id)+sum(comm_woi[,4]==id)]
  pcwoi <- table(cluster_woi)/length(cluster_woi)
  kn_woi[id] <- length(pcwoi)
  
  shannon_bray[id] <- (-sum(pcwi*log(pcwi)))-(-sum(pcwoi*log(pcwoi)))
  sil_val_bray[id] <- sil_bray$silhouette[sil_bray$k==k_bray[id]]
  
  j = vegan::vegdist(binmat, method="jaccard", binary=TRUE)
  j[is.na(j)] <- 0
  j_wi <- as.matrix(j)[1:sum(comm_wi[,4]==id),
                       1:sum(comm_wi[,4]==id)]
  mj_wi[id] <- mean(j_wi[upper.tri(j_wi)])
  mj_wi_woi[id] <- mean(as.matrix(j)[1:sum(comm_wi[,4]==id),
                                     sum(comm_wi[,4]==id)+1:sum(comm_woi[,4]==id)])
  mj_l2r[id] <- log2(mj_wi_woi[id]/mj_wi[id])
  
  sil_jac <- silhouette_k(D=j,k_max=rep-1)
  k_jac[id] <- sil_jac$k_opt
  cluster_wi_jac <- cutree(sil_jac$hc, k=k_jac[id])[1:sum(comm_wi[,4]==id)]
  pcwi_jac <- table(cluster_wi_jac)/length(cluster_wi_jac)
  kn_wi_jac[id] <- length(pcwi_jac)
  cluster_woi_jac <- cutree(sil_jac$hc, k=k_jac[id])[1:sum(comm_wi[,4]==id)+sum(comm_woi[,4]==id)]
  pcwoi_jac <- table(cluster_woi_jac)/length(cluster_woi_jac)
  kn_woi_jac[id] <- length(pcwoi_jac)
  
  shannon_jac[id] <- (-sum(pcwi_jac*log(pcwi_jac)))-(-sum(pcwoi_jac*log(pcwoi_jac)))
  sil_val_jac[id] <- sil_jac$silhouette[sil_jac$k==k_jac[id]]
}

init_depend_BC_wi <- foreach(i = 1:Nsp, .combine="c")%do% {
  mean(vegdist(comm_wi[keystonness$id == i,-c(1:4)],method="bray"))
}

init_depend_BC_woi <- foreach(i = 1:Nsp, .combine="c")%do% {
  mean(vegdist(comm_woi[keystonness$id == i,-c(1:4)],method="bray"))
}

init_depend_JC_wi <- foreach(i = 1:Nsp, .combine="c")%do% {
  mean(vegdist(ifelse(comm_wi[keystonness$id == i,-c(1:4)]>best_bin_th,1,0),
               method="jaccard", binary=TRUE))
}

init_depend_JC_woi <- foreach(i = 1:Nsp, .combine="c")%do% {
  mean(vegdist(ifelse(comm_woi[keystonness$id == i,-c(1:4)]>best_bin_th,1,0),
               method="jaccard", binary=TRUE))
}

mAB <- foreach(i = 1:Nsp, .combine="c")%do% {
  mean(keystonness[keystonness$id == i,"abundance"],na.rm=TRUE)
}

totalAB <- foreach(rp = 1:rep, .combine="c")%do% {
  sum(keystonness[keystonness$rep == rp,"abundance"],na.rm=TRUE)
}

# mBC <- foreach(i = 1:Nsp, .combine="c")%do% {
#   if(sum(keystonness[,"abundance"]>0 & keystonness$id == i)>0){
#     return(mean(keystonness[keystonness[,"abundance"]>0 & keystonness$id == i,"BC"],na.rm=TRUE))
#   }else{
#     return(NA)
#   }
# }

# mBCratio <- foreach(i = 1:Nsp, .combine="c")%do% {
#   if(sum(keystonness[,"abundance"]>0 & keystonness$id == i)>0){
#     k_pick <- keystonness[keystonness[,"abundance"]>0 & keystonness$id == i,]
#     return(mean(k_pick$BC/k_pick$BC_nonpert,na.rm=TRUE))
#   }else{
#     return(NA)
#   }
# }

mBC2 <- foreach(i = 1:Nsp, .combine="c")%do% {
  return(mean(keystonness[keystonness$id == i,
                          "BC_diffstart"],na.rm=TRUE))
}

mBC_l2r_perAB <- foreach(i = 1:Nsp, .combine="c")%do% {
  return(mdb_l2r[i]*mean(1-keystonness[keystonness$id == i,
                                       "abundance"]/totalAB))
}

# mJC <- foreach(i = 1:Nsp, .combine="c")%do% {
#   if(sum(keystonness[,"abundance"]>0 & keystonness$id == i)>0){
#     return(mean(keystonness[keystonness[,"abundance"]>0 & keystonness$id == i,"Jaccard"],na.rm=TRUE))
#   }else{
#     return(NA)
#   }
# }

mJC2 <- foreach(i = 1:Nsp, .combine="c")%do% {
  return(mean(keystonness[keystonness$id == i,
                          "Jaccard_diffstart"],na.rm=TRUE))
}

mJC_l2r_perAB <- foreach(i = 1:Nsp, .combine="c")%do% {
  return(mj_l2r[i]*mean(1-keystonness[keystonness$id == i,
                                      "abundance"]/totalAB))
}

# mBCperAB <- foreach(i = 1:Nsp, .combine="c")%do% {
#   if(sum(keystonness[,"abundance"]>0 & keystonness$id == i)>0){
#     k_pick <- keystonness[keystonness[,"abundance"]>0 & keystonness$id == i,]
#     totalAB_pick <- totalAB[k_pick$rep]
#     return(mean(k_pick[,"BC"]*(1-k_pick[,"abundance"]/totalAB_pick),na.rm=TRUE))
#   }else{
#     return(NA)
#   }
# }

# rmBCperAB <- foreach(i = 1:Nsp, .combine="c")%do% {
#   if(sum(keystonness[,"abundance"]>0 & keystonness$id == i)>0){
#     k_pick <- keystonness[keystonness[,"abundance"]>0 & keystonness$id == i,]
#     totalAB_pick <- totalAB[k_pick$rep]
#     return(mean((k_pick[,"BC"]/k_pick$BC_nonpert)*(1-k_pick[,"abundance"]/totalAB_pick),na.rm=TRUE))
#   }else{
#     return(NA)
#   }
# }

# mJCperAB <- foreach(i = 1:Nsp, .combine="c")%do% {
#   if(sum(keystonness[,"abundance"]>0 & keystonness$id == i)>0){
#     k_pick <- keystonness[keystonness[,"abundance"]>0 & keystonness$id == i,]
#     totalAB_pick <- totalAB[k_pick$rep]
#     return(mean(k_pick[,"Jaccard"]*(1-k_pick[,"abundance"]/totalAB_pick),na.rm=TRUE))
#   }else{
#     return(NA)
#   }
# }

multi_ind <- compute_multistability_metrics(Alist[[paste(set[1,1:2], collapse="_")]])

key_parm <- unique(keystonness[,-c(3,8:9)])

key_res <- data.frame(connectance=connectance,
                      init_depend_bray_wi=init_depend_BC_wi,
                      init_depend_bray_woi=init_depend_BC_woi,
                      l2r_id_BC=log2(init_depend_BC_wi/init_depend_BC_woi),
                      ncls_BC_wi=kn_wi,
                      ncls_BC_woi=kn_woi,
                      shannon_cls_bray=shannon_bray,
                      sil_bray=sil_val_bray,
                      init_depend_JC_wi=init_depend_JC_wi,
                      init_depend_JC_woi=init_depend_JC_woi,
                      l2r_id_JC=log2(init_depend_JC_wi/init_depend_JC_woi),
                      ncls_JC_wi=kn_wi_jac,
                      ncls_JC_woi=kn_woi_jac,
                      shannon_cls_jac=shannon_jac,
                      sil_jac=sil_val_jac,
                      #mBC=mBC,
                      mBC_diffstart=mBC2,
                      mBC_ds_l2r=mdb_l2r,
                      mBC_ds_l2r_perAB=mBC_l2r_perAB,
                      #mBCperAB=mBCperAB,
                      #mBCratio=mBCratio,
                      #mBCratio_perAB=rmBCperAB,
                      mJC_diffstart=mJC2,
                      #mJC=mJC,
                      #mJCperAB=mJCperAB,
                      mJC_ds_l2r=mj_l2r,
                      mJC_ds_l2r_perAB=mJC_l2r_perAB,
                      mAB=mAB,
                      multi_ind,
                      key_parm
)

saveRDS(key_res, sprintf("%s/seed%s/keystonness_summary_conn%s.rds",save.dir,seed,connectance))


g1 <- key_res|>
  pivot_longer(c(mAB,r,degree_out,
                 degree_out_abs, betweenness,   closenes,   pagerank,
                 delta_lambda_max,edge_asymmetry))|>
  ggplot(aes(x=value, y=mBC_ds_l2r))+
  geom_point(color="darkblue")+
  ggh4x::facet_grid2(N+C~name,scales="free", independent="all")+
  theme(aspect.ratio=1)
g12 <- key_res|>
  pivot_longer(c(mAB,r,degree_out,
                 degree_out_abs, betweenness,   closenes,   pagerank,
                 delta_lambda_max,edge_asymmetry))|>
  ggplot(aes(x=value, y=mBC_ds_l2r_perAB))+
  geom_point(color="blue")+
  ggh4x::facet_grid2(N+C~name,scales="free", independent="all")+
  theme(aspect.ratio=1)
g111 <- key_res|>
  pivot_longer(c(mAB,r,degree_out,
                 degree_out_abs, betweenness,   closenes,   pagerank,
                 delta_lambda_max,edge_asymmetry))|>
  ggplot(aes(x=value, y=mBC_diffstart))+
  geom_point(color="darkblue")+
  ggh4x::facet_grid2(N+C~name,scales="free", independent="all")+
  theme(aspect.ratio=1)
g121 <- key_res|>
  pivot_longer(c(mAB,r,degree_out,
                 degree_out_abs, betweenness,   closenes,   pagerank,
                 delta_lambda_max,edge_asymmetry))|>
  ggplot(aes(x=value, y=mJC_diffstart))+
  geom_point(color="purple")+
  ggh4x::facet_grid2(N+C~name,scales="free", independent="all")+
  theme(aspect.ratio=1)

g13 <- key_res|>
  pivot_longer(c(mAB,r,degree_out,
                 degree_out_abs, betweenness,   closenes,   pagerank,
                 delta_lambda_max,edge_asymmetry))|>
  ggplot(aes(x=value, y=mJC_ds_l2r))+
  geom_point(color="purple")+
  ggh4x::facet_grid2(N+C~name,scales="free", independent="all")+
  theme(aspect.ratio=1)
g14 <- key_res|>
  pivot_longer(c(mAB,r,degree_out,
                 degree_out_abs, betweenness,   closenes,   pagerank,
                 delta_lambda_max,edge_asymmetry))|>
  ggplot(aes(x=value, y=mJC_ds_l2r_perAB))+
  geom_point(color="darkorchid4")+
  ggh4x::facet_grid2(N+C~name,scales="free", independent="all")+
  theme(aspect.ratio=1)

g2 <- key_res|>
  pivot_longer(c(mAB,r,degree_out,
                 degree_out_abs, betweenness,closenes,pagerank,
                 delta_lambda_max,delta_modularity))|>
  ggplot(aes(x=value, y=shannon_cls_bray))+
  geom_point(aes(fill=sil_bray),shape=21)+
  ggh4x::facet_grid2(N+C~name,scales="free", independent="all")+
  scale_fill_viridis_c()+
  theme_bw()+
  theme(aspect.ratio=1)

g3 <- key_res|>
  pivot_longer(c(mAB,r,degree_out,
                 degree_out_abs, betweenness,closenes,pagerank,
                 delta_lambda_max,delta_modularity))|>
  ggplot(aes(x=value, y=shannon_cls_jac,fill=sil_jac))+
  geom_point(shape=21)+
  ggh4x::facet_grid2(N+C~name,scales="free", independent="all")+
  scale_fill_viridis_c()+
  theme_bw()+
  theme(aspect.ratio=1)

g22 <- key_res|>
  pivot_longer(c(mAB,r,degree_out,
                 degree_out_abs, betweenness,closenes,pagerank,
                 delta_lambda_max,delta_modularity))|>
  ggplot(aes(x=value, y=l2r_id_BC))+
  geom_point(aes(fill=sil_bray),shape=21)+
  ggh4x::facet_grid2(N+C~name,scales="free", independent="all")+
  scale_fill_viridis_c()+
  theme_bw()+
  theme(aspect.ratio=1)

g_merge <- g1+g12+g111+plot_layout(ncol=1)#+g13+g14+g2+g3+plot_layout(ncol=1)

ggsave(sprintf("%s/seed%s/keystonness_vs_commchange_bray_conn%s.png",save.dir,seed,connectance),
       g_merge,
       width=15, height=9)


g_merge <- g121+g13+g14+plot_layout(ncol=1)#+g2+g3+plot_layout(ncol=1)

ggsave(sprintf("%s/seed%s/keystonness_vs_commchange_jac_conn%s.png",save.dir,seed,connectance),
       g_merge,
       width=15, height=6)

g_merge <- g2+g22+g3+plot_layout(ncol=1)

ggsave(sprintf("%s/seed%s/keystonness_vs_commchange_even_conn%s.png",save.dir,seed,connectance),
       g_merge,
       width=15, height=10)

saveRDS(list(ab=mat0, bin_mat=mat2, bin_th=c(th=best_bin_th,
                                             cor=cr[max_cr[order(max_cr)[1]]])),
        sprintf("%s/seed%s/community_matrix_conn%s.rds",save.dir,seed,connectance))  

