
source("packages/fun_gLV_simulation_251229.R")

# %% loading library and functions
library(tidyverse)
library(ggplot2)
library(doParallel)
library(igraph)
library(deSolve)
library(vegan)
library(cluster)
library(patchwork)
library(rELA)
library(philentropy)

save.dir <- "06_01_gLV_samplesize_validation_260108"
dir.create(save.dir)








#setwd("/aptmp/mikihito/02_ELA_v060_241113_2")
#n.core <- 6
# Parameters
total_sample_seq <- c(20,50,100,150,200,300,500,800,1200)[nsamp]#1:8
nrep=10#number of reputation in evaluate variance & distance
time=100
Nsp_seq=c(20,30,40,50)
connectance=0.2
A_self = -1
#r0 = 1

# Main --------------------------------------------------------------------

dir.create(sprintf("%s/seed%s",save.dir,seed), showWarnings = FALSE)
set.seed(seed)


Nsp <- Nsp_seq[nsp] #nsp=1:4
# %% Make interaction network
set = expand.grid(Nsp, connectance)

if(nsamp==1){
  
Alist = foreach(n=1:nrow(set), .packages = c("igraph"))%do%{
  inf_ab <- TRUE
  count <- 1
  while(inf_ab){
    A = makeA(set[n,1], connectance = (set[n,2]), power=1, mu=-0.1,
              sigma = 0.3#ceiling(abs(A_self)/sqrt(set[n,1]*set[n,2])*100)/100
    )
    diag(A) = A_self
    eigA <- eigen(A, only.values = TRUE)$values
    # 最大実部固有値
    cat("Max Re eigenvalue:", max(Re(eigA)), "\n")
    test_mat = gen_dyn(50, dt=1000, A, runif(Nsp), runif(Nsp), type="gLV")
    
    if(!any(is.na(test_mat))&max(test_mat)<50){
      inf_ab <- FALSE
    }
   
    if(count>10000){
      inf_ab <- FALSE
    }
    count <- count +1
  }
  return(A)
}|> set_names(apply(set,1,paste, collapse="_"))

  
  saveRDS(Alist, sprintf("%s/seed%s/Amatrix_seed%s_Nsp%s.rds",save.dir,seed,connectance,Nsp))
  
  rlist = foreach(n=1:nrow(set))%do%{
    N=set[n,1]
    r <- runif(N, min=0, max=2)
    return(r)
  }|> set_names(apply(set,1,paste, collapse="_"))
  saveRDS(rlist, sprintf("%s/seed%s/rparameter_seed%s_Nsp%s.rds",save.dir,seed,connectance,Nsp))
  
}else{
  Alist <- readRDS( sprintf("%s/seed%s/Amatrix_seed%s_Nsp%s.rds",save.dir,seed,connectance,Nsp))
  rlist <- readRDS( sprintf("%s/seed%s/rparameter_seed%s_Nsp%s.rds",save.dir,seed,connectance,Nsp))
}

# Plot
# par(mfrow=c(2,3))
# map(Alist, ~plot( conv_to_graph(.x, weight=FALSE), edge.arrow.size=0.2 ))
# map(Alist, hist)


dist_bf <- c()
v_nss <- c()
nss <- NULL
for(nt in 1:length(total_sample_seq)){
  total_sample <- total_sample_seq[nt]
  dir.create(sprintf("%s/seed%s/samplesize_%s",save.dir,seed,total_sample), showWarnings = FALSE)
  
  rep = total_sample
  
  # %% Run simulation
  set = expand.grid(Nsp, connectance, 1:rep)|>
    set_names(c("N", "C", "rep"))|>
    as.matrix()
  
  each_nss <- c()
  bfreq <- list(NULL)
  for(rep_v in 1:nrep){
    
  cluster = makeCluster(n.core)
  registerDoParallel(cluster)  
  on.exit(stopCluster(cluster)) 
dynamics = foreach(n=1:nrow(set),
                   .packages = c("igraph", "tidyverse", "deSolve", "doParallel"))%dopar%{
  
  
  s = paste(set[n,1:2], collapse="_")
  A=Alist[[s]]
  N=set[n,1]
  r_set <- rlist[[s]]
  
   mat = gen_dyn(time, dt=1000, A, r_set, runif(N), type="gLV") #ode(y = N0(N), times = seq(1,500), func = lv_fun, parms = list(A=A, r=runif(N)))[,-1]
  
  return( list(dynamics=mat) )
  
}|> set_names(apply(set,1,paste, collapse="_"))

picks <- sample(length(dynamics),20)
pdf(sprintf("%s/seed%s/samplesize_%s/figs_dynamics_conn%s_rep%s.pdf",save.dir,seed,total_sample,connectance,rep_v), onefile = TRUE, w=6, h=6)
map2(dynamics[picks], names(dynamics)[picks], ~.x[[1]]|>plt_dynamics(title=.y) )
dev.off()

##create community matrix
non_pertabate = map(dynamics, ~.x[[1]])
saveRDS(non_pertabate, sprintf("%s/seed%s/samplesize_%s/without_perturbation_dyn_conn%s_rep%s.rds",save.dir,seed,total_sample,connectance,rep_v))

sim_set <- matrix(unlist(strsplit(names(non_pertabate),"_")),ncol=3,byrow=TRUE)

 
    cat("\nSample size:", total_sample, "rep:", rep_v, "\n")
    mat <- NULL
    for(i in sample(nrow(sim_set[sim_set[,2]==as.character(connectance),]),total_sample)){
      mat <- rbind(mat,non_pertabate[[which(sim_set[,2]==as.character(connectance))[i]]][sample(1:time,1),])
    }
    
     d_bray <- vegdist(mat,method="bray")
    
    qt <- quantile(mat, probs = c(0.1,0.9))
    bin_th <- seq(qt[1],qt[2], length.out = 30)
    
    cr <- c()
    for(i in 1:length(bin_th)){
      cat("\r", i, "/", length(bin_th))
      bth <- bin_th[i]
      mat2 <- mat
      mat2[mat2<=bth] <- 0
      mat2[mat2>0] <- 1
      d_jac <- vegdist(mat2,method="jaccard")
      cr[i] <- cor(d_bray, d_jac,method = "kendall")
    }
    
    max_cr <- which.max(cr)
    best_bin_th <- bin_th[max_cr[order(max_cr)[1]]]
    
    mat2 <- mat
    mat2[mat2<=best_bin_th] <- 0
    mat2[mat2>0] <- 1
    
    colnames(mat2) <- paste0("Sp",1:ncol(mat2))
    saveRDS(list(ab=mat, bin_mat=mat2, bin_th=c(th=best_bin_th,
                                                cor=cr[max_cr[order(max_cr)[1]]])),
            sprintf("%s/seed%s/samplesize_%s/community_matrix_conn%s_Nsp%s_nsamp%s_rep%s.rds",save.dir,seed,total_sample,connectance,Nsp,total_sample,rep_v))  
    ####
    #ELA parameters
    qth <- 10^-5 # do not change!!
    SS.itr <- 150000
    
    sa <- runSA(ocmat=as.matrix(mat2),
                qth=qth, rep=1280, threads=n.core)
    saveRDS(sa, file=sprintf("%s/seed%s/samplesize_%s/runSA_conn%s_Nsp%s_nsamp%s_rep%s.rds",save.dir,seed,total_sample,connectance,Nsp,total_sample,rep_v))
    
    param <- sa2params(sa)
    
    hge <- param[[4]]
    je <- param[[2]]
    
    ss <- SSestimate(hge, je)
    
    mst <- unique(ss[,-ncol(ss)])
    saveRDS(mst, file=sprintf("%s/seed%s/samplesize_%s/SScomp_conn%s_Nsp%s_nsamp%s_rep%s.rds",
                              save.dir,seed,total_sample,connectance,Nsp,total_sample,rep_v))
    
    nss <- rbind(nss,data.frame(nsample=total_sample,rep=rep_v,nbasin=nrow(mst)))
    each_nss[rep_v] <- nrow(mst)
    
    tb_bfreq <- table(apply(ss[,-ncol(ss)], 1, paste, collapse=''))
    bfreq[[rep_v]] <- as.numeric(tb_bfreq);names(bfreq[[rep_v]]) <- names(tb_bfreq)
  }
  v_nss[nt] <- sd(each_nss)

  all_basin <- unique(names(unlist(bfreq)))
  mat_bfreq <- do.call(rbind,map(bfreq, ~.x[all_basin]))/SS.itr
  colnames(mat_bfreq) <- all_basin
  mat_bfreq[is.na(mat_bfreq)] <- 0
  dist_bf[nt] <- mean(distance(mat_bfreq, method = "jensen-shannon"))
  cat("\nDistance between basin frequency distributions:", dist_bf[nt], "\n")
  saveRDS(cbind(nss,sd_nbasin=v_nss[nt],JSdist_bfreq=dist_bf[nt]),
          sprintf("%s/seed%s/samplesize_%s/summary_ELA_result_conn%s_Nsp%s.rds",
                  save.dir,seed,total_sample,connectance,Nsp))
}

# res <- merge(nss,data.frame(nsample=total_sample_seq,
#                SD_nb=v_nss,dist_bfreq=dist_bf),by="nsample")
# write_csv(res,sprintf("%s/seed%s/summary_Number_of_Basin_conn%s_Nsp%s.csv",
#                     save.dir,seed,connectance,Nsp))
# 
# g <- ggplot(res, aes(x=nsample, y=nbasin))+
#   geom_point()+
#   #geom_line()+
#   theme_bw()+
#   theme(text = element_text(size=14),
#         aspect.ratio = 0.5)+
#   labs(x="Number of samples", y="Number of detected basins")+
#   scale_x_continuous(breaks = total_sample_seq)
# 
# ggsave(sprintf("%s/seed%s/fig_Number_of_Basin_conn%s_Nsp%s.pdf",
#                save.dir,seed,connectance,Nsp), g,
#        width=8, height=4)
# 
# g <- ggplot(unique(res[,-c(2,3)]), aes(x=nsample, y=SD_nb))+
#   geom_point()+
#   geom_line()+
#   theme_bw()+
#   theme(text = element_text(size=14),
#         aspect.ratio = 0.5)+
#   labs(x="Number of samples", y="Standard deviation\n(Number of detected basins)")+
#   scale_x_continuous(breaks = total_sample_seq)
# 
# ggsave(sprintf("%s/seed%s/fig_SD_Number_of_Basin_conn%s_Nsp%s.pdf",
#                save.dir,seed,connectance,Nsp), g,
#        width=8, height=4)
# 
# g <- ggplot(unique(res[,-c(2,3)]), aes(x=nsample, y=dist_bf))+
#   geom_point()+
#   geom_line()+
#   theme_bw()+
#   theme(text = element_text(size=14),
#         aspect.ratio = 0.5)+
#   labs(x="Number of samples", y="Mean dissimilarity of basin propotions\n(Jensen-Shannon distance))")+
#   scale_x_continuous(breaks = total_sample_seq)
# 
# ggsave(sprintf("%s/seed%s/fig_mean_dist_basin_proportion_conn%s_Nsp%s.pdf",
#                save.dir,seed,connectance,Nsp), g,
#        width=8, height=4)