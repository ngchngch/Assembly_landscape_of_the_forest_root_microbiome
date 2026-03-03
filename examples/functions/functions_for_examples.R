# %% loading library and functions
library(igraph)
library(Matrix)
library(codyn)

# Functions ----------------------------------------------------------------

N0 = \(sp=4)(runif(sp, 0.1, 0.8))

makeA = function(N=100, connectance=0.3, power=1, mu=0, sigma=0.1, output_g=FALSE){
  
  g = igraph::sample_gnp(N, p=connectance, directed = TRUE)
  # g = igraph::sample_pa(N, power = power,
  #                       start.graph = g,  m=1)
  
  A = as_adjacency_matrix(g)
  link_posi = (A!=0)
  A[which(link_posi)] = rnorm(sum(link_posi), mu, sigma)#*rlnorm(sum(link_posi), 0, 1)
  
  if(output_g) return(g)
  return(as.matrix(A))
}

lv_fun = \(A, r, N) N * (r + A %*% N)+rnorm(length(N), 0, 0.01)*rlnorm(length(N), 0, 2)
ms_ricker = \(N, r, A) exp( r*(1+A%*%N) ) * N

gen_dyn = function(time, dt=1000, A, r, N0, type="ricker"){
  
  mat = matrix(0, ncol=length(N0), nrow=time)
  mat[1,] = N = N0
  fn = ifelse(type=="ricker", ms_ricker, lv_fun)
  
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

plt_dynamics = function(mat, title=""){
  plot(0, type="n", xlim=c(1,nrow(mat)), ylim=range(mat), main=title)
  for(i in 1:ncol(mat)) lines(mat[,i], col=i)
}

#####

SSchange <- function(state,sa,env_grad,env_cat=NA,steps=4,eq_steps=TRUE,
                     start=NA,
                     RA_label=NA,
                     range=NA,SS.itr=20000,threads=1,reporting=TRUE){
  #env_grad <- enmatf[,1];env_cat <- sp_info[1,-c(1,ncol(sp_info))]#Betula
  #range <- enmatf[,1]
  cluster = makeCluster(threads)
  registerDoParallel(cluster)
  on.exit(stopCluster(cluster))
  ##===============================================##
  ##GradELA
  s2 <- proc.time()[3]
  if(reporting){cat('Start SSestimation with gradient factor\n')}
  
  if(is.na(RA_label[1])){
    RA_label <- seq(1,step,1)
  }
  
  if(eq_steps){
    
    if(length(range)<2){
      de=range
    }else{
      if(is.na(range[1])){mi <- min(env_grad)}else{mi <- range[1]}
      if(is.na(range[1])){ma <- max(env_grad)}else{ma <- range[2]}
      
      dm <- (ma - mi)/(steps - 1)
      de <- seq(mi, ma, dm) 
    }
    
  }else{
    if(length(range)<steps){steps <- length(range)}
    ur <- ifelse(range>start,range,NA)
    RA_label <- RA_label[!is.na(ur)]
    
    ur2 <- ur[!is.na(ur)]
    if(length(ur2)>0){
      de <- c(start,ur2[order(ur2)])  
      RA_label <- RA_label[order(ur2)]
    }else{
      de <- NA
    }
  }
  
  if(!is.na(de[1])){
    ssprop <- foreach(i = 1:length(de),.packages = c("rELA","tidyr","doParallel","vegan"),
                      .combine="c")%do%{
                        #i <- 2
                        if(!is.na(env_cat[1])){
                          ee <- c(de[i],as.numeric(env_cat))  
                        }else{
                          ee <- de[i]
                        }
                        
                        param <- sa2params(sa,ee)
                        
                        hge <- param[[4]]
                        je <- param[[2]]
                        
                        ss <- SSestimate_given(hge, je, state)
                        
                        mst <- ss[,-ncol(ss)]
                        rownames(mst) <- rownames(state)
                        
                        id <- apply(mst, 1, paste, collapse='')
                        ssid <- apply(unique(mst),1,bin2id)
                        
                        umst <- data.frame(env=ee[1],id2=id[names(ssid)],as.numeric(table(id)),unique(mst))
                        dimnames(umst) <- list(ssid,c("env","id2","freq",names(hge)))
                        
                        names(ssid) <- apply(unique(mst), 1, paste, collapse='')
                        prop <- as.numeric(table(id)[names(ssid)])/SS.itr
                        return(list(list(df=data.frame(env_grad=de[i],ssid=ssid,
                                                       h1=vegan::diversity(prop),
                                                       prop=prop),
                                         simulation=id,
                                         ss_structure=umst)))
                      }
    
    sstr1 <- foreach(l=1:length(ssprop),.combine = "rbind")%do%{
      return(ssprop[[l]]$ss_structure)
    }
    sstr2 <- sstr1
    
    sstr <- unique(sstr1[,-c(1,3)])
    rownames(sstr) <- sstr[,"id2"]
    sdf <- foreach(l=1:length(ssprop),.combine = "rbind")%do%{
      return(ssprop[[l]]$df)
    }
    ssim <- foreach(l=1:length(ssprop),.combine = "cbind")%do%{
      return(ssprop[[l]]$simulation)
    }
    
    minss <- unique(ssim[,1])
    d_land <- c()
    for(i in 2:ncol(ssim)){#i <- 2
      ssp <- NULL
      for(j in 1:length(minss)){
        ssp_tb <- table(ssim[which(ssim[,1] %in% minss[j]),i])
        ssp <-  rbind(ssp,data.frame(ss1=minss[j],ss2=names(ssp_tb),
                                     count=as.numeric(ssp_tb)))
      }
      ssp_d <-foreach(k = 1:nrow(ssp),.packages = "vegan",.combine = "c")%dopar%{#k <- 1
        res <- ssp[k,3]*vegdist(rbind(sstr[ssp[k,1],-1],
                                      sstr[ssp[k,2],-1]),
                                method="jaccard")[1]/(SS.itr)
        if(is.na(res)){res <- 0}
        return(res)
      }
      d_land[i-1] <- sum(ssp_d)
    }
    
    res <- spread(sdf,key=ssid,value=prop)
    
    res[is.na(res)] <- 0
    
    #delta evenness
    d_even <- c(res$h1[2:length(res$h1)]-res$h1[1])
    
    if(reporting){cat(sprintf("Elapsed time %.2f sec\n", proc.time()[3] - s2))}
    
    return(list(skip=FALSE,
                SStable=cbind(rownames(sstr2),sstr2),
                df=res,
                result=data.frame(ra=RA_label,
                                  d_land=d_land,
                                  d_even=d_even)))
    
  }else{
    return(list(skip=TRUE,
                SStable=NA,
                df=NA,
                result=data.frame(ra=NA,
                                  d_land=NA,
                                  d_even=NA)))
  }
  
}


Quant_landshift <- function(bin, ex_var=NULL, abundance_focal = ab_mat[,sp],
                           qt_seq = c(0.25,0.5,0.75),
                            n.core=1,itr_fit=16,qth=10^-5,SS.itr=150000){
  enmat_ns <- abundance_focal
  
  if(is.null(ex_var[1])){
    enmatf <- matrix(scale(abundance_focal),ncol=1)
  }else{
    enmatf <- cbind(matrix(scale(abundance_focal),ncol=1),
                    ex_var)
  }
 
  dimnames(bin) <- list(1:nrow(bin), paste0("sp",1:ncol(bin)))
  sa <- runSA(ocmat = bin,enmat=enmatf,
              qth=qth, rep=itr_fit, threads=n.core)
  
  hg <- sa2params(sa)[[4]]
  state <- foreach(i=1:SS.itr,.combine="rbind")%do%{
    st <- runif(length(hg), 0, 2) |> as.integer()
  }
  rownames(state) <- sprintf("Start_%05d",1:SS.itr)
  colnames(state) <- names(hg)
    
  ra_perc <- quantile(enmatf[,1],qt_seq)
  ran_perc <- quantile(enmat_ns,qt_seq)
  
  sprop <- SSchange(state=state,
                    sa=sa,
                    steps=length(qt_seq),
                    RA_label=paste0("perc",qt_seq*100),
                    env_cat=ex_var,
                    reporting = FALSE,
                    start=min(enmatf),
                    range=ra_perc,
                    eq_steps = FALSE,
                    SS.itr=SS.itr,threads=n.core)
  
  md_sprop <- cbind(id=sp,ab=ran_perc,
                    scale_ab=ra_perc,
                    scale_ab0=c(min(enmatf),rep(NA,length(qt_seq)-1)),
                    sd_ab=c(sd(enmat_ns),
                            rep(NA,length(qt_seq)-1)),
                    mean_ab=c(mean(enmat_ns),
                              rep(NA,length(qt_seq)-1)),
                    sprop[["result"]])
  return(md_sprop)
  
}

##############
binarize_M2SD <- function(relf2=relf2){
  grelf <- gather(cbind(ID=rownames(relf2),as.data.frame(relf2)),
                  key=taxa,value="abundance",-1) 
  grelf2 <- grelf[grelf$abundance!=0,]
  grelf2$logRA <- log(grelf2$abundance)
  
  grelf2$noise <- "Pres"
  th_df <- NULL
  for(i in 1:length(unique(grelf2$taxa))){#i <- 1
    show.progress(i,1:length(unique(grelf2$taxa)))
    grf <- grelf2[which(grelf2$taxa==unique(grelf2$taxa)[i]),] 
    ra <- grf$logRA
    
    grf[which(grf$logRA < mean(ra)-2*sd(ra)),"noise"] <- "low (< M-2SD)"
    grelf2[which(grelf2$taxa==unique(grelf2$taxa)[i]),"noise"] <- grf$noise
    
    th_df <- rbind(th_df,
                   data.frame(unique(grelf2$taxa)[i],mean(ra)-2*sd(ra)))
  }
  
  dimnames(th_df) <- list(th_df[,1],c("taxa","th_2SD"))
  
  ocmat <- relf2
  for(i in 1:ncol(relf2)){
    ocmat[,i] <- ifelse(ocmat[,i]<exp(th_df[colnames(ocmat)[i],2]),0,1)
  }
  return(ocmat)
}

prioritize_adonis2 <- function(bin_mat,
                               ab_mat,
                               ex_var=NULL,
                               nperm=2,
                               n.core=1,
                               dist_method="bray"){
  if(!all(colnames(bin_mat)==colnames(ab_mat)&rownames(bin_mat)==rownames(ab_mat))){
    stop("ab_mat and bin_mat must have same rownames")
  }
  
  require(vegan)
  
  if(is.null(ex_var[1])){
    R2 <- NULL
    for(tx in 1:ncol(bin_mat)){#tx <- 1
      
      exp_df <- data.frame(Occ=bin_mat[,tx])
      perm <- adonis2(data=exp_df,
                      ab_mat~.,
                      by = "margin",
                      permutations=nperm,
                      parallel = n.core,
                      method = dist_method)
      
      R2 <- rbind(R2,data.frame(taxa=colnames(bin_mat)[tx],
                                R2=perm$R2[ncol(exp_df)],raw_p=perm$`Pr(>F)`[ncol(exp_df)]))
    }
    
  }else{
    if(!all(rownames(bin_mat)==rownames(ex_var))){
      stop("bin_mat and matrix of external variables must have same rownames")
    }  
    
    R2 <- NULL
    for(tx in 1:ncol(bin_mat)){#tx <- 1
      exp_df <- data.frame(ex_var,Occ=bin_mat[,tx])
      perm <- adonis2(data=exp_df,
                      ab_mat~.,
                      by = "margin",
                      permutations=nperm,
                      parallel = n.core,
                      method = dist_method)
      
      R2 <- rbind(R2,data.frame(taxa=colnames(bin_mat)[tx],
                                R2=perm$R2[ncol(exp_df)],raw_p=perm$`Pr(>F)`[ncol(exp_df)]))
    }
  }
  return(R2)
}


find_best_Spset <- function(ab_mat,
                            bin_mat,
                            priority,
                            min_nSp=20,
                            max_nSp=50,
                            n.core=1,
                            method_bin_dist="jaccard",
                            method_ab_dist="bray"){
  if(!all(colnames(bin_mat)==colnames(ab_mat)&rownames(bin_mat)==rownames(ab_mat))){
    stop("ab_mat and bin_mat must have same rownames")
  }
  require(vegan)
  require(doParallel)
  
  cluster = makeCluster(n.core)
  registerDoParallel(cluster)  
  on.exit(stopCluster(cluster)) 
  
  bc_dist <- as.matrix(vegdist(ab_mat,method = method_ab_dist))
  bc_dist[is.na(bc_dist)] <- 0
  diag(bc_dist) <- NA
  bc_dist[upper.tri(bc_dist)] <- NA
  
  Sp_range <- min_nSp:max_nSp
  mantel_r <- foreach(i = 1:length(Sp_range),.combine = "c",.packages = "vegan")%dopar%{
    nSp <- Sp_range[i]
    bin_mat2 <- bin_mat[rownames(ab_mat),priority[1:nSp]]
    jac_dist <- as.matrix(vegdist(bin_mat2,method = method_bin_dist))
    jac_dist[is.na(jac_dist)] <- 0
    diag(jac_dist) <- NA
    jac_dist[upper.tri(jac_dist)] <- NA
    
    col <- mantel(as.dist(bc_dist),as.dist(jac_dist),method="kendall",permutations=FALSE)
    return(col$statistic)
  }
  
  return(data.frame(nSp=Sp_range,tau=mantel_r))
}


showDG_mod <- function(ela, oc, label="",SS_colmat,
                       na.color="black",minor.color="gray50",fontsize=5,
                       annot_adj=c(0.75, 2.00)){
  if(length(ela[[1]])>1){
    s <- ncol(oc)
    grobj <- GraphObj(ela)
    
    DG_mod(grobj, s,ss_list=ela[[1]], DG_sample=label,
           annot_adj=c(annot_adj[1], annot_adj[2]),SS_colmat=SS_colmat,
           na.color=na.color,minor.color=minor.color,fontsize=fontsize)
  }else{
    return(cat("only one stable state found\n"))
  }}

DG_mod <- function(
    grobj, s,ss_list, DG_sample, annot_adj=c(0.75, 2.00),
    SS_colmat,na.color="black",minor.color="gray90",fontsize=5){
  # grobj <- GraphObj(ela);s <- ncol(ocmat[[fb]]);SS_colmat=gpr;label=sprintf("%s %s",fb,pl[i]);na.color="gray90"
  range_ <- c(min(grobj$energy), max(grobj$energy))
  grobj_pre <- grobj
  grobj_pre$len_cee <- lapply(
    grobj_pre$cee, function(x){length(x)}) %>% as.integer
  #d_list <- unlist(map(1: s, function(x){paste0('d', x)}))
  grobj_ <- cbind(
    grobj_pre,
    do.call(
      cbind,
      grobj_pre$point %>% map(function(x){id2bin(x, s)})) %>% t #%>%
    #'colnames<-'(d_list)
  )
  jun_pre <- grobj_[grobj_$len_cee == 1,]
  jun <- jun_pre[order(jun_pre$nodes2xposi),]
  jen_pre <- grobj_[grobj_$len_cee > 1,]
  jen <- jen_pre[order(jen_pre$nodes2xposi),]
  
  # grobj to plot
  grobj_to_plot <- data.frame()
  if(nrow(grobj_) != 1){
    for(i in 1: (nrow(grobj_))){
      aa <- grobj_[i,]
      aa$point_str <- SS_colmat[aa$point,"rename_SS"]
      bb_pre <- grobj_[grobj_$cee %>% map(
        function(x){all(unlist(aa$cee) %in% unlist(x))}) %>% unlist(),]
      bb <- bb_pre[bb_pre$ccc_str != aa$ccc_str,][1,]
      bb$point_str <- rep("", length(bb$point))
      # bb$point_str <- paste0('C', bb$point)
      between_aa_bb <- aa
      between_aa_bb$energy <- bb$energy
      between_aa_bb$point_str <- ''
      line_break <- between_aa_bb
      line_break$nodes2xposi <- NA
      line_break$energy <- NA
      grobj_to_plot_ <- rbind(aa, between_aa_bb, bb, line_break)
      grobj_to_plot_$id <- i
      grobj_to_plot <- grobj_to_plot %>% rbind(grobj_to_plot_)
    }
  }
  
  grobj_to_plot$point_str[grepl("c(",grobj_to_plot$ccc_str,fixed = TRUE)] <- ""
  grobj_to_plot$point_str[is.na(grobj_to_plot$point_str)] <- "otherstate"
  
  pick <- which(grobj_to_plot$point_str != "" )
  # pick <- intersect(which(grobj_to_plot$point_str != "" ),
  #                    grep("c(",grobj_to_plot$ccc_str,invert=TRUE,fixed=TRUE))
  # 
  #col[is.na(col)] <- na.color
  
  minor <- intersect(which(grobj_to_plot$point_str == "otherstate" ),pick)
  grobj_to_plot$shape <-  ""
  grobj_to_plot$shape[pick] <- "SS"
  grobj_to_plot$shape[minor] <- "MinorSS"
  
  grobj_to_plot$size <-  ""
  grobj_to_plot$size[pick] <- "SS"
  grobj_to_plot$size[minor] <- "MinorSS"
  
  grobj_to_plot$size <- factor(grobj_to_plot$size, levels=c("SS", "MinorSS", ""))
  grobj_to_plot$shape <- factor(grobj_to_plot$shape, levels=c("SS", "MinorSS", ""))
  
  grobj_to_plot$point_str[minor] <- " "
  grobj_to_plot$point_str[grobj_to_plot$point_str=="otherstate"] <- ""
  grobj_to_plot$point_str2 <- factor(grobj_to_plot$point_str,
                                     levels=c(setdiff(unique(grobj_to_plot[,"point_str"]),
                                                      c(" ", "")),
                                              " ", ""))
  
  col <- c(SS_colmat[match(setdiff(unique(grobj_to_plot[,"point_str"]),
                                   c(" ", "")),SS_colmat$rename_SS),"color"],minor.color,na.color)
  # scatter with annot and line
  g <- ggplot(
    grobj_to_plot,
    aes(x=nodes2xposi, y=energy))
  g <- g + geom_path() +
    geom_point(aes(color=point_str2,shape=shape,size=size),stroke=2,show.legend = F) +
    xlab("") +
    ylab("Energy")+
    #labs(title=DG_sample)+
    geom_text(hjust=annot_adj[1], vjust=annot_adj[2],size=fontsize,
              aes(fontface=2, label=point_str)) +
    theme_bw() +
    theme(axis.title = element_text(size=18),
          axis.text.y = element_text(size=12),
          aspect.ratio = 1,
          panel.background = element_rect(fill = "white", color = "black", size = 2))+
    scale_color_manual(values=col[which(levels(grobj_to_plot$point_str2) %in% unique(grobj_to_plot$point_str2))])+
    scale_shape_manual(values=c(16,16, 18)[which(levels(grobj_to_plot$shape) %in% unique(grobj_to_plot$shape))]) +
    scale_size_manual(values=c(5,3, 1)[which(levels(grobj_to_plot$size) %in% unique(grobj_to_plot$size))]) +
    coord_cartesian(ylim=c(min(grobj_to_plot$energy,na.rm = TRUE)-0.05*(max(grobj_to_plot$energy,na.rm = TRUE)-min(grobj_to_plot$energy,na.rm = TRUE)),NA))
  
  # pie
  #g <- g + geom_scatterpie(
  #  aes(x=nodes2xposi, y=energy),
  #  data=grobj_, cols=d_list)
  
  # plot
  g <- g + ggtitle(DG_sample)
  plot(g)
  return(NULL)}

Bin_2sd <- function(df){
  binmat <- NULL
  for(i in 1:ncol(df)){#i <- 1
    x <- df[,i]
    x2 <- log(x[x>0])
    lth <- exp(mean(x2)-2*sd(x2))
    binmat <- cbind(binmat,matrix(ifelse(x<lth,0,1),ncol=1))
  }
  return(binmat)
}

blockSample <- function(mat,el,name){
  block_rand <- c()
  names <- c()
  randnames <- c()
  for(i in 1:length(unique(el))){#i <- 1
    names <- c(names,name[which(el %in% unique(el)[i])])
    randnames <- c(randnames,sample(name[which(el %in% unique(el)[i])]))
  }
  if(is.vector(mat)){
    block_rand <- mat[names,]
    names(block_rand) <- randnames
    block_rand <- block_rand[names(mat),]
  }else{
    block_rand <- mat[names,]
    rownames(block_rand) <- randnames
    block_rand <- block_rand[rownames(mat),]
  }
  
  return(list(matrix=block_rand,rownames=randnames))
}


