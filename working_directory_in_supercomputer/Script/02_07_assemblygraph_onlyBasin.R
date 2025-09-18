











##########################################################################
#function
set.seed(1234)


MinpathonAG <- function(from_state_id,state_id,hge,je,tmax,threads){
  #########################################################################
  MinTippingPath_wComp <- function(s1, s2, hg, j, tmax){
    #s1 <- "EeW0000"; s2 <- "69W14G0"; hg <- hge; j <- je; tmax <- 10000
    n <- length(hg)
    ss1 <- id2bin(s1, n)
    ss2 <- id2bin(s2, n)
    dif <- ss2 - ss1
    pd <- which(abs(dif)==1)
    if(length(pd) > 1){
      seq <- sample(pd)
    }else{
      seq <- pd
    }
    
    #############
    
    if(s1 == s2){return(Energy(ss1, hg, j))}
    if(length(seq) <= 7){
      seq_perms <- gtools::permutations(length(seq), length(seq), seq)
      
      pathes <- apply(seq_perms, 1, function(x) {#x <- seq_perms[1,]
        samplepath <- foldList(replace_bin, ss1, x)
        sapply(samplepath, function(y) Energy(y, hg, j))
      })
      
      tote <- apply(pathes, 2, function(x) {
        dife <-  x[-1] - x[-length(x)]
        sum(dife[which(dife > 0)])})
      
      pe <- NULL
      for(mins in 1:sum(tote==min(tote))){
        pe <- rbind(pe,cbind(pathes[,which(tote==min(tote))[mins]],
                    min_total_e=min(tote),
                    matrix(unlist(foldList(replace_bin, ss1, seq_perms[which(tote==min(tote))[mins],])),
                           nrow=nrow(pathes),byrow=TRUE)))
      }
      
      
      
      minpath <- as.data.frame(pe)
      #if(is.null(ncol(pe))){minpath <- pe}else{minpath <- rowMeans(pe)}
    }
    else{
      tt <- 0
      pree <- Inf
      mintote <- Inf
      
      while (tt < tmax) {
        tt <- tt + 1
        rpl <- sample(length(seq), 2) 
        seqnew <- seq 
        seqnew[rpl] <- seqnew[c(rpl[2], rpl[1])]
        samplepath <- foldList(replace_bin, ss1, seqnew)
        pathe <- sapply(samplepath, function(y) Energy(y, hg, j))
        dife <- diff(pathe)
        tote <- sum(dife[dife > 0])
        
        if(tote == mintote){
          pe <- rbind(pe,cbind(as.numeric(pathe),
                               min_total_e=tote,
                               matrix(unlist(samplepath),
                                      nrow=length(pathe),byrow=TRUE)))
          minpath <- as.data.frame(pe)
          pree <- tote
        }else{
          if (tote < mintote) {
            pe <- cbind(as.numeric(pathe),
                        min_total_e=tote,
                        matrix(unlist(samplepath),
                               nrow=length(pathe),byrow=TRUE))
            minpath <- as.data.frame(pe)
            mintote <- tote
            pree <- tote
          }  
        }
        
        if (runif(1) < min(1, exp(pree) / exp(tote))) {
          seq <- seqnew
          pree <- tote
        }
      }}
    return(minpath)}
  #####
  #state_id <- ssid
  
  state_pairs_all <- matrix(t(combn(x=state_id,2)),ncol=2)
  if(nrow(state_pairs_all) > 1){
    state_pairs <- state_pairs_all[which(state_pairs_all[,1] %in% from_state_id | state_pairs_all[,2] %in% from_state_id),]
    
  }else{
    state_pairs <- matrix(state_pairs_all[which(state_pairs_all[,1] %in% from_state_id | state_pairs_all[,2] %in% from_state_id),],
                          ncol=2,byrow=TRUE)
    
  }
  cat(sprintf("Number of state pairs: %d\n", nrow(state_pairs)))
  
  cluster = makeCluster(threads)
  registerDoParallel(cluster)
  on.exit(stopCluster(cluster))
  
  node_en <- foreach(i=1:nrow(state_pairs),.packages = c("rELA","foreach","gtools"),
                     .combine = rbind) %dopar% {
                       mnTp_comp <- MinTippingPath_wComp(state_pairs[i,1], state_pairs[i,2], hge, je, tmax)
                       return(cbind(s1=state_pairs[i,1], s2=state_pairs[i,2],as.data.frame(mnTp_comp)))
                     }
  return(node_en)
}


as.AssemblyGraph <- function(nodes,node_id,start=1,end=NA,threads=8,return_edgelist=FALSE){
  #nodes <- nodes[,-1];node_id <- nodes[,1]
  if(is.na(end)){end <- nrow(nodes)}
  
  cluster = makeCluster(threads)
  registerDoParallel(cluster)
  on.exit(stopCluster(cluster))
  
  nodes2 <- unique(nodes)
  
  cat(sprintf("total states: %d\n",nrow(nodes2)))
  ls_edge <- NULL
  for(s1 in start:end){
    cat(sprintf("\rstates: %s",s1))
    ls_edge <- rbind(ls_edge,foreach(s2 = 1:nrow(nodes), .combine = rbind) %dopar% {
      if(s1>s2){
        if(sum(nodes[s1,]!=nodes[s2,])==1){return(c(s1,s2))}
      }
    })
  }
  
  if(return_edgelist){
    return(cbind(node_id[ls_edge[,1]], node_id[ls_edge[,2]]))
  }else{
    g_ass <- graph_from_edgelist(cbind(node_id[ls_edge[,1]], node_id[ls_edge[,2]]),
                                 directed = F)
    
    return(g_ass)
    
  }}
###############################################

####################
#setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA")
#install.packages("mgcv")
library(ggplot2)
library(ggtext)
library(ggnetwork)
library(parallel)
library('RColorBrewer')
#library('scatterpie')
library('cowplot')
library('patchwork')
library("rELA")
library('tidyverse')
library("doParallel")
library(foreach)
#library(renv)
library("Rcpp")
library("RcppArmadillo")
library(igraph)

#######














#########################################################################
save.dir <- "02_07_assemblygraph_250311"
dir.create(save.dir)

########################################################################
#read original functions
library(vegan)
source("packages/01_1_function.R")

########################################################################
tx_f <- readRDS("Base_data/Fungi/taxa_list_mod.rds")
tx_b <- as.data.frame(readRDS("Base_data/OTU_basedata_set/taxonomy_list_dada2.rds"))

####################################################

taxa <- "Family"

taxcol <- readRDS("color/cols_Family_02_Barplot_240720.rds")

#read data
rel_df <- readRDS(sprintf("%s/Df_RA_%s.rds",ELA_prep_dir,taxa))
ocmat <- list(Fungi=readRDS(sprintf("%s/ELA_input_ocmat_Fungi.rds",dir_02_05)),
              Prokaryota=readRDS(sprintf("%s/ELA_input_ocmat_Prokaryote.rds",dir_02_05)))


####################################################
tdf2 <- cbind(rel_df$Fungi[,colnames(ocmat$Fungi)],
              Others=rowSums(rel_df$Fungi[,setdiff(colnames(rel_df$Fungi),
                                                   c(colnames(ocmat$Fungi),"Unidentified"))]),
              Unidentified=rel_df$Fungi[,c("Unidentified")])

colnames(tdf2) <- paste0("Fng_",colnames(tdf2))
gtdf <- gather(cbind(ID=rownames(tdf2),as.data.frame(tdf2)),taxa,count,-c(1))

gtdf$taxa2 <- factor(gtdf$taxa,levels=c(setdiff(unique(gtdf$taxa),c("Fng_Others","Fng_Unidentified")),
                                        "Fng_Others","Fng_Unidentified"))

tdb2 <- cbind(rel_df$Prokaryota[,colnames(ocmat$Prokaryota)],
              Others=rowSums(rel_df$Prokaryota[,setdiff(colnames(rel_df$Prokaryota),
                                                        c(colnames(ocmat$Prokaryota),
                                                          "Unidentified"))]),
              Unidentified=rel_df$Prokaryota[,c("Unidentified")])

colnames(tdb2) <- paste0("Prok_",colnames(tdb2))

gtdb <- gather(cbind(ID=rownames(tdb2),as.data.frame(tdb2)),taxa,count,-c(1))

gtdb$taxa2 <- factor(gtdb$taxa,levels=c(setdiff(unique(gtdb$taxa),c("Prok_Others","Prok_Unidentified")),
                                        "Prok_Others","Prok_Unidentified"))

######################################
info <- readRDS(sprintf("%s/comp_sample_info_plant2.rds",ELA_prep_dir))

sp_info <- readRDS(sprintf("%s/ELA_input_plant.rds",ELA_prep_dir))
plmat <- unique(sp_info[,-c(1,ncol(sp_info))])


ss_plant <- list(Fungi=readRDS(sprintf("%s/Fungi/detected_SS_Fungi.rds",dir_02_06)),
                Prokaryota=readRDS(sprintf("%s/Prokaryota/detected_SS_Prokaryota.rds",dir_02_06)))

ss_info <- list(Fungi=readRDS(sprintf("%s/Fungi/graph_param_SStab_full_Fungi_Family.rds",dir_02_06)),
                 Prokaryota=readRDS(sprintf("%s/Prokaryota/graph_param_SStab_full_Prokaryota_Family.rds",dir_02_06)))
###########################################
sa_list <- list.files(dir_02_06,"runSA_ST_full",full.names = TRUE,recursive = TRUE)
ela_list <- list.files(dir_02_06,"^ELA_full",full.names = TRUE,recursive = TRUE)

fb <- "Fungi"
dir.create(sprintf("%s/%s",save.dir,fb))

#gpr <- readRDS(sprintf("%s/%s/graph_param_SStab_full_Fungi_Family_0.018_0.07.rds",sa_dir,fb))
saf <- sa_list[grep(fb,sa_list)]
elas <- ela_list[grep(fb,ela_list)]

sa <- readRDS(saf)
param_pl <- list(NULL)

pl_nam <- c()
for(pl in 1:nrow(plmat)){#pl <-1
  show.progress(pl,1:nrow(plmat))
  plmat2 <- plmat[pl,]
  if(sum(plmat2[1,]==1)==1){
    pl_nam[pl] <- colnames(plmat2)[which(plmat2[1,]==1)]
  }else{
    if(all(plmat2[1,]==0)){
      pl_nam[pl] <-colnames(sp_info)[ncol(sp_info)]
    }
  }
  
  sa_para <- sa2params(sa,as.numeric(plmat2))
  hge <- sa_para[[4]];je <- sa_para[[2]]
  
  ocm <- ocmat[[fb]][which(rownames(ocmat[[fb]]) %in% info[info$plant == pl_nam[pl],"ID"]),]
  ss_all <- readRDS(elas[grep(pl_nam[pl],elas)])[[1]]
  ss <- ss_all[ss_all %in% ss_plant[[fb]][ss_plant[[fb]][,"plant"] == pl_nam[pl],"SS"]]
  param_pl[[pl]] <- list(plant=pl_nam[pl],ssid=ss,hge=hge,je=je)
}

#ssid lidt all
ss_all <- unique(unlist(sapply(param_pl,function(x){x$ssid})))

minpath <- list(NULL)

for(p in 1:length(param_pl)){#p <-6
  show.progress(pl,1:length(param_pl))
  if(length(param_pl[[p]]$ssid)>1){
    mp <- MinpathonAG(from_state_id=param_pl[[p]]$ssid,
                      state_id=param_pl[[p]]$ssid,
                          param_pl[[p]]$hge,param_pl[[p]]$je,
                      10000,threads = n.core)
   
    minpath[[p]] <- cbind(sid=apply(mp[,-c(1:4)],1,bin2id),
                          energy=mp[,3],
                          as.data.frame(mp[,-c(1:4)]))
    colnames(minpath[[p]]) <- c("sid","energy",paste0("V",1:(ncol(minpath[[p]])-2)))
    
  }else{
    minpath[[p]] <- cbind(sid=param_pl[[p]]$ssid,
                          energy=cEnergy(id2bin(param_pl[[p]]$ssid,length(param_pl[[p]]$hge)),
                                         param_pl[[p]]$hge,param_pl[[p]]$je),
                          as.data.frame(matrix(id2bin(param_pl[[p]]$ssid,
                                               length(param_pl[[p]]$hge)),nrow=1)))
    colnames(minpath[[p]]) <- c("sid","energy",paste0("V",1:(ncol(minpath[[p]])-2)))
  }
  
   
  #e_scale[p] <- max(mp[,3]-min(mp[,3]))
}

names(minpath) <- pl_nam
saveRDS(minpath,sprintf("%s/%s/minpath_%s.rds",save.dir,fb,fb))


rmp <- do.call(rbind,minpath)

#color braek point for graphics

nodes <- unique(rmp[,-c(2)])

pE <- list(NULL)
max_ed <- 0 
for(p in 1:length(param_pl)){
  pee <- apply(nodes[,-1],1,function(x){cEnergy(as.numeric(x),
                                                param_pl[[p]]$hge,
                                                param_pl[[p]]$je)})
  pE[[p]] <- pee
  names(pE[[p]]) <- nodes[,1]
  if(max(pE[[p]])>max_ed){
    max_ed <- max(pE[[p]])
  }
  }

#saveRDS(nodes,sprintf("%s/%s/nodes_list_nodeNo%s.rds",save.dir,fb,nrow(nodes)))


ass_net <- as.AssemblyGraph(nodes[,-1],nodes[,1],threads=n.core)

saveRDS(ass_net,sprintf("%s/%s/AG_igarph.rds",save.dir,fb))

d_node <- vegdist(nodes[,-1],method = "jaccard")

cmd0 <- cmdscale(d_node,k=3,eig=TRUE)
cmd <- cmdscale(d_node,k=sum(cmd0$eig>0),eig=TRUE)
pcoa <- cmd$points
rownames(pcoa) <- nodes[,1]

gnet1 <- ggnetwork::ggnetwork(ass_net,layout = pcoa[,c(1,2)])

gnet <- merge(gnet1,cbind(x=gnet1$xend,
                          y=gnet1$yend,
                          name_end=gnet1$name),by=c("x","y"))

gnet$x <- pcoa[gnet$name,1]
gnet$y <- pcoa[gnet$name,2]
gnet$xend <- pcoa[gnet$name_end,1]
gnet$yend <- pcoa[gnet$name_end,2]

gnet$SSname <- as.data.frame(ss_info[[fb]])[gnet$name,"rename_SS"]

minmax <- data.frame(x=mean(gnet$x),y=mean(gnet$y),
                     energy_dif=c(0,
                                  ceiling(max_ed)))

saveRDS(gnet,sprintf("%s/all_AG_backborn_%s.rds",save.dir,fb))
saveRDS(minmax,sprintf("%s/minmax_inAG_%s.rds",save.dir,fb))

##
node_SS <- unique(gnet[,c("name","SSname")])
node_SS <- na.omit(node_SS)
saveRDS(list(composition=nodes,
             coord_SS=cbind(node_SS$SSname,as.data.frame(pcoa[node_SS$name,])),
             eig=cmd$eig,
             coord=pcoa),sprintf("%s/%s/coord_PCoA.rds",save.dir,fb))

##
gnsl <- list(NULL);g <- list(NULL)

for(p in 1:length(param_pl)){#p <-1
  show.progress(pl,1:length(param_pl))
  mp_sl <- minpath[[p]]
  gnet_sl1 <- gnet
  gnet_sl1$energy_dif <- pE[[p]][gnet$name]
  
  gnet_sl <- gnet_sl1[which(gnet_sl1$name %in% mp_sl[,1]),]
  
  nopath <- gnet_sl1[which(!gnet_sl1$name %in% mp_sl[,1]),c("x","y")]
  
  wpath <- apply(gnet_sl[,c("xend","yend")],1,
                 paste,collapse="_") %in% apply(nopath[,c("x","y")],1,
                                                paste,collapse="_")
  
  gnet_sl[wpath,c("xend","yend")] <- NA
  
  gnsl[[p]] <- cbind(gnet_sl,energy=mp_sl[match(gnet_sl$name,mp_sl$sid),"energy"])
  
  # nopath_st <- gnsl[[p]][is.na(gnsl[[p]]$energy),"name"] 
  # 
  # e_np_st <- sapply(nopath_st,function(x){
  #   cEnergy(id2bin(x,ncol(mp_sl)-2),param_pl[[p]]$hge,param_pl[[p]]$je)
  # })
  # 
  # gnsl[[p]][is.na(gnsl[[p]]$energy),"energy"] <- e_np_st
  # 
  gnsl[[p]]$node_info <- "Internal states"
  gnsl[[p]][gnsl[[p]]$name %in% param_pl[[p]]$ssid,
            "node_info"] <- "Basin of attraction"
 
  gnsl[[p]]$node_info <- factor(gnsl[[p]]$node_info,
                                levels=c("Basin of attraction",
                                         "Basin in other hosts",
                                         "Internal states"
                                         ))
  
  gnsl[[p]]$node_size <- ifelse(gnsl[[p]]$name %in% param_pl[[p]]$ssid,
                                5,4)
  
  g[[p]] <- ggplot(gnsl[[p]],aes(x=x,xend=xend,y=y,yend=yend))+
    geom_edges(data=gnet_sl1,aes(color=energy_dif),linewidth=0.3)+
    geom_nodes(data=gnet_sl1,aes(color=energy_dif),
               size=0.8,shape=16)+
    geom_edges(aes(color=energy_dif),linewidth=2)+
    geom_nodes(aes(shape=node_info,size=node_size,color=energy_dif))+
    geom_nodes(data=function(x){x[x$node_info != "Internal states",]},
               shape=15,size=6.5,color="white")+
    geom_nodes(data=function(x){x[x$node_info != "Internal states",]},
               aes(shape=node_info,size=node_size,color=energy_dif))+
    #geom_blank(data=minmax,aes(x=x,xend=x,y=y,yend=y,color=energy_dif))+
    theme_minimal()+
    labs(color="Instability (Energy)",
         shape="States",x="PCo1",y="PCo2")+
    guides(size="none",shape=guide_legend(override.aes = list(size=4)))+
    scale_size_continuous(range = c(4,6),breaks=c(4,6))+
    scale_color_gradientn(colors=c("black","purple4","blue","turquoise1",
                                   "green","gold","orange","red")
                          #,breaks=seq(min(minmax$energy_dif),max(minmax$energy_dif),length.out=5)
    )+
    scale_shape_manual(values=c("Basin of attraction"=15,
                                "Basin in other hosts"=12,
                                "Internal states"=16))+
    coord_fixed()
  ggsave(plot=g[[p]]+theme_bw()+theme(aspect.ratio=1,
                                      panel.background = element_rect(fill = "white", color = "black", size = 2),
                                      panel.grid.major = element_blank(),   #グリッドは入れない
                                      panel.grid.minor = element_blank()),
         filename=sprintf("%s/%s/AssemblyGraph_%s_%s.pdf",save.dir,fb,param_pl[[p]]$plant,fb),
         width=7,height=6)
  ggsave(plot=g[[p]]+
           geom_nodelabel_repel(data=function(x){x[x$node_info == "Basin of attraction",]},
                                aes(label=SSname),size=7,na.rm = FALSE,
                                fontface = "bold", box.padding = unit(1, "lines"))+
           theme_bw()+theme(aspect.ratio=1,
                            panel.background = element_rect(fill = "white", color = "black", size = 2),
                            panel.grid.major = element_blank(),   #グリッドは入れない
                            panel.grid.minor = element_blank()),
         filename=sprintf("%s/%s/AssemblyGraph_withLabel_%s_%s.pdf",save.dir,fb,param_pl[[p]]$plant,fb),
         width=7,height=6)
}
names(gnsl) <- pl_nam
saveRDS(gnsl,sprintf("%s/%s/AssemblyGraph_all_%s.rds",save.dir,fb,fb))

##overview
g_over <- ggplot(gnet,aes(x=x,xend=xend,y=y,yend=yend))+
  geom_edges(color="gray90",linewidth=0.3)+
  geom_nodes(color="gray50",shape=16,size=4)+
  geom_nodes(data=gnet[which(gnet$name %in% ss_all),],
             shape=15,size=5.5,color="white")+
  geom_nodes(data=gnet[which(gnet$name %in% ss_all),],
             shape=15,size=5,color="black")+
  theme_blank()+
  coord_fixed()

ggsave(plot=g_over,
       filename=sprintf("%s/%s/AssemblyGraph_overview_%s.png",save.dir,fb,fb),
       width=8,height=8)

##all minpath graph
g_all <- (g[[3]]+ggtitle("*Pinus*")+theme(aspect.ratio = 1,
                                        plot.title = element_markdown(size=25)))+
          (g[[1]]+ggtitle("*Betula*")+theme(aspect.ratio = 1,
                                            plot.title = element_markdown(size=25)))+
          (g[[4]]+ggtitle("*Acer*")+theme(aspect.ratio = 1,
                                          plot.title = element_markdown(size=25)))+
         plot_layout(guides = 'collect')
          

ggsave(plot=g_all,sprintf("%s/%s/AssemblyGraph_PBA_%s.pdf",save.dir,fb,fb),
       w=14,h=6)

g_all2 <- (g[[2]]+ggtitle("*Larix*")+theme(aspect.ratio = 1,
                                          plot.title = element_markdown(size=25)))+
  (g[[5]]+ggtitle("*Populus*")+theme(aspect.ratio = 1,
                                    plot.title = element_markdown(size=25)))+
  (g[[6]]+ggtitle("*Juglans*")+theme(aspect.ratio = 1,
                                  plot.title = element_markdown(size=25)))+
  plot_layout(guides = 'collect')


ggsave(plot=g_all2,sprintf("%s/%s/AssemblyGraph_LPJ_%s.pdf",save.dir,fb,fb),
       w=14,h=6)

######################################
fb <- "Prokaryota"
dir.create(sprintf("%s/%s",save.dir,fb))

#gpr <- readRDS(sprintf("%s/%s/graph_param_SStab_full_Fungi_Family_0.018_0.07.rds",sa_dir,fb))
saf <- sa_list[grep(fb,sa_list)]
elas <- ela_list[grep(fb,ela_list)]

sa <- readRDS(saf)
param_pl <- list(NULL)

pl_nam <- c()
for(pl in 1:nrow(plmat)){#pl <-1
  show.progress(pl,1:nrow(plmat))
  plmat2 <- plmat[pl,]
  if(sum(plmat2[1,]==1)==1){
    pl_nam[pl] <- colnames(plmat2)[which(plmat2[1,]==1)]
  }else{
    if(all(plmat2[1,]==0)){
      pl_nam[pl] <-colnames(sp_info)[ncol(sp_info)]
    }
  }
  
  sa_para <- sa2params(sa,as.numeric(plmat2))
  hge <- sa_para[[4]];je <- sa_para[[2]]
  
  ocm <- ocmat[[fb]][which(rownames(ocmat[[fb]]) %in% info[info$plant == pl_nam[pl],"ID"]),]
  ss_all <- readRDS(elas[grep(pl_nam[pl],elas)])[[1]]
  ss <- ss_all[ss_all %in% ss_plant[[fb]][ss_plant[[fb]][,"plant"] == pl_nam[pl],"SS"]]
  param_pl[[pl]] <- list(plant=pl_nam[pl],ssid=ss,hge=hge,je=je)
}

#ssid lidt all
ss_all <- unique(unlist(sapply(param_pl,function(x){x$ssid})))

minpath <- list(NULL)

for(p in 1:length(param_pl)){#p <-6
  show.progress(pl,1:length(param_pl))
  if(length(param_pl[[p]]$ssid)>1){
    mp <- MinpathonAG(from_state_id=param_pl[[p]]$ssid,
                      state_id=param_pl[[p]]$ssid,
                      param_pl[[p]]$hge,param_pl[[p]]$je,
                      10000,threads = n.core)
    
    minpath[[p]] <- cbind(sid=apply(mp[,-c(1:4)],1,bin2id),
                          energy=mp[,3],
                          as.data.frame(mp[,-c(1:4)]))
    colnames(minpath[[p]]) <- c("sid","energy",paste0("V",1:(ncol(minpath[[p]])-2)))
    
  }else{
    minpath[[p]] <- cbind(sid=param_pl[[p]]$ssid,
                          energy=cEnergy(id2bin(param_pl[[p]]$ssid,length(param_pl[[p]]$hge)),
                                         param_pl[[p]]$hge,param_pl[[p]]$je),
                          as.data.frame(matrix(id2bin(param_pl[[p]]$ssid,
                                                      length(param_pl[[p]]$hge)),nrow=1)))
    colnames(minpath[[p]]) <- c("sid","energy",paste0("V",1:(ncol(minpath[[p]])-2)))
  }
  
  
  #e_scale[p] <- max(mp[,3]-min(mp[,3]))
}

names(minpath) <- pl_nam
saveRDS(minpath,sprintf("%s/%s/minpath_%s.rds",save.dir,fb,fb))


rmp <- do.call(rbind,minpath)

#color braek point for graphics

nodes <- unique(rmp[,-c(2)])

pE <- list(NULL)
max_ed <- 0 
for(p in 1:length(param_pl)){
  pee <- apply(nodes[,-1],1,function(x){cEnergy(as.numeric(x),
                                                param_pl[[p]]$hge,
                                                param_pl[[p]]$je)})
  pE[[p]] <- pee 
  names(pE[[p]]) <- nodes[,1]
  if(max(pE[[p]])>max_ed){
    max_ed <- max(pE[[p]])
  }
}

#saveRDS(nodes,sprintf("%s/%s/nodes_list_nodeNo%s.rds",save.dir,fb,nrow(nodes)))


ass_net <- as.AssemblyGraph(nodes[,-1],nodes[,1],threads=n.core)

saveRDS(ass_net,sprintf("%s/%s/AG_igarph.rds",save.dir,fb))

d_node <- vegdist(nodes[,-1],method = "jaccard")
cmd0 <- cmdscale(d_node,k=3,eig=TRUE)
cmd <- cmdscale(d_node,k=sum(cmd0$eig>0),eig=TRUE)
pcoa <- cmd$points
rownames(pcoa) <- nodes[,1]

gnet1 <- ggnetwork::ggnetwork(ass_net,layout = pcoa[,c(1,2)])

gnet <- merge(gnet1,cbind(x=gnet1$xend,
                          y=gnet1$yend,
                          name_end=gnet1$name),by=c("x","y"))

gnet$x <- pcoa[gnet$name,1]
gnet$y <- pcoa[gnet$name,2]
gnet$xend <- pcoa[gnet$name_end,1]
gnet$yend <- pcoa[gnet$name_end,2]

gnet$SSname <- as.data.frame(ss_info[[fb]])[gnet$name,"rename_SS"]

minmax <- data.frame(x=mean(gnet$x),y=mean(gnet$y),
                     energy_dif=c(0,
                                  ceiling(max_ed)))

saveRDS(gnet,sprintf("%s/all_AG_backborn_%s.rds",save.dir,fb))
saveRDS(minmax,sprintf("%s/minmax_inAG_%s.rds",save.dir,fb))

##
node_SS <- unique(gnet[,c("name","SSname")])
node_SS <- na.omit(node_SS)
saveRDS(list(composition=nodes,
             coord_SS=cbind(node_SS$SSname,as.data.frame(pcoa[node_SS$name,])),
             eig=cmd$eig,
             coord=pcoa),sprintf("%s/%s/coord_PCoA.rds",save.dir,fb))


##

gnsl <- list(NULL);g <- list(NULL)

for(p in 1:length(param_pl)){#p <-1
  show.progress(pl,1:length(param_pl))
  mp_sl <- minpath[[p]]
  gnet_sl1 <- gnet
  gnet_sl1$energy_dif <- pE[[p]][gnet$name]
  
  gnet_sl <- gnet_sl1[which(gnet_sl1$name %in% mp_sl[,1]),]
  
  nopath <- gnet_sl1[which(!gnet_sl1$name %in% mp_sl[,1]),c("x","y")]
  
  wpath <- apply(gnet_sl[,c("xend","yend")],1,
                 paste,collapse="_") %in% apply(nopath[,c("x","y")],1,
                                                paste,collapse="_")
  
  gnet_sl[wpath,c("xend","yend")] <- NA
  
  gnsl[[p]] <- cbind(gnet_sl,energy=mp_sl[match(gnet_sl$name,mp_sl$sid),"energy"])
  
  # nopath_st <- gnsl[[p]][is.na(gnsl[[p]]$energy),"name"] 
  # 
  # e_np_st <- sapply(nopath_st,function(x){
  #   cEnergy(id2bin(x,ncol(mp_sl)-2),param_pl[[p]]$hge,param_pl[[p]]$je)
  # })
  # 
  # gnsl[[p]][is.na(gnsl[[p]]$energy),"energy"] <- e_np_st
  # 
  gnsl[[p]]$node_info <- "Internal states"
  gnsl[[p]][gnsl[[p]]$name %in% param_pl[[p]]$ssid,
            "node_info"] <- "Basin of attraction"
  
  gnsl[[p]]$node_info <- factor(gnsl[[p]]$node_info,
                                levels=c("Basin of attraction",
                                         "Basin in other hosts",
                                         "Internal states"
                                ))
  
  gnsl[[p]]$node_size <- ifelse(gnsl[[p]]$name %in% param_pl[[p]]$ssid,
                                5,4)
  
  g[[p]] <- ggplot(gnsl[[p]],aes(x=x,xend=xend,y=y,yend=yend))+
    geom_edges(data=gnet_sl1,aes(color=energy_dif),linewidth=0.3)+
    geom_nodes(data=gnet_sl1,aes(color=energy_dif),
               size=0.8,shape=16)+
    geom_edges(aes(color=energy_dif),linewidth=2)+
    geom_nodes(aes(shape=node_info,size=node_size,color=energy_dif))+
    geom_nodes(data=function(x){x[x$node_info != "Internal states",]},
               shape=15,size=6.5,color="white")+
    geom_nodes(data=function(x){x[x$node_info != "Internal states",]},
               aes(shape=node_info,size=node_size,color=energy_dif))+
    #geom_blank(data=minmax,aes(x=x,xend=x,y=y,yend=y,color=energy_dif))+
    theme_minimal()+
    labs(color="Instability (Energy)",
         shape="States",x="PCo1",y="PCo2")+
    guides(size="none",shape=guide_legend(override.aes = list(size=4)))+
    scale_size_continuous(range = c(4,6),breaks=c(4,6))+
    scale_color_gradientn(colors=c("black","purple4","blue","turquoise1",
                                   "green","gold","orange","red")
                          #,breaks=seq(min(minmax$energy_dif),max(minmax$energy_dif),length.out=5)
                          )+
    scale_shape_manual(values=c("Basin of attraction"=15,
                                "Basin in other hosts"=12,
                                "Internal states"=16))+
    coord_fixed()
  ggsave(plot=g[[p]]+theme_bw()+theme(aspect.ratio=1,
                                      panel.background = element_rect(fill = "white", color = "black", size = 2),
                                      panel.grid.major = element_blank(),   #グリッドは入れない
                                      panel.grid.minor = element_blank()),
         filename=sprintf("%s/%s/AssemblyGraph_%s_%s.pdf",save.dir,fb,param_pl[[p]]$plant,fb),
         width=7,height=6)
  ggsave(plot=g[[p]]+
           geom_nodelabel_repel(data=function(x){x[x$node_info == "Basin of attraction",]},
                                aes(label=SSname),size=7,na.rm = FALSE,
                                fontface = "bold", box.padding = unit(1, "lines"))+
           theme_bw()+theme(aspect.ratio=1,
                            panel.background = element_rect(fill = "white", color = "black", size = 2),
                            panel.grid.major = element_blank(),   #グリッドは入れない
                            panel.grid.minor = element_blank()),
         filename=sprintf("%s/%s/AssemblyGraph_withLabel_%s_%s.pdf",save.dir,fb,param_pl[[p]]$plant,fb),
         width=7,height=6)
}
names(gnsl) <- pl_nam
saveRDS(gnsl,sprintf("%s/%s/AssemblyGraph_all_%s.rds",save.dir,fb,fb))

##overview
g_over <- ggplot(gnet,aes(x=x,xend=xend,y=y,yend=yend))+
  geom_edges(color="gray90",linewidth=0.3)+
  geom_nodes(color="gray50",shape=16,size=4)+
  geom_nodes(data=gnet[which(gnet$name %in% ss_all),],
             shape=15,size=5.5,color="white")+
  geom_nodes(data=gnet[which(gnet$name %in% ss_all),],
             shape=15,size=5,color="black")+
  theme_blank()+
  coord_fixed()

ggsave(plot=g_over,
       filename=sprintf("%s/%s/AssemblyGraph_overview_%s.png",save.dir,fb,fb),
       width=8,height=8)

##all minpath graph
g_all <- (g[[3]]+ggtitle("*Pinus*")+theme(aspect.ratio = 1,
                                          plot.title = element_markdown(size=25)))+
  (g[[1]]+ggtitle("*Betula*")+theme(aspect.ratio = 1,
                                    plot.title = element_markdown(size=25)))+
  (g[[4]]+ggtitle("*Acer*")+theme(aspect.ratio = 1,
                                  plot.title = element_markdown(size=25)))+
  plot_layout(guides = 'collect')


ggsave(plot=g_all,sprintf("%s/%s/AssemblyGraph_PBA_%s.pdf",save.dir,fb,fb),
       w=14,h=6)

g_all2 <- (g[[2]]+ggtitle("*Larix*")+theme(aspect.ratio = 1,
                                           plot.title = element_markdown(size=25)))+
  (g[[5]]+ggtitle("*Populus*")+theme(aspect.ratio = 1,
                                     plot.title = element_markdown(size=25)))+
  (g[[6]]+ggtitle("*Juglans*")+theme(aspect.ratio = 1,
                                     plot.title = element_markdown(size=25)))+
  plot_layout(guides = 'collect')


ggsave(plot=g_all2,sprintf("%s/%s/AssemblyGraph_LPJ_%s.pdf",save.dir,fb,fb),
       w=14,h=6)
