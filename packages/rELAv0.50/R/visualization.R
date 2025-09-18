
# Arguments:
# ocmatrix, hest, jest: DataFrame
# export: Boolean
# Returns:
# None
#
#'PCplot: scatter plot of community compositions colored by basin membership
#' @useDynLib rELA, .registration=TRUE
#' @export
PCplot <- function(ocmat, sa=sa, env=NULL, ssrep=NULL, pruned=TRUE, export=FALSE){
  ######################
  sa2p <- sa2params(sa, env)
  hestp <- sa2p[[1]]
  jestp <- sa2p[[2]]
  gestp <- sa2p[[3]]
  hgestp <- sa2p[[4]] 
  ####################
    if(length(ssrep)!=0){
    ssrep_ <- ssrep %>% as.data.frame() %>% 'colnames<-'(
    c('stablestates', 'stablestates2'))
    stablestates <- ssrep_$stablestates
    stablestates2 <- ssrep_$stablestates2
  PCA_ocmatrix <- as.data.frame(prcomp(x=ocmat, scale=TRUE)[5]) %>%
    'colnames<-'(lapply(1: length(names(hgestp)), function(x){
      gsub(' ', '', paste('PC', as.character(x)), fixed=TRUE)}))
  states <- as.data.frame(t(as.data.frame(
    do.call(
      rbind,
      apply(ocmat, 1, function(x){
        t(Bi(x,
             hgestp, jestp))}))[, 1]))) %>% 'rownames<-'(
               rownames(ocmat)) %>% 'colnames<-'('stablestates')
  state_PC <- cbind(states, PCA_ocmatrix) %>% left_join(ssrep_, by='stablestates')
  if(pruned==FALSE){
  g <- ggplot(state_PC, aes(x=PC1, y=PC2, color=as.factor(stablestates)),
              width=8,height=8)}else{

        g <- ggplot(state_PC, aes(x=PC1, y=PC2, color=as.factor(stablestates2)),
                    width=8,height=8)
  }
  g <- g + geom_point() + theme(aspect.ratio = 1)
  plot(g)
  if(export){return(state_PC)}else{return(NULL)}
    }else{
        cat("Please set ssrep")
    }}


# Arguments:
# ww: float
# stablestates: list
# ssenergies: list
# tippingenergies: list of lists
# Returns:
# list of lists
#
#'calc_ccc: supporting function for DisconnectivityGraph
#' @useDynLib rELA, .registration=TRUE
calc_ccc <- function(ww, stablestates, ssenergies, tippingenergies){
  sta <- data.frame(ssenergies) %>% subset(ssenergies <= ww) %>% rownames() %>%
    as.integer
  coordinates <- as.data.frame(permutations(n=length(
    stablestates), r=2, v=1:length(stablestates), repeats.allowed = TRUE))
  tip <- cbind(1: length(tippingenergies), as.data.frame(tippingenergies)) %>%
    pivot_longer(cols=colnames(as.data.frame(tippingenergies))) %>% cbind(
      coordinates) %>% subset(value <= ww) %>% select(c('V1', 'V2'))
  rips <- rbind(as.data.frame(cbind(sta, sta)) %>% 'colnames<-'(c('V1', 'V2')), tip)
  tempmat <- cbind(rips, 1) %>% pivot_wider(names_from='V2', values_from='1') %>%
    column_to_rownames('V1')
  tempmat[is.na(tempmat)] <- 0
  ccc_pre1 <- as.vector(components(graph_from_adjacency_matrix(
    adjmatrix=as.matrix(tempmat) + t(tempmat), mode='undirected', diag=F))$membership)
  ccc_pre2 <- T(as.data.frame(rbind(ccc_pre1, 1: length(ccc_pre1)))) %>% 
    group_by(ccc_pre1) %>% group_map(as.list) %>% lapply(function(x){x$V2})
  ccc_pre2 <- ccc_pre2 %>% 'names<-'(1: length(ccc_pre2))
  ccc_pre2_length <- T(as.data.frame(lapply(
    ccc_pre2, function(x){length(x)}))) %>% 'rownames<-'(1: length(ccc_pre2))
  ccc_final1 <- ccc_pre2[ccc_pre2_length %>% subset(V1 >= 2) %>% rownames]
  ccc_final2_cond1 <- ccc_pre2_length %>% subset(V1 == 1) %>% rownames
  ccc_final2_cond2 <- cbind(ccc_pre2[ccc_final2_cond1], lapply(
    ccc_pre2[ccc_final2_cond1], function(x){length(intersect(x, sta))})) %>% 
    'colnames<-'(c('node', 'len')) %>% applymap(function(x){x[1]}) %>%
    as.data.frame
  ####################################################
  if(length(ccc_final2_cond2)>0){
      ccc_final2 <- ccc_pre2[rownames(ccc_final2_cond2['len'] == 1)]
      }else{ccc_final2 <- NULL}
  ####################################################
  ccc <- c(ccc_final1, ccc_final2)
  return(ccc)
}

# Arguments:
# parent: list
# daughter: list
# Returns:
# list
#
#'Assort: supporting function for DisconnectivityGraph
Assort <- function(parent, daughter){
  parent_ <- parent %>% unlist(recursive=FALSE)
  daughter_ <- daughter %>% unlist(recursive=FALSE)
  orderkey <- do.call(
    rbind, 
    daughter_  %>%
      lapply(function(x){list(
        match(TRUE, lapply(parent_, function(y){all(x %in% y)})), -length(x), x[1]
        )})) %>% as.data.frame %>% applymap(function(x){x[1]}) %>% as.data.frame %>%
    'colnames<-'(c('setloc', 'length', 'firstN')) %>% 
    cbind(index=1: length(daughter_))
  return(list(daughter_[order(orderkey$setloc, orderkey$length, orderkey$firstN)]))
}

# Arguments:
# hee: list
# Returns:
# cee: list
#
#'calc_cee: supporting function for DisconnectivityGraph
calc_cee <- function(hee){
  cee <- list()
  for(daughter in hee){
    daughter_ <- list(daughter)
    if(length(cee) == 0){
      cee <- append(cee, daughter_)
    }else{
      parent <- cee[length(cee)]
      cee <- append(cee, Assort(parent, daughter_))
    }
  }
  return(cee)
}

# Arguments:
# xpos_prev: list
# rcee_: list
# Returns:
# output: list
#
#'calc_xpos1: supporting function for DisconnectivityGraph
calc_xpos1 <- function(xpos_prev, rcee_){
  rcee__ <- rcee_
  xpos_prev_ <- unlist(xpos_prev, recursive=FALSE)
  if(length(rcee__) == length(unlist(xpos_prev_[1], recursive=FALSE))){
    output <- xpos_prev
  }else{
    xx <- xpos_prev_
    yy <- rcee__
    pp_pre <- lapply(yy, "length<-", max(lengths(yy))) %>% as.data.frame %>% T %>%
      'colnames<-'(1: max(unlist(lapply(yy, length)))) %>%
      cbind(index=1: length(yy)) %>%
      pivot_longer(cols=!index, names_to='columns', values_to='yy_') %>%
      na.omit %>% applymap(as.integer) %>% as.data.frame %>%
      column_to_rownames(var='yy_')
    pp <- pp_pre[as.character(unlist(xx[1])), 'index']
    pp2 <- cbind(pp, unlist(xx[2])) %>% as.data.frame %>% 'colnames<-'(c('i', 'j'))
    pp3 <- unique(pp) %>% map(function(x){
      transpose(transpose(pp2[pp2$i == x,]))})
    output <- list(list(
      yy %>% map(function(x){x[1]}),
      pp3 %>% lapply(function(x){
        sum(unlist(do.call(cbind, x)[, 'j'])) / length(unlist(x$j))})))
    }
  return(output)
}

# Arguments:
# rcee: list
# Returns:
# xpositions: list
#
#'calc_xpos: supporting function for DisconnectivityGraph
calc_xpos <- function(rcee){
  xpositions <- list()
  for(rcee_ in rcee){
    if(length(xpositions) == 0){
      xpositions <- append(
        xpositions, 
        list(list(rcee_ %>% lapply(function(x){x[1]}), 1: length(rcee_))))
    }else{
      xpos_prev <- xpositions[length(xpositions)]
      xpositions <- append(xpositions, calc_xpos1(xpos_prev, rcee_))
    }
  }
  return(xpositions)
}

# Arguments:
# eladata: tuple
# Returns:
# pd.DataFrame
#
#'GraphObj: Generating graph objects for DisconnectivityGraph
#' @useDynLib rELA, .registration=TRUE
#' @export
GraphObj <- function(ela){
  stablestates <- ela[[1]]
  ssenergies <- ela[[2]]
  tippingpoints <- ela[[3]]
  tippingenergies <- ela[[4]]
  energyall <- setdiff(union(ssenergies, unlist(tippingenergies)), Inf)
  rre_pre <- rbind(
    cbind(ssenergies, stablestates) %>% 'colnames<-'(c('energy', 'point')),
    cbind(
      tippingenergies %>% as.data.frame %>% 'rownames<-'(
        1: length(tippingenergies)) %>% 'colnames<-'(
        1: length(t(tippingenergies))) %>% pivot_longer(
          cols=1: length(t(tippingenergies))),
    tippingpoints %>% as.data.frame %>% 'rownames<-'(
      1: length(tippingpoints)) %>% 'colnames<-'(
        1: length(t(tippingpoints))) %>% pivot_longer(
          cols=1: length(t(tippingpoints)))) %>%
    'colnames<-'(c('V1.1', 'energy', 'V1.2', 'point')) %>%
    select(c('energy', 'point')))
  rre <- rre_pre[!duplicated(rre_pre),]
  eee <- cbind(energyall, lapply(energyall, function(ww){calc_ccc(
    ww, stablestates, ssenergies, tippingenergies)})) %>% 'colnames<-'(c('energy', 'ccc')) %>%
    as.data.frame
  eee$energy <- unlist(eee[,'energy'])
  eee <- eee[order(eee$energy, decreasing=TRUE),]
  ese <- eee[!duplicated(eee$ccc, fromLast=TRUE),]
  eps <- ese$ccc
  ori <- eps[1][[1]]$'1'
  hee <- apply(cbind(eps %>% map(function(x){t(t(x))}), eps %>% map(
    function(x){t(setdiff(ori, unlist(x)))})), 1, function(x){unlist(
      x, recursive=FALSE)})
  cee <- calc_cee(hee)
  rcee <- cee[length(cee): 1]
  xpositions <- do.call(rbind, calc_xpos(rcee))[, 2][length(rcee): 1]
  invalid_shape <- lapply(xpositions, function(x){dim(data.frame(x))[1]}) != 1
  xpositions[invalid_shape] <- lapply(
    xpositions[invalid_shape],
    function(x){lapply(unlist(x), function(y){y})})
  nodes2xposi_pre <- cbind(
    cee=unlist(cee, recursive=FALSE),
    xpositions=unlist(xpositions, recursive=FALSE)) %>% as.data.frame
  nodes2xposi_pre$xpositions <- nodes2xposi_pre$xpositions %>% as.double
  nodes2xposi <- nodes2xposi_pre[!duplicated(nodes2xposi_pre),]
  ff_pre <- do.call(
    rbind,
    ese %>% apply(1, function(x){cbind(ccc=x$ccc, energy=x$energy)})) %>%
    as.data.frame
  ff_pre2 <- cbind(
    ccc=ff_pre$ccc,
    ccc_str=ff_pre$ccc %>% as.character,
    energy=ff_pre$energy %>% as.double)
  ff <- ff_pre2 %>% as.data.frame %>% left_join(
    cbind(
      cee=nodes2xposi$cee,
      cee_str=nodes2xposi$cee %>% as.character,
      nodes2xposi=nodes2xposi$xpositions) %>% as.data.frame,
    by=c('ccc_str' = 'cee_str'))
  ff <- ff[, c('cee', 'ccc_str', 'energy', 'nodes2xposi')]
  ff$ccc_str <- ff$ccc_str %>% as.character
  ff[, c('energy', 'nodes2xposi')] <- ff[, c('energy', 'nodes2xposi')] %>%
    apply(1, as.double) %>% T
  graphinfo_pre <- ff[order(ff$ccc_str, ff$energy),]
  graphinfo_pre2 <- graphinfo_pre[!duplicated(graphinfo_pre$ccc_str),]
  graphinfo_pre2$len_cee <- lapply(
    graphinfo_pre2$cee, function(x){length(x)}) %>% as.integer
  graphinfo_pre2$sum_cee <- lapply(graphinfo_pre2$cee, sum) %>% as.integer
  graphinfo <- graphinfo_pre2[order(
    graphinfo_pre2$len_cee, graphinfo_pre2$sum_cee),]
  mm <- as.data.frame(as.character(graphinfo$energy))
  colnames(mm) <- "energy"
  mmm <- (mm  %>% left_join(as.data.frame(rre),by='energy'))
  point <- mmm$point
  graphinfo <- cbind(graphinfo, point)
  graphinfo_final <- graphinfo[, c(
    'cee', 'ccc_str', 'energy', 'nodes2xposi', 'point')]
  return(graphinfo_final)}


# Arguments:
# grobj: data.frame 
# s: int
# DG_sample: str
# annot_adj: tuple
# Returns:
# None
#'DisconnectivityGraph: generate DisconnectivityGraph
#' @useDynLib rELA, .registration=TRUE
#' @export
DisconnectivityGraph <- function(
    grobj, s, DG_sample='DG_sample', annot_adj=c(0.75, 2.00)){

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
      aa$point_str <- paste0('C', aa$point)
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

  # scatter with annot and line
  g <- ggplot(
    grobj_to_plot,
    aes(x=nodes2xposi, y=energy, label=point_str))
  g <- g + geom_point() + xlab("") + geom_text(
    hjust=annot_adj[1], vjust=annot_adj[2], aes(fontface=2)) + geom_path()

  # pie
  #g <- g + geom_scatterpie(
  #  aes(x=nodes2xposi, y=energy),
  #  data=grobj_, cols=d_list)

  # plot
  g <- g + ggtitle(DG_sample)
  plot(g)
  return(NULL)}


#'showDG: draw DisconnectivityGraph
#' @useDynLib rELA, .registration=TRUE
#' @export
showDG <- function(ela, oc, label=""){
  if(length(ela[[1]])>1){
  s <- ncol(oc)
  grobj <- GraphObj(ela)
  DisconnectivityGraph(grobj, s, label, annot_adj=c(0.75, 2.00))
  }else{
      return(cat("only one stable state found\n"))
  }}


# Arguments:
# elap: tuple
# sa: 
# annot_adj: tuple
# Returns:
# None
#
#'showIntrGraph: draw interaction matrix on a PCoA plot
#' @useDynLib rELA, .registration=TRUE
#' @export
showIntrGraph <- function(ela, sa, th=0., annot_adj=c(0.75, 2.00)){
  ######################
  sa2p <- sa2params(sa, env)
  hestp <- sa2p[[1]]
  jestp <- sa2p[[2]]
  gestp <- sa2p[[3]]
  hgestp <- sa2p[[4]] 
  ####################
  jest <- jestp
  s <- length(hestp)
  tss <- ela[[1]]
  tssen <- ela[[2]]
  tti <- ela[[3]]
  tien <- ela[[4]]
  #memberq <- do.call(cbind, tss %>% map(
  #  function(x){CIntegerDigits(ssid=x, n=s)})) %>%
  #  'colnames<-'(tss) %>%
  #  apply(2, function(x){names(x[x == 1]) %>% as.integer()}) %>%
  #  as.data.frame() %>% as.list() %>% 'names<-'(NULL)
  
  ############
  rn <- rownames(jest)
  colnames(jest) <- rn
  jj <- jest
  jest <- as.data.frame(jest)
  ############
  
  #fmq <- memberq %>% unlist() %>% unique() %>% sort()
  #tname <- colnames(je)[fmq]
  tname <- colnames(je)
  mat <- jest[tname, tname]
  mat$species1 <- rownames(mat)
  link_ <- mat %>% pivot_longer(cols=tname, names_to='species2') %>% as.data.frame()
  plink_pre <- subset(link_, value > th)[, c('species1', 'species2')]
  qlink_pre <- subset(link_, value < -th)[, c('species1', 'species2')]
  pcp <- as.data.frame(cmdscale(jj, k=2)) %>%
    'colnames<-'(lapply(1:2, function(x){
      gsub(' ', '', paste('PCo', as.character(x)), fixed=TRUE)}))
  lmap <- cbind(
    link_[, c('species1', 'species2')],
    lmap=link_$value %>% map(function(x){min(max(-1, x), 1)}) %>% unlist())
  tmap <- cbind(
    link_[, c('species1', 'species2')],
    tmap=link_$value %>% map(function(x){min(max(-1, x), 1)}) %>% unlist())
  ltmap <- cbind(lmap, tmap=tmap$tmap)
  plink <- plink_pre %>% left_join(ltmap, by=c('species1', 'species2'))
  qlink <- qlink_pre %>% left_join(ltmap, by=c('species1', 'species2'))
  
  # grobj to plot
  grobj_to_plot <- data.frame()
  if(nrow(plink) != 1){
    for(i in 1: nrow(plink)){
      grobj_ <- pcp[unlist(plink[i, c('species1', 'species2')]),]
      grobj_$lmap <- plink[i, 'lmap']
      grobj_$tmap <- plink[i, 'tmap']
      grobj_$id <- 'plink'
      rownames(grobj_) <- NULL
      line_break <- grobj_[1,]
      line_break$PCo1 <- NA
      line_break$PCo2 <- NA
      line_break$tmap <- NA
      line_break$lmap <- NA
      grobj_to_plot <- rbind(grobj_to_plot, grobj_, line_break)
    }
  }
  if(nrow(qlink) != 1){
    for(i in 1: nrow(qlink)){
      grobj_ <- pcp[unlist(qlink[i, c('species1', 'species2')]),]
      grobj_$lmap <- qlink[i, 'lmap']
      grobj_$tmap <- qlink[i, 'tmap']
      grobj_$id <- 'qlink'
      rownames(grobj_) <- NULL
      line_break <- grobj_[1,]
      line_break$PCo1 <- NA
      line_break$PCo2 <- NA
      line_break$tmap <- NA
      line_break$lmap <- NA
      grobj_to_plot <- rbind(grobj_to_plot, grobj_, line_break)
    }
  }
  grobj_to_plot$species <- NA
  for(i in 1: nrow(pcp)){
    pcp_ <- pcp[i,]
    grobj_ <- cbind(pcp_, lmap=Inf, tmap=1, id='species', species=rownames(pcp_))
    line_break <- grobj_ %>% sapply(function(x){NA})
    line_break$id <- 'species'
    line_break$species <- rownames(pcp_)
    grobj_to_plot <- rbind(grobj_to_plot, grobj_, line_break)
  }
  grobj_to_plot$lmap2 <- grobj_to_plot$lmap %>% sapply(function(x){pnorm(x)})
  grobj_to_plot$tmap2 <- grobj_to_plot$tmap %>% sapply(function(x){1 * abs(x)})
  
  g <- ggplot(
    grobj_to_plot,
    aes(x=PCo1, y=PCo2, label=species, color=id, alpha=lmap2))
  g <- g + geom_point() + geom_text(
    hjust=annot_adj[1], vjust=annot_adj[2], aes(fontface=2)) + geom_path() +
    scale_color_manual(values = c('blue', 'red', 'black'))
  g <- g + theme_classic()
  
  # plot
  plot(g)
  
  return(NULL)}


#'showSSD: draw stable state diagram
#' @useDynLib rELA, .registration=TRUE
#' @export
showSSD <- function(gela){
  mem <- foreach(i=gela[[1]]) %do% {
      i[[1]][[1]]
  }
  uqss <- unique(unlist(mem))

  memen <- foreach(i=gela[[1]]) %do% {
    i[[1]][[2]]
    }
  minen <- min(unlist(memen))
  maxen <- max(unlist(memen))

  cols <- rep(brewer.pal(8, "Set1"), times=100)

  yrange <- c(minen-abs(0.05*minen), maxen+abs(0.05*maxen))
  par(xpd=T)
  par(mar=c(6,6,3,6))

  cc=0
  foreach(j=uqss) %do% {
    cc=cc+1
    po <- foreach(i=gela[[1]]) %do% {
        pp <- which(i[[1]][[1]]==j)
        if(identical(pp,integer(0))){NA}else{i[[1]][[2]][[pp]]}
    }

    energy <- unlist(po)
    factor <- gela[[2]]
    par(new=T)

    if(cc==1){
    env = gela[[3]]
    plot(factor, energy, type="o", ylim=yrange, col=cols[cc], xlab=names(env[is.na(env)]), ylab="Energy")}else{
        plot(factor, energy, type="o", ylim=yrange, ann=F, col=cols[cc])
    }
  }

  legend(par()$usr[2] + 0.2, par()$usr[4], legend = uqss, col = cols, pch=1, lty=1)}


#'GELA3D: draw 3D plot of energy landscape
#' @useDynLib rELA, .registration=TRUE
#' @export
showGELA3D <- function(gelsobj) {
    # surface data
    surf3d <- gelsobj[[1]]
    # line data
    l3d <- gelsobj[[2]]
    GetMesh <- function (surf3d) {
        nx <- length(union(surf3d$x, c()))
        ny <- length(union(surf3d$y, c()))
        x <- matrix(0, nrow = nx, ncol = ny)
        y <- matrix(0, nrow = nx, ncol = ny)
        z <- matrix(0, nrow = nx, ncol = ny)
        for (i in 1:nx) {
            for (j in 1:ny) {
                x[i, j] <- surf3d[i + (j - 1) * nx, ]$x
                y[i, j] <- surf3d[i + (j - 1) * nx, ]$y
                z[i, j] <- surf3d[i + (j - 1) * nx, ]$z
            }
        }
        list(x, y, z)
    }
    x <- GetMesh(surf3d)[[1]]
    y <- GetMesh(surf3d)[[2]]
    z <- GetMesh(surf3d)[[3]]
    # surface plot
    persp3D(
        x, y, z, 
        theta = 25, phi = 50, expand = 0.5, facets=NA, ticktype = "detailed",
        alpha = 0.1, 
        main = "GradELA with plot3D",
        xlab = "X", ylab = "Y", zlab = "Energy"
    )
    # lines plot
    for (a_l3d in split(l3d, l3d$Index)) {
        ifelse(length(a_l3d$x) == 1, type <- "p", type <- "l")
        x <- unlist(a_l3d$x)
        y <- unlist(a_l3d$y)
        z <- unlist(a_l3d$z)
        scatter3D(
            x, y, z,
            type=type, add=TRUE, col='green', lwd=4
        )
        i <- a_l3d$Index[[1]]
        text3D(
            mean(x), mean(y), mean(z)*0.95, 
            labels = str_c("C", i), add=TRUE, cex=1.2
        )
    }
}

