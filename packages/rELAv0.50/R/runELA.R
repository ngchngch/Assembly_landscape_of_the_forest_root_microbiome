#'Energy landscape analysis in parallel mode
#'@description Identifying stable states and tipping points
#'
#'@param sa : output of runSA function.
#'@param SS.itr : number of iterations of the process to identify the stable states.
#'@param FindingTip.itr : number of iterations of the process to identify the tipping points.
#'@param threads : number of cores (threads) used for calculation. 
#'@param reporting : enable/disable the reporting function.
#'
#'@useDynLib rELA, .registration=TRUE
#'@importFrom Rcpp sourceCpp
#'@importFrom foreach foreach
#'@export
ELA <- function(sa=sa, env=NULL, SS.itr=20000, FindingTip.itr=10000,
                 threads=1, reporting=TRUE, pp=TRUE){
  ######################
  sa2p <- sa2params(sa, env)
  hestp <- sa2p[[1]]
  jestp <- sa2p[[2]]
  gestp <- sa2p[[3]]
  hgestp <- sa2p[[4]] 
  
  ####################
  ## ||||||||||||||||||||||||||||||||||||| ##
  if(reporting){cat('Start ELA:\n')}
  start <- proc.time()[3]
	speciesName <- rownames(jestp)

	## ||||||||||||||||||||||||||||||||||||| ##
	## -- Stable state estimatin
	minsets <- as.data.frame(SSestimate(hgestp, jestp, itr = SS.itr))
    uqm <- unique(minsets[,-length(minsets[1,])])
    energy <- apply(uqm,1,function(x){Energy(x,hgestp,jestp)})
    minsets <- cbind(uqm, energy)
    minsets <- as.matrix(minsets[order(minsets[,ncol(minsets)]),])
    colnames(minsets) <- c(speciesName, "energy")
    if(reporting){cat(nrow(minsets),"stable states were found.\n")}
	sampleSS <- data.frame(minsets)
	ssid <- sprintf("SS_%s", formatC(1:nrow(minsets), width = nchar(ncol(minsets)),
	                                 flag = "0"))
	rownames(minsets) <- ssid

	## ||||||||||||||||||||||||||||||||||||| ##
	## -- Tipping point estimation

	if(length(ssid)>1){
	comb <- expand.grid(1:nrow(minsets), 1:nrow(minsets))
	comb <- comb[comb[, 1] <comb[, 2], 2:1]

	tpnodeID <- sprintf("TPnode_%s", formatC(1:nrow(comb), width = nchar(nrow(comb)),
	                             flag = "0"))
	if(pp){cluster = makeCluster(threads)
	    registerDoParallel(cluster)
	    on.exit(stopCluster(cluster))}
  if(reporting){cat('Checking', nrow(comb), 'tipping points.\n')}
	FindTip <- foreach (k = 1:nrow(comb), .combine=rbind, .packages="rELA") %dopar% {

	    minsetsub <- minsets[as.integer(comb[k, ]), ]
	    ss1 = minsetsub[1, -ncol(minsets)]
	    ss2 = minsetsub[2, -ncol(minsets)]

	    tippoint.tmp = FindingTippingpoint_cpp(s1 = ss1,  s2 = ss2,
	                                           alpha = hgestp, jj = jestp,
	                                           tmax = FindingTip.itr)

	    colnames(tippoint.tmp) <- c(speciesName, "energy")

	    ss <- data.frame(SS1=rownames(minsetsub)[1], SS1.energy=minsetsub[1,ncol(minsetsub)],
	               SS2=rownames(minsetsub)[2], SS2.energy=minsetsub[2,ncol(minsetsub)])
	    tp <- data.frame(TP=tpnodeID[k], tippoint.tmp)
	    return( cbind(ss, tp))
	}
	}

	## ||||||||||||||||||||||||||||||||||||| ##
	## -- Summarize
	if(length(ssid)>1){
	tpState <- (FindTip[,-c(1:5)]);
    uqt <- unique(tpState[,-length(tpState[1,])])
    energy <- apply(uqt,1,function(x){Energy(x,hgestp,jestp)})
    tpStateUni <- cbind(uqt, energy)
	tpid <- sprintf("TP_%s", formatC(1:nrow(tpStateUni), width = nchar(nrow(tpStateUni)),
	                                 flag = "0"))
	rownames(tpStateUni) <- tpid

	stateInfo <- rbind(data.frame(stateID=rownames(minsets), state=rep('Stable state', nrow(minsets)), minsets),
	                   data.frame(stateID=rownames(tpStateUni), state=rep('Tipping point', nrow(tpStateUni)), tpStateUni))
	                   }
	                   else{
	                   stateInfo <- data.frame(stateID=rownames(minsets), state=rep('Stable state', nrow(minsets)), minsets)
	                   }

	# ||||||||||||||||||||||||||||||||||||| ##
	if(length(ssid)>1){
	Stable <- FindTip[,c(1:4)]

	tmp <- rownames(tpStateUni); names(tmp) <-apply(tpStateUni[,-ncol(tpStateUni)], 1, paste, collapse='')
	ordered <- apply(tpState[,-ncol(tpState)], 1, paste, collapse='')
	tpidVec <- tmp[ordered]

	network <- data.frame(Stable, TP=tpidVec, Tp.energy=FindTip[,ncol(FindTip)])
	elasummary <- list(stateInfo = stateInfo, network =network)
	}else{
	elasummary <- list(stateInfo = stateInfo, network =NULL)
	}

	# ||||||||||||||||||||||||||||||||||||| ##
	if(reporting==TRUE){cat("converting...\n")}
  elaconv <- elaconvert(elasummary)
  end <- proc.time()[3]
	if(reporting==TRUE){cat(sprintf("ELA: elapsed time %.2f sec\n", end - start))}
	return(elaconv)}


#'ELA convert
#'@description Converting output of the ELA function
#'
#'@param ela : output of the ELA function.
#'
#'@export
elaconvert <- function(ela){
  stablestates <- apply(ela$stateInfo[ela$stateInfo$state=='Stable state', 3:(length(colnames(ela$stateInfo))-1)],
                        1,
                        bin2id) |> as.vector()

  ssenergy <- ela$stateInfo[ela$stateInfo$state=='Stable state','energy']

  # number of stable state
  num_ss <- nrow(ela$stateInfo[grepl('SS', ela$stateInfo$stateID),])

  # ela$stateInfo$stateID ‚ÌŒ…”. ELAparallel(runELA.R) ‚Ì ssid ‰Šú‰»‚ðŽQÆ.

  num_d <- nchar(ncol(ela$stateInfo)-2)
  get_tp <- function(i, j) {
      ss1 <- sprintf("SS_%s", formatC(i, width = num_d,flag="0"))
      ss2 <- sprintf("SS_%s", formatC(j, width = num_d,flag="0"))
      return(ela$network[ela$network$SS1==ss1 & ela$network$SS2==ss2,
                        "TP"])}

  # tippingpoints
  if(num_ss>1){
  tippingpoints <- lapply(matrix(1:num_ss), function(i)
  apply(matrix(1:num_ss), 1, function(j) {
      tp <- get_tp(i, j)
      if(length(tp)==0){return("Inf")}
      bin <- ela$stateInfo[ela$stateInfo$stateID==tp, 3:(length(colnames(ela$stateInfo))-1)]
                                          return(apply(bin,1, bin2id))
                                          })) %>% 'names<-'(sapply(1:num_ss, \(i) paste(i, sep="")))
  }else{
      tippingpoints <- lapply(matrix(1:1),function(i) apply(matrix(1:1),1, function(j) {
      return("Inf")
      })) %>% 'names<-'(sapply(1:1, \(i) paste(i, sep="")))
  }

  ## tpenergy
  if(num_ss>1){
  tpenergy <- lapply(matrix(1:num_ss),function(i) apply(matrix(1:num_ss),1, function(j) {
      tp <- get_tp(i, j)
      if(length(tp)==0){return(Inf)}
      return(ela$stateInfo[ela$stateInfo$stateID==tp, 'energy'])
      })) %>% 'names<-'(sapply(1:num_ss, \(i) paste(i, sep="")))
  }else{
      tpenergy <- lapply(matrix(1:1),function(i) apply(matrix(1:1),1, function(j) {
      return(Inf)
      })) %>% 'names<-'(sapply(1:1, \(i) paste(i, sep="")))
  }

  return(list(stablestates, ssenergy, tippingpoints, tpenergy))}


#'ELPruning
#'@description Pruning shallow basins of energy landscape
#'
#'@param ela : output of ELA function.
#'@param th : ratio of the depth of the basin to be pruned with respect to the deepest basin.
#'@param threads : number of cores (threads) used for calculation. 
#'@param reporting : enable/disable the reporting function.
#'
#'@importFrom foreach foreach
#'@export
ELPruning <- function(ela, th=0.05, threads=1, reporting=TRUE, pp=TRUE){
  start <- proc.time()[3]
  if(reporting){cat('Start pruning:\n')}
  
  if(pp){cluster = makeCluster(threads)
  registerDoParallel(cluster)
  on.exit(stopCluster(cluster))}
  
  tss <- ela[[1]]
  tssen <- ela[[2]]
  tti <- ela[[3]]
  ttien <- ela[[4]]
  if(length(tss)>1){
      rat <- as.data.frame(permutations(
          n=length(tss), r=2, v=1:length(tss), repeats.allowed = TRUE)) %>% subset(
            V1 < V2)
      if(reporting==TRUE){cat("*")}
      cmpax <- calc_mist_pax(rat, ttien, tssen)
      if(reporting==TRUE){cat(".")}
      mist <- as.numeric(unlist(cmpax[,(dim(cmpax)[2]-1):(dim(cmpax)[2])]))
      pax <- cmpax[,-dim(cmpax)[2]:-(dim(cmpax)[2]-1)]
      paxmax <- max(mist)
      ssrep <- tss
      tss0 <- tss
      tt <- 0
  while(min(pax['val1']) < th * paxmax && length(tss) > 1){
      tt <- tt + 1
      rat <- as.data.frame(permutations(
          n=length(tss), r=2, v=1:length(tss), repeats.allowed = TRUE)) %>% subset(
            V1 < V2)
      cmpax <- calc_mist_pax(rat, ttien, tssen)
      if(reporting==TRUE){if(tt%%10 == 0){cat("*")}else{cat(".")}}
      mist <- as.numeric(unlist(cmpax[,(dim(cmpax)[2]-1):(dim(cmpax)[2])]))
      pax <- cmpax[,-dim(cmpax)[2]:-(dim(cmpax)[2]-1)]
      paxmax <- max(mist)
    if(min(pax['val1']) >= th * paxmax){
      break
      }else{
        pp <- pax[order(pax$val1),][1, c('val2_1', 'val2_2')] %>% as.list() %>% unlist()
        ssrep <- replace_value(tss[pp[2]], ssrep, tss[pp[1]])
        tss <- tss[tss != tss[pp[2]]]
        tssen <- tssen[tssen != tssen[pp[2]]]
        tti_pre <- tti %>% as.data.frame()
        tti <- tti_pre[rownames(tti_pre) != rownames(tti_pre)[pp[2]],
                       colnames(tti_pre) != colnames(tti_pre)[pp[2]]] %>% as.list()
        ttien_pre <- ttien %>% as.data.frame()
        ttien <- ttien_pre[rownames(ttien_pre) != rownames(ttien_pre)[pp[2]],
                           colnames(ttien_pre) != colnames(ttien_pre)[pp[2]]] %>% as.list()
      }
    if(length(ssrep) == 1){break}
  }
  tti <- tti %>% as.data.frame() %>% 'colnames<-'(1: length(tti)) %>% as.list()
  ttien <- ttien %>% as.data.frame() %>% 'colnames<-'(1: length(ttien)) %>% as.list()
  end <- proc.time()[3]
  if(reporting==TRUE){cat(sprintf("\nELPruning: elapsed time %.2f sec\n", end - start))}
  states <- as.data.frame(cbind(tss0, ssrep))
  colnames(states) <- c("ss.before.pruning", "ss.after.pruning")
  return(list(list(tss, tssen, tti, ttien), states))
  }else{
  end <- proc.time()[3]
  if(reporting==TRUE){cat(sprintf("\nELPruning: elapsed time %.2f sec\n", end - start))}
  states <- as.data.frame(cbind(tss, tss))
  colnames(states) <- c("ss.before.pruning", "ss.after.pruning")
  return(list(list(tss, tssen, tti, ttien), states))
  }}


#'Calc mist pax
#'@description supporting function for ELpruning
#'@importFrom foreach foreach
#'
calc_mist_pax <- function(rat, ttien, tssen){
  pax <- rat
  npax <- max(pax['V2'])
  ttien.m <- as.matrix(ttien)
  
  a <- foreach(i = pax[,'V1'], j = pax[,'V2'], .combine="rbind", .packages=c("tidyverse","rELA")) %dopar% {
    row_ <- as.character(npax*(i-1) + j)
    te <- ttien.m[[j]][i]
    se <- c(tssen[i], tssen[j])
    pax[row_, 'val1'] <- min(te - se[1], te - se[2])
    val2_pre <- pax[row_, c('V1', 'V2')] %>% T %>% cbind(se)
    val2_pre <- val2_pre[order(val2_pre$se)][,row_]
    pax[row_, 'val2_1'] <- val2_pre[1]
    pax[row_, 'val2_2'] <- val2_pre[2]
    pax[row_, 'mist1'] <- te - se[1]
    pax[row_, 'mist2'] <- te - se[2]
    pax[row_,]
    }
  return(a)}

#'Replace value
#'@description supporting function for ELpruning
#'
replace_value <- function(value_before, list_before, value_after){
  list_after <- list_before
  for(i in 1: length(list_after)){
    if(list_after[i] == value_before){list_after[i] <- value_after}
  }
  return(list_after)
}


#'SteepestDescent
#' @useDynLib rELA, .registration=TRUE
#' @export
SteepestDescent <- function(state, alpha, beta){
	res <- SteepestDescent_cpp(state, alpha, beta)
	return(res)
}

#'Bi
#' @useDynLib rELA, .registration=TRUE
#' @export
Bi <- function(xi, h, j){
  sd <- SteepestDescent(xi,h,j)
  stateyen <- list(as.vector(sd[-length(sd)]), sd[length(sd)])
  state <- stateyen[[1]]
  energy <- stateyen[[2]]
  return(
    list(bin2id(state),
      energy))}

#'cEnergy
#' @export
cEnergy <- function(x, h, j){
  return(
    (- unlist(x) %*% unlist(h) - as.vector(unlist(as.vector(x))) %*%
       as.vector((unlist(x) %*% j) / 2))[1, 1]
    )
}

#'Energy
#' @useDynLib rELA, .registration=TRUE
#' @export
Energy <- function(state, alpha, beta){
	res <- cEnergy(state, alpha, beta)
	return(res)
}


#'GradELA
#'@description Estimate energy landscape across environmental gradients
#'
#'@param sa : output of runSA function.
#'@param eid : position or name of the environmental factor to be changed.
#'@param enmat : array of the vectors of environmental factors for each sample. Normalization required.
#'@param env : vector of environmental factors for specifying a single environmental condition. The environmental factor at the location specified by eid is changed. If not specified, the mean vector of all samples is used.
#'@param range : range over which the environmental factor specified by eid is changed. If not specified, the minimum and maximum observed values are used..
#'@param steps : number of steps in the range.
#'@param th : ratio of the depth of the basin to be pruned with respect to the deepest basin.
#'@param threads : number of cores (threads) used for calculation. 
#'@param SS.itr : number of iterations of the process to identify the stable states.
#'@param FindingTip.itr : number of iterations of the process to identify the tipping points.
#'
#'@useDynLib rELA, .registration=TRUE
#'@importFrom foreach foreach
#'@export
GradELA <- function(sa=sa, eid=NULL, enmat=enmat, env=NULL, range=NULL, steps=16, th=0.05, threads=1, SS.itr=20000, FindingTip.itr=10000){
  start <- proc.time()[3]
  
  cluster = makeCluster(threads)
  registerDoParallel(cluster)
  on.exit(stopCluster(cluster))
  
  if(is.null(eid)){return(cat("eid not specified prease set eid=position/colname\n"))}
  if(is.character(eid)){ei <- which(colnames(enmat)==eid)}else{ei<-eid}
  if(is.null(env)){refe <- apply(enmat, 2, mean)
  cat("env not specified, the mean of enmat is used\n")
  }else{refe <- env}
  if(is.null(range)){mi <- min(enmat[,ei])}else{mi <- range[1]}
  if(is.null(range)){ma <- max(enmat[,ei])}else{ma <- range[2]}
  dm <- (ma - mi)/(steps - 1)
  de <- seq(mi, ma, dm)
  cat("processing(")
  cat(length(de))
  cat(") |")
  els <- foreach(i=de) %do% {
      cat("=")
      ee <- replace(refe,ei,i)
      elanp <- ELA(sa, ee, SS.itr=SS.itr, FindingTip.itr=FindingTip.itr, threads=threads, reporting=FALSE, pp=FALSE)
      ELPruning(elanp, th, threads=threads, reporting=FALSE, pp=FALSE)
  }
  cat("|\n")
  end <- proc.time()[3]
  refe[eid] <- NA
  cat(sprintf("Elapsed time %.2f sec\n", end - start))
  return(list(els, de, refe))
}


#'SSentropy
#' @useDynLib rELA, .registration=TRUE
#' @export
SSentropy <- function(state, ss,
          				alpha, beta, 
          				seitr=1000, convTime=10000){
	res <- SSentropy_cpp(uoc= state, ss= ss,
          				alpha= alpha, beta=beta, 
          				seitr=seitr, convTime=convTime)
	return(res)
}


#'stability
#'@description Calculate energy gap and stable state entropy as stability indices and returns a dataframe.
#'
#'@param sa : output of runSA function.
#'@param ocmat : array of the vectors of presence/absence status for each sample. Row is sample, column is species.
#'@param enmat : array of the vectors of environmental factors for each sample. Normalization required.
#'@param threads : number of cores (threads) used for calculation. 
#'
#'@importFrom foreach foreach
#' @useDynLib rELA, .registration=TRUE
#' @export
Stability <- function(sa, ocmat, enmat=NULL, threads=1, reporting=TRUE){
  if(reporting){cat('Start Stability:\n')}
  start <- proc.time()[3]
  sadim <- dim(sa[[1]])
  if(reporting && sadim[2] > sadim[1]+1 && is.null(enmat)){cat("sa obtained with enmat, but enmat is not given.\nassigned a zero vector as the environmental condition for all data points.\n")}
  if(reporting && sadim[2] == sadim[1]+1 && !is.null(enmat)){cat("sa obtained without enmat, but enmat is given. enmat is omitted.\n")}
  
  no <- ncol(ocmat)
  if(!is.null(enmat) && sadim[2] > sadim[1]+1){
      ne <- ncol(enmat)
      cluster = makeCluster(threads)
      registerDoParallel(cluster)
      on.exit(stopCluster(cluster))
      stab <- foreach(x=t(cbind(as.data.frame(ocmat), enmat)), .combine="rbind", .packages=c("rELA", "tidyverse", "foreach")) %dopar%{
      env <- x[(no+1):(no+ne)]
      sa2p <- sa2params(sa, env)
      je <- sa2p[[2]]
      hge <- sa2p[[4]]
      state <- bin2id(x[1:no]) 
      # basin #
      bi <- Bi(x[1:no], hge, je)
      ssid <- bi[[1]]
      e.stable <- bi[[2]] 
      # energy gap #
      e.realize <- cEnergy(x[1:no], hge, je)
      energy.gap <- as.numeric(e.realize - e.stable)
      # ss entropy #
      minsets <- SSestimate(hge, je, itr = 10000)
      minsets <- unique(minsets[,-ncol(minsets)])
      sse <- SSentropy(
          as.matrix(t(as.data.frame(x[1:no]))),
                       ss=minsets, alpha= hge, beta=je, seitr=1000, convTime=10000)
      ss.entropy <- sse[,1]
      cbind(as.data.frame(t(as.matrix(c(energy.gap, ss.entropy, e.realize, e.stable)))),
            as.data.frame(t(as.matrix(c(state, ssid))))) %>% "colnames<-" (c("energy.gap", "ss.entropy","e.realize","e.stable","state.id","stable.state.id"))
      }
     rownames(stab) <- rownames(ocmat)
     }else{
      sa2p <- sa2params(sa)
      je <- sa2p[[2]]
      hge <- sa2p[[4]]
      sampleSS <- t(apply(ocmat, 1, SteepestDescent, alpha=hge, beta=je))
      sse <- SSentropy(ocmat, ss=unique(sampleSS[,-ncol(sampleSS)]), alpha= hge, beta=je, seitr=1000, convTime=10000)
      if(sum(sse[,3])>1){cat("convTime may be too small, try convTime>10000.\n")}
      ss.entropy <- sse[,1]
      names(ss.entropy) <- rownames(ocmat)
      e.stable <- sampleSS[, ncol(sampleSS)]
      e.realize <- apply(ocmat, 1, cEnergy, hge, je)
      energy.gap <- e.realize - e.stable
      state.id <- apply(ocmat, 1, bin2id)
      stable.state.id <- apply(sampleSS[, -ncol(sampleSS)], 1, bin2id)
      stab <- cbind(as.data.frame(cbind(energy.gap, ss.entropy, e.realize, e.stable)),
      as.data.frame(cbind(state.id, stable.state.id)))}
  end <- proc.time()[3]
  if(reporting==TRUE){cat(sprintf("Stability: elapsed time %.2f sec\n\n", end - start))}
  return(stab)}


#'gstability
#'@description Calculate energy gap, energy barrier and stable state entropy as stability indices and returns a list of 2 elements: the dataframe for pruned/non-pruned energy landscape, respectively.
#'
#'@param sa : output of runSA function.
#'@param ocmat : array of the vectors of presence/absence status for each sample. Row is sample, column is species.
#'@param enmat : array of the vectors of environmental factors for each sample. Normalization required.
#'@param th : ratio of the depth of the basin to be pruned with respect to the deepest basin.
#'@param threads : number of cores (threads) used for calculation. 
#'
#' @importFrom foreach foreach
#' @useDynLib rELA, .registration=TRUE
#' @export
gStability <- function(sa, ocmat, enmat=NULL, th=0., threads=1, reporting=TRUE){
  if(reporting){cat('Start gStability:\n')}
  start <- proc.time()[3]
  sadim <- dim(sa[[1]])
  if(reporting && sadim[2] > sadim[1]+1 && is.null(enmat)){cat("sa obtained with enmat, but enmat is not given.\nassigned a zero vector as the environmental condition for all data points.\n")}
  if(reporting && sadim[2] == sadim[1]+1 && !is.null(enmat)){cat("sa obtained without enmat, but enmat is given. enmat is omitted.\n")}

  cluster = makeCluster(threads)
  registerDoParallel(cluster)
  on.exit(stopCluster(cluster))

  if(!is.null(enmat) && sadim[2] > sadim[1]+1){
    no <- ncol(ocmat)
    mo <- nrow(ocmat)
    ne <- ncol(enmat)
    ### ela
    elall <- foreach(x=t(cbind(as.data.frame(ocmat), enmat)), .combine="rbind", .packages=c("rELA", "tidyverse", "foreach", "gtools")) %dopar%{
        env <- x[(no+1):(no+ne)]
        elanp <- ELA(sa, env=env, threads=threads, reporting=FALSE, pp=FALSE)
        ela <- ELPruning(elanp, th=th, reporting=FALSE, pp=FALSE)
        list(elanp, ela)}
     
    ### stability evaluation (1)
    stab <- foreach(x=t(cbind(as.data.frame(ocmat), enmat)), ela=elall[(mo+1):(2*mo)], elanp=elall[1:mo], .combine="rbind") %do% {
        env <- x[(no+1):(no+ne)]
        sa2p <- sa2params(sa, env)
        je <- sa2p[[2]]
        hge <- sa2p[[4]]
        state <- bin2id(x[1:no])
        #####
        stablestates <- ela[[1]][[1]]
        stablen <- ela[[1]][[2]]
        tippingpoints <- ela[[1]][[3]]
        tippingen <- ela[[1]][[4]]
        stablestates.np <- elanp[[1]]
        stablen.np <- elanp[[2]]
        tippingpoints.np <- elanp[[3]]
        tippingen.np <- elanp[[4]]
        # basin #
        bi <- Bi(x[1:no], hge, je)
        ssid.np <- bi[[1]]
        ssid <- ela[[2]][,2][ela[[2]][,1]==bi[[1]]]
        e.stable <- stablen[stablestates==ssid]
        e.stable.np <- bi[[2]]
        # energy gap #
        e.realize <- cEnergy(x[1:no], hge, je)
        energy.gap <- as.numeric(e.realize - e.stable)
        energy.gap.np <- as.numeric(e.realize - e.stable.np)
        # energy barrier #
        if(length(stablestates) > 1){
            e.tipping <- min(unlist(c(as.data.frame(tippingen)[stablestates==ssid,],as.data.frame(tippingen)[,stablestates==ssid])))
            energy.barrier <- e.tipping - e.stable}else{
            e.tipping <- Inf
            energy.barrier <- Inf}
        if(length(stablestates.np) > 1){
            e.tipping.np <- min(unlist(c(as.data.frame(tippingen.np)[stablestates.np==ssid.np,],as.data.frame(tippingen.np)[,stablestates.np==ssid.np])))
            energy.barrier.np <- e.tipping.np - e.stable.np}else{
            e.tipping.np <- Inf
            energy.barrier.np <- Inf}
        # ss entropy #
        minsets <- as.matrix(t(sapply(stablestates, function(x){id2bin(x, no)})))
        sse <- SSentropy(as.matrix(t(as.data.frame(x[1:no]))), ss=minsets, alpha= hge, beta=je, seitr=1000, convTime=10000)
        ss.entropy <- sse[,1]
        minsets.np <- as.matrix(t(sapply(stablestates.np, function(x){rELA::id2bin(x, no)})))
        sse.np <- SSentropy(as.matrix(t(as.data.frame(x[1:no]))), ss=minsets.np, alpha= hge, beta=je, seitr=1000, convTime=10000)
        ss.entropy.np <- sse.np[,1]
        cbind(as.data.frame(t(as.matrix(c(energy.gap, ss.entropy, energy.barrier, e.realize, e.stable, e.tipping)))),
              as.data.frame(t(as.matrix(c(state, ssid)))),
              as.data.frame(t(as.matrix(c(energy.gap.np, ss.entropy.np, energy.barrier.np, e.realize, e.stable.np, e.tipping.np)))),
              as.data.frame(t(as.matrix(c(state, ssid.np))))) %>% "colnames<-" (
                  c("energy.gap", "ss.entropy", "energy.barrier", "e.realize", "e.stable", "e.tipping", "state.id", "stable.state.id",
                    "energy.gap.np", "ss.entropy.np", "energy.barrier.np", "e.realize", "e.stable.np", "e.tipping.np", "state.id.np", "stable.state.id.np"))}
     rownames(stab) <- rownames(ocmat)

     ### Stability evaluation (2)
     stab2 <- foreach(x=t(cbind(as.data.frame(ocmat), enmat)), ela=elall[(mo+1):(2*mo)], elanp=elall[1:mo]) %do% {
        env <- x[(no+1):(no+ne)]
        sa2p <- sa2params(sa, env)
        stablestates <- ela[[1]][[1]]
        stablen <- ela[[1]][[2]]
        bin <- as.list(lapply(stablestates, function(x){id2bin(x, ncol(ocmat))}))
        sstable <- as.data.frame(cbind(stablestates,
                                       stablen, t(as.data.frame(bin)))) %>% 'colnames<-'(c('ID', 'Energy', colnames(ocmat))) %>% 'rownames<-'(1: length(stablestates))
        list(list(sa2p, sstable), list(ocmat, env, sa, ela, elanp))}
     stab2 <- list(lapply(stab2, function(x) x[[1]]), lapply(stab2, function(x) x[[2]]))
     }else{
      no <- ncol(ocmat)
      sa2p <- sa2params(sa)
      je <- sa2p[[2]]
      hge <- sa2p[[4]]
      ### ela
      elanp <- ELA(sa, env=NULL, threads=threads, reporting=FALSE, pp=FALSE)
      ela <- ELPruning(elanp, th=th, reporting=FALSE, pp=FALSE)

      ### Stability evaluation (1)
      stablestates <- ela[[1]][[1]]
      stablen <- ela[[1]][[2]]
      tippingpoints <- ela[[1]][[3]]
      tippingen <- ela[[1]][[4]]
      stablestates.np <- elanp[[1]]
      stablen.np <- elanp[[2]]
      tippingpoints.np <- elanp[[3]]
      tippingen.np <- elanp[[4]]
      stab <- foreach(x=t(as.data.frame(ocmat)), .combine="rbind", .packages=c("rELA", "tidyverse", "foreach")) %dopar%{
        state <- bin2id(x[1:no])
        # basin #
        bi <- Bi(x[1:no], hge, je)
        ssid.np <- bi[[1]]
        ssid <- ela[[2]][,2][ela[[2]][,1]==bi[[1]]]
        e.stable <- stablen[stablestates==ssid]
        e.stable.np <- bi[[2]]
        # energy gap #
        e.realize <- cEnergy(x[1:no], hge, je)
        energy.gap <- as.numeric(e.realize - e.stable)
        energy.gap.np <- as.numeric(e.realize - e.stable.np)
        # energy barrier #
        if(length(stablestates) > 1){
            e.tipping <- min(unlist(c(as.data.frame(tippingen)[stablestates==ssid,],as.data.frame(tippingen)[,stablestates==ssid])))
            energy.barrier <- e.tipping - e.stable}else{
            e.tipping <- Inf
            energy.barrier <- Inf}
        if(length(stablestates.np) > 1){
            e.tipping.np <- min(unlist(c(as.data.frame(tippingen.np)[stablestates.np==ssid.np,],as.data.frame(tippingen.np)[,stablestates.np==ssid.np])))
            energy.barrier.np <- e.tipping.np - e.stable.np}else{
            e.tipping.np <- Inf
            energy.barrier.np <- Inf}
        # ss entropy #
        minsets <- as.matrix(t(sapply(stablestates, function(x){id2bin(x, no)})))
        sse <- SSentropy(as.matrix(t(as.data.frame(x[1:no]))), ss=minsets, alpha= hge, beta=je, seitr=1000, convTime=10000)
        ss.entropy <- sse[,1]
        minsets.np <- as.matrix(t(sapply(stablestates.np, function(x){id2bin(x, no)})))
        sse.np <- SSentropy(as.matrix(t(as.data.frame(x[1:no]))), ss=minsets.np, alpha= hge, beta=je, seitr=1000, convTime=10000)
        ss.entropy.np <- sse.np[,1]
        cbind(as.data.frame(t(as.matrix(c(energy.gap, ss.entropy, energy.barrier, e.realize, e.stable, e.tipping)))),
              as.data.frame(t(as.matrix(c(state, ssid)))),
              as.data.frame(t(as.matrix(c(energy.gap.np, ss.entropy.np, energy.barrier.np, e.realize, e.stable.np, e.tipping.np)))),
              as.data.frame(t(as.matrix(c(state, ssid.np))))) %>% "colnames<-" (
                  c("energy.gap", "ss.entropy", "energy.barrier", "e.realize", "e.stable", "e.tipping", "state.id", "stable.state.id",
                    "energy.gap.np", "ss.entropy.np", "energy.barrier.np", "e.realize", "e.stable.np", "e.tipping.np", "state.id.np", "stable.state.id.np"))}
      rownames(stab) <- rownames(ocmat)

      ### Stability evaluation (2)
      bin <- as.list(lapply(stablestates, function(x){id2bin(x, ncol(ocmat))}))
      sstable <- as.data.frame(cbind(stablestates,
                                      stablen, t(as.data.frame(bin)))) %>% 'colnames<-'(c('ID', 'Energy', colnames(ocmat))) %>% 'rownames<-'(1: length(stablestates))
      stab2 <- list(list(sa2p, sstable), list(ocmat, NULL, sa, ela, elanp))
     }
  end <- proc.time()[3]
  if(reporting==TRUE){cat(sprintf("gStability: elapsed time %.2f sec\n\n", end - start))}
  return(list(as.data.frame(stab)[1:8], as.data.frame(stab)[9:16], stab2[[1]], stab2[[2]]))
  }

#'sstable
#'@description Summary of stable states.
#'
#'@param ela: output of ELA or the first return value of ELPruning.
#'
#' @useDynLib rELA, .registration=TRUE
#' @export
sstable <- function(ela){
  stablestates  <- ela[[1]]
  stablen <- ela[[2]]
  tippingpoints <- ela[[3]]
  tippingen <- ela[[4]]
  bin <- as.list(lapply(stablestates, function(x){id2bin(x, ncol(ocmat))}))
  names(bin) <- stablestates
  sstable <- as.data.frame(cbind(stablestates, stablen, t(as.data.frame(bin)))) %>%
    'colnames<-'(c('ID', 'Energy', colnames(ocmat))) %>%
    'rownames<-'(1: length(stablestates))
  return(sstable)}

#'tptable
#'@description Summary of tipping points.
#'
#'@param ela: output of ELA or the first return value of ELPruning.
#'
#' @useDynLib rELA, .registration=TRUE
#' @export
tptable <- function(ela){
  stablestates  <- ela[[1]]
  stablen <- ela[[2]]
  tippingpoints <- ela[[3]]
  tippingen <- ela[[4]]
  tpid.df <- as.data.frame(tippingpoints)
  tpid <- unlist(tpid.df[upper.tri(tpid.df, diag = FALSE)], use.names = FALSE)
  tpen.df <- as.data.frame(tippingen)
  tpen <- unlist(tpen.df[upper.tri(tpen.df, diag = FALSE)], use.names = FALSE)
  tpbin <- as.list(lapply(tpid, function(x){id2bin(x, ncol(ocmat))}))
  tp.ss1 <- stablestates[unlist(lapply(seq_along(tippingpoints), function(i) {
      non_inf_count <- sum(tippingpoints[[i]] != "Inf")
      return(rep(i, non_inf_count))}))]
  tp.ss2 <-stablestates[unlist(lapply(tippingpoints, function(x) which(x != "Inf")))]

  tptable <- as.data.frame(cbind(tpid, tp.ss1, tp.ss2, tpen, t(as.data.frame(tpbin)))) %>%
    'colnames<-'(c('TP','SS1','SS2', 'Energy', colnames(ocmat))) %>%
    'rownames<-'(1: length(tpid))
  return(tptable)
}

####################################
####################################

#'Calculation stability indices.
#'@description Claculating stable state energy, difference between sample's energy and stable state energy and stable state entropy.
#'
#'@param data : data is a matrix.
#'@param alpha : alpha is a vector which explicit/implicit preference. Use output of runSA function.
#'@param J : J is a matrix which shows species preference. Use output of runSA function.
#'@param 
#'
#' @useDynLib rELA, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @export
calcStability <- function(data=NULL, alpha, J,
						  plun=0,
						  seitr=1000, convTime=10000){
						  	
    start <- proc.time()[3]
    sampleSS <- t(apply(data, 1, SteepestDescent_cpp, alpha=alpha, beta=J))
   	
   	uniqueSS <- data.frame(unique(sampleSS[,-ncol(sampleSS)]))
	uniqueSS$energy <- NA
	
	for(i in 1:nrow(uniqueSS)){
	    
	    state <- uniqueSS[i,-ncol(uniqueSS)]
	    detect <- apply(sampleSS[,-ncol(sampleSS)], 1, function(x){ all(x==state) })
	    
	    uniqueSS[i, ncol(uniqueSS)] <- mean(sampleSS[detect, ncol(sampleSS)])
	}
	sampleSS <- uniqueSS[order(uniqueSS$energy),]
	
    ## ================================ ##
	if(plun>0){
		
		## |||||||||||||||||||||||||||| ##
		## plunned landscape
	    elasummary <- ELAparallel(alpha=res[,1], J=res[,-1], threads=8)
		redela <- ELplunning(elasummary, th=0.1)
		
		stateInfo <- elasummary[[1]]
		ssInfo <- stateInfo[stateInfo$state=='Stable state',]
		
		## |||||||||||||||||||||||||||| ##
		## -- Convert to deeper landscape
		
		sampleSSpluned <- as.matrix(sampleSS)
		for(i in 1:nrow(ssInfo)){
		    
		    state <- ssInfo[i,-c(1, 2, ncol(ssInfo))]
		    detect <- apply(sampleSS[,-ncol(sampleSS)], 1, function(x){ all(x==state) })
		    
		    if(sum(detect)>1){
		        id <- ssInfo[i,1]
		        
		        if( !any(redela$log[,2]==id)){
		            log <- redela$log
		            change <- log[id,2]
		            
		            sampleSSpluned[detect,] <- matrix( unlist(ssInfo[change, -c(1,2)]), 
		            								   ncol=ncol(sampleSStmp), nrow=sum(detect), 
		            								   byrow=TRUE)
		        }
		        
		    }
		
		}
		sampleSS <- sampleSSpluned
        ## |||||||||||||||||||||||||||| ##
	}
    
    stablestate=t(apply(sampleSS, 1, 
                    function(x){ 
                        sep <- round(length(x)/20)
                        line <- c()
                        for(i in 1:sep){
                            y <- x[-length(x)]
                            bi <- na.omit(y[((1+20*(i-1)):(20*i))])
                            line <- c(line, str2i( paste0(bi, collapse='')))
                            
                        }
                        
                        ssid=paste(line, collapse='-')
                        c(ssid, x[length(x)]) }))
                  
    ssenergy <- sampleSS[, ncol(sampleSS)]
    energy <- apply(data, 1, cEnergy, alpha= alpha, beta=J)
    
    energy.gap <- energy-ssenergy
    ssent <- SSentropy_cpp(uoc=data, ss=unique(sampleSS[,-ncol(sampleSS)]),
          				alpha= alpha, beta=J, 
          				seitr=seitr, convTime=convTime)
          	  
    return( data.frame(StableState.id=stablestate[,1], SSenergy=ssenergy,
                       energy.gap=energy.gap, SSentropy= ssent[,1]))
                       
    ## ================================ ##
    end <- proc.time()[3]
    cat(sprintf('Stability calculation done.\nElapsed time %.2f sec\n\n', end-start))                   
}

#'ELAall 
#'@description Run all ELA pipeline
#'@param data : data is a binary matrix.
#'@param alpha : alpha is a vector which explicit/implicit preference. Use output of runSA function.
#'@param J : J is a matrix which shows species preference. Use output of runSA function.
#'@param 
#'
#' @useDynLib rELA, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @export
ELAall <- function( data=NULL, alpha, J,
				    SS.itr=20000, FindingTip.itr=10000,
					plun=0,
					seitr=1000, convTime=10000,
					threads=1){
						  	
    start <- proc.time()[3]
    
   	elasummary <- ELAparallel(alpha=res[,1], J=res[,-1], 
   							  SS.itr= SS.itr, FindingTip.itr=10000, threads=threads)
   	
    ## ================================ ##    
    sampleSS <- t(apply(data, 1, SteepestDescent_cpp, alpha=alpha, beta=J))
   	
	if(plun>0){
		
		## |||||||||||||||||||||||||||| ##
		## plunned landscape
	    
		redela <- ELplunning(elasummary, th= plun)
		
		stateInfo <- elasummary[[1]]
		ssInfo <- stateInfo[stateInfo$state=='Stable state',]
		
		## |||||||||||||||||||||||||||| ##
		## -- Convert to deeper landscape
		
		sampleSSpluned <- as.matrix(sampleSS)
		ssid <- rep(NA, nrow(sampleSS))
		for(i in 1:nrow(ssInfo)){
		    
		    state <- ssInfo[i,-c(1, 2, ncol(ssInfo))]
		    detect <- apply(sampleSS[,-ncol(sampleSS)], 1, function(x){ all(x==state) })
		    ssid[detect] <- rownames(ssInfo)[i]
        
            if(sum(detect)>1){
		        id <- ssInfo[i,1]
		        
		        if( !any(redela$log[,2]==id)){
		            log <- redela$log
		            change <- log[id,2]
		            
		            sampleSSpluned[detect,] <- matrix( unlist(ssInfo[change, -c(1,2)]), 
		            								   ncol=ncol(sampleSSpluned), nrow=sum(detect), 
		            								   byrow=TRUE)
		        }
		        
		    }
		
		}
		sampleSS <- sampleSSpluned
        ## |||||||||||||||||||||||||||| ##
	}
                      
    ssenergy <- sampleSS[, ncol(sampleSS)]
    energy <- apply(data, 1, cEnergy, alpha= alpha, beta=J)
    
    energy.gap <- energy-ssenergy
    ssent <- SSentropy_cpp(uoc=data, ss=unique(sampleSS[,-ncol(sampleSS)]),
          				alpha= alpha, beta=J, 
          				seitr=seitr, convTime=convTime)
          	  
    stability=data.frame(StableState.id=ssid, SSenergy=ssenergy,
                     energy.gap=energy.gap, SSentropy= ssent[,1])


	return( list(elaResult=elasummary, elaPlunning=redela, stability=stability ))
                    
    ## ================================ ##
    end <- proc.time()[3]
    cat(sprintf('Stability calculation done.\nElapsed time %.2f sec\n\n', end-start))                   
}