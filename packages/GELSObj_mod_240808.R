
GELSObj_mod <- function(gela, sa, threads = 1){
  de <- gela[[2]]
  s <- length(sa[1][[1]][,1])
  refenv <- gela[[3]]
  #derange <- range((sa[[2]])[,is.na(refenv)])
  #de0 <- seq(derange[[1]], derange[[2]], length.out = length(de))
  
  #elsurf <- foreach(x=gela[[1]], .packages=c("tidyverse","rELA","gtools","igraph")) %dopar% {GraphObj(x[[1]])}
  elsurf <- foreach(x=gela[[1]]) %do% {GraphObj(x[[1]])}
  ssp <- lapply(elsurf, function(x) {
    x[which(sapply(x$ccc_str, function(y) length(eval(parse(text = y))) == 1)), ]
  })
  sssp <- foreach(x=ssp, y=de) %do% {
    mat <- matrix(c(rep(y, nrow(x)), x[,3], x[,5]), nrow(x))
    split(mat, seq_len(nrow(mat)))}
  sssp <- unlist(sssp, recursive = FALSE)
  sslist <- unique(unlist(lapply(ssp, function(df) df[, 5])))
  if(length(sslist)>1){
    ss_combi <- combn(sslist, 2)
    
    #neibor states in assembly graph
    mm <- apply(ss_combi, 2, function(pair) {
      dist <- hamming_distance(id2bin(pair[1], s), id2bin(pair[2], s))
      if (dist == 1) return(pair)})
    mm <- Filter(function(x) !is.null(x), mm)
    
  }else{
    mm <- NULL
  }
  
  concompo <- gen_concompo(mm)
  
  grp <- c(concompo, sapply(setdiff(sslist, unlist(concompo)), function(x) list(x)))
  ssspa <- lapply(grp, function(gg) {
    lapply(sssp, function(ss) {
      ss[grepl(paste(gg, collapse = "|"), ss[3])]
    })})
  
  ssspa <- lapply(ssspa, function(x) x[sapply(x, function(y) length(y) > 0)]) #��̗v�f�̍폜
  
  if(length(grp)==1){
    ssspa <- lapply(ssspa, function(x) x[order(as.numeric(sapply(x, '[[', 1)))]) #�e�v�f���̕��ёւ�  
  }else{
    ssspa <- sapply(ssspa, function(x) x[order(as.numeric(sapply(x, '[[', 1)))]) #�e�v�f���̕��ёւ�
  }
  
  ssspa <- ssspa[order(sapply(ssspa, length))] #�������̕��ёւ�
  ssspa <- ssspa[order(sapply(ssspa, function(x){as.numeric(x[[1]][[1]])}))] #�n�_�ɂ����ёւ�
  ssspa_unlist <- unlist(ssspa, recursive = FALSE)
  spl <- con_split(ssspa_unlist, sapply(ssspa_unlist, '[[', 3))
  tssspastb <- lapply(spl, function(x) ssextend(x, de))
  ssid <- sapply(tssspastb, function(x){x[[2]][[1]]})
  
  if(length(ssid) > 1){
    hh <- apply(matrix(stats::embed(ssid, 2)[, 2:1],ncol=2),
                1, function(x) hamming_distance(id2bin(x[1], s), id2bin(x[2], s)))
    ml <- max(7, max(hh))
    
    posit <- Reduce(function(x, y) x + y, hh, accumulate = TRUE, init = ml + 1)
  }else{
    ml <- 7
    posit <- 8
  }
  
  cluster = makeCluster(threads)
  registerDoParallel(cluster)
  on.exit(stopCluster(cluster))
  
  filledsurf <- foreach(i=seq(length(de)), .packages=c("rELA", "tidyverse", "purrr", "foreach", "gtools"),
                        .export = c("tssspastb","refenv","sa","posit","ml", "MinTippingPath", "MarginePath", "foldList", "replace_bin")) %dopar% {
                          v <- lapply(tssspastb, function(x) x[i,])
                          v <- do.call(rbind, lapply(v, function(x) as.matrix(x[,])))
                          xde <- as.numeric(v[1,1])
                          basenv <- refenv
                          basenv[is.na(basenv)] <- xde
                          basenv <- t(basenv)
                          sa2p <- sa2params(sa, c(basenv))
                          he <- sa2p[[1]]
                          je <- sa2p[[2]]
                          ge <- sa2p[[3]]
                          hge <- sa2p[[4]]
                          ###
                          pstb <- as.numeric(v[,3])
                          sid <- ssid[as.numeric(v[,3])==1]
                          zx <- posit[as.numeric(v[,3])==1]
                          ###
                          scape <- unlist(map(seq_along(sid)[-length(sid)], function(j) {
                            MinTippingPath(sid[j], sid[j + 1], hge, je, 10000)}))
                          if(length(scape) == 0){
                            scaledscp <- data.frame(z = zx,  energy = MinTippingPath(sid[[1]], sid[[1]], hge, je, 10000))
                          }else{
                            scaledscp <- data.frame(z = seq(zx[1], zx[length(zx)], length.out = length(scape) - 1), energy = scape[-length(scape)])}
                          
                          leftmargine <- MarginePath(ssid[1], zx[1] - 1, hge, je, 10000)
                          leftmargine <- data.frame(
                            z = seq_along(leftmargine[-1]),
                            energy = rev(leftmargine[-1]))
                          rightmargine <- MarginePath(ssid[length(ssid)], ml + (posit[length(posit)] - zx[length(zx)]), hge, je, 10000)
                          rightmargine <- data.frame(
                            z = seq(zx[length(zx)] + 1, posit[length(posit)] + 1 + ml),
                            energy = rightmargine)
                          rbind(leftmargine, scaledscp, rightmargine)
                        }
  
  cc <- 0
  intpol <- lapply(filledsurf, function(df) {
    cc <<- cc + 1
    spline_fit <- smooth.spline(df$z, df$energy, spar = 0.2)  # �����\�ȕ����p�����[�^ spar ��ݒ�
    pol <- predict(spline_fit, seq(df$z[1], df$z[length(df$z)], length.out = 100))$y
    dx <- seq(1, df$z[length(df$z)], length.out = 100)
    dex <- rep(de[cc], length(pol))#dex <- rep(de0[cc], length(pol))
    data.frame(dx = dx, dex = dex, pol = pol)
  })
  surf3d <- lapply(intpol, function(df) {
    spline_fit <- smooth.spline(df$dx, df$pol, spar = 0.2)
    smoothed_pol <- predict(spline_fit, df$dx)$y
    filtered_pol <- stats::filter(smoothed_pol, filter = rep(1/5, 5), sides = 2, circular = TRUE)
    data.frame(x = df$dx, y = df$dex, z = filtered_pol)
  })
  surf3dx <- dplyr::bind_rows(surf3d)
  
  l3d <- lapply(tssspastb, function(y){
    si <- which(y$stb==1)
    xe <- posit[which(ssid==y[1,2])]
    xval <- surf3d[[1]][,1]
    nearest_x <- which.min(abs(xval - xe))
    data.frame(t(sapply(si, function(x) surf3d[[x]][nearest_x,])))})
  l3dx <- dplyr::bind_rows(lapply(seq_along(l3d), function(i) {
    df <- l3d[[i]]
    df$Index <- i
    df$ssid <- sapply(tssspastb,function(x){unique(x[,2])})[i]
    df <- df[,c(4,5,1,2,3)]
    return(df)}))
  return(list(surf3dx, l3dx))
}
