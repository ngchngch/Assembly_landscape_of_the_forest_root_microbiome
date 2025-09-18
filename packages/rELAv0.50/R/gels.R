# Arguments:
#
#'hamming_distance: 
#' @useDynLib rELA, .registration=TRUE
hamming_distance <- function(x, y) {
  sum(x != y)
}

# Arguments:
#
#'gen_concompo: 
#' @useDynLib rELA, .registration=TRUE
gen_concompo <- function(mm){
  concompo <<- list()
  # mmの各要素に関数を適用 ## concompoへの要素追加 (連続する安定状態のグループ化判定)
  lapply(mm, function(mmm) {
  # concompoの中にmmの要素を含むか判定
  matching_indices <- sapply(seq_along(concompo), function(i) {
    any(sapply(mmm, function(x) x %in% unlist(concompo[[i]])))
  })

  # concompoに新しい要素を加えるか、既存の要素に文字列を追加する
  if (!any(matching_indices)) {
    concompo <- c(concompo, list(mmm))
  } else {
    matching_elements <- concompo[matching_indices]
    concompo <- c(concompo[!matching_indices], list(union(unlist(matching_elements), mmm)))
  }})
  return(concompo)
  }

# Arguments:
#
#'con_split: 
#' @useDynLib rELA, .registration=TRUE
con_split <- function(data, factor_variable){group_flag <- rle(factor_variable)$lengths
                        grouped_data <- split(data, rep(seq_along(group_flag), group_flag))
                        return(grouped_data)}

# Arguments:
#
#'ssextend: 
#' @useDynLib rELA, .registration=TRUE
ssextend <- function(spli, de) {
  aaa <- sapply(spli, function(x) x[[1]])
  stb <- ifelse(de %in% aaa, 1, 0)
  data.frame(de, rep(spli[[1]][[3]], length(de)), stb)
}

# Arguments:
#
#'replace_bin: 
#' @useDynLib rELA, .registration=TRUE
#' @export
replace_bin <- function(vec, index) {
  vec[index] <- abs(vec[index] - 1)
  return(vec)
}

# Arguments:
#
#'foldList: 
#' @useDynLib rELA, .registration=TRUE
#' @export
foldList <- function(func, init, vec) {
  result <- list(init)
  for (i in seq_along(vec)) {
    result <- c(result, list(func(result[[length(result)]], vec[i])))
  }
  return(result)
}

# Arguments:
#
#'MinTippingPath: 
#' @useDynLib rELA, .registration=TRUE
#' @export
MinTippingPath <- function(s1, s2, hg, j, tmax){
    n <- length(hg)
    ss1 <- id2bin(s1, n)
    ss2 <- id2bin(s2, n)
    dif <- ss2 - ss1
    pd <- which(abs(dif)==1)
    seq <- sample(pd)
    #############

    if(s1 == s2){return(Energy(ss1, hg, j))}
    if(length(seq) <= 7){
      seq_perms <- permutations(length(seq), length(seq), seq)

      pathes <- apply(seq_perms, 1, function(x) {
          samplepath <- foldList(replace_bin, ss1, x)
          sapply(samplepath, function(y) Energy(y, hg, j))
          })

      tote <- apply(pathes, 2, function(x) {
          dife <-  x[-1] - x[-length(x)]
          sum(dife[which(dife > 0)])})
      pe <- pathes[,which(tote==min(tote))]
      if(is.null(ncol(pe))){minpath <- pe}else{minpath <- rowMeans(pe)}
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
      if (tote < mintote) {
        minpath <- pathe
        mintote <- tote
        pree <- tote
      }
      if (runif(1) < min(1, exp(pree) / exp(tote))) {
        seq <- seqnew
        pree <- tote
      }
    }}
    return(minpath)}


# Arguments:
#
#'MarginePath: 
#' @useDynLib rELA, .registration=TRUE
#' @export
MarginePath <- function(s1, l, hg, j, tmax) {
  n <- length(hg)
  ss1 <- id2bin(s1, n)
  
  pathes <- replicate(min(1024, 2^l), {
    samplepath <- replicate(l, {
        red <-  Reduce(function(vec, i) replace_bin(vec, sample(length(vec), 1)), seq_len(l), init = ss1, accumulate = TRUE)
        sapply(red, function(x) Energy(x, hg, j))
    })
    samplepath
  })
  
  mean_path <- apply(pathes, 1, mean)
  return(mean_path)}

# Arguments:
#
#'GELSObj: 
#' @useDynLib rELA, .registration=TRUE
#' @export
GELSObj <- function(gela, sa, threads = 1){
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
  ss_combi <- combn(sslist, 2)
  mm <- apply(ss_combi, 2, function(pair) {
    dist <- hamming_distance(id2bin(pair[1], s), id2bin(pair[2], s))
    if (dist == 1) return(pair)})
  mm <- Filter(function(x) !is.null(x), mm)
  concompo <- gen_concompo(mm)
  grp <- c(concompo, sapply(setdiff(sslist, unlist(concompo)), function(x) list(x)))
  ssspa <- lapply(grp, function(gg) {
    lapply(sssp, function(ss) {
      ss[grepl(paste(gg, collapse = "|"), ss[3])]
    })})
  ssspa <- lapply(ssspa, function(x) x[sapply(x, function(y) length(y) > 0)]) #空の要素の削除
  ssspa <- sapply(ssspa, function(x) x[order(as.numeric(sapply(x, '[[', 1)))]) #各要素内の並び替え
  ssspa <- ssspa[order(sapply(ssspa, length))] #長さ順の並び替え
  ssspa <- ssspa[order(sapply(ssspa, function(x){as.numeric(x[[1]][[1]])}))] #始点による並び替え
  ssspa_unlist <- unlist(ssspa, recursive = FALSE)
  spl <- con_split(ssspa_unlist, sapply(ssspa_unlist, '[[', 3))
  tssspastb <- lapply(spl, function(x) ssextend(x, de))
  ssid <- sapply(tssspastb, function(x){x[[2]][[1]]})
  hh <- apply(embed(ssid, 2)[, 2:1],1, function(x) hamming_distance(id2bin(x[1], s), id2bin(x[2], s)))
  ml <- max(7, max(hh))
  posit <- Reduce(function(x, y) x + y, hh, accumulate = TRUE, init = ml + 1)

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
          scaledscp <- data.frame(z = zx,  energy = MinTippingPath(sid[[1]], sid[[1]], hge, je, 10000))}
      else{
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
    spline_fit <- smooth.spline(df$z, df$energy, spar = 0.2)  # 調整可能な平滑パラメータ spar を設定
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
    df <- df[,c(4,1,2,3)]
    return(df)}))
  return(list(surf3dx, l3dx))
}