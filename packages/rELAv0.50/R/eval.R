### データ生成用関数 ###

#' @export
logmon <- function(y, h, J) {
  # 行列のサイズを取得
  n <- nrow(y)
  m <- ncol(y)

  # 内部の計算を実行
  term <- -h - y %*% J
  result <- 1 / (1 + exp(term))

  return(result)
}

#' @export
OnestepHBS <- function(y, logmo) {
  n <- nrow(y)
  m <- ncol(y)

  # Function to update one row of y
  update_row <- function(row_y, row_logmo) {
    position <- sample(1:m, 1)  # Randomly choose a position
    ne <- ifelse(runif(1) < row_logmo[position], 1, 0)  # Compare with logmo
    row_y[position] <- ne  # Update the chosen position
    return(row_y)
  }

  # Apply the update_row function to each row
  updated_y <- t(mapply(update_row, as.data.frame(t(y)), as.data.frame(t(logmo))))

  return(updated_y)
}

#' @export
HeatBath <- function(steps, nproc, h, J) {
  # 初期状態の設定
  nn <- length(h)
  y <- array(0, dim = c(nproc, nn))

  # 各ステップでの更新
  for (i in 1:steps) {
    logmo <- apply(y, 1, logmon, h = h, J = J)
    y <- OnestepHBS(y, t(logmo))
  }
  colnames(y) <- seq(length(h))
  return(y)
}

### パラメータ生成 ###
#' @export
hb.paramgen <- function(nn){
  h.act <<- runif(nn, min = -2, max = 2)

  rnp <- function() {
    tuples <- t(combn(nn, 2))
    rnp_list <- lapply(seq_len(nrow(tuples)), function(i) {
      if (runif(1) < 0.5) {
        return(rev(tuples[i, ]))
      } else {
        return(tuples[i, ])
      }
    })
    return(do.call(rbind, rnp_list))
  }

  generate_matrix <- function(rnp_list) {
    ma <- matrix(0, nn, nn)
    for (i in seq_len(nrow(rnp_list))) {
      index <- rnp_list[i, ]
      ma[index[1], index[2]] <- ifelse(runif(1) < 1, runif(1, min = -2, max = 2), 0)
    }
    return(ma)
  }

  rnp_list <- rnp()
  ma <- generate_matrix(rnp_list)
  j.act <<- ma + t(ma)
}

### Gradient descent ###
#' @export
GradientDescent <- function(data, maxit) {
  units <- ncol(data)
  tdata <- t(data)
  e <- 10^-8
  datalength <- nrow(data)
  tuples <- as.matrix(expand.grid(rep(list(0:1), units)))
  lt <- nrow(tuples)
  
  semp <- colMeans(data)
  ssemp <- t(data) %*% data / datalength * abs(diag(units) - 1)
  htt <- rep(0, units)
  jtt <- matrix(0, units, units)
  ttt <- 0
  while (ttt < maxit) {
    ht <- htt
    jt <- jtt
    e0 <- apply(tuples, 1, Energy, htt, jtt)
    en <- exp(-e0)
    etot <- sum(en)
    prb <- en / etot
    pt <- prb * tuples
    smod <- colSums(pt)
    ssmod <- (t(pt) %*% tuples) * abs(diag(units) - 1)
    dh <- 0.1 * ((semp + e) - (smod + e))
    dj <- 0.1 * ((ssemp + e) - (ssmod + e))
    htt <- ht + dh
    jtt <- jt + dj
    ttt <- ttt + 1
  }
  
  list(h.est = htt, j.est = jtt)
}