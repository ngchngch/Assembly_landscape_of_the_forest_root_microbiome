
showGELA3D_mod <- function(gelsobj,title,xlab,ylab,zlab="Energy",new_X=NA
) {
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
  par(mar = c(5, 5, 5, 5))
  persp3D(
    x, y, z, 
    theta = 25, phi = 50, expand = 0.5, facets=NA, ticktype = "detailed",
    alpha = 0.1, 
    main = title,
    xlab = xlab, ylab = ylab, zlab = zlab
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