library(cluster)
library(vegan)




##################
save.dir <- "03_02_00_prep_spatial_block"
dir.create(save.dir, showWarnings = FALSE)

info <- readRDS(sprintf("%s/comp_sample_info_plant2.rds",ELA_prep_dir))

d <- dist(info[,c("x_m","y_m")])

k_range　<- 2:15
sil_width <- numeric(length(k_range))

set.seed(123)

for (i in seq_along(k_range)) {
  k <- k_range[i]
  
  km <- kmeans(info[,c("x_m","y_m")], centers = k, nstart = 50)
  
  sil <- silhouette(km$cluster, d)
  
  sil_width[i] <- mean(sil[, "sil_width"])
}

pdf(sprintf("%s/silhouette_width_plot.pdf", save.dir), width = 6, height = 5)
plot(k_range, sil_width, type = "b", xlab = "Number of clusters K", ylab = "Average silhouette width")
dev.off()

km_best <- kmeans(info[,c("x_m","y_m")], centers = k_range[which.max(sil_width)], nstart = 50)

print(table(km_best$cluster))

info$sp_block <- km_best$cluster

saveRDS(info, sprintf("%s/comp_sample_info_plant2_spblock.rds", save.dir))