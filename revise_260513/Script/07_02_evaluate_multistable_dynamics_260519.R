source("packages/fun_gLV_simulation_revise_260514.R")
library(Rcpp)
library(RcppArmadillo)
# =============================================================================
# 0. Define Rcpp code
# =============================================================================
cpp_code <- '
#include <RcppArmadillo.h>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec glv_step_rcpp(arma::vec N, arma::vec r, arma::mat A,
                        double dt_inner, double ext_th,
                        double noise_sd,
                        double lnorm_meanlog,
                        double lnorm_sdlog,
                        double immi_p,
                        double immi_s) {

    int n_sp = N.n_elem;
    int n_iter = (int)dt_inner;

    arma::vec dN(n_sp);
    arma::vec noise(n_sp);

    for (int i = 0; i < n_iter; ++i) {
        dN = N % (r + A * N);

        // Generate additive noise that does not depend on biomass
        // However, do not add noise to absent species with N == 0
        arma::vec normal_noise = randn(n_sp) * noise_sd;
        arma::vec lognormal_noise = exp(lnorm_meanlog + lnorm_sdlog * randn(n_sp));

        // Mask to apply noise only to present species
        arma::vec present_mask = conv_to<arma::vec>::from(N > 0.0);

        // Do not scale the noise by biomass N
        noise = (normal_noise % lognormal_noise) % present_mask;

        N += (dN + noise) / dt_inner;

        N.elem(find(N < ext_th)).fill(0.0);

        N = max(N, arma::zeros(n_sp));
    }

    // Immigration
    arma::uvec immi_idx = find(randu(n_sp) < immi_p);

    if (immi_idx.n_elem > 0) {
        N.elem(immi_idx) += immi_s;
    }

    return N;
}

// [[Rcpp::export]]
arma::mat run_screen_rcpp(
    arma::mat A,
    arma::vec r,
    int n_init,
    int scr_time,
    double dt_inner,
    double ext_th,
    double noise_sd,
    double lnorm_meanlog,
    double lnorm_sdlog,
    double immi_p,
    double immi_s,
    int endpoint_window
) {

    int n_sp = A.n_cols;

    endpoint_window = std::min(endpoint_window, scr_time);

    arma::mat end_states(n_init, n_sp);

    for (int i = 0; i < n_init; ++i) {

        // Draw initial abundances independently from U(0, 1) for all species
        arma::vec N0 = randu(n_sp);

        // Randomly retain 70% of species at initialization
        arma::vec mask = conv_to<arma::vec>::from(randu(n_sp) < 0.7);

        N0 = N0 % mask;
        arma::vec N = N0;

        // Accumulator for averaging the final endpoint_window time points
        arma::vec N_sum(n_sp, arma::fill::zeros);
        int n_recorded = 0;

        for (int t = 0; t < scr_time; ++t) {

            N = glv_step_rcpp(
                N, r, A,
                dt_inner,
                ext_th,
                noise_sd,
                lnorm_meanlog,
                lnorm_sdlog,
                immi_p,
                immi_s
            );

            // Record only the final endpoint_window time points
            if (t >= scr_time - endpoint_window) {
                N_sum += N;
                n_recorded++;
            }
        }

        // Store the mean abundance over the final endpoint_window time points
        end_states.row(i) = (N_sum / n_recorded).t();
    }

    return end_states;
}
'

sourceCpp(code = cpp_code)

cat("Rcpp compilation completed.\n")


glv_step <- function(N, r, A,
                     dt_inner = dt,
                     ext_th = extinction_th,
                     noise_sd = noise_sigma,
                     lnorm_meanlog = noise_lnorm_meanlog,
                     lnorm_sdlog = noise_lnorm_sdlog,
                     immi_p = immi_rate,
                     immi_s = immi_strength) {
  result <- glv_step_rcpp(
    N,
    r,
    A,
    dt_inner,
    ext_th,
    noise_sd,
    lnorm_meanlog,
    lnorm_sdlog,
    immi_p,
    immi_s
  )
  
  as.numeric(result)
}

run_long_simulation <- function(A, r, N0,
                                total_t = total_time,
                                burn = burn_in) {
  mat <- matrix(0, nrow = total_t, ncol = length(N0))
  mat[1, ] <- N0
  
  for (t in seq(2, total_t)) {
    mat[t, ] <- glv_step(mat[t - 1, ], r, A)
  }
  
  mat_obs <- mat[(burn + 1):total_t, , drop = FALSE]
  
  return(list(
    full = mat,
    observed = mat_obs
  ))
}

#####################
# Set the working directory
setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_260513")

# Load required packages
library(tidyverse)
library(ggplot2)
library(patchwork)
library(vegan)
library(cluster)
library(parallel)
library(doParallel)
library(foreach)

# Create the output directory
save.dir <- "Analysis/07_02_evaluate_multistable_dynamics_260519"
dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)

set.seed(12345678)

####
# Noise parameters
noise_sigma          <- 0.01
noise_lnorm_meanlog <- 0
noise_lnorm_sdlog   <- 2

# Extinction and immigration parameters
extinction_th    <- 0.01
immi_rate        <- 0.002
immi_strength    <- 0.02

# Time settings
dt               <- 1000

n_cores <- max(1, parallel::detectCores() - 1)
#####
dir_07_01 <- "Analysis/07_01_generate_multistable_dynamics_lognormal_260518_top3"
lf <- list.files(dir_07_01, pattern = "dynamics_rank", full.names = TRUE,recursive = TRUE) 
lf2 <- lf[grep("rds",lf)]

paral <- list.files(dir_07_01, pattern = "best_result", full.names = TRUE,recursive = TRUE) 

lnam <- sapply(lf2, function(x) strsplit(x, "/")[[1]][3])

para_list <- lapply(paral,readRDS)
mat_list <- lapply(lf2,readRDS)
  
for(lm in 1:length(mat_list)){
  cat("\r", lm, "/", length(mat_list))
  cat("binarization\n")
  dir.create(paste0(save.dir, "/", lnam[lm]), showWarnings = FALSE, recursive = TRUE)
  mat0 <- mat_list[[lm]][sample(1:nrow(mat_list[[lm]]),200),]
  
d_bray <- vegdist(mat0,method="bray")

qt <- quantile(mat0, probs = c(0.1,0.9))
bin_th <- seq(qt[1],qt[2], length.out = 30)

cr <- c()
for(i in 1:length(bin_th)){
  cat("\r", i, "/", length(bin_th))
  bth <- bin_th[i]
  mat2 <- mat0
  mat2[mat2<=bth] <- 0
  mat2[mat2>0] <- 1
  d_jac <- vegdist(mat2,method="jaccard")
  cr[i] <- cor(d_bray, d_jac,method = "kendall")
}
rm(i)

max_cr <- which.max(cr)
best_bin_th <- bin_th[max_cr[order(max_cr)[1]]]

mat2 <- mat0
mat2[mat2<=best_bin_th] <- 0
mat2[mat2>0] <- 1

cat(sum(colMeans(mat2)>0.1&colMeans(mat2)<0.9),"\n")
saveRDS(list(ab=mat0, bin_mat=mat2, bin_th=c(th=best_bin_th,
                                             cor=cr[max_cr[order(max_cr)[1]]])),
        sprintf("%s/%s/community_matrix_%s.rds",save.dir,lnam[lm],lnam[lm]))  

cat("keystoneness\n")

set <- para_list[[lm]]
A <- set$A
r_set <- set$r

res = c()
for(rep in 1:200){
  cat("\r", rep, "/200")
    results <- mclapply(1:Nsp,function(i){
      init <- c(runif(Nsp-1),0.5)[order(c(c(1:Nsp)[-i],i))]
      mat_woi = run_long_simulation(total_t=100, A=A[-i,-i], r=r_set[-i],
                                    N0=init[-i], burn=0)$full
      mat_wi = run_long_simulation(total_t=100, A=A, r=r_set,
                                   N0=init, burn=0)$full
      
      abruptness_diffstart = calc_aburptness(before=cbind(mat_woi,0)[,order(c(c(1:Nsp)[-i],i))],
                                             after=mat_wi,
                                             ab_th = best_bin_th)
      names(abruptness_diffstart) <- paste0(names(abruptness_diffstart),"_diffstart")
      
      diff = c(id=i, rep=rep,
                           abruptness_diffstart[1:2])
      return(diff)
    }, mc.cores = n_cores)
    
    res <- rbind(res,do.call(rbind, results))
}

keystoneness = 
  data.frame(res,r=r_set) |>
  left_join(
    calc_centrality(A),
    by="id"
  )

saveRDS(keystoneness , sprintf("%s/%s/keystonness_vs_index.rds",save.dir,lnam[lm]))

###

mBC2 <- foreach(i = 1:Nsp, .combine="c")%do% {
  return(mean(keystoneness[keystoneness$id == i,
                          "BC_diffstart"],na.rm=TRUE))
}

mJC2 <- foreach(i = 1:Nsp, .combine="c")%do% {
  return(mean(keystoneness[keystoneness$id == i,
                          "Jaccard_diffstart"],na.rm=TRUE))
}

key_parm <- unique(keystoneness[,-c(2:4)])

key_res <- data.frame(connectance=0.4,
                      mBC_diffstart=mBC2,
                      mJC_diffstart=mJC2,
                      key_parm
)
saveRDS(key_res, sprintf("%s/%s/keystoneness_summary.rds",save.dir,lnam[lm]))

}
