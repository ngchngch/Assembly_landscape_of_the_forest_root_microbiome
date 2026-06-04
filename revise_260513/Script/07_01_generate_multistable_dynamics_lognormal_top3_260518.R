# Set the working directory
setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_260513")

# Load required packages
library(tidyverse)
library(igraph)
library(ggplot2)
library(patchwork)
library(Rcpp)
library(RcppArmadillo)
library(vegan)
library(cluster)
library(parallel)

# Create the output directory
save.dir <- "Analysis/07_01_generate_multistable_dynamics_lognormal_260518_top3"
dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)

set.seed(123)

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

# =============================================================================
# 1. Parameter settings
# =============================================================================
Nsp              <- 80
connectance <- 0.20
A_self      <- -6
mu_A        <- 0
sigma_A     <- 0.1

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
total_time       <- 1000
burn_in          <- 100

# Screening settings
n_init_screen    <- 200
screen_time      <- 200

# Number of final time points used to define endpoint states
endpoint_window_screen <- 6

# Minimum mean silhouette width required for accepting multistability
silhouette_th <- 0.30

# Threshold for detecting state shifts
shift_th         <- 0.25

# effective species number for evalustion 
n_esp_th <- 35

# Window size for state-shift detection
shift_window     <- 20

# Number of candidates to be collected
n_required_candidates <- 3

# Required number of shifts for accepting a candidate
n_required_shifts <- 1

# Number of top candidates used for final visualization and storage
n_top_candidates <- 3

# Parallel computing settings
n_cores <- max(1, parallel::detectCores() - 1)

# Number of trials evaluated in one parallel batch
parallel_batch_size <- n_cores

# =============================================================================
# 2. Helper functions
# =============================================================================
safe_max <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  
  max(x, na.rm = TRUE)
}

# =============================================================================
# 3. Function definitions
# =============================================================================
makeA_multistable <- function(N = Nsp,
                              conn = connectance,
                              mu = mu_A,
                              sigma = sigma_A,
                              a_self = A_self) {
  repeat {
    g <- igraph::sample_gnp(N, p = conn, directed = TRUE)
    A <- as.matrix(igraph::as_adjacency_matrix(g))
    
    link_pos <- A != 0
    A[link_pos] <- rnorm(sum(link_pos), mu, sigma)*rlnorm(sum(link_pos), 0, 2)
    
    diag(A) <- a_self
    
    eig_max <- max(Re(eigen(A, only.values = TRUE)$values))
    
    if (eig_max < 0.2) {
      break
    }
  }
  
  return(A)
}

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

screen_multistability <- function(A, r,
                                  n_init = n_init_screen,
                                  scr_time = screen_time,
                                  ext_th = extinction_th,
                                  noise_sd = noise_sigma,
                                  lnorm_meanlog = noise_lnorm_meanlog,
                                  lnorm_sdlog = noise_lnorm_sdlog,
                                  immi_p = immi_rate,
                                  immi_s = immi_strength,
                                  endpoint_window = endpoint_window_screen,
                                  silhouette_threshold = silhouette_th) {
  
  end_states_arma <- run_screen_rcpp(
    A = A,
    r = r,
    n_init = n_init,
    scr_time = scr_time,
    dt_inner = dt,
    ext_th = ext_th,
    noise_sd = noise_sd,
    lnorm_meanlog = lnorm_meanlog,
    lnorm_sdlog = lnorm_sdlog,
    immi_p = immi_p,
    immi_s = immi_s,
    endpoint_window = endpoint_window
  )
  
  end_states <- as.matrix(end_states_arma)
  
  end_states_nz <- end_states[rowSums(end_states) > 0, , drop = FALSE]
  
  if (nrow(end_states_nz) < 3) {
    return(list(
      score = 0,
      n_clusters = 1,
      end_states = end_states,
      cluster = rep(1, nrow(end_states)),
      mean_between_bc = 0,
      min_between_bc = 0,
      mean_within_bc = NA_real_,
      separation_ratio = 0,
      silhouette_threshold = silhouette_threshold,
      passed_silhouette = FALSE,
      is_multistable = FALSE
    ))
  }
  
  d <- vegan::vegdist(end_states_nz, method = "bray")
  
  if (length(d) == 0 || all(is.na(d))) {
    return(list(
      score = 0,
      n_clusters = 1,
      end_states = end_states_nz,
      cluster = rep(1, nrow(end_states_nz)),
      mean_between_bc = 0,
      min_between_bc = 0,
      mean_within_bc = NA_real_,
      separation_ratio = 0,
      silhouette_threshold = silhouette_threshold,
      passed_silhouette = FALSE,
      is_multistable = FALSE
    ))
  }
  
  hc <- hclust(d, method = "ward.D2")
  
  max_k <- min(10, nrow(end_states_nz) - 1)
  
  if (max_k < 2) {
    return(list(
      score = 0,
      n_clusters = 1,
      end_states = end_states_nz,
      cluster = rep(1, nrow(end_states_nz)),
      mean_between_bc = 0,
      min_between_bc = 0,
      mean_within_bc = NA_real_,
      separation_ratio = 0,
      silhouette_threshold = silhouette_threshold,
      passed_silhouette = FALSE,
      is_multistable = FALSE
    ))
  }
  
  sil_scores <- sapply(2:max_k, function(k) {
    cl <- cutree(hc, k = k)
    
    if (length(unique(cl)) < 2) {
      return(0)
    }
    
    sil <- cluster::silhouette(cl, d)
    mean(sil[, 3], na.rm = TRUE)
  })
  
  best_k <- which.max(sil_scores) + 1
  best_sil <- sil_scores[best_k - 1]
  best_cl <- cutree(hc, k = best_k)
  
  d_mat <- as.matrix(d)
  pair_idx <- which(upper.tri(d_mat), arr.ind = TRUE)
  
  same_cluster <- best_cl[pair_idx[, 1]] == best_cl[pair_idx[, 2]]
  
  between_vals <- d_mat[pair_idx[!same_cluster, , drop = FALSE]]
  within_vals  <- d_mat[pair_idx[same_cluster, , drop = FALSE]]
  
  mean_between_bc <- if (length(between_vals) > 0) {
    mean(between_vals, na.rm = TRUE)
  } else {
    0
  }
  
  min_between_bc <- if (length(between_vals) > 0) {
    min(between_vals, na.rm = TRUE)
  } else {
    0
  }
  
  mean_within_bc <- if (length(within_vals) > 0) {
    mean(within_vals, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  separation_ratio <- if (!is.na(mean_within_bc) && mean_within_bc > 0) {
    mean_between_bc / mean_within_bc
  } else {
    mean_between_bc
  }
  
  passed_silhouette <- is.finite(best_sil) && best_sil >= silhouette_threshold
  is_multistable <- best_k >= 2 && passed_silhouette
  
  return(list(
    score = best_sil,
    n_clusters = best_k,
    end_states = end_states_nz,
    cluster = best_cl,
    mean_between_bc = mean_between_bc,
    min_between_bc = min_between_bc,
    mean_within_bc = mean_within_bc,
    separation_ratio = separation_ratio,
    silhouette_threshold = silhouette_threshold,
    passed_silhouette = passed_silhouette,
    is_multistable = is_multistable
  ))
}

detect_state_shifts <- function(mat_obs,
                                window = shift_window,
                                shift_th = shift_th) {
  n_time <- nrow(mat_obs)
  n_windows <- floor(n_time / window)
  
  if (n_windows < 2) {
    return(list(
      bc_dist = numeric(0),
      max_bc_dist = NA_real_,
      shift_times = numeric(0),
      n_shifts = 0,
      window_states = matrix(NA_real_, nrow = n_windows, ncol = ncol(mat_obs)),
      shift_th = shift_th
    ))
  }
  
  window_states <- matrix(0, nrow = n_windows, ncol = ncol(mat_obs))
  
  for (w in seq_len(n_windows)) {
    idx <- ((w - 1) * window + 1):(w * window)
    window_states[w, ] <- colMeans(mat_obs[idx, , drop = FALSE])
  }
  
  bc_dist <- numeric(n_windows - 1)
  shift_detected <- logical(n_windows - 1)
  
  for (w in seq_len(n_windows - 1)) {
    d <- vegan::vegdist(
      rbind(window_states[w, ], window_states[w + 1, ]),
      method = "bray"
    )[1]
    
    bc_dist[w] <- d
    shift_detected[w] <- is.finite(d) && d > shift_th
  }
  
  shift_times <- which(shift_detected) * window
  
  return(list(
    bc_dist = bc_dist,
    max_bc_dist = safe_max(bc_dist),
    shift_times = shift_times,
    n_shifts = length(shift_times),
    window_states = window_states,
    shift_th = shift_th
  ))
}

# =============================================================================
# 4. Parallel candidate evaluation
# =============================================================================
evaluate_candidate <- function(trial_id,
                               n_required_shifts = n_required_shifts,
                               n_init_screen = n_init_screen,
                               scr_time = screen_time,
                               ext_th = extinction_th,
                               shift_th = shift_th,
                               n_esp_th = 40,
                               window = shift_window,
                               total_t = total_time,
                               burn = burn_in,
                               endpoint_window = endpoint_window_screen,
                               silhouette_threshold = silhouette_th) {
  
  A_cand <- makeA_multistable()
  r_cand <- runif(Nsp)
  
  scr <- screen_multistability(
    A = A_cand,
    r = r_cand,
    n_init = n_init_screen,
    scr_time = scr_time,
    ext_th = ext_th,
    # Screening without noise and immigration
    noise_sd = 0,
    lnorm_meanlog = 0,
    lnorm_sdlog = 2,
    immi_p = 0,
    immi_s = 0,
    endpoint_window = endpoint_window,
    silhouette_threshold = silhouette_threshold
  )
  
  if (!scr$is_multistable) {
    return(list(
      accepted = FALSE,
      trial_id = trial_id,
      reason = "multistability criterion was not met",
      n_clusters = scr$n_clusters,
      screen_score = scr$score,
      silhouette_threshold = scr$silhouette_threshold,
      mean_between_bc = scr$mean_between_bc,
      n_shifts = NA_integer_,
      max_bc_dist = NA_real_
    ))
  }
  
  # Run a long-term simulation
  N0_cand <- runif(Nsp)
  
  sim_cand <- run_long_simulation(
    A = A_cand,
    r = r_cand,
    N0 = N0_cand,
    total_t = total_t,
    burn = burn
  )
  
  mat_obs_cand <- sim_cand$observed
  
  # ---- effective species criterion ----
  n_effective_species <- sum(
    colMeans(mat_obs_cand > 0) > 0.05 &
      colMeans(mat_obs_cand == 0) > 0.05
  )
  
  accepted_by_effective_species <- n_effective_species > n_esp_th
  
  if (!accepted_by_effective_species) {
    return(list(
      accepted = FALSE,
      trial_id = trial_id,
      reason = "effective species criterion was not met",
      n_effective_species = n_effective_species,
      n_clusters = scr$n_clusters,
      screen_score = scr$score,
      silhouette_threshold = scr$silhouette_threshold,
      mean_between_bc = scr$mean_between_bc,
      n_shifts = NA_integer_,
      max_bc_dist = NA_real_
    ))
  }
  
  shifts <- detect_state_shifts(
    mat_obs = mat_obs_cand,
    window = window,
    shift_th = shift_th
  )
  
  accepted_by_shift <- shifts$n_shifts >= n_required_shifts
  
  if (!accepted_by_shift) {
    return(list(
      accepted = FALSE,
      trial_id = trial_id,
      reason = "number of shifts did not meet the criterion",
      n_effective_species = n_effective_species,
      n_clusters = scr$n_clusters,
      screen_score = scr$score,
      silhouette_threshold = scr$silhouette_threshold,
      mean_between_bc = scr$mean_between_bc,
      n_shifts = shifts$n_shifts,
      max_bc_dist = shifts$max_bc_dist
    ))
  }
  
  result <- list(
    accepted = TRUE,
    trial_id = trial_id,
    A = A_cand,
    r = r_cand,
    N0_screen_long = N0_cand,
    sim_screen_long = sim_cand,
    mat_obs_screen_long = mat_obs_cand,
    n_effective_species = n_effective_species,
    end_states = scr$end_states,
    cluster = scr$cluster,
    n_clusters = scr$n_clusters,
    screen_score = scr$score,
    silhouette_threshold = scr$silhouette_threshold,
    passed_silhouette = scr$passed_silhouette,
    is_multistable = scr$is_multistable,
    mean_between_bc = scr$mean_between_bc,
    min_between_bc = scr$min_between_bc,
    mean_within_bc = scr$mean_within_bc,
    separation_ratio = scr$separation_ratio,
    n_shifts = shifts$n_shifts,
    shift_times = shifts$shift_times,
    bc_dist = shifts$bc_dist,
    max_bc_dist = shifts$max_bc_dist,
    window_states = shifts$window_states,
    shift_th = shift_th
  )
  
  return(result)
}
# =============================================================================
# 5. Parameter screening with parallel processing
# =============================================================================
cat("=== Starting parameter screening ===\n")
cat(sprintf(
  "Target number of candidates: %d, number of initial conditions per candidate: %d, dt = %d\n",
  n_required_candidates,
  n_init_screen,
  dt
))
cat(sprintf(
  "Multistability silhouette threshold = %.3f\n",
  silhouette_th
))
cat(sprintf(
  "State-shift threshold shift_th = %.3f, window = %d\n",
  shift_th,
  shift_window
))
cat(sprintf(
  "Parallel settings: n_cores = %d, parallel_batch_size = %d\n\n",
  n_cores,
  parallel_batch_size
))


screen_until_shifts <- function(n_required_shifts = n_required_shifts,
                                n_required_candidates = n_required_candidates,
                                n_init_screen = n_init_screen,
                                scr_time = screen_time,
                                ext_th = extinction_th,
                                shift_th = shift_th,
                                n_esp_th = n_esp_th,
                                window = shift_window,
                                total_t = total_time,
                                burn = burn_in,
                                endpoint_window = endpoint_window_screen,
                                silhouette_threshold = silhouette_th,
                                n_cores = n_cores,
                                batch_size = parallel_batch_size,
                                max_trials = Inf) {
  
  found_candidates <- list()
  trial_id <- 0
  
  cat("Parallel screening settings:\n")
  cat(sprintf("  Number of cores: %d\n", n_cores))
  cat(sprintf("  Parallel batch size: %d\n", batch_size))
  cat(sprintf("  Required accepted candidates: %d\n", n_required_candidates))
  cat("\n")
  
  while (length(found_candidates) < n_required_candidates && trial_id < max_trials) {
    
    remaining_trials <- max_trials - trial_id
    
    current_batch_size <- if (is.finite(max_trials)) {
      min(batch_size, remaining_trials)
    } else {
      batch_size
    }
    
    if (current_batch_size <= 0) {
      break
    }
    
    trial_ids <- seq.int(
      from = trial_id + 1,
      length.out = current_batch_size
    )
    
    trial_id <- max(trial_ids)
    
    cat("------------------------------------------------------------\n")
    cat(sprintf(
      "Running trials %d to %d in parallel | Accepted candidates: %d / %d\n",
      min(trial_ids),
      max(trial_ids),
      length(found_candidates),
      n_required_candidates
    ))
    
    batch_results <- parallel::mclapply(
      X = trial_ids,
      FUN = function(id) {
        evaluate_candidate(
          trial_id = id,
          n_required_shifts = n_required_shifts,
          n_init_screen = n_init_screen,
          scr_time = scr_time,
          ext_th = ext_th,
          shift_th = shift_th,
          n_esp_th = n_esp_th,
          window = window,
          total_t = total_t,
          burn = burn,
          endpoint_window = endpoint_window,
          silhouette_threshold = silhouette_threshold
        )
      },
      mc.cores = n_cores,
      mc.set.seed = TRUE
    )
    
    for (res in batch_results) {
      
      bc_txt <- if (
        !is.null(res$max_bc_dist) &&
        is.finite(res$max_bc_dist)
      ) {
        sprintf("%.3f", res$max_bc_dist)
      } else {
        "NA"
      }
      
      shift_txt <- if (
        !is.null(res$n_shifts) &&
        is.finite(res$n_shifts)
      ) {
        sprintf("%d", res$n_shifts)
      } else {
        "NA"
      }
      
      effsp_txt <- if (
        !is.null(res$n_effective_species) &&
        is.finite(res$n_effective_species)
      ) {
        sprintf("%d", res$n_effective_species)
      } else {
        "NA"
      }
      
      sil_txt <- if (
        !is.null(res$screen_score) &&
        is.finite(res$screen_score)
      ) {
        sprintf("%.3f", res$screen_score)
      } else {
        "NA"
      }
      
      cluster_txt <- if (
        !is.null(res$n_clusters) &&
        is.finite(res$n_clusters)
      ) {
        sprintf("%d", res$n_clusters)
      } else {
        "NA"
      }
      
      if (isTRUE(res$accepted)) {
        
        found_candidates[[length(found_candidates) + 1]] <- res
        
        cat(sprintf(
          paste0(
            "[ACCEPT] ",
            "trial=%d | ",
            "clusters=%s | ",
            "sil=%.3f | ",
            "maxBC=%s | ",
            "shifts=%s | ",
            "effSp=%s | ",
            "accepted=%d/%d\n"
          ),
          res$trial_id,
          cluster_txt,
          res$screen_score,
          bc_txt,
          shift_txt,
          effsp_txt,
          length(found_candidates),
          n_required_candidates
        ))
        
        if (length(found_candidates) >= n_required_candidates) {
          break
        }
        
      } else {
        
        cat(sprintf(
          paste0(
            "[REJECT] ",
            "trial=%d | ",
            "reason=%s | ",
            "clusters=%s | ",
            "sil=%s | ",
            "maxBC=%s | ",
            "shifts=%s | ",
            "effSp=%s\n"
          ),
          res$trial_id,
          res$reason,
          cluster_txt,
          sil_txt,
          bc_txt,
          shift_txt,
          effsp_txt
        ))
      }
    }
  }
  
  if (length(found_candidates) < n_required_candidates) {
    warning(sprintf(
      "max_trials was reached. The number of accepted candidates is %d / %d.",
      length(found_candidates),
      n_required_candidates
    ))
  }
  
  found_candidates <- lapply(found_candidates, function(x) {
    x$accepted <- NULL
    x
  })
  
  return(found_candidates)
}

screened_candidates <- screen_until_shifts(
  n_required_shifts = n_required_shifts,
  n_required_candidates = n_required_candidates,
  n_init_screen = n_init_screen,
  scr_time = screen_time,
  ext_th = extinction_th,
  shift_th = shift_th,
  n_esp_th = n_esp_th,
  window = shift_window,
  total_t = total_time,
  burn = burn_in,
  endpoint_window = endpoint_window_screen,
  silhouette_threshold = silhouette_th,
  n_cores = n_cores,
  batch_size = parallel_batch_size
)

# Check whether any candidates were detected
if (length(screened_candidates) == 0) {
  stop("No parameter sets with detected shifts were found.")
}

# =============================================================================
# 6. Select the top candidates
# =============================================================================
# Extract the maximum number of effective species for each candidate
max_n_esp <- sapply(screened_candidates, function(x) x$n_effective_species)

# Select the top candidates in descending order of max_bc_dist
top_idxes <- order(max_n_esp, decreasing = TRUE)[
  seq_len(min(n_top_candidates, length(screened_candidates)))
]

cat("\n============================================================\n")
cat(sprintf("Processing the top %d candidates in order.\n", length(top_idxes)))
cat(sprintf("Selected candidates: %s\n", paste(top_idxes, collapse = ", ")))
cat("============================================================\n\n")

# List for storing the top candidate results
top_results <- vector("list", length(top_idxes))

# =============================================================================
# 7. Process the top candidates one by one
# =============================================================================
set.seed(12345678)
for (rank_id in seq_along(top_idxes)) {
  
  best_idx <- top_idxes[rank_id]
  best_result <- screened_candidates[[best_idx]]
  
  A_best <- best_result$A
  r_best <- best_result$r
  
  cat("\n============================================================\n")
  cat(sprintf("Rank %d / Top %d\n", rank_id, length(top_idxes)))
  cat(sprintf("Selected candidate: candidate %d\n", best_idx))
  cat(sprintf("n_clusters = %d\n", best_result$n_clusters))
  cat(sprintf("screen_score = %.3f\n", best_result$screen_score))
  cat(sprintf("silhouette_threshold = %.3f\n", best_result$silhouette_threshold))
  cat(sprintf("screen mean_between_bc = %.3f\n", best_result$mean_between_bc))
  cat(sprintf("shift max_bc_dist = %.3f\n", best_result$max_bc_dist))
  cat(sprintf("shift n_shifts = %d\n", best_result$n_shifts))
  cat("============================================================\n\n")
  
  # Create an output directory specific to this candidate
  save.dir.rank <- sprintf("%s/top%d_candidate%d", save.dir, rank_id, best_idx)
  
  if (!dir.exists(save.dir.rank)) {
    dir.create(save.dir.rank, recursive = TRUE)
  }
  
  # =============================================================================
  # 8. Run a long-term time-series simulation
  # =============================================================================
  cat(sprintf("Rank %d: Running a long-term time-series simulation.\n", rank_id))
  
  # Set the initial condition
  N0_long <- runif(Nsp)
  
  # Run the simulation
  sim_long <- run_long_simulation(
    A = A_best,
    r = r_best,
    N0 = N0_long,
    total_t = 3000,
    burn = burn_in
  )
  
  # Extract the observed data after burn-in
  mat_obs <- sim_long$observed
  
  saveRDS(
    mat_obs,
    sprintf("%s/dynamics_rank%d_candidate%d.rds", save.dir.rank, rank_id, best_idx)
  )
  
  cat(sprintf("Length of the observed time series: %d steps\n", nrow(mat_obs)))
  
  # =============================================================================
  # 9. Detect state shifts
  # =============================================================================
  shifts <- detect_state_shifts(
    mat_obs = mat_obs,
    window = shift_window,
    shift_th = shift_th
  )
  
  cat("\n--- State-shift detection results ---\n")
  cat(sprintf("Number of detected shifts: %d\n", shifts$n_shifts))
  cat(sprintf("Maximum Bray-Curtis dissimilarity: %.3f\n", shifts$max_bc_dist))
  cat(sprintf("Shift threshold: %.3f\n", shift_th))
  
  if (shifts$n_shifts > 0) {
    cat(sprintf("Shift time steps: %s\n", paste(shifts$shift_times, collapse = ", ")))
  }
  
  # =============================================================================
  # 10. Visualization
  # =============================================================================
  obs_steps <- seq_len(nrow(mat_obs))
  
  ts_df <- as.data.frame(mat_obs) |>
    setNames(paste0("sp", seq_len(Nsp))) |>
    mutate(time = obs_steps) |>
    pivot_longer(
      -time,
      names_to = "species",
      values_to = "abundance"
    )
  
  # =============================================================================
  # 10-1. Visualize time series of all species
  # =============================================================================
  p1 <- ggplot(ts_df, aes(x = time, y = abundance, color = species)) +
    geom_line(alpha = 0.3, linewidth = 0.3) +
    {
      if (shifts$n_shifts > 0) {
        geom_vline(
          xintercept = shifts$shift_times,
          linetype = "dashed",
          color = "red",
          alpha = 0.6
        )
      }
    } +
    scale_color_viridis_d(option = "plasma") +
    labs(
      title = sprintf(
        "Long-term gLV Dynamics Rank %d Candidate %d",
        rank_id,
        best_idx
      ),
      subtitle = sprintf(
        "Nsp=%d, Shifts=%d, max BC=%.3f, threshold=%.3f",
        Nsp,
        shifts$n_shifts,
        shifts$max_bc_dist,
        shift_th
      ),
      x = "Time step after burn-in",
      y = "Abundance",
      color = "Species"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none")
  
  ggsave(
    filename = sprintf("%s/dynamics_all_rank%d_candidate%d.png", save.dir.rank, rank_id, best_idx),
    plot = p1,
    height = 6,
    width = 8,
    dpi = 150
  )
  
  # =============================================================================
  # 10-2. Visualize the presence/absence heatmap
  # =============================================================================
  sample_sp <- unique(ts_df$species)[1:min(20, Nsp)]
  
  pa_df_sample <- ts_df |>
    filter(species %in% sample_sp) |>
    mutate(present = as.integer(abundance > extinction_th))
  
  p2 <- ggplot(
    pa_df_sample |>
      mutate(present = factor(present, levels = c(0, 1))),
    aes(x = time, y = species, fill = present)
  ) +
    geom_tile() +
    {
      if (shifts$n_shifts > 0) {
        geom_vline(
          xintercept = shifts$shift_times,
          linetype = "dashed",
          color = "red",
          alpha = 0.8
        )
      }
    } +
    scale_fill_manual(
      values = c("0" = "white", "1" = "steelblue"),
      labels = c("Absent", "Present"),
      na.value = "gray80"
    ) +
    labs(
      title = sprintf(
        "Presence/Absence Heatmap Rank %d Candidate %d",
        rank_id,
        best_idx
      ),
      x = "Time step after burn-in",
      y = "Species",
      fill = "State"
    ) +
    theme_bw(base_size = 10)
  
  ggsave(
    filename = sprintf("%s/presence_absence_rank%d_candidate%d.png", save.dir.rank, rank_id, best_idx),
    plot = p2,
    height = 6,
    width = 8,
    dpi = 150
  )
  
  # =============================================================================
  # 10-3. Visualize between-window Bray-Curtis distances
  # =============================================================================
  bc_df <- data.frame(
    window = seq_along(shifts$bc_dist) * shift_window,
    bc_dist = shifts$bc_dist
  )
  
  p3 <- ggplot(bc_df, aes(x = window, y = bc_dist)) +
    geom_line(color = "gray40") +
    geom_point(aes(color = bc_dist > shift_th), size = 1.5) +
    geom_hline(
      yintercept = shift_th,
      linetype = "dashed",
      color = "red"
    ) +
    scale_color_manual(
      values = c("FALSE" = "steelblue", "TRUE" = "red"),
      labels = c("No shift", "Shift")
    ) +
    labs(
      title = sprintf(
        "Between-window Bray-Curtis Distance Rank %d Candidate %d",
        rank_id,
        best_idx
      ),
      subtitle = sprintf("max BC=%.3f, threshold=%.3f", shifts$max_bc_dist, shift_th),
      x = "Time step",
      y = "Bray-Curtis Distance",
      color = "Shift event"
    ) +
    theme_bw(base_size = 10)
  
  ggsave(
    filename = sprintf("%s/timewindow_BC_rank%d_candidate%d.png", save.dir.rank, rank_id, best_idx),
    plot = p3,
    height = 6,
    width = 8,
    dpi = 150
  )
  
  # =============================================================================
  # 10-4. Visualize screening endpoint states
  # =============================================================================
  end_states_df <- as.data.frame(best_result$end_states) |>
    setNames(paste0("sp", seq_len(Nsp)))
  
  end_df <- end_states_df |>
    mutate(
      init_id = seq_len(nrow(end_states_df)),
      cluster = factor(best_result$cluster)
    ) |>
    pivot_longer(
      cols = starts_with("sp"),
      names_to = "species",
      values_to = "abundance"
    ) |>
    filter(grepl("^sp([1-9]|1[0-9]|20)$", species))
  
  p4 <- ggplot(end_df, aes(x = init_id, y = species, fill = abundance)) +
    geom_tile() +
    scale_fill_viridis_c(option = "plasma") +
    labs(
      title = sprintf(
        "Endpoint States Rank %d Candidate %d",
        rank_id,
        best_idx
      ),
      subtitle = sprintf(
        "N_clusters=%d, silhouette=%.3f, threshold=%.3f",
        best_result$n_clusters,
        best_result$screen_score,
        best_result$silhouette_threshold
      ),
      x = "Initial condition ID",
      y = "Species",
      fill = "Abundance"
    ) +
    theme_bw(base_size = 10)
  
  ggsave(
    filename = sprintf("%s/endpoints_rank%d_candidate%d.png", save.dir.rank, rank_id, best_idx),
    plot = p4,
    height = 6,
    width = 8,
    dpi = 150
  )
  
  # =============================================================================
  # 10-5. Combine and save the figures
  # =============================================================================
  g_all <- (p1 / p2 / p3 / p4) +
    plot_layout(heights = c(2, 1, 1, 1))
  
  ggsave(
    filename = sprintf(
      "%s/gLV_multistable_dynamics_rank%d_candidate%d.png",
      save.dir.rank,
      rank_id,
      best_idx
    ),
    plot = g_all,
    width = 6,
    height = 12,
    dpi = 150
  )
  
  # =============================================================================
  # 10-6. Calculate PCoA
  # =============================================================================
  bc_dist_matrix <- vegan::vegdist(mat_obs, method = "bray")
  
  pcoa_result <- cmdscale(bc_dist_matrix, k = 2, eig = TRUE)
  
  eigenvalues <- pcoa_result$eig
  eigenvalues_positive <- pmax(eigenvalues, 0)
  var_explained <- eigenvalues_positive / sum(eigenvalues_positive) * 100
  
  pcoa_df <- data.frame(
    PC1 = pcoa_result$points[, 1],
    PC2 = pcoa_result$points[, 2],
    time = seq_len(nrow(mat_obs))
  )
  
  cat(sprintf("PCoA axis 1 explained variance: %.1f%%\n", var_explained[1]))
  cat(sprintf("PCoA axis 2 explained variance: %.1f%%\n", var_explained[2]))
  
  p_pcoa <- ggplot(pcoa_df, aes(x = PC1, y = PC2, colour = time)) +
    geom_path(linewidth = 0.4, alpha = 0.7) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_colour_viridis_c(
      option = "plasma",
      name = "Time step"
    ) +
    labs(
      title = sprintf(
        "PCoA Rank %d Candidate %d",
        rank_id,
        best_idx
      ),
      subtitle = sprintf(
        "PC1: %.1f%%, PC2: %.1f%% | Detected shifts: %d",
        var_explained[1],
        var_explained[2],
        shifts$n_shifts
      ),
      x = sprintf("PC1 %.1f%%", var_explained[1]),
      y = sprintf("PC2 %.1f%%", var_explained[2])
    ) +
    theme_classic(base_size = 12) +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  print(p_pcoa)
  
  ggsave(
    filename = sprintf("%s/pcoa_bray_curtis_rank%d_candidate%d.png", save.dir.rank, rank_id, best_idx),
    plot = p_pcoa,
    width = 6,
    height = 5,
    dpi = 150
  )
  
  # =============================================================================
  # 11. Save results
  # =============================================================================
  top_results[[rank_id]] <- list(
    rank = rank_id,
    candidate_index = best_idx,
    best_result = best_result,
    A = A_best,
    r = r_best,
    N0_long = N0_long,
    sim_long = sim_long,
    mat_obs = mat_obs,
    shifts = shifts,
    bc_df = bc_df,
    pcoa_df = pcoa_df,
    var_explained = var_explained
  )
  
  saveRDS(
    object = top_results[[rank_id]],
    file = sprintf("%s/result_rank%d_candidate%d.rds", save.dir.rank, rank_id, best_idx)
  )
  
  saveRDS(
    object = best_result,
    file = sprintf("%s/best_result_rank%d_candidate%d.rds", save.dir.rank, rank_id, best_idx)
  )
  
  write.csv(
    bc_df,
    file = sprintf("%s/bray_curtis_timewindow_rank%d_candidate%d.csv", save.dir.rank, rank_id, best_idx),
    row.names = FALSE
  )
  
  write.csv(
    pcoa_df,
    file = sprintf("%s/pcoa_scores_rank%d_candidate%d.csv", save.dir.rank, rank_id, best_idx),
    row.names = FALSE
  )
  
  cat("\nVisualization results were saved.\n")
  cat(sprintf("Output directory for rank %d, candidate %d: %s\n", rank_id, best_idx, save.dir.rank))
}

# =============================================================================
# 12. Save overall results
# =============================================================================
saveRDS(
  object = screened_candidates,
  file = sprintf("%s/screened_candidates.rds", save.dir)
)

saveRDS(
  object = top_results,
  file = sprintf("%s/top_results.rds", save.dir)
)

# Summary table for the top candidates
top_summary <- data.frame(
  rank = seq_along(top_idxes),
  candidate_index = top_idxes,
  screening_max_bc_dist = sapply(top_idxes, function(i) screened_candidates[[i]]$max_bc_dist),
  screening_n_shifts = sapply(top_idxes, function(i) screened_candidates[[i]]$n_shifts),
  n_clusters = sapply(top_idxes, function(i) screened_candidates[[i]]$n_clusters),
  n_effective_species = sapply(top_idxes, function(i) screened_candidates[[i]]$n_effective_species),
  screen_score = sapply(top_idxes, function(i) screened_candidates[[i]]$screen_score),
  silhouette_threshold = sapply(top_idxes, function(i) screened_candidates[[i]]$silhouette_threshold),
  passed_silhouette = sapply(top_idxes, function(i) screened_candidates[[i]]$passed_silhouette),
  is_multistable = sapply(top_idxes, function(i) screened_candidates[[i]]$is_multistable),
  mean_between_bc = sapply(top_idxes, function(i) screened_candidates[[i]]$mean_between_bc),
  min_between_bc = sapply(top_idxes, function(i) screened_candidates[[i]]$min_between_bc),
  mean_within_bc = sapply(top_idxes, function(i) screened_candidates[[i]]$mean_within_bc),
  separation_ratio = sapply(top_idxes, function(i) screened_candidates[[i]]$separation_ratio)
)

write.csv(
  top_summary,
  file = sprintf("%s/top_summary.csv", save.dir),
  row.names = FALSE
)

cat("\n============================================================\n")
cat(sprintf("Processing of all top %d candidates has been completed.\n", length(top_idxes)))
cat(sprintf("Output directory: %s\n", save.dir))
cat("============================================================\n")
