# Simulation Monte-Carlo : comparaison MC vs DQ pour lois alpha-stables
# Koami AMOU - Septembre 2026


library(stabledist)
library(ggplot2)
library(dplyr)
library(tidyr)
library(writexl)
library(parallel)

set.seed(12345)

DOSSIER_SORTIE <- "resultats_30"
dir.create(DOSSIER_SORTIE, showWarnings = FALSE)
dir.create(file.path(DOSSIER_SORTIE, "figures"), showWarnings = FALSE)
dir.create(file.path(DOSSIER_SORTIE, "tableaux"), showWarnings = FALSE)

N_FIXE <- 10


#  Grille qstable 
cat("Pre-calcul grille qstable\n")
t0 <- Sys.time()

alpha_grid <- seq(1.005, 1.995, by = 0.005)
tau_grid   <- seq(0.50, 0.995, by = 0.005)
n_alpha    <- length(alpha_grid)
n_tau      <- length(tau_grid)

Q_grid <- matrix(NA_real_, n_alpha, n_tau)
for (i in 1:n_alpha) {
  for (j in 1:n_tau) {
    Q_grid[i, j] <- qstable(tau_grid[j], alpha = alpha_grid[i],
                            beta = 0, gamma = 1, delta = 0, pm = 1)
  }
}
cat("  fait en", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")

qstable_fast <- function(tau, alpha) {
  i <- max(1, min(n_alpha - 1, findInterval(alpha, alpha_grid)))
  j <- max(1, min(n_tau - 1, findInterval(tau, tau_grid)))
  a1 <- alpha_grid[i]; a2 <- alpha_grid[i + 1]
  t1 <- tau_grid[j];   t2 <- tau_grid[j + 1]
  wa <- (alpha - a1) / (a2 - a1)
  wt <- (tau - t1) / (t2 - t1)
  (1 - wa) * (1 - wt) * Q_grid[i, j] + (1 - wa) * wt * Q_grid[i, j + 1] +
    wa * (1 - wt) * Q_grid[i + 1, j] + wa * wt * Q_grid[i + 1, j + 1]
}


#  Coefficients b_k 
cache_bk <- new.env(parent = emptyenv())

b_k <- function(k, alpha) {
  key <- paste0(round(alpha, 4), "_", k)
  if (!exists(key, envir = cache_bk)) {
    val <- (-1)^(k - 1) / factorial(k) * gamma(alpha * k) * sin(k * pi * alpha / 2)
    assign(key, val, envir = cache_bk)
  }
  get(key, envir = cache_bk)
}

get_bk_vec <- function(alpha) {
  key <- paste0("vec_", round(alpha, 4))
  if (!exists(key, envir = cache_bk)) {
    assign(key, sapply(1:N_FIXE, function(k) b_k(k, alpha)), envir = cache_bk)
  }
  get(key, envir = cache_bk)
}

P_N_std <- function(q, alpha, N = N_FIXE) {
  if (q <= 0) return(NA_real_)
  bk <- get_bk_vec(alpha)
  sum(bk * q^(-alpha * (1:N))) / pi
}

# Version robuste : balayage log-espacé pour trouver un bracket
q_truncated_std <- function(tau, alpha, N = N_FIXE) {
  target <- 1 - tau
  f <- function(q) P_N_std(q, alpha, N) - target
  
  qs <- 10^seq(-2, 6, by = 0.2)
  fs <- sapply(qs, f)
  
  # On garde les intervalles ou f change de signe
  sign_change <- which(diff(sign(fs)) != 0)
  if (length(sign_change) == 0) return(NA_real_)
  
  # On prend le premier changement de signe (f decroissante, on veut f=0)
  i <- sign_change[1]
  tryCatch(
    uniroot(f, lower = qs[i], upper = qs[i + 1], tol = 1e-8)$root,
    error = function(e) NA_real_
  )
}

S_alpha_N_raw <- function(tau, alpha) {
  q <- q_truncated_std(tau, alpha)
  if (is.na(q)) return(NA_real_)
  bk    <- get_bk_vec(alpha)
  k_all <- 1:N_FIXE
  terms <- bk * q^(1 - alpha * k_all) / (alpha * k_all - 1)
  q + sum(terms) / (pi * (1 - tau))
}


#  NOUVEAUX NIVEAUX TAU 
taus_MC  <- c(0.96, 0.97, 0.98, 0.99)
taus_DQ  <- c(0.99, 0.98, 0.97, 0.96)
taus_all <- sort(unique(c(taus_MC, taus_DQ)))
idx_MC   <- match(taus_MC, taus_all)
idx_DQ   <- match(taus_DQ, taus_all)

alpha_fine   <- seq(1.05, 1.98, by = 0.002)
n_alpha_fine <- length(alpha_fine)

cat("Pre-calcul S_alpha_N...\n")
t0 <- Sys.time()
S_grid <- matrix(NA_real_, n_alpha_fine, length(taus_all))
for (j in 1:n_alpha_fine) {
  for (i in seq_along(taus_all)) {
    S_grid[j, i] <- S_alpha_N_raw(taus_all[i], alpha_fine[j])
  }
}
cat("  fait en", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")

n_na <- sum(is.na(S_grid))
cat("NA dans S_grid :", n_na, "/", length(S_grid), "\n")
if (n_na > 0) {
  na_idx <- which(is.na(S_grid), arr.ind = TRUE)
  cat("  alpha concernes :", range(alpha_fine[na_idx[, 1]]), "\n")
  cat("  taus concernes  :", range(taus_all[na_idx[, 2]]), "\n")
}

S_alpha_N_fast <- function(tau_idx, alpha) {
  if (alpha < alpha_fine[1] || alpha > alpha_fine[n_alpha_fine]) return(NA_real_)
  j <- findInterval(alpha, alpha_fine)
  j <- max(1, min(n_alpha_fine - 1, j))
  w <- (alpha - alpha_fine[j]) / (alpha_fine[j + 1] - alpha_fine[j])
  (1 - w) * S_grid[j, tau_idx] + w * S_grid[j + 1, tau_idx]
}


# Generation alpha-stable 
r_alpha_stable_sym <- function(n, alpha) {
  V <- runif(n, -pi / 2, pi / 2)
  W <- rexp(n, rate = 1)
  sin(alpha * V) / (cos(V))^(1 / alpha) *
    (cos((1 - alpha) * V) / W)^((1 - alpha) / alpha)
}

r_alpha_stable <- function(n, alpha, sigma, mu) {
  mu + sigma * r_alpha_stable_sym(n, alpha)
}


#  Superquantile empirique 
SQ_empirique_tri <- function(X, tau) {
  n <- length(X)
  Xs <- sort(X)
  k <- floor(n * tau)
  lambda <- n * tau - k
  if (k >= n) return(Xs[n])
  if (lambda == 0) {
    return(mean(Xs[(k + 1):n]))
  } else {
    queue <- if (k + 2 <= n) sum(Xs[(k + 2):n]) else 0
    ((1 - lambda) * Xs[k + 1] + queue) / (n * (1 - tau))
  }
}


#  Methode MC 
estimer_MC <- function(X, taus, theta0 = c(1.5, 0, 0)) {
  SQ_emp <- sapply(taus, function(tau) SQ_empirique_tri(X, tau))
  if (any(is.na(SQ_emp))) return(c(NA, NA, NA))
  idx <- match(taus, taus_all)
  
  objectif <- function(theta) {
    alpha <- theta[1]; eta <- theta[2]; mu <- theta[3]
    sigma <- exp(eta)
    if (alpha < 1.05 || alpha > 1.98) return(1e10)
    S_vals <- sapply(idx, function(j) S_alpha_N_fast(j, alpha))
    if (any(is.na(S_vals))) return(1e10)
    sum((SQ_emp - (mu + sigma * S_vals))^2)
  }
  
  res <- tryCatch(
    optim(theta0, objectif, method = "L-BFGS-B",
          lower = c(1.05, -Inf, -Inf),
          upper = c(1.98, Inf, Inf),
          control = list(maxit = 200, factr = 1e7)),
    error = function(e) NULL
  )
  if (is.null(res)) return(c(NA, NA, NA))
  res$par
}


#  Methode DQ 
estimer_DQ <- function(X, taus = c(0.99, 0.98, 0.97, 0.96)) {
  SQ_emp <- sapply(taus, function(tau) SQ_empirique_tri(X, tau))
  if (any(is.na(SQ_emp))) return(c(NA, NA, NA))
  
  den <- SQ_emp[3] - SQ_emp[4]
  if (is.na(den) || abs(den) < 1e-10) return(c(NA, NA, NA))
  
  R_hat <- (SQ_emp[1] - SQ_emp[2]) / den
  idx <- match(taus, taus_all)
  
  R_theo <- function(alpha) {
    if (alpha < 1.05 || alpha > 1.98) return(NA_real_)
    S <- sapply(idx, function(j) S_alpha_N_fast(j, alpha))
    if (any(is.na(S))) return(NA_real_)
    num <- S[1] - S[2]; d <- S[3] - S[4]
    if (abs(d) < 1e-10) return(NA_real_)
    num / d
  }
  
  f <- function(alpha) {
    v <- R_theo(alpha)
    if (is.na(v)) return(NA_real_)
    v - R_hat
  }
  
  alpha_hat <- tryCatch(
    uniroot(f, lower = 1.05 + 1e-6, upper = 1.98 - 1e-6, tol = 1e-4)$root,
    error = function(e) NA_real_
  )
  if (is.na(alpha_hat)) return(c(NA, NA, NA))
  
  S_hat <- sapply(idx, function(j) S_alpha_N_fast(j, alpha_hat))
  if (any(is.na(S_hat))) return(c(NA, NA, NA))
  
  den_sigma <- S_hat[1] - S_hat[2]
  if (abs(den_sigma) < 1e-10) return(c(NA, NA, NA))
  
  sigma_hat <- (SQ_emp[1] - SQ_emp[2]) / den_sigma
  if (is.na(sigma_hat) || sigma_hat <= 0) return(c(NA, NA, NA))
  
  mu_hat <- SQ_emp[1] - sigma_hat * S_hat[1]
  c(alpha_hat, log(sigma_hat), mu_hat)
}


#  Scenarios 
scenarios <- data.frame(
  id     = 1:6,
  n      = c(100, 500, 1000, 500, 500, 500),
  alpha0 = c(1.5, 1.5, 1.5, 1.2, 1.8, 1.5),
  sigma0 = c(1,   1,   1,   1,   1,   2),
  mu0    = c(0,   0,   0,   0,   0,   log(2))
)
M <- 500


simuler_scenario <- function(scenario, M = 500) {
  n      <- scenario$n
  alpha0 <- scenario$alpha0
  sigma0 <- scenario$sigma0
  mu0    <- scenario$mu0
  
  res_MC <- matrix(NA_real_, M, 3)
  res_DQ <- matrix(NA_real_, M, 3)
  t_MC   <- numeric(M)
  t_DQ   <- numeric(M)
  
  for (m in 1:M) {
    X <- r_alpha_stable(n, alpha0, sigma0, mu0)
    t_MC[m] <- system.time({ res_MC[m, ] <- estimer_MC(X, taus_MC) })[3]
    t_DQ[m] <- system.time({ res_DQ[m, ] <- estimer_DQ(X, taus_DQ) })[3]
  }
  
  theta_vrai <- c(alpha0, log(sigma0), mu0)
  
  metriques <- function(res, temps, methode) {
    data.frame(
      Methode   = methode,
      Parametre = c("alpha", "eta", "mu"),
      Biais     = colMeans(res, na.rm = TRUE) - theta_vrai,
      Variance  = apply(res, 2, var, na.rm = TRUE),
      EQM       = (colMeans(res, na.rm = TRUE) - theta_vrai)^2 +
        apply(res, 2, var, na.rm = TRUE),
      Temps     = mean(temps, na.rm = TRUE),
      Taux_NA   = colMeans(is.na(res))
    )
  }
  
  tab <- rbind(metriques(res_MC, t_MC, "MC"),
               metriques(res_DQ, t_DQ, "DQ"))
  
  est <- rbind(
    data.frame(Scenario = scenario$id, Methode = "MC", Rep = 1:M,
               alpha = res_MC[, 1], eta = res_MC[, 2], mu = res_MC[, 3]),
    data.frame(Scenario = scenario$id, Methode = "DQ", Rep = 1:M,
               alpha = res_DQ[, 1], eta = res_DQ[, 2], mu = res_DQ[, 3])
  )
  
  list(tableau = tab, estimations = est)
}


#  Lancement 
n_cores <- max(1, detectCores() - 1)
cat("Simulations sur", n_cores, "coeurs\n")
t_debut <- Sys.time()

if (.Platform$OS.type == "windows") {
  cl <- makeCluster(n_cores)
  clusterSetRNGStream(cl, 12345)
  clusterExport(cl, c(
    "scenarios", "simuler_scenario", "M",
    "r_alpha_stable", "r_alpha_stable_sym", "SQ_empirique_tri",
    "b_k", "get_bk_vec", "cache_bk", "P_N_std", "q_truncated_std",
    "S_alpha_N_fast", "S_grid", "alpha_fine", "n_alpha_fine",
    "qstable_fast", "alpha_grid", "tau_grid", "n_alpha", "n_tau", "Q_grid",
    "N_FIXE", "taus_MC", "taus_DQ", "taus_all", "idx_MC", "idx_DQ",
    "estimer_MC", "estimer_DQ"
  ), envir = environment())
  clusterEvalQ(cl, library(stabledist))
  resultats <- parLapply(cl, 1:nrow(scenarios), function(i) {
    cat("scenario", i, "\n")
    simuler_scenario(scenarios[i, ], M = M)
  })
  stopCluster(cl)
} else {
  resultats <- mclapply(1:nrow(scenarios), function(i) {
    cat("scenario", i, "\n")
    simuler_scenario(scenarios[i, ], M = M)
  }, mc.cores = n_cores)
}

cat("Termine en", round(as.numeric(difftime(Sys.time(), t_debut, units = "mins")), 2), "min\n")


#  Sauvegarde
cat("Sauvegarde des resultats\n")

synthese <- do.call(rbind, lapply(1:nrow(scenarios), function(i) {
  tab <- resultats[[i]]$tableau
  tab$Scenario <- i
  tab
}))

write.csv(synthese,
          file.path(DOSSIER_SORTIE, "tableaux", "synthese_globale.csv"),
          row.names = FALSE)
write_xlsx(synthese,
           file.path(DOSSIER_SORTIE, "tableaux", "synthese_globale.xlsx"))

for (i in 1:nrow(scenarios)) {
  write.csv(synthese[synthese$Scenario == i, ],
            file.path(DOSSIER_SORTIE, "tableaux", sprintf("scenario_%d.csv", i)),
            row.names = FALSE)
}

estimations <- do.call(rbind, lapply(resultats, function(r) r$estimations))
write.csv(estimations,
          file.path(DOSSIER_SORTIE, "tableaux", "estimations_brutes.csv"),
          row.names = FALSE)

temps_tab <- synthese %>%
  group_by(Scenario, Methode) %>%
  summarise(Temps = mean(Temps), .groups = "drop") %>%
  pivot_wider(names_from = Methode, values_from = Temps)

write.csv(temps_tab,
          file.path(DOSSIER_SORTIE, "tableaux", "temps_calcul.csv"),
          row.names = FALSE)


#  Figures 
cat("Sauvegarde des figures\n")

plot_box <- function(est_df, y_var, y_lab, vraie, file_name, scenario_id) {
  p <- ggplot(est_df, aes(x = Methode, y = .data[[y_var]], fill = Methode)) +
    geom_boxplot() +
    geom_hline(yintercept = vraie, color = "red", linetype = "dashed") +
    labs(title = paste0("Estimation de ", y_lab, " - scenario ", scenario_id),
         y = y_lab, x = "Methode") +
    theme_minimal()
  ggsave(file_name, plot = p, width = 6, height = 5, dpi = 300)
}

# Boucle sur les 6 scenarios
for (i in 1:nrow(scenarios)) {
  est_i  <- resultats[[i]]$estimations
  vrai_i <- c(scenarios$alpha0[i],
              log(scenarios$sigma0[i]),
              scenarios$mu0[i])
  
  plot_box(est_i, "alpha", "alpha", vrai_i[1],
           file.path(DOSSIER_SORTIE, "figures",
                     sprintf("boxplots_alpha_scenario_%d.png", i)), i)
  
  plot_box(est_i, "eta", "eta", vrai_i[2],
           file.path(DOSSIER_SORTIE, "figures",
                     sprintf("boxplots_eta_scenario_%d.png", i)), i)
  
  plot_box(est_i, "mu", "mu", vrai_i[3],
           file.path(DOSSIER_SORTIE, "figures",
                     sprintf("boxplots_mu_scenario_%d.png", i)), i)
}

# Figures EQM pour alpha, eta et mu 
# Boucle generique pour eviter la duplication de code
parametres_eqm <- list(
  list(nom = "alpha", label = expression(alpha), fichier = "eqm_alpha.png"),
  list(nom = "eta",   label = expression(eta == log(sigma)), fichier = "eqm_eta.png"),
  list(nom = "mu",    label = expression(mu),  fichier = "eqm_mu.png")
)

for (p in parametres_eqm) {
  df_p <- synthese[synthese$Parametre == p$nom, ]
  g <- ggplot(df_p, aes(x = factor(Scenario), y = EQM, fill = Methode)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = paste0("EQM de ", p$nom, " par scenario"),
         x = "Scenario", y = "EQM") +
    theme_minimal()
  ggsave(file.path(DOSSIER_SORTIE, "figures", p$fichier),
         plot = g, width = 8, height = 5, dpi = 300)
}

# Temps de calcul (inchange)
ggsave(file.path(DOSSIER_SORTIE, "figures", "temps_calcul.png"),
       ggplot(temps_tab %>% pivot_longer(-Scenario, names_to = "Methode",
                                         values_to = "Temps"),
              aes(x = factor(Scenario), y = Temps, fill = Methode)) +
         geom_bar(stat = "identity", position = "dodge") +
         labs(title = "Temps de calcul par scenario",
              x = "Scenario", y = "Temps (s)") +
         theme_minimal(),
       width = 8, height = 5, dpi = 300)