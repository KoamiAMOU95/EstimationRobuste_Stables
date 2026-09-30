
# Comparaison des trois nouvelles statistiques de test
# avec les statistiques N1, N2 et N3 de Pitera et al. (2021)

# Ce script evalue la puissance de six statistiques pour tester
# l'adequation d'un echantillon a une loi alpha-stable symetrique :
#   - N1, N2, N3 : statistiques de Pitera (QCV sur intervalles bornes)
#   - N^SQ        : differences de superquantiles
#   - N^BSQ       : differences de superquantiles de Bregman
#   - N^SQCV      : ratios de variance conditionnelle sur la queue
#
# Protocole : H_0 : alpha = 1.5, H_1 : alpha varie de 1.1 a 2.0,
# test bilateral au seuil 5%, coefficients des nouvelles statistiques
# calibres sous H_0 par Monte-Carlo, seuils critiques empiriques.
#
# Remarque : les superquantiles empiriques sont calcules avec la
# ponderation fractionnaire introduite a l'equation (2.30) du chapitre 2.


# Packages necessaires 

suppressPackageStartupMessages({
  library(stabledist)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

set.seed(2024)
dir.create("comparaison_pitera", showWarnings = FALSE)
dir.create("comparaison_pitera/figures", showWarnings = FALSE)
dir.create("comparaison_pitera/tableaux", showWarnings = FALSE)

#  1- Generation de lois alpha-stables symetriques
# Methode de Chambers-Mallows-Stuck (1976)
r_stable_sym <- function(n, alpha) {
  V <- runif(n, -pi/2, pi/2)
  W <- rexp(n, rate = 1)
  sin(alpha * V) / (cos(V))^(1/alpha) *
    (cos((1 - alpha) * V) / W)^((1 - alpha) / alpha)
}
r_stable <- function(n, alpha, sigma = 1, mu = 0) {
  mu + sigma * r_stable_sym(n, alpha)
}

# 2 Fonctionnelles empiriques 

# Variance conditionnelle sur un intervalle quantile (a, b)
QCV_emp <- function(X, a, b) {
  n <- length(X)
  Xs <- sort(X)
  ra <- floor(n * a)
  rb <- floor(n * b)
  if (rb - ra < 2) return(NA_real_)
  var(Xs[(ra + 1):rb])
}

# Superquantile empirique avec ponderation fractionnaire (eq. 2.30)
# SQ_hat = [(1-lambda) X_(k+1) + sum_{j=k+2}^n X_(j)] / [n(1-tau)]
SQ_emp <- function(X, tau) {
  n <- length(X)
  Xs <- sort(X)
  k <- floor(n * tau)
  lambda <- n * tau - k
  if (k >= n) return(NA_real_)
  if (lambda == 0) {
    return(mean(Xs[(k + 1):n]))
  }
  queue <- if (k + 2 <= n) sum(Xs[(k + 2):n]) else 0
  ((1 - lambda) * Xs[k + 1] + queue) / (n * (1 - tau))
}

# Superquantile de Bregman empirique avec gamma(x) = |x|^p / p
# Ponderation fractionnaire appliquee a l'observation situee au seuil.
BSQ_emp <- function(X, tau, p = 1.25) {
  n <- length(X)
  Xs <- sort(X)
  k <- floor(n * tau)
  lambda <- n * tau - k
  if (k >= n) return(NA_real_)
  tail_vals <- Xs[(k + 1):n]
  m <- length(tail_vals)
  if (m < 1) return(NA_real_)
  w <- c(1 - lambda, rep(1, m - 1))
  w <- w / sum(w)
  phi_vals <- sign(tail_vals) * abs(tail_vals)^(p - 1)
  m_phi <- sum(w * phi_vals)
  sign(m_phi) * abs(m_phi)^(1 / (p - 1))
}

# SQCV empirique : variance conditionnelle autour du superquantile.
# Definition formelle introduite au chapitre 4 :
#   SQCV_tau(X) = E[(X - SQ_tau(X))^2 | X >= q_tau].
SQCV_emp <- function(X, tau) {
  n <- length(X)
  Xs <- sort(X)
  r <- floor(n * tau)
  if (r >= n - 1) return(NA_real_)
  tail <- Xs[(r + 1):n]
  if (length(tail) < 2) return(NA_real_)
  sq <- SQ_emp(X, tau)
  if (is.na(sq)) return(NA_real_)
  mean((tail - sq)^2)
}

#  3- Statistiques de Pitera (coefficients publies dans l'article
N1_stat <- function(X) {
  n <- length(X)
  q1 <- QCV_emp(X, 0.05, 0.25)
  q2 <- QCV_emp(X, 0.25, 0.75)
  q3 <- QCV_emp(X, 0.75, 0.95)
  qf <- QCV_emp(X, 0.05, 0.95)
  if (is.na(qf) || qf == 0) return(NA_real_)
  sqrt(n) * (1.00 * q1 - 1.01 * q2 + 1.00 * q3) / qf
}

N2_stat <- function(X) {
  n <- length(X)
  q1 <- QCV_emp(X, 0.005, 0.25)
  q2 <- QCV_emp(X, 0.25, 0.75)
  q3 <- QCV_emp(X, 0.75, 0.995)
  qf <- QCV_emp(X, 0.005, 0.995)
  if (is.na(qf) || qf == 0) return(NA_real_)
  sqrt(n) * (0.60 * q1 - 1.61 * q2 + 0.60 * q3) / qf
}

N3_stat <- function(X) {
  n <- length(X)
  q1 <- QCV_emp(X, 0.005, 0.04)
  q2 <- QCV_emp(X, 0.04, 0.96)
  q3 <- QCV_emp(X, 0.96, 0.995)
  qf <- QCV_emp(X, 0.005, 0.995)
  if (is.na(qf) || qf == 0) return(NA_real_)
  sqrt(n) * (1.15 * q1 - 0.17 * q2 + 1.15 * q3) / qf
}

#  4- Nouvelles statistiques : composantes brutes
# Chaque fonction renvoie un vecteur de trois ratios, invariants par
# translation et par changement d'echelle. Les coefficients d seront
# calibres ensuite.

NSQ_raw <- function(X) {
  s1 <- SQ_emp(X, 0.01); s2 <- SQ_emp(X, 0.05)
  s3 <- SQ_emp(X, 0.25); s4 <- SQ_emp(X, 0.75); s5 <- SQ_emp(X, 0.95)
  den <- s4 - s5
  if (is.na(den) || den == 0) return(c(NA, NA, NA))
  c((s1 - s2) / den, (s2 - s3) / den, (s3 - s4) / den)
}

NBSQ_raw <- function(X, p = 1.25) {
  b1 <- BSQ_emp(X, 0.01, p); b2 <- BSQ_emp(X, 0.05, p)
  b3 <- BSQ_emp(X, 0.25, p); b4 <- BSQ_emp(X, 0.75, p); b5 <- BSQ_emp(X, 0.95, p)
  den <- b4 - b5
  if (is.na(den) || den == 0) return(c(NA, NA, NA))
  c((b1 - b2) / den, (b2 - b3) / den, (b3 - b4) / den)
}

NSQCV_raw <- function(X) {
  s_ref <- SQCV_emp(X, 0.70)
  s1 <- SQCV_emp(X, 0.80)
  s2 <- SQCV_emp(X, 0.90)
  s3 <- SQCV_emp(X, 0.95)
  if (is.na(s_ref) || s_ref == 0) return(c(NA, NA, NA))
  c(s1 / s_ref, s2 / s_ref, s3 / s_ref)
}

 #  5- Calibration des coefficients d sous H0
# On cherche d tel que d^T Delta_0 = 0 et d^T Sigma_0 d = 1, ou
# Delta_0 et Sigma_0 sont la moyenne et la covariance asymptotiques
# du vecteur de ratios bruts sous H0. Convention : d3 = -d1.

calibrer <- function(raw_fn, alpha0 = 1.5, n_cal = 500, M_cal = 1000) {
  V <- matrix(NA_real_, M_cal, 3)
  for (m in 1:M_cal) {
    X <- r_stable(n_cal, alpha0, 1, 0)
    V[m, ] <- raw_fn(X)
  }
  V <- V[complete.cases(V), , drop = FALSE]
  if (nrow(V) < 100) {
    warning("Trop peu d'echantillons valides pour la calibration")
    return(c(1, 0, -1))
  }
  
  mu <- colMeans(V)
  Sigma <- n_cal * cov(V)
  
  if (abs(mu[2]) < 1e-10) return(c(1, 0, -1))
  ratio <- -(mu[1] - mu[3]) / mu[2]
  v <- c(1, ratio, -1)
  var_v <- as.numeric(t(v) %*% Sigma %*% v)
  if (var_v <= 0) return(c(1, 0, -1))
  a <- 1 / sqrt(var_v)
  c(a, a * ratio, -a)
}

# Statistiques completes a partir des composantes brutes et de d
NSQ_stat   <- function(X, d) { raw <- NSQ_raw(X);   if (any(is.na(raw))) NA_real_ else sqrt(length(X)) * sum(d * raw) }
NBSQ_stat  <- function(X, d) { raw <- NBSQ_raw(X);  if (any(is.na(raw))) NA_real_ else sqrt(length(X)) * sum(d * raw) }
NSQCV_stat <- function(X, d) { raw <- NSQCV_raw(X); if (any(is.na(raw))) NA_real_ else sqrt(length(X)) * sum(d * raw) }

#  6- Calibration effective sous H_0
cat(" Calibration des coefficients sous H_0\n")
d_NSQ   <- calibrer(NSQ_raw,   1.5, 500, 1000)
d_NBSQ  <- calibrer(NBSQ_raw,  1.5, 500, 1000)
d_NSQCV <- calibrer(NSQCV_raw, 1.5, 500, 1000)
cat(sprintf("  d[N^SQ]   = (%.4f, %.4f, %.4f)\n", d_NSQ[1],   d_NSQ[2],   d_NSQ[3]))
cat(sprintf("  d[N^BSQ]  = (%.4f, %.4f, %.4f)\n", d_NBSQ[1],  d_NBSQ[2],  d_NBSQ[3]))
cat(sprintf("  d[N^SQCV] = (%.4f, %.4f, %.4f)\n", d_NSQCV[1], d_NSQCV[2], d_NSQCV[3]))
cat("\n")

# Fonction qui calcule les six statistiques sur un echantillon
calculer_stats <- function(X, d_NSQ, d_NBSQ, d_NSQCV) {
  c(N1    = N1_stat(X),
    N2    = N2_stat(X),
    N3    = N3_stat(X),
    NSQ   = NSQ_stat(X, d_NSQ),
    NBSQ  = NBSQ_stat(X, d_NBSQ),
    NSQCV = NSQCV_stat(X, d_NSQCV))
}

# Seuils critiques sous H0 : quantile a 95% de |N| sous H0
seuils_sous_H0 <- function(n, M0 = 500) {
  S <- matrix(NA_real_, M0, 6)
  colnames(S) <- c("N1", "N2", "N3", "NSQ", "NBSQ", "NSQCV")
  for (m in 1:M0) {
    X <- r_stable(n, 1.5, 1, 0)
    S[m, ] <- calculer_stats(X, d_NSQ, d_NBSQ, d_NSQCV)
  }
  sapply(1:6, function(j) quantile(abs(S[, j]), 0.95, na.rm = TRUE))
}

 # 7- Etude de puissance 
n_vec <- c(50, 100, 200, 500)
alpha1_vec <- seq(1.1, 2.0, by = 0.1)
M1 <- 500

cat(" Etude de puissancen")
results <- data.frame()

for (n in n_vec) {
  cat("  Calcul des seuils pour n =", n, "...\n")
  seuils <- seuils_sous_H0(n, M0 = 500)
  
  for (a1 in alpha1_vec) {
    S <- matrix(NA_real_, M1, 6)
    colnames(S) <- c("N1", "N2", "N3", "NSQ", "NBSQ", "NSQCV")
    for (m in 1:M1) {
      X <- r_stable(n, a1, 1, 0)
      S[m, ] <- calculer_stats(X, d_NSQ, d_NBSQ, d_NSQCV)
    }
    for (j in 1:6) {
      pw <- mean(abs(S[, j]) > seuils[j], na.rm = TRUE)
      results <- rbind(results, data.frame(
        n = n, alpha1 = a1, Stat = colnames(S)[j], Puissance = pw))
    }
  }
}

#  8- Sauvegarde des resultats
write.csv(results, "comparaison_pitera/tableaux/puissance_globale.csv",
          row.names = FALSE)

for (n in n_vec) {
  tab <- results %>% filter(n == !!n) %>% select(-n) %>%
    pivot_wider(names_from = Stat, values_from = Puissance) %>% arrange(alpha1)
  write.csv(tab, sprintf("comparaison_pitera/tableaux/puissance_n%d.csv", n),
            row.names = FALSE)
}

#  9- Figures
lvls <- c("N1", "N2", "N3", "NSQ", "NBSQ", "NSQCV")
labs <- c("N1", "N2", "N3", "N^SQ", "N^BSQ", "N^SQCV")
pal  <- c("N1" = "blue", "N2" = "cyan3", "N3" = "deepskyblue4",
          "NSQ" = "red", "NBSQ" = "darkgreen", "NSQCV" = "purple")
results$Stat <- factor(results$Stat, levels = lvls)

# Figure 1 : quatre panneaux (un par taille d'echantillon)
p1 <- ggplot(results, aes(x = alpha1, y = Puissance, color = Stat)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  facet_wrap(~ n, ncol = 2,
             labeller = labeller(n = function(x) paste0("n = ", x))) +
  labs(title = "Comparaison des six statistiques de test",
       subtitle = expression(paste("H"[0], " : ", alpha, " = 1.5, seuil bilateral 5%")),
       x = expression(alpha[1]), y = "Puissance empirique") +
  scale_color_manual(values = pal, labels = labs) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("comparaison_pitera/figures/puissance_par_n.png", p1,
       width = 10, height = 8, dpi = 300)

# Figure 2 : zoom sur n = 500
p2 <- results %>% filter(n == 500) %>%
  ggplot(aes(x = alpha1, y = Puissance, color = Stat)) +
  geom_line(linewidth = 1.3) + geom_point(size = 3) +
  labs(title = "Comparaison des statistiques (n = 500)",
       subtitle = expression(paste("H"[0], " : ", alpha, " = 1.5, seuil 5%")),
       x = expression(alpha[1]), y = "Puissance empirique") +
  scale_color_manual(values = pal, labels = labs) +
  ylim(0, 1) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("comparaison_pitera/figures/puissance_n500.png", p2,
       width = 9, height = 6, dpi = 300)

#  10- Synthese textuelle 
cat("\nSynthese : meilleure statistique par configuration \n")
for (n in n_vec) {
  cat("\nPour n =", n, ":\n")
  tab <- results %>% filter(n == !!n)
  for (a1 in alpha1_vec) {
    best <- tab %>% filter(alpha1 == a1) %>% slice_max(Puissance, n = 1)
    cat(sprintf("  alpha1 = %.1f  ->  %-6s  (puissance = %.3f)\n",
                a1, as.character(best$Stat), best$Puissance))
  }
}

