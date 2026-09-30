
# REPRODUCTION DES FIGURES DE PITERA ET AL. (2021)
# Goodness-of-fit pour les lois alpha-stables basé sur la QCV
# Auteur : Koami AMOU
# Date : Septembre 2026


# 1 Chargement des bibliothèques
library(stabledist)
library(ggplot2)
library(gridExtra)
library(dplyr)

# 2 Fonction QCV empirique 
qcv_empirique <- function(x, a, b) {
  q_a <- quantile(x, probs = a, na.rm = TRUE)
  q_b <- quantile(x, probs = b, na.rm = TRUE)
  subset <- x[x >= q_a & x <= q_b]
  return(var(subset, na.rm = TRUE))
}

# 3 Fonctions pour les statistiques N1, N2, N3 
calcul_N1 <- function(x) {
  a1 <- 0.05; a2 <- 0.25; a3 <- 0.75; a4 <- 0.95
  sigma12 <- qcv_empirique(x, a1, a2)
  sigma23 <- qcv_empirique(x, a2, a3)
  sigma34 <- qcv_empirique(x, a3, a4)
  sigma14 <- qcv_empirique(x, a1, a4)
  d1 <- 1.00; d2 <- -1.01; d3 <- 1.00
  return((d1 * sigma12 + d2 * sigma23 + d3 * sigma34) / sigma14)
}

calcul_N2 <- function(x) {
  a1 <- 0.005; a2 <- 0.25; a3 <- 0.75; a4 <- 0.995
  sigma12 <- qcv_empirique(x, a1, a2)
  sigma23 <- qcv_empirique(x, a2, a3)
  sigma34 <- qcv_empirique(x, a3, a4)
  sigma14 <- qcv_empirique(x, a1, a4)
  d1 <- 0.60; d2 <- -1.61; d3 <- 0.60
  return((d1 * sigma12 + d2 * sigma23 + d3 * sigma34) / sigma14)
}

calcul_N3 <- function(x) {
  a1 <- 0.005; a2 <- 0.04; a3 <- 0.96; a4 <- 0.995
  sigma12 <- qcv_empirique(x, a1, a2)
  sigma23 <- qcv_empirique(x, a2, a3)
  sigma34 <- qcv_empirique(x, a3, a4)
  sigma14 <- qcv_empirique(x, a1, a4)
  d1 <- 1.15; d2 <- -0.17; d3 <- 1.15
  return((d1 * sigma12 + d2 * sigma23 + d3 * sigma34) / sigma14)
}

# FIGURE 1 : DYNAMIQUE DE LA QCV


set.seed(12345)
n_sim <- 100000



# Panneau gauche : QCV en fonction de alpha

alpha_seq <- seq(1.1, 2, length.out = 30)
beta_seq <- c(-0.8, -0.4, 0, 0.4, 0.8)
a <- 0.2
b <- 0.8

resultats_gauche <- data.frame()
for (beta in beta_seq) {
  for (alpha in alpha_seq) {
    x <- rstable(n_sim, alpha, beta, 1, 0)
    qcv <- qcv_empirique(x, a, b)
    resultats_gauche <- rbind(resultats_gauche,
                              data.frame(alpha = alpha, beta = beta, qcv = qcv))
  }
}

p_gauche <- ggplot(resultats_gauche, aes(x = alpha, y = qcv, color = as.factor(beta))) +
  geom_line(size = 1.2) +
  labs(title = "QCV en fonction de alpha",
       x = expression(alpha), y = "QCV", color = expression(beta)) +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))

# Panneau milieu : QCV en fonction de beta

beta_seq2 <- seq(-1, 1, length.out = 30)
alpha_values <- c(1.2, 1.5, 1.8)

resultats_milieu <- data.frame()
for (alpha in alpha_values) {
  for (beta in beta_seq2) {
    x <- rstable(n_sim, alpha, beta, 1, 0)
    qcv <- qcv_empirique(x, a, b)
    resultats_milieu <- rbind(resultats_milieu,
                              data.frame(alpha = alpha, beta = beta, qcv = qcv))
  }
}

p_milieu <- ggplot(resultats_milieu, aes(x = beta, y = qcv, color = as.factor(alpha))) +
  geom_line(size = 1.2) +
  labs(title = "QCV en fonction de beta",
       x = expression(beta), y = "QCV", color = expression(alpha)) +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))

# Panneau droite : QCV en fonction de l'intervalle (a, 1-a)

a_seq <- seq(0.01, 0.4, length.out = 20)
alpha_fixe <- 1.5
beta_fixe <- 0

resultats_droite <- data.frame()
for (a_val in a_seq) {
  b_val <- 1 - a_val
  x <- rstable(n_sim, alpha_fixe, beta_fixe, 1, 0)
  qcv <- qcv_empirique(x, a_val, b_val)
  resultats_droite <- rbind(resultats_droite,
                            data.frame(a = a_val, qcv = qcv))
}

p_droite <- ggplot(resultats_droite, aes(x = a, y = qcv)) +
  geom_line(size = 1.2, color = "blue") +
  labs(title = "QCV en fonction de (a, 1-a)", x = "a", y = "QCV") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Assemblage et sauvegarde
pdf("figure1_qcv_dynamique.pdf", width = 12, height = 4)
grid.arrange(p_gauche, p_milieu, p_droite, ncol = 3,
             top = "Figure 1 - Dynamique de la QCV")
dev.off()
cat("Figure 1 enregistrée : figure1_qcv_dynamique.pdf\n")


# FIGURE 2 : VALEUR DE N(alpha) EN FONCTION DE alpha


cat("Génération de la Figure 2...\n")

n_sim_alpha <- 50000
alpha_seq2 <- seq(1.0, 2.0, by = 0.05)
resultats_N <- data.frame()

for (alpha in alpha_seq2) {
  x <- rstable(n_sim_alpha, alpha, 0, 1, 0)
  N1 <- calcul_N1(x)
  N2 <- calcul_N2(x)
  N3 <- calcul_N3(x)
  resultats_N <- rbind(resultats_N,
                       data.frame(alpha = alpha, N1 = N1, N2 = N2, N3 = N3))
}

p_fig2 <- ggplot(resultats_N, aes(x = alpha)) +
  geom_line(aes(y = N1, color = "N1"), size = 1.2) +
  geom_line(aes(y = N2, color = "N2"), size = 1.2) +
  geom_line(aes(y = N3, color = "N3"), size = 1.2) +
  scale_color_manual(values = c("N1" = "blue", "N2" = "red", "N3" = "green")) +
  labs(title = "Figure 2 - Valeur de N(alpha) en fonction de alpha",
       x = expression(alpha), y = "N(alpha)", color = "Statistique") +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("figure2_n_alpha.pdf", p_fig2, width = 8, height = 6)
cat("Figure 2 enregistrée : figure2_n_alpha.pdf\n")


# FIGURE 3 : PUISSANCE DES TESTS POUR n = 500




# Fonctions de puissance

calcul_puissance_N1 <- function(n, alpha1, M0, M1) {
  stats_H0 <- numeric(M0)
  for (i in 1:M0) {
    x <- rstable(n, 1.5, 0, 1, 0)
    stats_H0[i] <- calcul_N1(x)
  }
  seuil <- quantile(stats_H0, probs = 0.95, na.rm = TRUE)
  rejet <- 0
  for (i in 1:M1) {
    x <- rstable(n, alpha1, 0, 1, 0)
    stat <- calcul_N1(x)
    if (stat > seuil) rejet <- rejet + 1
  }
  return(rejet / M1)
}

calcul_puissance_N2 <- function(n, alpha1, M0, M1) {
  stats_H0 <- numeric(M0)
  for (i in 1:M0) {
    x <- rstable(n, 1.5, 0, 1, 0)
    stats_H0[i] <- calcul_N2(x)
  }
  seuil <- quantile(stats_H0, probs = 0.95, na.rm = TRUE)
  rejet <- 0
  for (i in 1:M1) {
    x <- rstable(n, alpha1, 0, 1, 0)
    stat <- calcul_N2(x)
    if (stat > seuil) rejet <- rejet + 1
  }
  return(rejet / M1)
}

calcul_puissance_N3 <- function(n, alpha1, M0, M1) {
  stats_H0 <- numeric(M0)
  for (i in 1:M0) {
    x <- rstable(n, 1.5, 0, 1, 0)
    stats_H0[i] <- calcul_N3(x)
  }
  seuil <- quantile(stats_H0, probs = 0.95, na.rm = TRUE)
  rejet <- 0
  for (i in 1:M1) {
    x <- rstable(n, alpha1, 0, 1, 0)
    stat <- calcul_N3(x)
    if (stat > seuil) rejet <- rejet + 1
  }
  return(rejet / M1)
}

# Calcul pour n = 500

n_exemple <- 500
M0 <- 1000
M1 <- 1000
alpha1_seq <- seq(1.1, 2.0, by = 0.1)

resultats_puissance <- data.frame()
for (alpha1 in alpha1_seq) {
  puiss_N1 <- calcul_puissance_N1(n_exemple, alpha1, M0, M1)
  puiss_N2 <- calcul_puissance_N2(n_exemple, alpha1, M0, M1)
  puiss_N3 <- calcul_puissance_N3(n_exemple, alpha1, M0, M1)
  resultats_puissance <- rbind(resultats_puissance,
                               data.frame(alpha1 = alpha1,
                                          N1 = puiss_N1,
                                          N2 = puiss_N2,
                                          N3 = puiss_N3))
}

p_fig3 <- ggplot(resultats_puissance, aes(x = alpha1)) +
  geom_line(aes(y = N1, color = "N1"), size = 1.2) +
  geom_line(aes(y = N2, color = "N2"), size = 1.2) +
  geom_line(aes(y = N3, color = "N3"), size = 1.2) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("N1" = "blue", "N2" = "red", "N3" = "green")) +
  labs(title = paste("Figure 3 - Puissance des tests pour n =", n_exemple),
       x = expression(alpha[1]), y = "Puissance", color = "Statistique") +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("figure3_puissance.pdf", p_fig3, width = 8, height = 6)
cat("Figure 3 enregistrée : figure3_puissance.pdf\n")


# FIGURE 4 : PUISSANCE POUR ALPHA PROCHE DE 2 (n = 2000)




n_large <- 2000
M0_large <- 1000
M1_large <- 1000
alpha_proche <- seq(1.80, 2.00, by = 0.02)

resultats_proche <- data.frame()
for (alpha1 in alpha_proche) {
  puiss_N2 <- calcul_puissance_N2(n_large, alpha1, M0_large, M1_large)
  puiss_N3 <- calcul_puissance_N3(n_large, alpha1, M0_large, M1_large)
  resultats_proche <- rbind(resultats_proche,
                            data.frame(alpha1 = alpha1,
                                       N2 = puiss_N2,
                                       N3 = puiss_N3))
}

p_fig4 <- ggplot(resultats_proche, aes(x = alpha1)) +
  geom_line(aes(y = N2, color = "N2"), size = 1.2) +
  geom_line(aes(y = N3, color = "N3"), size = 1.2) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("N2" = "red", "N3" = "green")) +
  labs(title = paste("Figure 4 - Puissance des tests pour alpha proche de 2, n =", n_large),
       x = expression(alpha[1]), y = "Puissance", color = "Statistique") +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("figure4_puissance_proche2.pdf", p_fig4, width = 8, height = 6)
