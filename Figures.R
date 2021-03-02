#
#
# FIGURES
# Visualize and analyze the results.
#
#



# =============================================================================
# ----- load data -----

results <- read.csv("Results.csv")



# =============================================================================
# ----- create figures -----

par(mfrow = c(1,3))

# Nestedness vs. SD
plot(nestedness ~ SD, results, las = 1,
              ylab = "Nestedness Temperature", xlab = "Abundance SD",
              pch = 16, col = rgb(0, 0, 0, 0.2))

# Generality vs. SD
plot(generality ~ SD, results, las = 1,
              ylab = "Association Generality", xlab = "Abundance SD",
              pch = 16, col = rgb(0, 0, 0, 0.2))

# Connectance vs. SD
plot(connectance ~ SD, results, las = 1,
              ylab = "Network Connectance", xlab = "Abundance SD",
              pch = 16, col = rgb(0, 0, 0, 0.2))



# =============================================================================
# ----- determine fits -----

# Power law fit for network connectance vs. SD.
lm(log(generality) ~ SD, results) 
# intercept = 4.5798, slope = -0.7002
exp(4.5798) 
# = 97.49489
# Thus, generality = 97.49489 * SD ^ -0.7002
cor(log(results$generality), results$SD)^2
# R^2 = 0.6878211

lm(log(connectance) ~ SD, results)
# intercept = -0.9252, slope = -0.9089
exp(-0.9252)
# = 0.3964521
# Thus, generality = 0.3964521 * SD ^ -0.9089
cor(log(results$connectance), results$SD)^2
# R^2 = 0.865125