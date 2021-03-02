#
#
# SIMULATION
# Simulate neutral interaction networks to create a null model. Record the
# resulting network characteristics while varying the skewedness of the
# species abundance distribution.
#
#



# =============================================================================
# ----- load libraries -----

# Install packages
# install.packages("bipartite")

# Load packages
library(bipartite)



# =============================================================================
# ----- set parameters -----

# Number of species of each type, hosts (I) and visitors (J).
I <- 200
J <- 200

# Population of each type, hosts and visitors.
pop <- 10000

# Run conditions
sd.min <- 1
sd.max <- 10
sd.step <- 0.1
replications <- 20


# =============================================================================
# ----- simulate data -----

# Create a data frame to store outputs.
results <- data.frame("SD" = numeric(),
                      "evenness" = numeric(),
                      "nestedness" = numeric(),
                      "generality" = numeric(),
                      "connectance" = numeric())

# Loop through SD values to manipulate the skewedness of the species abundance
# distribution.
for (SD in seq(from = sd.min, to = sd.max, by = sd.step)) {
  # Repeat each condition
  for (r in 1:replications) {
    # Generate a vector of host species abundances.Then convert abundances to
    # relative abundances by dividing by the sum. Finally, multiply by the
    # population size.
    I.vec <- rlnorm(I, 0, SD)
    I.vec <- round(I.vec/sum(I.vec) * pop)
    
    # Do the same for visitor species.
    J.vec <- rlnorm(J, 0, SD)
    J.vec <- round(J.vec/sum(J.vec) * pop)
    
    # Before we construct the network, we must ensure that the total number of
    # host individuals equals the total number of visitor individuals. To do
    # this, add individuals randomly to species of the smaller of the two types
    # until the types have the same totals.
    while (sum(I.vec) != sum(J.vec)) {
      if(sum(I.vec) < sum(J.vec)) {  # Add a random individual to I.vec
        # Select a random index in I.vec
        index <- sample(1:I, 1)
        I.vec[index] <- I.vec[index] + 1
      } else {  # Add a random individual to J.vec
        # Select a random index in J.vec
        index <- sample(1:J, 1)
        J.vec[index] <- J.vec[index] + 1
      }
    }
    
    # Generate a random network from the species abundances.
    network <- r2dtable(1, I.vec, J.vec)
    network <- as.data.frame(network)
    
    # Calculate host species evenness
    H <- diversity(c(I.vec + J.vec))
    evenness <- H / log(specnumber(c(I.vec + J.vec)))
    
    # Calculate connectance.
    connectance <- sum(sign(network)) / (I * J)
    
    # Calculate nestedness temperature.
    nw.matrix <- as.matrix(network)
    nestedness <- nestedness(nw.matrix, null.models = F)
    generality <- grouplevel(nw.matrix, index = "vulnerability", nrep = 1,
                             level = "lower")
    
    # Record the results
    results[nrow(results) + 1,] <-
      data.frame(SD, evenness, nestedness$temperature, generality, connectance)
  }
}



# =============================================================================
# ----- download data -----

write.csv(results, "Results.csv", row.names = F)