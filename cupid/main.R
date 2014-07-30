# Copyright 2013 Panos Toulis, Donald B. Rubin
# Author: Panos Toulis(ptoulis@fas.harvard.edu)
#
# Simple examples demonstrating the functionality of the code.

source("population.R")
# creates a population of 10 units.
pop = new.population(50000, singles.pct = 0.2)

# rerandomization.
pop = population.rerandomize(pop)

# compute treatment effect for units in treatment-units control
y.obs = population.Y(pop, obs=T)

# control males w/ control spouses.
males.00 = population.filter(pop, is.type = "m", has.treatment = 0, match.treatment = 0,
                             obs=F)
# control males w/ treated spouses.
males.01 = population.filter(pop, is.type="m", has.treatment = 0, match.treatment = 1,
                             obs=F)

print(sprintf("Cupid effect=%.2f", mean(y.obs[males.01]) - mean(y.obs[males.00])))

# Primary females effect
females.1 = population.filter(pop, is.type="f", has.treatment=1)
females.0 = population.filter(pop, is.type="f", has.treatment = 0, obs=F)
print(sprintf("Primary female effect=%.2f", mean(y.obs[females.1]) - mean(y.obs[females.0])))
