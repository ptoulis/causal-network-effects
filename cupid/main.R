# Copyright 2013 Panos Toulis, Donald B. Rubin
# Author: Panos Toulis(ptoulis@fas.harvard.edu)
#
# Simple examples demonstrating the functionality of the code.

source("population.R")
# creates a population of 10 units.
pop = new.population(10000, 0.2)
pop.obs = population.obs(pop)
Y = sample.Y(pop)
# compute treatment effect for units in treatment-units control
y.obs = population.Yobs(pop, Y)

# control males w/ control spouses.
warning("Using full population to compute cupid effect.")
males.00 = population.filter(pop, is.type = "M", has.treatment = 0, match.treatment = 0)
# control males w/ treated spouses.
males.01 = population.filter(pop, is.type="M", has.treatment = 0, match.treatment = 1)
print(sprintf("Cupid effect=%.2f", mean(y.obs[males.01]) - mean(y.obs[males.00])))

# Primary females effect
females.1 = population.filter(pop, is.type="F", has.treatment= 1)
females.0 = population.filter(pop, is.type="F", has.treatment = 0)

print(sprintf("Primary female effect=%.2f", mean(y.obs[females.1]) - mean(y.obs[females.0])))
