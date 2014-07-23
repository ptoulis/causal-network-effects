# Simple examples demonstrating the functionality of the code.

source("population.R")
# creates a population of 10 units.
pop = new.population(50000, singles.pct = 0.2)

# rerandomization.
pop = population.rerandomize(pop)

# compute treatment effect for units in treatment-units control
y.obs = population.Y.obs(pop)
z = population.treatment(pop)
treated = which(z==1)
control = which(z==0)
print(sprintf("Average effect %.2f", mean(y.obs[treated]) - mean(y.obs[control])))

# get inteference effect
warning("Interference effect: Using unobserved data for that.")
all.types = population.types.com(pop)  # complete types (unobserved)
# get the males in control i.e., Y(0, *)
control.males = intersect(which(z==0), types.males(all.types))
# get their matched females.
match.females = types.female.match(all.types, control.males)
control.match.females = intersect(which(z==0), match.females)
# Get control males with control matched females
# i.e. Y(0, 0)
males.00 = types.male.match(all.types, control.match.females)
males.01 = types.male.match(all.types, setdiff(match.females, control.match.females))

print(sprintf("Interference effect=%.2f", mean(y.obs[males.01]) - mean(y.obs[males.00])))

