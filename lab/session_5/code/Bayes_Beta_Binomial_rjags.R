#
# Getting beta-binomial via rjags 
#
 
library(rjags)

# BUGS Modell
binom.model <- "model{
   y~dbin(p,N)
   p ~ dbeta(1,1) # Prior
}"

write(binom.model, "Bayes_Binom_Beta.bug")

# Data 
jags.data <- list(y = 60, N = 100)

# Running JAGS
jags.reg <- jags.model(file="Bayes_Binom_Beta.bug",
          data=jags.data, n.chains=3)

update(jags.reg, 1000)
parameters <- c("p")
jags.out <- coda.samples(jags.reg,
         variable.name= parameters,n.iter=1000, thin=1)


# Simple description of posterior
summary(jags.out)
plot(jags.out)

# Which percentage of posterior p>0.5 ?
p <- unlist(jags.out)

table(p>0.5)

hist(p)
plot(density(p))




