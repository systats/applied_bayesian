#------------------------------------------------------------------------------#
#------------------           Bayesian Statistics   ---------------------------#
#------------------       Deviance and DIC                   ------------------#
#------------------------------------------------------------------------------#

# Set Working Directory
setwd("c:/R-Programme_sync/Bayes_Kurs/Daten")

# packages
library(rjags)
library(MASS)

# Loading data
load("Bayes_Student_Survey.RData")
# A reduced dataset of Student Survey
#          during the Lecture in Introduction to Political Methodology
#          Winter term 2016/2017
#
# poleff: Political Efficacy (Likert Score based on 7 items)
#             A larger value = A higher level of efficacy
# friend: Number of alteri in friendship network
# poldisc: Number of alteri in political discussion network
# lr.self: Ideological orientation (left right self-placement)
#             1: Left <- -> 11: Right
# univ.election: Vote intention at the next university election
#             1: Yes; 0: other (No and DK)


#
# Task: Calculate Deviance Information Criterion of the following Poisson Model
#

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Poisson regression
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# ------------------------------------------------------------------------------
#    Calculatng deviance I
#

# Calculate LL in JAGS
pois.model.loglike <- "model{
   for (i in 1:N){
       y[i] ~ dpois(lambda[i])
       log(lambda[i]) <- inprod(X[i,],beta[])
       # Log-Likelihood of individual response
       loglike[i] <- y[i]*log(lambda[i]) - lambda[i] - logfact(y[i])
   }

  for (j in 1:J){
    beta[j] ~ dnorm(0,0.001)
  }
  
  LL <- sum(loglike[]) 
}"

write(pois.model.loglike,
      "Bayes_Poisson_loglike_Student_Survey.bug")

y <- dat$poldisc
X <- cbind(1,dat$lr.self)
X <- na.omit(X)                   # Omit missings in X
y <- y[-attributes(X)$na.action]  # Omit the correspnding obs in y

N <- length(y)
J <- dim(X)[2]

jags.data <- list(y=y,N=N,J=J,X=X)
jags.pois <- jags.model(file=
                        "Bayes_Poisson_loglike_Student_Survey.bug",
                      data=jags.data, n.chains=3)

update(jags.pois, 2000)

jags.pois.out <- coda.samples(jags.pois,
                            variable.names=c("beta","LL"),
                            n.iter=5000, thin=1)


LL <- unlist(jags.pois.out[,"LL"])

beta.posterior <- as.matrix(jags.pois.out[,c("beta[1]","beta[2]")])

dev.LL <- LL * -2

dev.bar <- mean(dev.LL)

# calculating deviance hat
mean.posterior <- apply(beta.posterior,2,mean)
lambda.hat <- X %*% mean.posterior
lambda.hat <- exp(lambda.hat)
loglike.hat <- y*log(lambda.hat) - lambda.hat - lfactorial(y)
dev.hat <- sum(loglike.hat * -2)


pD <- dev.bar - dev.hat
dev.bar
pD
dev.bar + pD



# ------------------------------------------------------------------------------
#    Calculatng deviance II (using dpois())
#
beta.posterior <- as.matrix(jags.pois.out[,c("beta[1]","beta[2]")])

# Calculating deviance bar
lambda <- X %*% t(beta.posterior)
lambda <- exp(lambda)
loglike <- lambda
for (i in 1:ncol(lambda)){
     loglike[,i] <- dpois(y,lambda=lambda[,i],log=TRUE)
}
dev.LL2 <- loglike * -2
dev.LL2 <- apply(dev.LL2,2,sum,na.rm=T)

plot(dev.LL,dev.LL2)  # check with deviance

dev.bar <- mean(dev.LL2)

# calculating deviance hat
mean.posterior <- apply(beta.posterior,2,mean)
lambda.hat <- X %*% mean.posterior
lambda.hat <- exp(lambda.hat)
loglike.hat <- dpois(y,lambda=lambda.hat,log=TRUE)
dev.hat <- sum(loglike.hat * -2)


pD <- dev.bar - dev.hat
dev.bar
pD
dev.bar + pD


# ------------------------------------------------------------------------------
#    Calculatng deviance III  (using dic-module)
#
load.module('dic')

jags.pois.out <- coda.samples(jags.pois,
                            variable.names=c("beta","LL","deviance"),
                            n.iter=5000, thin=1)

# Likelihood (y|beta) and deviance
posterior.dev <- unlist(jags.pois.out[,"deviance"])
beta.posterior <- as.matrix(jags.pois.out[,c("beta[1]","beta[2]")])

dev.bar <- mean(posterior.dev)

# calculating deviance hat
mean.posterior <- apply(beta.posterior,2,mean)
lambda.hat <- X %*% mean.posterior
lambda.hat <- exp(lambda.hat)
loglike.hat <- dpois(y,lambda=lambda.hat,log=TRUE)
dev.hat <- sum(loglike.hat * -2)


pD <- dev.bar - dev.hat
dev.bar
pD
dev.bar + pD


# ------------------------------------------------------------------------------
#    Calculatng deviance IV   (using dic.samples())
#
pois.dic.out <- dic.samples(jags.pois,n.iter=2000, thin=1)
pois.dic.out


