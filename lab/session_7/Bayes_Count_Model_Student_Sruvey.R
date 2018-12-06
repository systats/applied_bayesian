#------------------------------------------------------------------------------#
#------------------           Bayesian Statistics   ---------------------------#
#------------------       Regression for Count Variables     ------------------#
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
# male: Dummy Variable for male students
# poleff: Political Efficacy (Likert Score based on 7 items)
#             A larger value = A higher level of efficacy
# friend: Number of alteri in friendship network
# poldisc: Number of alteri in political discussion network
# lr.self: Ideological orientation (left right self-placement)
#             1: Left <- -> 11: Right
# univ.election: Vote intention at the next university election
#             1: Yes; 0: other (No and DK)
# polint: interest at university politics
#             1: not interested at all <- -> 5 strongly interested
# tuition: opinion on the general tuition fee for German universities
#             1: support; 2: reject; 3: indifferent
# acceptable: acceptable level of the tuition fee (in Euro per Semester)
#             (Only those who support the tuition fee or indifferent)
# protest1 - protest6: willingness to participate a protest action
#                      against the general tuition fee
#                      1: yes; 0: no
#     protest1: demonstration in Konstanz 
#     protest2: demonstration in Stuttgart 
#     protest3: giving signature at petitions 
#     protest4: strike 
#     protest5: occupation of university buildings 
#     protest6: legal dispute at courts


#
# Task: Estimate the effect of the left-right ideology on
#                                        the political discussion network
#

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Poisson regression
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#


# Simple plot
hist(dat$poldisc)
hist(dat$lr.self)

plot(dat$poldisc ~ dat$lr.self)
pois.out <- glm(poldisc ~ lr.self, family=poisson,data=dat)
summary(pois.out)

# ----------------------------------------------------------------------------
# From here Bayesian estimation
# ----------------------------------------------------------------------------
# First trial
# ----------------------------------------------------------------------------

pois.model <- "model{
   for (i in 1:N){
       y[i] ~ dpois(lambda[i])
       log(lambda[i]) <- inprod(X[i,],beta[])
   }

  for (j in 1:J){
    beta[j] ~ dnorm(0,0.001)
  }
}"

write(pois.model,
      "Bayes_Poisson_Student_Survey.bug")

y <- dat$poldisc
X <- cbind(1,dat$lr.self)

N <- length(y)
J <- dim(X)[2]

jags.data <- list(y=y,N=N,J=J,X=X)
jags.pois <- jags.model(file=
                        "Bayes_Poisson_Student_Survey.bug",
                      data=jags.data, n.chains=3)

# This does not work....
#   due to missings in X!

# ----------------------------------------------------------------------------
# Second trial: Delete observations with missings!
# ----------------------------------------------------------------------------

pois.model <- "model{
   for (i in 1:N){
       y[i] ~ dpois(lambda[i])
       log(lambda[i]) <- inprod(X[i,],beta[])
   }

  for (j in 1:J){
    beta[j] ~ dnorm(0,0.001)
  }
}"

write(pois.model,
      "Bayes_Poisson_Student_Survey.bug")

y <- dat$poldisc
X <- cbind(1,dat$lr.self)
X <- na.omit(X)                   # Omit missings in X
y <- y[-attributes(X)$na.action]  # Omit the correspnding obs in y

N <- length(y)
J <- dim(X)[2]

jags.data <- list(y=y,N=N,J=J,X=X)
jags.pois <- jags.model(file=
                        "Bayes_Poisson_Student_Survey.bug",
                      data=jags.data, n.chains=3)




update(jags.pois, 2000)

jags.pois.out <- coda.samples(jags.pois,
                            variable.names=c("beta"),
                            n.iter=5000, thin=1)
plot(jags.pois.out)
gelman.plot(jags.pois.out)
summary(jags.pois.out)

# ----------------------------------------------------------------------------
# Third trial: Impute missing values!
# ----------------------------------------------------------------------------

pois.model <- "model{
   for (i in 1:N){
       y[i] ~ dpois(lambda[i])
       log(lambda[i]) <- inprod(X[i,],beta[])
       X[i,2] ~ dunif(1,11)   # Impute by using a uniform distribution
   }

  for (j in 1:J){
    beta[j] ~ dnorm(0,0.001)
  }
}"

write(pois.model,
      "Bayes_Poisson_Student_Survey.bug")

y <- dat$poldisc
X <- cbind(1,dat$lr.self)

N <- length(y)
J <- dim(X)[2]

jags.data <- list(y=y,N=N,J=J,X=X)
jags.pois <- jags.model(file=
                        "Bayes_Poisson_Student_Survey.bug",
                      data=jags.data, n.chains=3)

update(jags.pois, 2000)

jags.pois.out <- coda.samples(jags.pois,
                            variable.names=c("beta"),
                            n.iter=20000, thin=20)
plot(jags.pois.out)
gelman.plot(jags.pois.out)
summary(jags.pois.out)

# For comparison: A single impuration with the mid value
summary(glm(poldisc ~ lr.self, family=poisson,data=dat))

lr.self.imp <- ifelse(is.na(dat$lr.self),6,dat$lr.self)
cbind(lr.self.imp,dat$lr.self)     # check!
summary(glm(poldisc ~ lr.self.imp, family=poisson,data=dat))


#
# Plotting the predicted size of network
#

beta <- as.matrix(jags.pois.out)

# Prediction by using quantiles
x.scale <- range(dat$lr.self,na.rm=T)
x.scale <- seq(x.scale[1],x.scale[2],length.out=100)
X0 <- cbind(1,x.scale)
predicted <- beta %*% t(X0)
predicted <- exp(predicted)
predicted <- apply(predicted,2,quantile,c(.05,.5,.95))


# Making graphics
#pdf("C:/localtexmf/docs/Bayes_Kurs_Count_Student_Survey_Effect%01d.pdf",
 #                                      width=6,height=6,pointsize=15,onefile=F)
for (i in 1:2){
plot(dat$lr.self,jitter(dat$poldis),xlab="Left Right Ideology",
                                                ylab="Size of pol.disc.network")
if (i==2){
   lines(x.scale,predicted[2,],col="red")
   lines(x.scale,predicted[1,],col="red",lty=2)
   lines(x.scale,predicted[3,],col="red",lty=2)
   }
}
#dev.off()



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Negative binomial regression
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

nb.out <- glm.nb(poldisc ~ lr.self,data=dat)
summary(nb.out)

# for comparison
pois.out <- glm(poldisc ~ lr.self, family=poisson,data=dat)
summary(pois.out)


# ------------------------------------------------------------------------------
# JAGS
# ------------------------------------------------------------------------------

nb.model <- "model{
   for (i in 1:N){
       y[i] ~ dpois(lambda2[i])
       lambda2[i] <- lambda[i]* rho[i]
       rho[i] ~ dgamma(r,r)
       log(lambda[i]) <- inprod(X[i,],beta[])
       X[i,2] ~ dunif(1,11)   # Impute by using a uniform distribution
   }

  for (j in 1:J){
    beta[j] ~ dnorm(0,0.001)
  }

 r ~ dgamma(0.0001,0.0001)
}"

write(nb.model,
      "Bayes_NB_Student_Survey_JAGS.bug")

y <- dat$poldis
N <- length(y)
X <- cbind(1,dat$lr.self)
J <- dim(X)[2]
jags.data <- list(y=y,N=N,J=J,X=X)
jags.nb <- jags.model(file=
                        "Bayes_NB_Student_Survey_JAGS.bug",
                      data=jags.data, n.chains=3)

parameters <- c("beta","r")

update(jags.nb, 2000)

jags.nb.out <- coda.samples(jags.nb,
                            variable.names=parameters,
                            n.iter=20000, thin=10)
plot(jags.nb.out)
gelman.plot(jags.nb.out)

summary(jags.nb.out)

# ----------------------------------------------------------------------------
# Negative binomial in another specification
# ----------------------------------------------------------------------------

nb2.model <- "model{
 for (i in 1:N) {
    y[i] ~ dnegbin(p[i],r)
    p[i] <- r/(r+lambda[i])
    log(lambda[i]) <- inprod(X[i,],beta[])
    X[i,2] ~ dunif(1,11)   # Impute by using a uniform distribution
 }

 for (j in 1:J){
    beta[j] ~ dnorm(0,0.001)
 }

 r ~ dgamma(0.001,0.001)
}"

write(nb2.model,"Bayes_NB2_Student_Survey_JAGS.bug")

jags.nb2 <- jags.model(file=
                        "Bayes_NB2_Student_Survey_JAGS.bug",
                      data=jags.data, n.chains=3)

parameters <- c("beta","r")

update(jags.nb2, 2000)

jags.nb2.out <- coda.samples(jags.nb2,
                            variable.names=parameters,
                            n.iter=2000, thin=1)
summary(jags.nb2.out)

# to compare
summary(jags.nb.out)









