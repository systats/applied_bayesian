#------------------------------------------------------------------------------#
#------------------           Bayesian Statistics   ---------------------------#
#------------------           Binary Logit Model    ---------------------------#
#------------------------------------------------------------------------------#

# Set Working Directory
setwd("c:/R-Programme_Sync/Bayes_Kurs/Daten")

# packages
library(rjags)

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
# Task: Estimate the effect of pol. efficacy on vote intention
#


# JAGS Modell
logit.model <- "model{
  for (i in 1:N){
    y[i] ~ dbern(p[i])
    logit(p[i]) <- ystar[i]
    ystar[i] <- alpha + beta * x[i]
  }

  alpha ~ dnorm(0,0.0001)
  beta ~ dnorm(0,0.0001)

}"

write(logit.model, "Bayes_Binary_Logit_Student_Survey.bug")

y <- dat$univ.election
x <- dat$poleff
N <- length(y)

jags.data <- list(y=y,x=x,N=N)

jags.logit <- jags.model(file="Bayes_Binary_Logit_Student_Survey.bug",
                       data=jags.data, n.chains=3)

jags.logit.out <- coda.samples(jags.logit,
                             variable.names=c("alpha","beta"),
                             n.iter=2000, thin=1)

# Checking convergence and chain auto-correlation
gelman.plot(jags.logit.out)
autocorr.plot(jags.logit.out)

summary(jags.logit.out)
plot(jags.logit.out)


# for comparison with ML results
ml.logit.out <- glm(univ.election ~ poleff,
                    family=binomial(link="logit"),data=dat)
summary(ml.logit.out)


# Making a prediction with interval for different x values 

# extract coefficients and build a matrix
alpha.vec <- unlist(jags.logit.out[,"alpha"])
beta.vec <- unlist(jags.logit.out[,"beta"])
coef.mat <- cbind(alpha.vec,beta.vec)

# Set up a matrix for x 
range(dat$poleff) # check the range
x.var <- seq(min(dat$poleff),max(dat$poleff),by=0.1)
x.mat <- cbind(1,x.var)

# Prediction
predict.mat <- coef.mat%*% t(x.mat)
predict.mat <- exp(predict.mat)/(1+exp(predict.mat))

# Calculating the percentiles and mean of predictions
predict.ci <- apply(predict.mat,2,quantile,pr=c(0.025,0.5,0.975))
predict.mean <- apply(predict.mat,2,mean)

# Plot everything
plot(x.var,predict.ci[1,],type="l",ylim=c(0,1))
par(new=T)
plot(x.var,predict.ci[2,],type="l",ylim=c(0,1))
par(new=T)
plot(x.var,predict.ci[3,],type="l",ylim=c(0,1))
par(new=T)
plot(x.var,predict.mean,type="l",ylim=c(0,1),col="red")


