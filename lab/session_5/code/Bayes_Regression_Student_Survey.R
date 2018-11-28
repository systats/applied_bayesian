#------------------------------------------------------------------------------#
#------------------           Bayesian Statistics   ---------------------------#
#------------------           Bivariate Regression  ---------------------------#
#------------------------------------------------------------------------------#

# packages
library(rjags)

# Loading data
load("Data/Bayes_Student_Survey.RData")
# 
# A reduced dataset of Student Panel Survey 
#          during the Lecture in Introduction to Political Methodology
#          Winter term 2016/2017
#          at the University of Konstanz 
#
# poleff: Political Efficacy (Likert Score based on 7 items)
#             A larger value = A higher level of efficacy
# friend: Number of alteri in friendship network
# poldisc: Number of alteri in political discussion network
# lr.self: Ideological orientation (left right self-placement)
#             1: Left <- -> 11: Right
# lr.self.2: Ideological orientation (left right self-placement, second measurement)
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
# Task: Estimate the effect of the size of friendship network on pol. efficacy
#

#------------------------------------------------------------------------------#
# OLS regression
#------------------------------------------------------------------------------#


# Simple plot
par(mfrow=c(1,2))
hist(dat$poleff)
hist(dat$friend)

plot(log(dat$friend+1)  , dat$poleff  , 
     type="p",main="Student Survey",
     xlab="No of friends",ylab="Pol. Efficacy")

plot(jitter(log(dat$friend+1))  , jitter(dat$poleff)  , 
     type="p",main="Student Survey",
     xlab="No of friends",ylab="Pol. Efficacy")

# OLS Regression
ols.out <- lm(poleff ~ log(friend+1), data=dat)
summary(ols.out)

abline(reg=ols.out, col="red")


#------------------------------------------------------------------------------#
# Regression via rjags
#------------------------------------------------------------------------------#

# JAGS Modell
reg.model <- "model{
  for (i in 1:N){
    y[i] ~ dnorm(mu[i],tau)
    mu[i] <- beta0 + beta1 * x[i]
  }
  
  beta0 ~ dnorm(0,0.0001)
  beta1 ~ dnorm(0,0.0001) 
  
  tau ~ dgamma(0.001,0.001) 
  sigma <- 1/sqrt(tau)
}"

write(reg.model, "Bayes_Bivariate_Reg_Student_Survey.bug")

y <- dat$poleff
x <- log(dat$friend+1)
N <- length(y)

jags.data <- list(y=y,x=x,N=N)
# three different intial values for beta1
jags.inits.1 <- list(beta1=323)
jags.inits.2 <- list(beta1=5000)
jags.inits.3 <- list(beta1=-10)

jags.reg <- jags.model(file="Bayes_Bivariate_Reg_Student_Survey.bug",
                       inits=list(jags.inits.1,jags.inits.2,jags.inits.3),
                       data=jags.data, n.chains=3)

update(jags.reg, 2000)
jags.reg.out <- coda.samples(jags.reg,
                             variable.names=c("beta0","beta1","sigma"),
                             n.iter=2000, thin=1)
summary(jags.reg.out)
plot(jags.reg.out)

# for comparison with OLS results
summary(ols.out)

# Plotting posterior
plot(jags.reg.out)


# Covariance of beta0 and beta1??
beta0.vec <- unlist(jags.reg.out[,"beta0"])
beta1.vec <- unlist(jags.reg.out[,"beta1"])

plot(beta0.vec,beta1.vec,type="p")
cor(beta0.vec,beta1.vec)

summary(ols.out,correlation = TRUE)$correlation  # correlation matrix


# Checking convergence and chain auto-correlation
gelman.plot(jags.reg.out)
autocorr.plot(jags.reg.out)


# DIC
jags.reg.dic.out <- dic.samples(jags.reg,
                             variable.names=parameters,
                             n.iter=2000, thin=1)
jags.reg.dic.out



