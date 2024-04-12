if (interactive()) {
Y = matrix(rnorm(100*4), ncol=4)

A = matrix(rnorm(16), ncol=4)
Y = Y %*% A + rnorm(4)
Ybar = apply(Y, 2, mean)

Sigma = cov(Y)

# this gives recursion but results are stored in reverse order :-(
S = Sigma[4:1, 4:1]


U = chol(solve(S))
U.inv = solve(U)

t(U) %*% U  # is Sigma^-1
solve(S)
(U.inv) %*% t(U.inv) # Sigma
S

betas = diag(rep(1,4))-diag((1/diag(U))) %*% U

alpha = Ybar[4:1] - betas %*% Ybar[4:1]
# Conditional mean 
# E[Y_j | Y <j] = E[Y_j] + S_{j<j} S{<j<j}^{-1}(Y_{<j} - E[Y_{<j}])
#               = E[Y_j] + Beta (Y_{<j} - E[Y_{<j}])
#     alpha     = E[Y_j] - Beta E[Y_{<j}] + Beta Y_{<j}
cbind(alpha, betas)
# verified that is correct for slopes and intercepts
coef(lm(Y[,4] ~ Y[,3:1]))
coef(lm(Y[,3] ~ Y[,2:1]))
coef(lm(Y[,2] ~ Y[,1]))


# Sigma = L^T L
# Sigma^-1 = (L-1 L^-T) = L U   


S = Sigma  # do not permute!  
U.inv = solve(chol(S))
betas = diag(rep(1,4)) - diag(1/diag(U.inv)) %*% t(U.inv)
alpha = Ybar - betas %*% Ybar
# E[Y_j] = Ybar_j + betas(Y_<j - Ybar_<j)

cbind(alpha, betas)
# verified that is correct for slopes and intercepts

coef(lm(Y[,2] ~ Y[,1]))
coef(lm(Y[,3] ~ Y[,1:2]))
coef(lm(Y[,4] ~ Y[,1:3]))


#note this gives the right solution too and in the correct order!!!

# C/FORTRAN Check
#SS: 200 iterations
#0  0.955000  191.000000 159.000000 14.000000 46.000000 
#1  0.805000  0.000000 161.000000 11.000000 16.000000 
#2  0.080000  0.000000 0.000000 16.000000 6.000000 
#3  0.265000  0.000000 0.000000 0.000000 53.000000 

#if (interactive()) {
mean.cov = scan()
0  0.950000  190.000000 137.000000 55.000000 73.000000 
1  0.715000  0.000000 143.000000 39.000000 25.000000 
2  0.315000  0.000000 0.000000 63.000000 34.000000 
3  0.410000  0.000000 0.000000 0.000000 82.000000 

mean.cov = matrix(mean.cov, ncol=6, byrow = T)
ybar = mean.cov[,2]
SS = mean.cov[1:4, 3:6]
m = 200; lambda = 1
Cov = (SS + t(SS))
diag(Cov) = diag(SS)
Cov = (Cov - m*ybar %*% t(ybar))/m + diag(lambda,4)
U = chol(Cov)
U.inv = solve(U)
diag(1, 4) - diag(1/diag(U.inv))%*%t(U.inv)
}


data(Hald)
set.seed(42)
hald.amcmc = bas.lm(Y ~ ., prior = "ZS-null", modelprior = uniform(),
                       data = Hald, method = "AMCMC", burnin.iteration = 5000, 
                       MCMC.iterations = 5000, thin = 1, delta = .01, renormalize = TRUE)

amc.sample = hald.amcmc$sampleprobs > 0
plot(hald.amcmc$sampleprobs[amc.sample], hald.amcmc$postprobs.RN[amc.sample])
abline(0,1)
probs.HT = function(obj, M) {
  pi.inc = 1 - (1 - obj$sampleprobs)^M
  wt = exp(obj$logmarg + log(obj$priorprobs) - log(pi.inc)); 
  wt[!amc.sample] = 0.0
  return(wt/sum(wt))
}
plot(hald.amcmc$postprobs, probs.HT(hald.amcmc, 10000), col=amc.sample+1); abline(0,1)
hist(wt)
ESS = 1/sum(wt^2)


data(UScrime, package="MASS")
#UScrime[,-2] = log(UScrime[,-2])
set.seed(42)
crime.amcmc =  bas.lm(log(y) ~ log(M) + So + log(Ed) + log(Po1) + log(Po2)
                    + log(LF) + log(M.F) + log(Pop) + log(NW) +
                      log(U1) + log(U2) + log(GDP) + log(Ineq) + log(Prob)+
                      log(Time), 
                    data=UScrime, n.models=2^15, prior="BIC", 
                    method = "AMCMC", burnin.iteration = 5000, 
                    MCMC.iterations = 5000, thin = 15, delta = .01, importance.sampling = FALSE, 
                    renormalize = FALSE)

plot(crime.amcmc$sampleprobs, crime.amcmc$postprobs.MCMC)
abline(0,1)

set.seed(42)
crime.ais =  bas.lm(log(y) ~ log(M) + So + log(Ed) + log(Po1) + log(Po2)
                      + log(LF) + log(M.F) + log(Pop) + log(NW) +
                        log(U1) + log(U2) + log(GDP) + log(Ineq) + log(Prob)+
                        log(Time), 
                      data=UScrime, n.models=2^15, prior="BIC", 
                      method = "AMCMC", burnin.iteration = 5000, 
                      MCMC.iterations = 5000, thin = 15, delta = .01, importance.sampling = TRUE)

plot(crime.ais$sampleprobs, crime.ais$postprobs.RN)
abline(0,1)

set.seed(42)
crime.mcmc = bas.lm(log(y) ~ log(M) + So + log(Ed) + log(Po1) + log(Po2)
                    + log(LF) + log(M.F) + log(Pop) + log(NW) +
                      log(U1) + log(U2) + log(GDP) + log(Ineq) + log(Prob)+
                      log(Time), 
                    data=UScrime, n.models=2^15, prior="BIC", 
                    method = "MCMC", burnin.iteration = 5000, 
                    MCMC.iterations = 5000, thin = 15, delta = .01, renormalize = FALSE)

plot(crime.mcmc$postprobs.MCMC, crime.mcmc$postprobs.RN)
abline(0,1)

crime.det = bas.lm(log(y) ~ log(M) + So + log(Ed) + log(Po1) + log(Po2)
                   + log(LF) + log(M.F) + log(Pop) + log(NW) +
                     log(U1) + log(U2) + log(GDP) + log(Ineq) + log(Prob)+
                     log(Time), 
                   data=UScrime, n.models=2^15, prior="BIC", 
                   method = "deterministic")
 
df = data.frame(TRUTH=crime.det$probne0.RN,AIS = crime.ais$probne0.RN, AMCMC=crime.ais$probne0.MCMC, 
                MCMC=crime.mcmc$probne0.MCMC)                  
plot(df)

cor(df)


match.model = function(pop, methods, df.m) {
  m = length(methods)
  for (i in 1:m) {
    sample = methods[[i]] 
    marg.pop = pop$logmarg + log(pop$priorprobs)
    marg.sample = sample$logmarg + log(sample$priorprobs)
    loc = match( marg.sample, marg.pop)
    df.m[loc, i+1] = sample$postprobs
  }
  return(df.m)
}

df.m = data.frame(TRUTH = crime.det$postprobs.RN, AIS=NA, AMCMC=NA, MCMC=NA)
df.mod = match.model(crime.det, list(crime.ais, crime.amcmc, crime.mcmc), df.m)
plot(df.mod)
## tecator  

data(tecator, package="FuncNN")
# Extract data and target
X <- tecator$absorp.fdata$data[1:172,]
fat <- tecator$y$Fat[1:172] 
data = data.frame(fat, X)

subsamp = seq(1,100, by=4)
p = length(subsamp)
b.it = 2000000
mc.it = 2000000
set.seed(42)
tecator.IS = bas.lm(fat ~ ., data=data[, c(subsamp,101)], prior="g-prior", alpha = 100, 
                     method = "AMCMC", 
                     burnin.iterations = b.it, thin=100,
                     MCMC.iterations = mc.it,  lambda = 102, # renormalize = TRUE,
                     delta = 0.01, n.models = 2^19, importance.sampling = TRUE)
plot(tecator.IS$sampleprobs[tecator.IS$sampleprobs>0], tecator.IS$postprobs.RN[tecator.IS$sampleprobs>0])

set.seed(42)
tecator.amc = bas.lm(fat ~ ., data=data[, c(subsamp,101)], prior="g-prior", alpha = 100, 
                     method = "AMCMC", 
                     burnin.iterations = 2000000, thin = 100,
                     MCMC.iterations = 2000000, lambda = 102,
                     delta = 0.01, n.models = 2^19, importance.sampling = FALSE)
plot(tecator.amc$sampleprobs, tecator.amc$postprobs.MCMC)
plot(tecator.amc$postprobs.RN, tecator.amc$postprobs.MCMC)
abline(0,1)

plot(tecator.IS$probne0, tecator.amc$probne0)
set.seed(42)
tecator.mcmc = bas.lm(fat ~ ., data=data[, c(subsamp,101)], prior="g-prior", alpha = 100, method = "MCMC", 
                     burnin.iterations = b.it, thin = 100,
                     MCMC.iterations = mc.it,
                     delta = 0.01, n.models = 2^19)
plot(tecator.mcmc$postprobs.RN, tecator.mcmc$postprobs.MCMC)
abline(0,1)
plot(tecator.mcmc$probne0.MCMC, tecator.IS$probne0)
tecator.det = bas.lm(fat ~ ., data=data[, c(subsamp,101)], prior="g-prior", alpha = 100,
                     method = "deterministic", thin = 100,
                     burnin.iterations = b.it,
                      MCMC.iterations = mc.it,
                      delta = 0.01, n.models = 2^p)
df = data.frame(IS = tecator.IS$probne0, AMC = tecator.amc$probne0.MCMC, 
                MCMC=tecator.mcmc$probne0.MCMC, BAS=tecator.bas$probne0)
plot(df)
cor(df)

df.m = data.frame(TRUTH = tecator.det$postprobs.RN, AIS=NA, AMCMC=NA, MCMC=NA)
df.mod = match.model(tecator.det, list(tecator.IS, tecator.amc, tecator.mcmc), df.m)
plot(df.mod)
cor(df.mod, use = "pair")

bias = function(yhat, ytrue) {
  diff = abs(yhat -ytrue)
  mean(diff, na.rm = TRUE)
}

mse= function(yhat, ytrue) {
  diff = (yhat -ytrue)^2
  mean(diff, na.rm = TRUE)
}
c(bias(df.mod$TRUTH, df.mod$AIS), bias(df.mod$TRUTH, df.mod$AMCMC), bias(df.mod$TRUTH, df.mod$MCMC))
c(mse(df.mod$TRUTH, df.mod$AIS), mse(df.mod$TRUTH, df.mod$AMCMC), mse(df.mod$TRUTH, df.mod$MCMC))
