######################################################
## Correlated particle marginal Metropolis-Hastings ##
######################################################
#     Stochastic-Volatility model:
# 
#       X_t|X_{t-1} ~ N(phi*X_{t-1}, sigma^2)
#       Y_t|X_t     ~ N(0, exp(2(gamma+X_t)))
#
#     Code is an adapted version of PMMH code
#     by Paul Fearnhead: http://www.maths.lancs.ac.uk/~fearnhea/GTP/
#
#     Random-walk on theta=(logit(phi), log(tau), gamma)
#     where tau^2=sigma^2/(1-phi^2)
#
source('CPM_fun.R')


L = 100 # length of HMM

data=SVsim(L,c(0.9, sqrt(1-0.9^2), 1)) # (phi,sigma,gamma); $x - latent , $y - observed


p.m=c(3,0,1);p.sd=c(10,10,10) # prior distribution parameters - uninformative normals

N = 50 # number of particles
M = 5000 # length of MCMC

U = matrix(rnorm(2*N*L), ncol = L) # initial auxiliary variable pool


rho = 0.95 # correlation of U's - setting rho=0 gives PMMH
prop.sd=c(0.5,0.2,0.1) # sd's for proposal distribution (Gaussian)
init.prior = c(0, 1)

chain = matrix(0, nrow=M, ncol = 4)
phi = 0.9
tau = 1
gam = 1
sig = tau*sqrt(1-phi^2)
init.par = c(0,tau)
par = c(phi, sig, gam)
set.seed(100)
chain.cur=c(log(phi/(1-phi)), log(tau), gam, 0)
chain.cur[4] = corr_filterSV(N,data$y, U, par, init.par)$l
for(i in 1:M){
  U.prop = rho*U+sqrt(1-rho^2)*matrix(rnorm(2*N*L), ncol = L) # U'|U --------------- change depending on which vars are unknown
  #chain.prop = c(rnorm(1,chain.cur[1],prop.sd[1]), chain.cur[2], chain.cur[3]) # theta'|theta
  #chain.prop = c(rnorm(1,chain.cur[1],prop.sd[1]), rnorm(1,chain.cur[2],prop.sd[2]), rnorm(1,chain.cur[3], prop.sd[3]))
  chain.prop = c(chain.cur[1], rnorm(1,chain.cur[2],prop.sd[2]), chain.cur[3])
  par.prop = c(exp(chain.prop[1])/(1+exp(chain.prop[1])),exp(chain.prop[2])*sqrt(1-(exp(chain.prop[1])/(1+exp(chain.prop[1])))^2), chain.prop[3])
  init.prop = c(0, exp(2*chain.prop[2])) # X_1 ~ N(0, exp(2*tau))
  Z.prop = corr_filterSV(N, data$y, U.prop, par.prop, init.prop)$l
  
  l.acc = Z.prop-chain.cur[4] + sum(dnorm(chain.prop, p.m, p.sd, log = T))-sum(dnorm(chain.cur[1:3], p.m, p.sd, log = T))
  if(runif(1)<exp(l.acc)){
    chain.cur = c(chain.prop, Z.prop)
    U = U.prop
  }
  chain[i,] = chain.cur
  print(i)
}

par(mfrow =c(2,2))
plot(chain[,1], type = 'l', xlab = 'logit(phi)')
plot(chain[,2], type = 'l', xlab = 'log(tau)')
plot(chain[,3], type = 'l', xlab = 'gamma')
plot(chain[,4], type = 'l', xlab = 'Log-Lik')

par(mfrow =c(1,3)) #transform log(tau) to sigma^2
plot(exp(2*chain[,2])*(1-0.9^2), type = 'l', xlab = 'iteration', ylab = 'sigma^2')
plot(acf(chain[,2], plot = F), type = 'l')
hist(exp(2*chain[,2])*(1-0.9^2))


