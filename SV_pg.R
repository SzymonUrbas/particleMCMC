########################################################
## Particle Gibbs with and without ancestral sampling ##
########################################################

#     Stochastic-Volatility model:
# 
#       X_t|X_{t-1} ~ N(phi*X_{t-1}, sigma^2)
#       Y_t|X_t     ~ N(0, exp(2(gamma+X_t)))
#
#       - phi and gamma are assumed known
#       - Inverse-Gamma conjugate prior for sigma^2


L = 100 # length of HMM

data=SVsim(L,c(0.9, sqrt(1-0.9^2), 1)) # (phi,sigma,gamma)

plot(data$y, type='l', col = 'blue')
lines(data$x, type='l', col = 'red')

N = 50 # number of particles
M = 1000 # length of MCMC

alpha = 0.01; beta = 0.01 # hyperparameters for IG prior
prior = c(alpha,beta)
y = data$y
## run PG and PGAS
par(mfrow = c(2,3))
samp1 = pg(M, y,prior,N, 0.2)
plot(samp1$sig2, type='l')
acf(samp1$sig2)
hist(samp1$sig2)
samp2 = pgas(M, y,prior,N, 0.2)
plot(samp2$sig2, type='l')
acf(samp2$sig2)
hist(samp2$sig2)
summary(samp1$sig2)
summary(samp2$sig2)


plot(acf(samp2$sig2, plot = F), type = 'l')
#____________________________________________
# Mixing - Path degeneracy illustration
#____________________________________________
par(mfrow = c(1,2))

N = 40
sig2 =0.19 # assume known
X = matrix(0,2,L)

## PGAS ##
init.state = cpf_as(y, sig2, N, cond = F)
particles = init.state$x;w=init.state$w[,L]
J = sample(1:N, prob = w, size = 1)
X[1,] = particles[J,]
plot(X[1,], type = 'l') # reference

alpha = prior[1];beta = prior[2]

cur.state = cpf_as(y, sig2[1], N,x.c = X[1,], cond = T)
particles = cur.state$x;w=cur.state$w[,L]
J = sample(1:N, prob = w, size = 1)
X[2,] = particles[J,]
lines(X[2,], type = 'l', col = 'blue') # sampled




## PG ##


plot(X[1,], type = 'l') # reference

cur.state = cpf(y, sig2, N,x.c = X[1,], cond = T)
particles = cur.state$x;w=cur.state$w[,L]
J = sample(1:N, prob = w, size = 1)
X[2,] = particles[J,]
lines(X[2,], type = 'l', col = 'blue') # sampled
#____________________________________________________________
# Comparing changes to sampled paths
#____________________________________________________________

n = 100
pg_change = matrix(0,5, n)
pgas_change = matrix(0,5, n)
NN = c(10, 50, 80, 100, 150)

init.state = cpf(y, sig2, N, cond = F);particles = init.state$x;w=init.state$w[,L] # reference path
J = sample(1:N, prob = w, size = 1)
x.ref = particles[J,]

for(i in 1:n){
  for(j in 1:5){
    cur.state = cpf(y, sig2[1], NN[j],x.c = x.ref, cond = T)
    particles = cur.state$x;w=cur.state$w[,L]
    J = sample(1:NN[j], prob = w, size = 1)
    x.pg = particles[J,]
    pg_change[j,] = pg_change[j,] + (x.pg!=x.ref)
    
    cur.state = cpf_as(y, sig2[1], NN[j],x.c = x.ref, cond = T)
    particles = cur.state$x;w=cur.state$w[,L]
    J = sample(1:NN[j], prob = w, size = 1)
    x.pgas = particles[J,]
    pgas_change[j,] = pgas_change[j,] + (x.pgas!=x.ref)
  }
  
  
}

pg_change = pg_change/n
pgas_change = pgas_change/n

cls =c('red','yellow','purple', 'blue', 'brown', 'grey')

par(mfrow = c(1,2))
plot(pg_change[1,], type = 'l', col = cls[1], lwd = 2, ylim = c(0,1), xlab = 't', ylab = 'Proportion changed')

lines(pg_change[2,], type = 'l', col = cls[2], lwd = 2)
lines(pg_change[3,], type = 'l', col = cls[3], lwd = 2)
lines(pg_change[4,], type = 'l', col = cls[4], lwd = 2)
lines(pg_change[5,], type = 'l', col = cls[5], lwd = 2)

abline(h=(NN[1]-1)/NN[1], lty = 'dotted', col = cls[1])
abline(h=(NN[2]-1)/NN[2], lty = 'dotted', col = cls[2])
abline(h=(NN[3]-1)/NN[3], lty = 'dotted', col = cls[3])
abline(h=(NN[4]-1)/NN[4], lty = 'dotted', col = cls[4])
abline(h=(NN[5]-1)/NN[5], lty = 'dotted', col = cls[5])
#legend('topleft',col = cls, legend=c('N = 10','N = 50','N = 80','N = 100','N = 150'), lty=1)

plot(pgas_change[1,], type = 'l', col = cls[1], lwd = 2, ylim = c(0,1), xlab = 't', ylab = '')
lines(pgas_change[2,], type = 'l', col = cls[2], lwd = 2)
lines(pgas_change[3,], type = 'l', col = cls[3], lwd = 2)
lines(pgas_change[4,], type = 'l', col = cls[4], lwd = 2)
lines(pgas_change[5,], type = 'l', col = cls[5], lwd = 2)

abline(h=(NN[1]-1)/NN[1], lty = 'dotted', col = cls[1])
abline(h=(NN[2]-1)/NN[2], lty = 'dotted', col = cls[2])
abline(h=(NN[3]-1)/NN[3], lty = 'dotted', col = cls[3])
abline(h=(NN[4]-1)/NN[4], lty = 'dotted', col = cls[4])
abline(h=(NN[5]-1)/NN[5], lty = 'dotted', col = cls[5])
legend('bottomright',  col = cls, legend=c('N = 10','N = 50','N = 80','N = 100','N = 150'), lty=1)
















