# Cedric vocal sampling problem

library(rethinking)

# synthetic data

N <- 20 # number of individuals
M <- 10 # number of unique vocalizations
p <- sort(runif(M)) # prob of each vocalization, ordered rare to common, among individuals
L <- rlnorm(M,0.25,0.5)/10 # rate of each vocalization, when individual uses it

# simulate individual repertoires
# individuals in rows, tokens in columns
x <- t( replicate( N , rbern(M,prob=p) ) )

# simulate vocal samples with varying sample sizes for each individual
max_tokens <- 20 # max number of tokens per recording
n_recordings <- rpois(N,10) + 1
id <- rep(1:N,times=n_recordings) # individual making each recording
y <- matrix(NA,sum(n_recordings),max_tokens)
d <- 0.2 * rpois(sum(n_recordings),5) + 1 # durations

if ( FALSE ) {
# use Gillespie algorithm to simulate vocalizations
# this gives intervals between tokens, but we don't need that right now
for ( i in 1:sum(n_recordings) ) {    
    t0 <- 0 # time in recording
    j <- 1
    while( (t0 < d[i]) & (j <= max_tokens) ) {
        r <- runif(2)
        R <- sum(L*x[id[i],]) # total rate of vocalizations
        t1 <- (1/R)*log(1/r[1]) # time to next vocalization
        # which vocalization?
        k <- L*x[id[i],] # rates for this individual
        y[i,j] <- sample(1:M,size=1,prob=k) # draw 1 proportional to relative rates
        # update
        t0 <- t1
        j <- j + 1
        # check for max_tokens breech and truncate recording at this point
        if ( j > max_tokens ) d[i] <- t1
    }#t0 < d
}#i
}

# count number of each token in each recording
XX <- matrix(NA,sum(n_recordings),M)
for ( i in 1:sum(n_recordings) ) {
    for ( j in 1:M )
        XX[i,] <- rpois( M , L*x[id[i],]*d[i] )
}#i

# estimate

dat <- list(
    N=N,
    M=M,
    J=sum(n_recordings),
    K=max_tokens,
    d=d,
    id=id,
    Y=XX )

m0 <- cstan(file="model0.stan",data=dat,chains=1)

precis(m0,2)
post <- extract.samples(m0)

# blank2(w=2)
par(mfrow=c(1,2))

plot( L , apply(post$L,2,mean) , xlab="true rate" , ylab="posterior mean rate" , col=2 , lwd=2 )
abline(a=0,b=1)

plot( p , apply(post$p,2,mean) , xlab="true prob" , ylab="posterior mean prob" , col=2 , lwd=2 )
abline(a=0,b=1)

# individual repertoires

phi <- apply(post$phi,2:3,mean)
# round(phi,2)

par(mfrow=c(1,2))
image(x)
image(phi)

# compare estimated repertoire size to observed
Robs = rep(0,N)
Rtrue = rep(0,N)
for ( i in 1:N ) {
    for ( m in 1:M ) {
        # does individual i ever do m?
        ii <- which(id==i)
        if ( any( XX[ii,m] > 0 ) )
            Robs[i] = Robs[i] + 1
    }#m
    Rtrue[i] = sum(x[i,])
}#i

par(mfrow=c(1,2))
plot(Rtrue,Robs,xlab="true repertoire size",ylab="observed size",col=2,lwd=2)
abline(a=0,b=1)

plot(Rtrue,apply(post$R,2,mean),xlab="true repertoire size",ylab="posterior mean size",col=2,lwd=2)
abline(a=0,b=1)
