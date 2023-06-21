# Cedric vocal sampling problem

library(rethinking)

# synthetic data

sim_repertoire <- function(
    N = 50, # number of individuals
    M = 50, # number of unique vocalizations
    NR = 20, # average recordings per individual
    age_model = "NULL",
    beta = 2,
    age,
    p,L,
    X
) {

    a <- NA
    
    if ( missing(p) )
        p <- sort(rbeta(M,5,2)) # prob of each vocalization, ordered rare to common, among individuals
    if ( missing(L) )
        L <- rlnorm(M,0.25,0.5)/10 # rate of each vocalization, when individual uses it

    # simulate individual repertoires
    # individuals in rows, tokens in columns
    x <- matrix(NA,nrow=N,ncol=M)
    if ( age_model=="NULL" )
        x <- t( replicate( N , rbern(M,prob=p) ) )

    # age as individual covariate
    if ( age_model=="monotonic" ) {
        # continuous age
        if ( missing(age) )
            a <- runif(N)
        else
            a <- age
        alpha <- p # asymptotic probabilities of each 
        # curve( 0.6*(1-exp(-beta*x)) , from=0 , to=1 , ylim=c(0,1) )
        for ( i in 1:N ) {
            x[i,] <- rbern( M , prob=alpha*(1-exp(-beta*a[i])) )
        }#i
    }

    if ( age_model=="cat" ) {
        # categorical age
        if ( missing(age) )
            a <- runif(N)
        else
            a <- age
        for ( i in 1:N ) {
            x[i,] <- rbern( M , prob=p*(1-exp(-beta*a[i])) )
        }#i
        a_cat <- ifelse( a < 0.25 , 1 , a )
        a_cat <- ifelse( a < 0.5 & a > 0.25 , 2 , a_cat )
        a_cat <- ifelse( a > 0.5 , 3 , a_cat )
        a <- a_cat
    }

    # general categorical covariates, table X supplied for individual covariates, p matrix supplied as argument in which X values already used to generate p values for each individual (set of covariates)
    if ( age_model=="X" ) {
        for ( i in 1:N ) {
            x[i,] <- rbern( M , prob=p[i,] )
        }#i
    }

    # simulate vocal samples with varying sample sizes for each individual
    max_tokens <- 20 # max number of tokens per recording
    n_recordings <- rpois(N,NR) + 1
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

    dat <- list(
        N=N, # number of individuals
        M=M, # number of unique vocalizations
        J=sum(n_recordings), # number of recordings (all individuals)
        d=d, # durations of each recording
        id=id, # individual IDs, length J
        Y=XX, # matrix with recordings on rows and counts of vocalizations in columns
        a=a, # age covariate
        L=L, # rates
        p=p, # base rates
        Rtrue=x ) # true repertoires 

    return(dat)
}

################################################################################
# estimate - basic no covariate model

dat <- sim_repertoire(N=50,M=50,NR=10)
dat$a <- NULL

m0 <- cstan(file="model0.stan",data=dat,chains=1)

precis(m0,2)
post <- extract.samples(m0)

# blank2(w=2)
par(mfrow=c(1,2))

plot( dat$L , apply(post$L,2,mean) , xlab="true rate" , ylab="posterior mean rate" , col=2 , lwd=2 )
abline(a=0,b=1)

plot( dat$p , apply(post$p,2,mean) , xlab="true prob" , ylab="posterior mean prob" , col=2 , lwd=2 )
abline(a=0,b=1)

# individual repertoires

phi <- apply(post$phi,2:3,mean)
# round(phi,2)

par(mfrow=c(1,2))
image(dat$Rtrue)
image(phi)

# compare estimated repertoire size to observed
Robs = rep(0,dat$N)
Rtrue = rep(0,dat$N)
for ( i in 1:dat$N ) {
    for ( m in 1:dat$M ) {
        # does individual i ever do m?
        ii <- which(dat$id==i)
        if ( any( dat$Y[ii,m] > 0 ) )
            Robs[i] = Robs[i] + 1
    }#m
    Rtrue[i] = sum(dat$Rtrue[i,])
}#i

par(mfrow=c(1,2))
plot(Rtrue,Robs,xlab="true repertoire size",ylab="observed size",col=2,lwd=2)
abline(a=0,b=1)

plot(Rtrue,apply(post$R,2,mean),xlab="true repertoire size",ylab="posterior mean size",col=2,lwd=2)
abline(a=0,b=1)

# ulam version with post-sampling computaiton of repertoires
# convert Y matrix to long form with m column
YY <- rep(NA,length(dat$Y))
MM <- YY
ID <- YY
JJ <- YY
k <- 1
for ( j in 1:nrow(dat$Y) ) for ( m in 1:ncol(dat$Y) ) {
    YY[k] <- dat$Y[j,m]
    MM[k] <- m
    ID[k] <- dat$id[j]
    JJ[k] <- j
    k <- k + 1
}

dat_list <- list(
    Y = YY,    # behavior counts
    ID = ID,   # individual ID
    BID = MM,  # behavior ID
    J = JJ,    # recording ID
    D = dat$d, # duration of recording,
    M = max(MM)
)

m0u <- ulam(
    alist(
        # Y > 0
        Y|Y>0 ~ custom( log(p[BID]) + poisson_log_lpmf(Y|log_lambda1) ),
        # Y = 0
        Y|Y==0 ~ custom( log_sum_exp( 
            log(p[BID]) + poisson_log_lpmf(0|log_lambda1) , 
            log1m(p[BID])
         ) ),
        # log-linear model for rate of behavior
        log_lambda1 <- log(L[BID]) + log(D[J]),
        # priors
        vector[M]:L ~ exponential(1),
        vector[M]:p ~ beta(2,2)
    ), data=dat_list , chains=1 , sample=TRUE )

postu <- extract.samples(m0u)

# function to compute repertoires
Rcalc <- function(fit) {
    N <- max(fit@data$ID)
    M <- fit@data$M
    post <- extract.samples(fit)
    S <- dim(post$p)[1]
    R <- array(NA,c(N,M,S))
    for ( n in 1:N ) {
        for ( m in 1:M ) {
            # sum all m observations for individual n
            dm <- 0 # sum duration
            ym <- 0 # sum count
            for ( i in 1:length(fit@data$Y) ) {
                if ( fit@data$ID[i]==n & fit@data$BID[i]==m ) {
                    dm <- dm + fit@data$D[ fit@data$J[i] ]
                    ym <- ym + fit@data$Y[i]
                }
            }#i
            # now invert probability
            if ( ym > 0 ) {
                # observed
                R[n,m,] <- 1
            } else {
                # bayes
                pp1 <- log(post$p[,m]) + dpois(0,post$L[,m]*dm,log=TRUE);
                pp2 <- log1p(-post$p[,m]);
                Z <- rep(NA,S)
                for ( s in 1:S )
                    Z[s] <- log_sum_exp( c( pp1[s] , pp2[s] ) );
                R[n,m,] <- exp( pp1 - Z );
            }
        }#m
    }#n
    return(R)
}

phiu <- Rcalc(m0u)
phi <- post$phi
plot( apply(phi,2:3,mean) , apply(phiu,1:2,mean) , xlab="stan" , ylab="ulam" , col=2 , lwd=2 )
abline(a=0,b=1)

######################################################################################
# categorical covariate model

dat <- sim_repertoire(N=50,M=50,NR=10,age_model="cat",beta=2)
dat$Na <- 3

m0 <- cstan(file="model1_cat.stan",data=dat,chains=1)

precis(m0,2)
post <- extract.samples(m0)

# blank2(w=2)

plot( dat$L , apply(post$L,2,mean) , xlab="true rate" , ylab="posterior mean rate" , col=2 , lwd=2 )
abline(a=0,b=1)

# individual repertoires
phi <- apply(post$phi,2:3,mean)
# round(phi,2)

par(mfrow=c(1,2))
image(dat$Rtrue)
image(phi)

# compare estimated repertoire size to observed
Robs = rep(0,dat$N)
Rtrue = rep(0,dat$N)
for ( i in 1:dat$N ) {
    for ( m in 1:dat$M ) {
        # does individual i ever do m?
        ii <- which(dat$id==i)
        if ( any( dat$Y[ii,m] > 0 ) )
            Robs[i] = Robs[i] + 1
    }#m
    Rtrue[i] = sum(dat$Rtrue[i,])
}#i
Rest <- apply(post$R,2,mean)

par(mfrow=c(1,2))
plot(Rtrue,Robs,xlab="true repertoire size",ylab="observed size",col=2,lwd=2)
abline(a=0,b=1)

plot(Rtrue,Rest,xlab="true repertoire size",ylab="posterior mean size",col=2,lwd=2)
abline(a=0,b=1)

# age against repertoire
yr <- max(abs(Robs - Rtrue))
plot( dat$a , Robs - Rtrue , xlim=c(0,4) , xlab="age" , ylab="repertoire size" , lwd=2 , col=2 , ylim=c(-yr,yr) )
points( dat$a , Rest - Rtrue , lwd=2 , col=1 )
pi <- apply(post$R,2,PI,0.89)
for ( i in 1:dat$N ) lines( rep(dat$a[i],2) , pi[,i] - Rtrue[i] , col=1 )
abline(h=0,lty=2)


######################################################################################
# covariate model

dat <- sim_repertoire(N=50,M=50,NR=10,age_model="monotonic",beta=2)

m0 <- cstan(file="model1.stan",data=dat,chains=1)

precis(m0,2)
post <- extract.samples(m0)

# blank2(w=2)
par(mfrow=c(1,3))

plot( dat$L , apply(post$L,2,mean) , xlab="true rate" , ylab="posterior mean rate" , col=2 , lwd=2 )
abline(a=0,b=1)

plot( dat$p , apply(post$alpha,2,mean) , xlab="true prob" , ylab="posterior mean prob" , col=2 , lwd=2 )
abline(a=0,b=1)

dens(post$beta)
abline(v=2,col=2)

# individual repertoires
phi <- apply(post$phi,2:3,mean)
# round(phi,2)

par(mfrow=c(1,2))
image(dat$Rtrue)
image(phi)

# compare estimated repertoire size to observed
Robs = rep(0,dat$N)
Rtrue = rep(0,dat$N)
for ( i in 1:dat$N ) {
    for ( m in 1:dat$M ) {
        # does individual i ever do m?
        ii <- which(dat$id==i)
        if ( any( dat$Y[ii,m] > 0 ) )
            Robs[i] = Robs[i] + 1
    }#m
    Rtrue[i] = sum(dat$Rtrue[i,])
}#i
Rest <- apply(post$R,2,mean)

par(mfrow=c(1,2))
plot(Rtrue,Robs,xlab="true repertoire size",ylab="observed size",col=2,lwd=2)
abline(a=0,b=1)

plot(Rtrue,Rest,xlab="true repertoire size",ylab="posterior mean size",col=2,lwd=2)
abline(a=0,b=1)

# age against repertoire
yr <- max(abs(Robs - Rtrue))
plot( dat$a , Robs - Rtrue , xlim=c(0,1) , xlab="age" , ylab="repertoire size" , lwd=2 , col=2 , ylim=c(-yr,yr) )
points( dat$a , Rest - Rtrue , lwd=2 , col=1 )
pi <- apply(post$R,2,PI,0.89)
for ( i in 1:dat$N ) lines( rep(dat$a[i],2) , pi[,i] - Rtrue[i] , col=1 )
abline(h=0,lty=2)

# functions
plot(NULL,xlim=c(0,1),ylim=c(0,1),xlab="age",ylab="p")
xseq <- seq(from=0,to=1,len=30)
for ( m in dat$M:dat$M ) {
    for ( s in 1:100 ) {
        pm <- sapply( xseq , function(a) mean( post$alpha[s,m] * ( 1 - exp( -post$beta[s]*a ) ) ) )
        lines( xseq , pm , lwd=1 )
    }
    pm_true <- sapply( xseq , function(a) ( dat$p[m] * ( 1 - exp( -2*a ) ) ) )
    lines( xseq , pm_true , lwd=6 , col="white" )
    lines( xseq , pm_true , lwd=2 , col=2 )
}

##########
# prior predictive for age function

S <- 20
prior <- list()
prior$alpha <- rbeta(S,2,2)
prior$beta <- rexp(S,0.25)
prior$gamma <- rexp(S,2) + 1

# blank2(w=2)

plot(NULL,xlim=c(0,1),ylim=c(0,1),xlab="age",ylab="prob in repertoire")
for ( i in 1:S ) with( prior , curve( alpha[i]*(1-exp(-beta[i]*x))^gamma[i] , add=TRUE , lwd=2 ) )
