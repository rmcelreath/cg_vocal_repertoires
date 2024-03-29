# general covariate model - sim example

# function to compute repertoires post-sampling
Rcalc <- function(fit,S=500) {
    N <- max(fit@data$ID)
    M <- fit@data$M
    post <- extract.samples(fit)
    phi <- array(NA,c(S,N,M))
    R <- array(NA,c(S,N))
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
                phi[,n,m] <- rep(1,S)
            } else {
                # bayes
                # find a row i for this individual and behavior, so we fetch right p for its covariates
                idx <- which(fit@data$ID==n & fit@data$BID==m)[1]
                # compute
                pp1 <- log(post$p[,idx]) + dpois(0,post$L[,m]*dm,log=TRUE)
                pp2 <- log1p(-post$p[,idx])
                Z <- rep(NA,S)
                for ( s in 1:S )
                    Z[s] <- log_sum_exp( c( pp1[s] , pp2[s] ) )
                phi[,n,m] <- exp( pp1 - Z )
            }
        }#m
        for ( s in 1:S ) R[s,n] <- sum( rbern(M,prob=phi[s,n,]) )
    }#n
    return(list(phi=phi,R=R))
}

# covariates for each individual
N <- 50
M <- 50
age <- sample(1:3,size=N,replace=TRUE)
sex <- sample(1:2,size=N,replace=TRUE)
vlength <- standardize(log(1:M)) # lengths of each vocalization
#p_base <- rbeta(M,2,2) # model without length
p_base <- inv_logit( 0 - vlength )
p_sim <- matrix(NA,nrow=N,ncol=M)
for ( n in 1:N ) {
    p_sim[n,] <- inv_logit(
        logit(p_base) + 0.5*(age[n]-1) + 0.5*(sex[n]-1)
    )
}

dat <- sim_repertoire(N=N,M=M,NR=10,age_model="X",p=p_sim)

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
    M = max(MM),
    age=age,
    sex=sex,
    vlen=vlength
)

# random group label
dat_list$group <- sample(1:2,size=N,replace=TRUE)

m0u <- ulam(
    alist(
        # Y > 0
        Y|Y>0 ~ custom( log(p) + poisson_log_lpmf(Y|log_lambda1) ),
        # Y = 0
        Y|Y==0 ~ custom( log_sum_exp( 
            log(p) + poisson_log_lpmf(0|log_lambda1) , 
            log1m(p)
         ) ),
        # models for rate of behavior and prob possess behavior
        log_lambda1 <- log(L[BID]) + log(D[J]),
        save> logit(p) <- 
            # vocialization intercepts for each group
            a_mean[BID] + a[BID,group[ID]]*sigma + 
            # age-sex offsets for each group
            B_mean[age[ID],sex[ID]] + B[age[ID],sex[ID],group[ID]]*tau + 
            # association of vocalization length with intercept
            bvl*vlen[BID],
        # priors
        bvl ~ normal(0,1),
        vector[M]:L ~ exponential(1),
        vector[M]:a_mean ~ normal(0,1),
        matrix[M,2]:a ~ normal(0,1),
        matrix[3,2]:B_mean ~ normal(0,1),
        real["3,2,2"]:B ~ normal(0,1),
        c(sigma,tau) ~ exponential(1)
    ), data=dat_list , chains=4 , cores=4 , sample=TRUE , iter=500 )

precis(m0u,3,pars=c("a_mean","B","sigma","tau","bvl"))


# compute probability of each vocalization for each age,sex combination
post <- extract.samples(m0u)

S <- 1000 # number of samples
gid <- 1 # which group you want to plot

p_age_sex <- array(NA,c(S,dat_list$M,3,2)) # sample,vocalization,age,sex
for ( age in 1:3 ) for ( sex in 1:2 ) for ( v in 1:dat_list$M ) {
    p_age_sex[,v,age,sex] <- sapply( 1:S , function(s) with(post,
        inv_logit( 
            a_mean[s,v] + a[s,v,gid]*sigma[s] + 
            B_mean[s,age,sex] + B[s,age,sex,gid]*tau[s] +
            bvl[s]*vlength[v]
        )
    ))
}

# plot age groups for specific sex
plot(NULL,xlim=c(1,dat_list$M),ylim=c(0,1),xlab="vocalization",ylab="probability")
p_mean <- apply(p_age_sex,2:4,mean)
p_ci <- apply(p_age_sex,2:4,PI,prob=0.89)
the_sex <- 1
for ( a in 1:3 ) points( 1:dat_list$M , p_mean[,a,the_sex] , col=a , pch=16 )
# intervals for age group 1 as example
the_age <- 1
# for ( v in 1:dat_list$M ) lines( c(v,v) , p_ci[,v,the_age,the_sex] , col=the_age )

# repertoire sizes for each age,sex
rep_age_sex <- array(NA,c(S,3,2))
for ( age in 1:3 ) for ( sex in 1:2 ) for ( s in 1:S ) {
    rep_age_sex[s,age,sex] <- sum( p_age_sex[s,,age,sex] )
}
apply( rep_age_sex , 2:3 , mean )

# calculate distribution of longest vocalization (sequence) for each age,sex
# vector of sequence lengths
vlen <- 1:dat_list$M # replace with real data
maxlen_age_sex <- array(NA,c(S,3,2))
for ( age in 1:3 ) for ( sex in 1:2 ) for ( s in 1:S ) {
    x <- rbern( dat_list$M , p_age_sex[s,,age,sex] ) # random realization of repertoire
    maxlen_age_sex[s,age,sex] <- max( vlen * x )
}
apply( maxlen_age_sex , 2:3 , mean )

# compare estimated repertoire size to observed - for actual individuals in sample
r <- Rcalc(m0u,S=1000)
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

plot(Rtrue,apply(r$R,2,mean),xlab="true repertoire size",ylab="posterior mean size",col=2,lwd=2)
abline(a=0,b=1)
