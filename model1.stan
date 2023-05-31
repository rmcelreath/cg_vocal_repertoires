// model with age covariate
data{
    int N; // individuals
    int M; // tokens
    int J; // number of recordings
    vector[J] d; // durations of each recording
    array[J] int id; // individual
    array[J,M] int Y; // observed token counts
    vector[N] a; // age of each individual
}
parameters {
    vector<lower=0>[M] L; // rates of each token, conditional on individual possessing it
    vector<lower=0,upper=1>[M] alpha; // asymptotic prob of each vocalization
    real<lower=0> beta;
}
model{
    L ~ exponential(1);
    alpha ~ beta(2,2);
    beta ~ normal(2,0.5);
    // loop over recordings
    for ( j in 1:J ) {
        for ( m in 1:M ) {
            // compute vocalization probs for individual with this a[id[j]]
            real log_pa = log(alpha[m]) + log( 1-exp( -beta*a[id[j]] ) );
            if ( Y[j,m] >0 ) {
                // observed non-zero count, estimate rate
                Y[j,m] ~ poisson( L[m] * d[j] );
                target += log_pa;
            } else {
                // observed zero so account for mixture
                target += log_sum_exp( 
                    log_pa + poisson_lpmf(0|L[m]*d[j]) , // has it, but produced zero
                    log1m_exp(log_pa) // doesn't have it
                );
            }
        }//m
    }//j
}
generated quantities {
    // calculate posterior prob of each token for each individual in sample
    // this allows estimation of undercounting due to sampling and rate variation
    matrix[N,M] phi;
    vector[N] R; // repertiore size
    for ( i in 1:N ) {
        vector[2] pp;
        real Z;
        for ( m in 1:M ) {
            // Pr(p[m]==1|Y) = Pr(Y|p[m]==1) p[m] / Z
            // Z = Pr(Y|p[m]==1) p[m] + Pr(Y|p[m]==0) (1-p[m])
            // need to sum Y across all recordings for individual i
            // need sum duration of these recordings as well
            real dd = 0;
            real yy = 0;
            for ( j in 1:J ) {
                if ( id[j]==i ) {
                    dd += d[j];
                    yy += Y[j,m];
                }
            }//j

            if ( yy>0 )
                // individual has it for sure, so no calculation needed
                phi[i,m] = 1;
            else {
                // bayes it up
                real log_pa = log(alpha[m]) + log(1-exp(-beta*a[i]));
                pp[1] = log_pa + poisson_lpmf(0|L[m]*dd);
                pp[2] = log1m_exp(log_pa);
                Z = log_sum_exp( pp[1] , pp[2] );
                phi[i,m] = exp( pp[1] - Z );
            }
        }//m
        R[i] = sum( bernoulli_rng( phi[i,] ) );
    }//i
}
