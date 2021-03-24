functions {
    real forward(real pi, real p, int M, int[] A, int[] R) {
        real even = binomial_lpmf(A[1]|A[1] + R[1], p);
        real odd = log(0);
        for(i in 2:M) {
            real next_even = log_sum_exp(even + log(1-pi), odd + log(pi)) + binomial_lpmf(A[i]|A[i] + R[i], p);
            real next_odd = log_sum_exp(even + log(pi), odd + log(1-pi)) + binomial_lpmf(R[i]|A[i] + R[i], p);
            even = next_even;
            odd = next_odd;
        }
        return log_sum_exp(even, odd);
    }
}

data {
    int<lower=0> M; // num het sites in this gene
    int<lower=0> A[M]; // alt read counts for the M sites
    int<lower=0> R[M]; // ref read counts for the M sites
    real<lower=0> sigma; // variance parameter for prior on theta
}

parameters {
    real<lower=0> theta; // amount of ASE
    real<lower=0,upper=1> pi;
}

transformed parameters {
    real p = theta / (1 + theta);
}

model {
    pi ~ beta(1, 10);
    log2(theta) ~ normal(0, sigma);
    target += -log(theta * log(2));
    target += forward(pi, p, M, A, R);
}
