// =======================================================================
// This is version 4 of iBEASTIE that work with het-1 gene and added deepest function
// =======================================================================

functions {
  real forward(array[] real pi, real p, int M, array[] int A, array[] int R, array[] real missingPi, int fix) {
    real even = binomial_lpmf(A[fix] | A[fix] + R[fix], p);
    real odd = log(0);
    for (i in (fix + 1):M) {
      real PI = pi[i - 1]; 
      if (PI < 0) { PI = missingPi[i - 1]; }
      real next_even = log_sum_exp(even + log(1 - PI), odd + log(PI)) + binomial_lpmf(A[i] | A[i] + R[i], p);
      real next_odd = log_sum_exp(even + log(PI), odd + log(1 - PI)) + binomial_lpmf(R[i] | A[i] + R[i], p);
      even = next_even; 
      odd = next_odd;
    }
    for (j in 1:(fix - 1)) {
      int reverse_i = fix - j;
      real reverse_PI = pi[reverse_i]; 
      if (reverse_PI < 0) { reverse_PI = missingPi[reverse_i]; }
      real next_even = log_sum_exp(even + log(1 - reverse_PI), odd + log(reverse_PI)) + binomial_lpmf(A[reverse_i] | A[reverse_i] + R[reverse_i], p);
      real next_odd = log_sum_exp(even + log(reverse_PI), odd + log(1 - reverse_PI)) + binomial_lpmf(R[reverse_i] | A[reverse_i] + R[reverse_i], p);
      even = next_even; 
      odd = next_odd;
    }
    return log_sum_exp(even, odd);
  }

  int deepest(int M, array[] int A, array[] int R) {
    int fix = 1;
    for (i in 2:M) {
      if (A[i] + R[i] > A[fix] + R[fix]) {
        fix = i;
      }
    }
    return fix;
  }
}

data {
  int<lower=1> M; // num het sites in this gene
  array[M] int A; // alt read counts for the M sites
  array[M] int R; // ref read counts for the M sites
  real<lower=0> sigma; // variance parameter for prior on theta
  int<lower=0> N_MISSING_PI; // number of missing values in pi array
  array[N_MISSING_PI] real<lower=-1,upper=1> pi; // -1 if missing, empty if N_MISSING_PI is 0
}
transformed data {
  int fix = deepest(M, A, R);
  int<lower=0> N_PI = max(M - 1, 0); // ensure this is not negative
}
parameters {
  real<lower=0> theta; // amount of ASE
  array[N_MISSING_PI] real<lower=0,upper=1> missingPi; // latent PIs for missing data
}
transformed parameters {
  real<lower=0,upper=1> p = theta / (1 + theta);
}
model {
  // Prior for theta
  log2(theta) ~ normal(0, sigma);
  target += -log(theta * log(2)); // Jacobian adjustment for log transform

  // Direct calculation for single-site genes
  if (M == 1) {
    target += binomial_lpmf(A[1] | A[1] + R[1], p);
  } else {
    // Priors for missing PI values
    for (i in 1:N_MISSING_PI) {
      missingPi[i] ~ beta(1, 10);
    }
    // Use the forward function for genes with more than 1 het site
    target += forward(pi, p, M, A, R, missingPi, fix);
  }
}
