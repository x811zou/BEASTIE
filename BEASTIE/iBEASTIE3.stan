// =======================================================================
// This is version 3 of iBEASTIE that accommodates sites with no information
// about phasing error rate (due to missing LD information)
// =======================================================================

functions {
real forward(real[] pi,real p,int M,int[] A,int[] R,real[] missingPi) {
   real even=binomial_lpmf(A[1]|A[1]+R[1],p);
   real odd=log(0);
   int j=1;
   for(i in 2:M) {
      real PI=pi[i-1]; real next_even; real next_odd;
      if(PI<0) { PI=missingPi[j]; j+=1; }
      next_even = log_sum_exp(even+log(1-PI),odd+log(PI)) +
         binomial_lpmf(A[i]|A[i]+R[i],p);
      next_odd = log_sum_exp(even+log(PI),odd+log(1-PI)) +
         binomial_lpmf(R[i]|A[i]+R[i],p);
      even=next_even; odd=next_odd;
   }
   return log_sum_exp(even,odd);
}}
data {
   int<lower=0> M; // num het sites in this gene
   int<lower=0> A[M]; // alt read counts for the M sites
   int<lower=0> R[M]; // ref read counts for the M sites
   real<lower=0> sigma; // variance parameter for prior on theta
   int N_MISSING_PI; // number of missing values in pi array
   real<lower=-1,upper=1> pi[M-1]; // -1 if missing
}
parameters {
   real<lower=0> theta; // amount of ASE
   real<lower=0,upper=1> missingPi[N_MISSING_PI]; // latent PIs
}
transformed parameters {
  real<lower=0,upper=1> p=theta/(1+theta);
}
model {
   for(i in 1:N_MISSING_PI) missingPi[i] ~ beta(1,10);
   log2(theta) ~ normal(0,sigma);
   target += -log(theta*log(2));
   target+=forward(pi,p,M,A,R,missingPi);
}