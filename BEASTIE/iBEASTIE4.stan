// =======================================================================
// This is version 4 of iBEASTIE that work with het-1 gene and added deepest function
// =======================================================================
 
functions {
real forward(real[] pi,real p,int M,int[] A,int[] R,real[] missingPi,int fix) {
   real even=binomial_lpmf(A[fix]|A[fix]+R[fix],p);
   real odd=log(0);
   if(M>1){
      for(i in (fix+1):M) {
         real PI=pi[i-1]; 
         if(PI<0) { PI=missingPi[i-1]; }
         real next_even = log_sum_exp(even+log(1-PI),odd+log(PI)) + binomial_lpmf(A[i]|A[i]+R[i],p);
         real next_odd = log_sum_exp(even+log(PI),odd+log(1-PI)) + binomial_lpmf(R[i]|A[i]+R[i],p);    
         even=next_even; odd=next_odd;
      }
      for(j in 1:(fix-1)) {
         int reverse_i=fix-j;
         real reverse_PI=pi[reverse_i]; 
         if(reverse_PI<0) { reverse_PI=missingPi[reverse_i]; }
         real next_even = log_sum_exp(even+log(1-reverse_PI),odd+log(reverse_PI)) + binomial_lpmf(A[reverse_i]|A[reverse_i]+R[reverse_i],p);
         real next_odd = log_sum_exp(even+log(reverse_PI),odd+log(1-reverse_PI)) + binomial_lpmf(R[reverse_i]|A[reverse_i]+R[reverse_i],p);    
         even=next_even; odd=next_odd;
      }
   }
   return log_sum_exp(even,odd);
}
int deepest(int M,int[] A,int[] R) {
   int fix=1;
   for(i in 2:M) {
      if(A[i]+R[i]>A[fix]+R[fix]) { fix=i; }
   }
   return fix;
}
}
data {
   int<lower=1> M; // num het sites in this gene
   int<lower=0> A[M]; // alt read counts for the M sites
   int<lower=0> R[M]; // ref read counts for the M sites
   real<lower=0> sigma; // variance parameter for prior on theta
   int N_MISSING_PI; // number of missing values in pi array
   real<lower=-1,upper=1> pi[M-1]; // -1 if missing
}
transformed data {
   int fix=deepest(M,A,R);
   int<lower=0> N_PI=M-1;
}
parameters {
   real<lower=0> theta; // amount of ASE
   real<lower=0,upper=1> missingPi[M-1]; // latent PIs
}
transformed parameters {
  real<lower=0,upper=1> p=theta/(1+theta);
}
model {
   for(i in 1:N_PI) missingPi[i] ~ beta(1,10);
   log2(theta) ~ normal(0,sigma);
   target += -log(theta*log(2));
   target+=forward(pi,p,M,A,R,missingPi,fix);
}
