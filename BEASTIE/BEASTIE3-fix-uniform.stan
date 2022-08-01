functions {
real forward(real pi,real p,int M,int[] A,int[] R,int fix) {
    real even=binomial_lpmf(A[fix]|A[fix]+R[fix],p);
    real odd=log(0);
    for(i in (fix+1):M) {
        real next_even = log_sum_exp(even+log(1-pi),odd+log(pi)) + binomial_lpmf(A[i]|A[i]+R[i],p);
        real next_odd = log_sum_exp(even+log(pi),odd+log(1-pi)) + binomial_lpmf(R[i]|A[i]+R[i],p);    
        even=next_even; odd=next_odd;
    }
    even=log_sum_exp(even,odd);
    odd=log(0);
    for(j in 1:(fix-1)) {
        int i=fix-j;
        real next_even = log_sum_exp(even+log(1-pi),odd+log(pi)) + binomial_lpmf(A[i]|A[i]+R[i],p);
        real next_odd = log_sum_exp(even+log(pi),odd+log(1-pi)) + binomial_lpmf(R[i]|A[i]+R[i],p);    
        even=next_even; odd=next_odd;
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
   int<lower=0> M; // num het sites in this gene
   int<lower=0> A[M]; // alt read counts for the M sites
   int<lower=0> R[M]; // ref read counts for the M sites
   real<lower=0> sigma; // variance parameter for prior on theta
}
 
transformed data {
   int fix=deepest(M,A,R);
}
 
parameters {
   real<lower=0> theta; // amount of ASE
}
 
transformed parameters {
  real p=theta/(1+theta);
}
 
model {
   log2(theta) ~ normal(0,sigma);
   target += -log(theta*log(2));
   target+=forward(0.5,p,M,A,R,fix);
}
