data {
  int<lower=0> N; // Number of observations
  int<lower=0> P; // Number of plots
  matrix[N, P] M; // Model matrix for plot.year
  
  vector[N] n2o; // Daily nitrous oxide fluxes, g/ha/day
  vector[N] no3; // NO3 concentrations, ug/g soil
  vector[N] nh4; // NH4 concentrations, ug/g soil
  vector[N] wfps; // Water-filled pore space, %
  vector[N] soilT; // Soil temperature, degrees C
}
parameters {
  vector<lower=1, upper=10000>[P] PDR; // Potential denitrification rate
  real<lower=0.1, upper=10> DWa; // Multiplier for WFPS relationship
  real<lower=5, upper=120> DKm; // half-saturation constant (ug/g soil)
  real<lower=1, upper=2.5> DQ10; // Q10
  
  vector<lower=1, upper=1000>[P] PNR; // Potential nitrification rate
  real<lower=1, upper=50> NKm; // Half saturation for ammonium (ug/g soil)
  real<lower=7, upper=15> NWa; // Alpha parameter for moisture
  real<lower=1.01, upper=12> NWb; // Beta parameter for moisture
  real<lower=4.5, upper=13> NQ10; // Q10
  
  real<lower=0, upper=30> sig; // Variance

}
transformed parameters {
  // Denitrification constraints
  vector[N] Dw; // Soil moisture
  vector[N] Dn; // Nitrate concentration
  vector[N] Dt; // Temperature
  // Nitrification constraints
  vector[N] Nw1; // First step for soil moisture
  vector[N] Nw; // Soil moisture
  vector[N] Nn; // Ammonium concentration
  vector[N] Nt; // Temperature
  
  vector[N] est; //Estimated flux
  
  Dw <- -(1-wfps)*DWa;
  Dt <- (soilT - 20) * log(DQ10) / 10;
  Dn <- log(no3) - log(no3 + DKm);
  
  for(n in 1:N){
    Nw1[n] <- beta_log(wfps[n], NWa, NWb);
  }
  Nw <- Nw1 - beta_log((NWa-1) / (NWa+NWb-2), NWa, NWb);
  Nt <- (soilT - 20) * log(NQ10) / 10;
  Nn <- log(nh4) - log(nh4 + NKm);
  
  est <- exp(log(M * PDR) + Dw + Dt + Dn) + exp(log(M * PNR) + Nw + Nt + Nn);
}
model {
  n2o ~ normal(est, sig);
	  
}