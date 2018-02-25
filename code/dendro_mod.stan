data {
  int<lower=0> n;
  int<lower=0> nsp;
  int<lower=0> ntree;
  int<lower=0> nplanting;
  int<lower=0> nyear;
  int<lower=0> nrem;
  int<lower=0> nsoil;
  vector[n] y;
  int<lower=0> id[n];
  int<lower=0> year[n];
  int<lower=0> sp[n];
  int<lower=0> planting[n];
  int<lower=0> soil_tex[n];
  int<lower=0> remnant[n];
  vector[n] rain;
  vector[n] maxt;
  vector[n] uthor;
  vector[n] soilk;
  vector[n] site_age;
  vector[n] sp_rich;
  vector[n] tree_size;
  vector[n] nfix_present;
  vector[n] tree_dens;
  vector[n] basal_area;
}

parameters {
  real alpha;
  vector[nsp] alpha_sp;
  vector[nsp] beta_rain;
  vector[nsp] beta_maxt;
  vector[nsp] beta_uthor;
  vector[nsp] beta_soilk;
  vector[nsp] beta_rich;
  vector[nsp] beta_age;
  vector[nsp] beta_size;
  vector[nsp] beta_nfix;
  vector[nsp] beta_dens;
  vector[nsp] beta_basal;
  vector[nsp] beta_soil[nsoil];
  vector[ntree] gamma_tree;
  vector[nplanting] gamma_planting;
  vector[nyear] gamma_year;
  vector[nrem] gamma_remnant;
  real<lower=0> sigma_rain;
  real<lower=0> sigma_maxt;
  real<lower=0> sigma_uthor;
  real<lower=0> sigma_soilk;
  real<lower=0> sigma_rich;
  real<lower=0> sigma_age;
  real<lower=0> sigma_size;
  real<lower=0> sigma_dens;
  real<lower=0> sigma_basal;
  real<lower=0> sigma_soil[nsoil];
  real<lower=0> sigma_nfix;
  real<lower=0> sigma_tree;
  real<lower=0> sigma_planting;
  real<lower=0> sigma_year;
  real<lower=0> sigma_remnant;
  real<lower=0> sigma_main;
  real rain_mean;
  real maxt_mean;
  real uthor_mean;
  real soilk_mean;
  real rich_mean;
  real age_mean;
  real size_mean;
  real nfix_mean;
  real dens_mean;
  real basal_mean;
  real soil_mean[nsoil];
}

transformed parameters {
  vector[n] mu;
  
  for (i in 1:n) {
    mu[i] = alpha +
      alpha_sp[sp[i]] +
      beta_rain[sp[i]] * rain[i] +
      beta_maxt[sp[i]] * maxt[i] +
      beta_uthor[sp[i]] * uthor[i] +
      beta_soilk[sp[i]] * soilk[i] +
      beta_rich[sp[i]] * sp_rich[i] +
      beta_age[sp[i]] * site_age[i] +
      beta_size[sp[i]] * tree_size[i] +
      beta_nfix[sp[i]] * nfix_present[i] +
      beta_dens[sp[i]] * tree_dens[i] +
      beta_basal[sp[i]] * basal_area[i] +
      beta_soil[soil_tex[i]][sp[i]] +
      gamma_tree[id[i]] +
      gamma_planting[planting[i]] +
      gamma_remnant[remnant[i]] +
      gamma_year[year[i]];
  }
}

model {
  // likelihood
  y ~ normal(mu, sigma_main);
  
  // priors
  alpha ~ normal(0.0, 1.0);
  alpha_sp[1] ~ normal(0.0, 0.0001);
  for (i in 2:nsp)
    alpha_sp[i] ~ normal(0.0, 1.0);
  beta_rain ~ normal(rain_mean, sigma_rain);
  beta_maxt ~ normal(maxt_mean, sigma_maxt);
  beta_uthor ~ normal(uthor_mean, sigma_uthor);
  beta_soilk ~ normal(soilk_mean, sigma_soilk);
  beta_rich ~ normal(rich_mean, sigma_rich);
  beta_age ~ normal(age_mean, sigma_age);
  beta_size ~ normal(size_mean, sigma_size);
  beta_nfix ~ normal(nfix_mean, sigma_nfix);
  beta_dens ~ normal(dens_mean, sigma_dens);
  beta_basal ~ normal(basal_mean, sigma_basal);
  beta_soil[1] ~ normal(soil_mean[1], sigma_soil[1]);
  for (i in 2:nsoil)
    beta_soil[i] ~ normal(soil_mean[i], sigma_soil[i]);
  
  rain_mean ~ normal(0.0, 2.0);
  maxt_mean ~ normal(0.0, 2.0);
  uthor_mean ~ normal(0.0, 2.0);
  soilk_mean ~ normal(0.0, 2.0);
  rich_mean ~ normal(0.0, 2.0);
  age_mean ~ normal(0.0, 2.0);
  size_mean ~ normal(0.0, 2.0);
  nfix_mean ~ normal(0.0, 2.0);
  dens_mean ~ normal(0.0, 2.0);
  basal_mean ~ normal(0.0, 2.0);
  soil_mean[1] ~ normal(0.0, 0.0001);
  for (i in 2:nsoil)
    soil_mean[i] ~ normal(0.0, 2.0);

  // exchangeable priors
  gamma_tree ~ normal(0.0, sigma_tree);
  gamma_planting ~ normal(0.0, sigma_planting);
  gamma_remnant ~ normal(0.0, sigma_remnant);
  gamma_year ~ normal(0.0, sigma_year);
  
  // overall variance
  sigma_rain ~ normal(0.0, 2.0);
  sigma_maxt ~ normal(0.0, 2.0);
  sigma_uthor ~ normal(0.0, 2.0);
  sigma_soilk ~ normal(0.0, 2.0);
  sigma_rich ~ normal(0.0, 2.0);
  sigma_age ~ normal(0.0, 2.0);
  sigma_size ~ normal(0.0, 2.0);
  sigma_dens ~ normal(0.0, 2.0);
  sigma_basal ~ normal(0.0, 2.0);
  sigma_nfix ~ normal(0.0, 2.0);
  sigma_soil[1] ~ normal(0.0, 0.0001);
  for (i in 2:nsoil)
    sigma_soil[i] ~ normal(0.0, 2.0);
  sigma_tree ~ normal(0.0, 2.0);
  sigma_planting ~ normal(0.0, 2.0);
  sigma_year ~ normal(0.0, 2.0);
  sigma_remnant ~ normal(0.0, 2.0);
  sigma_main ~ normal(0.0, 2.0);
}

generated quantities {
  vector[n] loglik;
  
  for (i in 1:n)
    loglik[i] = normal_lpdf(y[i] | mu[i], sigma_main);
}
