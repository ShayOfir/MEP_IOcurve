data {
  int<lower=1> N;
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> K;

  array[N] int<lower=1,upper=I> subj;
  array[N] int<lower=1,upper=J> coil;
  array[N] int<lower=1,upper=K> side;

  vector[N] Intensity;
  vector[N] Y;

  vector[J] MSO_limit;

  // Prior hyperparameters
  real prior_mu_theta_mean;
  real prior_mu_theta_sd;

  real prior_mu_slope_mean;
  real prior_mu_slope_sd;

  real prior_mu_A_mean;
  real prior_mu_A_sd;

  real prior_sigma_theta_loc;
  real prior_sigma_theta_scale;

  real prior_sigma_slope_loc;
  real prior_sigma_slope_scale;

  real prior_sigma_A_loc;
  real prior_sigma_A_scale;

  real prior_sigma_loc;
  real prior_sigma_scale;

  // Coil/side effects on θ
  real prior_delta_theta_coil_loc;
  real prior_delta_theta_coil_scale;

  real prior_delta_theta_side_loc;
  real prior_delta_theta_side_scale;

  // Coil/side effects on slope
  real prior_delta_slope_coil_loc;
  real prior_delta_slope_coil_scale;

  real prior_delta_slope_side_loc;
  real prior_delta_slope_side_scale;

  // Coil/side effects on logA
  real prior_delta_logA_coil_loc;
  real prior_delta_logA_coil_scale;

  real prior_delta_logA_side_loc;
  real prior_delta_logA_side_scale;

  real sigma_ubound;
  real logA_cap;

  int<lower=0,upper=1> use_likelihood;
}

parameters {
  // Subject-level non-centered random effects
  vector[I] theta_z;
  vector[I] slope_z;
  vector[I] logA_z;

  // Population means
  real mu_theta;
  real mu_slope;
  real mu_logA;

  // Population scales
  real<lower=0, upper=sigma_ubound> sigma_theta;
  real<lower=0, upper=sigma_ubound> sigma_slope;
  real<lower=0, upper=sigma_ubound> sigma_logA;
  real<lower=0, upper=sigma_ubound> sigma;

  // RAW coil/side fixed effects (centered later)
  vector[J] delta_theta_coil_raw;
  vector[K] delta_theta_side_raw;

  vector[J] delta_slope_coil_raw;
  vector[K] delta_slope_side_raw;

  vector[J] delta_logA_coil_raw;
  vector[K] delta_logA_side_raw;
}

transformed parameters {
  // Centered fixed effects
  vector[J] delta_theta_coil =
    delta_theta_coil_raw - mean(delta_theta_coil_raw);
  vector[K] delta_theta_side =
    delta_theta_side_raw - mean(delta_theta_side_raw);

  vector[J] delta_slope_coil =
    delta_slope_coil_raw - mean(delta_slope_coil_raw);
  vector[K] delta_slope_side =
    delta_slope_side_raw - mean(delta_slope_side_raw);

  vector[J] delta_logA_coil =
    delta_logA_coil_raw - mean(delta_logA_coil_raw);
  vector[K] delta_logA_side =
    delta_logA_side_raw - mean(delta_logA_side_raw);

  array[I, J, K] real theta_raw;
  array[I, J, K] real slope_raw;
  array[I, J, K] real logA_raw;
  array[I, J, K] real A_raw;

  vector[N] mu_curve;
  vector[N] mu_log;

  // Build IO-curve parameters
  for (i in 1:I)
    for (j in 1:J)
      for (k in 1:K) {

        // logA: subject RE + coil/side FE
        logA_raw[i,j,k] =
          mu_logA +
          delta_logA_coil[j] +
          delta_logA_side[k] +
          sigma_logA * logA_z[i];

        // cap to avoid overflow
        real logA_capped = fmin(logA_raw[i,j,k], logA_cap);
        A_raw[i,j,k] = exp(logA_capped);

        // θ: subject RE + coil/side FE
        theta_raw[i,j,k] =
          mu_theta +
          delta_theta_coil[j] +
          delta_theta_side[k] +
          sigma_theta * theta_z[i];

        // slope: subject RE + coil/side FE
        slope_raw[i,j,k] =
          mu_slope +
          delta_slope_coil[j] +
          delta_slope_side[k] +
          sigma_slope * slope_z[i];
      }

  // Observation-level means
  for (n in 1:N) {
    int i = subj[n];
    int j = coil[n];
    int k = side[n];

    mu_curve[n] =
      A_raw[i,j,k] /
      (1 + exp(-slope_raw[i,j,k] *
               (Intensity[n] - theta_raw[i,j,k])));

    mu_log[n] = log(mu_curve[n] + 1e-9);
  }
}

model {
  // Hyperpriors
  mu_theta ~ normal(prior_mu_theta_mean, prior_mu_theta_sd);
  mu_slope ~ normal(prior_mu_slope_mean, prior_mu_slope_sd);
  mu_logA  ~ normal(log(prior_mu_A_mean), prior_mu_A_sd);

  sigma_theta ~ student_t(3, prior_sigma_theta_loc, prior_sigma_theta_scale);
  sigma_slope ~ student_t(3, prior_sigma_slope_loc, prior_sigma_slope_scale);
  sigma_logA  ~ student_t(3, prior_sigma_A_loc,    prior_sigma_A_scale);
  sigma       ~ student_t(3, prior_sigma_loc,      prior_sigma_scale);

  // Subject random effects
  theta_z ~ normal(0,1);
  slope_z ~ normal(0,1);
  logA_z  ~ normal(0,1);

  // RAW coil/side fixed effects
  delta_theta_coil_raw ~ normal(prior_delta_theta_coil_loc,
                                prior_delta_theta_coil_scale);
  delta_theta_side_raw ~ normal(prior_delta_theta_side_loc,
                                prior_delta_theta_side_scale);

  delta_slope_coil_raw ~ normal(prior_delta_slope_coil_loc,
                                prior_delta_slope_coil_scale);
  delta_slope_side_raw ~ normal(prior_delta_slope_side_loc,
                                prior_delta_slope_side_scale);

  delta_logA_coil_raw ~ normal(prior_delta_logA_coil_loc,
                               prior_delta_logA_coil_scale);
  delta_logA_side_raw ~ normal(prior_delta_logA_side_loc,
                               prior_delta_logA_side_scale);

  // Likelihood
  if (use_likelihood == 1)
    for (n in 1:N) {
      int j = coil[n];
      if (Intensity[n] <= MSO_limit[j])
        Y[n] ~ lognormal(mu_log[n], sigma);
    }
}

generated quantities {
  vector[N] Y_rep;
  vector[N] Y_prior;

  // -----------------------------
  // POSTERIOR PREDICTIVE
  // -----------------------------
  // Save full IO-curve parameters for plotting
  array[I, J, K] real theta_out;
  array[I, J, K] real slope_out;
  array[I, J, K] real A_out;

  for (i in 1:I)
    for (j in 1:J)
      for (k in 1:K) {
        theta_out[i,j,k] = theta_raw[i,j,k];
        slope_out[i,j,k] = slope_raw[i,j,k];
        A_out[i,j,k]     = A_raw[i,j,k];
      }


  for (n in 1:N) {
    int j = coil[n];
    if (Intensity[n] <= MSO_limit[j])
      Y_rep[n] = lognormal_rng(mu_log[n], sigma);
    else
      Y_rep[n] = -1;
  }

  // -----------------------------
  // PRIOR PREDICTIVE
  // -----------------------------
  {
    // Draw hyperparameters from priors
    real mu_theta_prior = normal_rng(prior_mu_theta_mean, prior_mu_theta_sd);
    real mu_slope_prior = normal_rng(prior_mu_slope_mean, prior_mu_slope_sd);
    real mu_logA_prior  = normal_rng(log(prior_mu_A_mean), prior_mu_A_sd);

    real sigma_theta_prior =
      fmin(sigma_ubound, abs(student_t_rng(3, prior_sigma_theta_loc, prior_sigma_theta_scale)));
    real sigma_slope_prior =
      fmin(sigma_ubound, abs(student_t_rng(3, prior_sigma_slope_loc, prior_sigma_slope_scale)));
    real sigma_logA_prior =
      fmin(sigma_ubound, abs(student_t_rng(3, prior_sigma_A_loc, prior_sigma_A_scale)));
    real sigma_prior =
      fmin(sigma_ubound, abs(student_t_rng(3, prior_sigma_loc, prior_sigma_scale)));

    // Raw coil/side effects from priors
    vector[J] delta_theta_coil_prior_raw;
    vector[K] delta_theta_side_prior_raw;

    vector[J] delta_slope_coil_prior_raw;
    vector[K] delta_slope_side_prior_raw;

    vector[J] delta_logA_coil_prior_raw;
    vector[K] delta_logA_side_prior_raw;

    for (j in 1:J) {
      delta_theta_coil_prior_raw[j] = normal_rng(prior_delta_theta_coil_loc,
                                                 prior_delta_theta_coil_scale);
      delta_slope_coil_prior_raw[j] = normal_rng(prior_delta_slope_coil_loc,
                                                 prior_delta_slope_coil_scale);
      delta_logA_coil_prior_raw[j]  = normal_rng(prior_delta_logA_coil_loc,
                                                 prior_delta_logA_coil_scale);
    }

    for (k in 1:K) {
      delta_theta_side_prior_raw[k] = normal_rng(prior_delta_theta_side_loc,
                                                 prior_delta_theta_side_scale);
      delta_slope_side_prior_raw[k] = normal_rng(prior_delta_slope_side_loc,
                                                 prior_delta_slope_side_scale);
      delta_logA_side_prior_raw[k]  = normal_rng(prior_delta_logA_side_loc,
                                                 prior_delta_logA_side_scale);
    }

    // Center them (sum-to-zero)
    vector[J] delta_theta_coil_prior =
      delta_theta_coil_prior_raw - mean(delta_theta_coil_prior_raw);
    vector[K] delta_theta_side_prior =
      delta_theta_side_prior_raw - mean(delta_theta_side_prior_raw);

    vector[J] delta_slope_coil_prior =
      delta_slope_coil_prior_raw - mean(delta_slope_coil_prior_raw);
    vector[K] delta_slope_side_prior =
      delta_slope_side_prior_raw - mean(delta_slope_side_prior_raw);

    vector[J] delta_logA_coil_prior =
      delta_logA_coil_prior_raw - mean(delta_logA_coil_prior_raw);
    vector[K] delta_logA_side_prior =
      delta_logA_side_prior_raw - mean(delta_logA_side_prior_raw);

    // Subject-level prior z's
    vector[I] theta_z_prior;
    vector[I] slope_z_prior;
    vector[I] logA_z_prior;

    for (i in 1:I) {
      theta_z_prior[i] = normal_rng(0,1);
      slope_z_prior[i] = normal_rng(0,1);
      logA_z_prior[i]  = normal_rng(0,1);
    }

    // Build prior IO-curve parameters
    array[I, J, K] real theta_prior;
    array[I, J, K] real slope_prior;
    array[I, J, K] real logA_prior;
    array[I, J, K] real A_prior;

    for (i in 1:I)
      for (j in 1:J)
        for (k in 1:K) {

          logA_prior[i,j,k] =
            mu_logA_prior +
            delta_logA_coil_prior[j] +
            delta_logA_side_prior[k] +
            sigma_logA_prior * logA_z_prior[i];

          // cap to avoid overflow
          real logA_prior_capped = fmin(logA_prior[i,j,k], logA_cap);
          A_prior[i,j,k] = exp(logA_prior_capped);

          theta_prior[i,j,k] =
            mu_theta_prior +
            delta_theta_coil_prior[j] +
            delta_theta_side_prior[k] +
            sigma_theta_prior * theta_z_prior[i];

          slope_prior[i,j,k] =
            mu_slope_prior +
            delta_slope_coil_prior[j] +
            delta_slope_side_prior[k] +
            sigma_slope_prior * slope_z_prior[i];
        }

    // PRIOR predictive
    for (n in 1:N) {
      int i = subj[n];
      int j = coil[n];
      int k = side[n];

      if (Intensity[n] <= MSO_limit[j]) {

        real mu_curve_prior =
          A_prior[i,j,k] /
          (1 + exp(-slope_prior[i,j,k] *
                   (Intensity[n] - theta_prior[i,j,k])));

        real mu_log_prior = log(mu_curve_prior + 1e-9);

        Y_prior[n] = lognormal_rng(mu_log_prior, sigma_prior);

      } else {
        Y_prior[n] = -1;
      }
    }
  }
}