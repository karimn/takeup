hyper_mu ~ normal(0, 2);

raw_strata_mu ~ std_normal();
strata_mu_tau ~ normal(0, 1);

sigma ~ normal(0, 1); 

{ 
  int wtp_stratum_pos = 1;

  for (stratum_index in 1:num_strata) {
    int curr_wtp_stratum_size = wtp_strata_sizes[stratum_index];
    int wtp_stratum_end = wtp_stratum_pos + curr_wtp_stratum_size - 1;

    for (wtp_obs_index in wtp_stratum_pos:wtp_stratum_end) {
      if (wtp_response[wtp_obs_index] == -1) { // Decided to KEEP
        if (gift_choice[wtp_obs_index] == -1) { // Initial choice: BRACELET
          target += log(Phi_approx((- scaled_wtp_offer[wtp_obs_index] - strata_mu[stratum_index]) / sigma));
        } else { // Initial choice: CALENDAR
          target += log(1 - Phi_approx((scaled_wtp_offer[wtp_obs_index] - strata_mu[stratum_index]) / sigma));
        }
      } else { // Decided to SWITCH
        target += log(gift_choice[wtp_obs_index] *
          (Phi_approx((gift_choice[wtp_obs_index] * scaled_wtp_offer[wtp_obs_index] - strata_mu[stratum_index]) / sigma)
           - Phi_approx(- strata_mu[stratum_index] / sigma)));
      }
    }

    wtp_stratum_pos = wtp_stratum_end + 1;
  }
}