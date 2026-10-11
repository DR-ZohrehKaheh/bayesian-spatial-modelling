// Project supervisor: Dr. Zohreh Kaheh.
// Undergraduate project collaborator: Jingxuan Hu.
// Initial student coding was developed with supervisory teaching and improvement guidance.
#include "bd_secr_core.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace {

double mean(const std::vector<int>& values) {
  return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

}  // namespace

int main() {
  bdsecr::SamplerConfig config;
  config.n_iter = 300000;
  config.burnin = 30000;
  config.thin = 3;
  config.x_traps = {0.25};
  config.y_traps = {0.25};
  config.y_obs = {0};
  config.habitat = {0};
  config.effort = {0.0};
  config.mask_x = {0.0, 0.0, 1.0, 1.0};
  config.mask_y = {0.0, 1.0, 0.0, 1.0};
  config.init_N = 0;
  config.init_sigma = 0.2;
  config.init_beta0 = 0.0;
  config.init_beta1 = 0.0;
  config.prior_mean_N = 7.5;
  config.prior_log_sigma_mean = std::log(0.2);
  config.prior_log_sigma_sd = 0.5;
  config.seed = 7319;
  config.n_location_updates = 0;
  config.store_density = false;

  const bdsecr::SamplerResult fit = bdsecr::fit(config);
  require(static_cast<int>(fit.N.size()) ==
              1 + (config.n_iter - 1 - config.burnin) / config.thin,
          "retained-draw count is incorrect");
  require(std::abs(mean(fit.N) - 7.5) < 0.06,
          "boundary-sensitive prior-only mean check failed");
  require(fit.cache_max_abs_error < 1e-10,
          "cached intensity drift exceeded tolerance");

  bdsecr::SamplerConfig invalid = config;
  invalid.habitat = {2};
  bool rejected = false;
  try {
    (void)bdsecr::fit(invalid);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  require(rejected, "non-binary habitat input was not rejected");

  for (int field=0;field<3;++field) {
    for (double bad : {std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::infinity(),
                       -std::numeric_limits<double>::infinity()}) {
      invalid=config;
      if(field==0) invalid.prior_log_sigma_mean=bad;
      if(field==1) invalid.prior_beta0_mean=bad;
      if(field==2) invalid.prior_beta1_mean=bad;
      rejected=false;
      try { (void)bdsecr::fit(invalid); }
      catch(const std::invalid_argument&) { rejected=true; }
      require(rejected,"non-finite prior mean was not rejected");
    }
  }
  std::cout << "core_compile_and_smoke_test=PASS\n";
  std::cout << "retained_draws=" << fit.N.size() << "\n";
  std::cout << "prior_only_target_mean=7.5\n";
  std::cout << "prior_only_sample_mean=" << mean(fit.N) << "\n";
  std::cout << "cache_max_abs_error=" << fit.cache_max_abs_error << "\n";
  return 0;
}
