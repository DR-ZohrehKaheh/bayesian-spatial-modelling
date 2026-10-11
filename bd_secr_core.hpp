// Project supervisor: Dr. Zohreh Kaheh.
// Undergraduate project collaborator: Jingxuan Hu.
// Initial student coding was developed with supervisory teaching and improvement guidance.
#ifndef BD_SECR_CORE_HPP
#define BD_SECR_CORE_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace bdsecr {

inline double sqr(double value) { return value * value; }

inline double squared_distance(
    double x1, double y1, double x2, double y2) {
  return sqr(x1 - x2) + sqr(y1 - y2);
}

inline double normal_log_density(double x, double mean, double sd) {
  const double z = (x - mean) / sd;
  return -0.5 * z * z - std::log(sd) - 0.5 * std::log(2.0 * M_PI);
}

inline double log_factorial(int value) {
  return std::lgamma(static_cast<double>(value) + 1.0);
}

struct SamplerConfig {
  int n_iter = 20000;
  int burnin = 5000;
  int thin = 5;

  std::vector<double> x_traps;
  std::vector<double> y_traps;
  std::vector<int> y_obs;
  std::vector<int> habitat;
  std::vector<double> effort;
  std::vector<double> mask_x;
  std::vector<double> mask_y;

  double res_primary = 1.0;
  double res_secondary = 1.0;
  int init_N = 10;
  double init_sigma = 1.0;
  double init_beta0 = 0.0;
  double init_beta1 = 0.0;

  double prior_mean_N = 20.0;
  double prior_log_sigma_mean = std::log(2.0);
  double prior_log_sigma_sd = 1.0;
  double prior_beta0_mean = 0.0;
  double prior_beta1_mean = 0.0;
  double prior_beta0_sd = 3.16;
  double prior_beta1_sd = 3.16;

  std::uint64_t seed = 1;
  int n_location_updates = 1;
  double log_sigma_step = 0.08;
  double beta0_step = 0.20;
  double beta1_step = 0.20;
  double location_step = 0.50;
  int rebuild_every = 250;
  bool store_density = false;
  bool store_ppc = false;
};

struct MoveCounts {
  std::int64_t sigma = 0;
  std::int64_t beta0 = 0;
  std::int64_t beta1 = 0;
  std::int64_t location = 0;
  std::int64_t birth = 0;
  std::int64_t death = 0;
};

struct SamplerResult {
  std::vector<int> N;
  std::vector<double> sigma;
  std::vector<double> beta0;
  std::vector<double> beta1;
  std::vector<double> log_likelihood;
  std::vector<double> density;
  std::vector<double> ppc_total, ppc_zeros, ppc_pearson_observed, ppc_pearson_replicate;
  MoveCounts proposed;
  MoveCounts accepted;
  int density_columns = 0;
  int final_N = 0;
  double final_sigma = 0.0;
  double final_beta0 = 0.0;
  double final_beta1 = 0.0;
  double final_log_likelihood = -std::numeric_limits<double>::infinity();
  double cache_max_abs_error = 0.0;
  double xmin = 0.0;
  double xmax = 0.0;
  double ymin = 0.0;
  double ymax = 0.0;
};

inline void validate_config(const SamplerConfig& config) {
  const std::size_t detectors = config.x_traps.size();
  if (config.n_iter <= 0) throw std::invalid_argument("n_iter must be positive");
  if (config.burnin < 0 || config.burnin >= config.n_iter) {
    throw std::invalid_argument("burnin must satisfy 0 <= burnin < n_iter");
  }
  if (config.thin <= 0) throw std::invalid_argument("thin must be positive");
  if (detectors == 0) throw std::invalid_argument("at least one detector is required");
  if (config.y_traps.size() != detectors || config.y_obs.size() != detectors ||
      config.habitat.size() != detectors || config.effort.size() != detectors) {
    throw std::invalid_argument(
        "detector coordinates, counts, habitat and effort must have equal lengths");
  }
  if (config.mask_x.empty() || config.mask_x.size() != config.mask_y.size()) {
    throw std::invalid_argument("mask coordinates must be non-empty and paired");
  }
  if (config.init_N < 0) throw std::invalid_argument("init_N must be non-negative");
  if (!std::isfinite(config.init_sigma) || config.init_sigma <= 0.0) {
    throw std::invalid_argument("init_sigma must be finite and positive");
  }
  if (!std::isfinite(config.init_beta0) || !std::isfinite(config.init_beta1)) {
    throw std::invalid_argument("initial beta values must be finite");
  }
  if (!std::isfinite(config.prior_mean_N) || config.prior_mean_N <= 0.0) {
    throw std::invalid_argument("prior_mean_N must be finite and positive");
  }
  if (!std::isfinite(config.prior_log_sigma_mean) ||
      !std::isfinite(config.prior_beta0_mean) ||
      !std::isfinite(config.prior_beta1_mean)) {
    throw std::invalid_argument("normal-prior means must be finite");
  }
  if (!std::isfinite(config.prior_log_sigma_sd) ||
      config.prior_log_sigma_sd <= 0.0 ||
      !std::isfinite(config.prior_beta0_sd) || config.prior_beta0_sd <= 0.0 ||
      !std::isfinite(config.prior_beta1_sd) || config.prior_beta1_sd <= 0.0) {
    throw std::invalid_argument("prior standard deviations must be finite and positive");
  }
  if (!std::isfinite(config.res_primary) || config.res_primary <= 0.0 ||
      !std::isfinite(config.res_secondary) || config.res_secondary <= 0.0) {
    throw std::invalid_argument("distance multipliers must be finite and positive");
  }
  if (config.n_location_updates < 0) {
    throw std::invalid_argument("n_location_updates must be non-negative");
  }
  if (!std::isfinite(config.log_sigma_step) || config.log_sigma_step <= 0.0 ||
      !std::isfinite(config.beta0_step) || config.beta0_step <= 0.0 ||
      !std::isfinite(config.beta1_step) || config.beta1_step <= 0.0 ||
      !std::isfinite(config.location_step) || config.location_step <= 0.0) {
    throw std::invalid_argument("proposal step sizes must be finite and positive");
  }
  if (config.rebuild_every < 0) {
    throw std::invalid_argument("rebuild_every must be non-negative");
  }

  for (std::size_t j = 0; j < detectors; ++j) {
    if (!std::isfinite(config.x_traps[j]) || !std::isfinite(config.y_traps[j])) {
      throw std::invalid_argument("detector coordinates must be finite");
    }
    if (config.y_obs[j] < 0) {
      throw std::invalid_argument("counts must be non-negative integers");
    }
    if (config.habitat[j] != 0 && config.habitat[j] != 1) {
      throw std::invalid_argument(
          "habitat must be binary 0/1; recode raw categories explicitly");
    }
    if (!std::isfinite(config.effort[j]) || config.effort[j] < 0.0) {
      throw std::invalid_argument("effort must be finite and non-negative");
    }
    if (config.effort[j] == 0.0 && config.y_obs[j] > 0) {
      throw std::invalid_argument("positive counts cannot have zero effort");
    }
  }
  for (std::size_t k = 0; k < config.mask_x.size(); ++k) {
    if (!std::isfinite(config.mask_x[k]) || !std::isfinite(config.mask_y[k])) {
      throw std::invalid_argument("mask coordinates must be finite");
    }
  }
}

inline SamplerResult fit(const SamplerConfig& config) {
  validate_config(config);
  const int detectors = static_cast<int>(config.x_traps.size());
  const int mask_size = static_cast<int>(config.mask_x.size());
  const double xmin = *std::min_element(config.mask_x.begin(), config.mask_x.end());
  const double xmax = *std::max_element(config.mask_x.begin(), config.mask_x.end());
  const double ymin = *std::min_element(config.mask_y.begin(), config.mask_y.end());
  const double ymax = *std::max_element(config.mask_y.begin(), config.mask_y.end());
  if (!(xmax > xmin) || !(ymax > ymin)) {
    throw std::invalid_argument("the rectangular state space must have positive area");
  }

  const int n_save = 1 + (config.n_iter - 1 - config.burnin) / config.thin;
  SamplerResult result;
  result.N.resize(n_save);
  result.sigma.resize(n_save);
  result.beta0.resize(n_save);
  result.beta1.resize(n_save);
  result.log_likelihood.resize(n_save);
  result.xmin = xmin;
  result.xmax = xmax;
  result.ymin = ymin;
  result.ymax = ymax;
  if (config.store_density) {
    result.density_columns = mask_size;
    result.density.resize(static_cast<std::size_t>(n_save) * mask_size);
  }

  std::mt19937_64 generator(config.seed);
  std::mt19937_64 predictive_generator(config.seed ^ 0x9e3779b97f4a7c15ULL);
  std::uniform_real_distribution<double> uniform01(0.0, 1.0);
  std::uniform_real_distribution<double> uniform_x(xmin, xmax);
  std::uniform_real_distribution<double> uniform_y(ymin, ymax);
  std::normal_distribution<double> standard_normal(0.0, 1.0);

  auto uniform_open = [&]() {
    return std::max(uniform01(generator), std::numeric_limits<double>::min());
  };
  auto accept_log_ratio = [&](double log_ratio) {
    // A finite proposal against an impossible current state has log ratio +Inf
    // and must be accepted; only NaN denotes an invalid comparison.
    return !std::isnan(log_ratio) &&
           (log_ratio >= 0.0 || std::log(uniform_open()) < log_ratio);
  };
  auto birth_probability = [](int abundance) {
    return abundance == 0 ? 1.0 : 0.5;
  };
  auto death_probability = [&](int abundance) {
    return 1.0 - birth_probability(abundance);
  };

  int abundance = config.init_N;
  double sigma = config.init_sigma;
  double beta0 = config.init_beta0;
  double beta1 = config.init_beta1;
  std::vector<double> centres_x(abundance);
  std::vector<double> centres_y(abundance);
  for (int individual = 0; individual < abundance; ++individual) {
    centres_x[individual] = uniform_x(generator);
    centres_y[individual] = uniform_y(generator);
  }

  auto contribution = [&](double x, double y, int detector, double scale,
                          double intercept, double habitat_effect) {
    if (config.effort[detector] == 0.0) return 0.0;
    const double multiplier = config.habitat[detector] == 1
                                  ? config.res_primary
                                  : config.res_secondary;
    const double log_value = std::log(config.effort[detector]) + intercept +
                             habitat_effect * config.habitat[detector] -
                             multiplier * multiplier *
                                 squared_distance(x, y, config.x_traps[detector],
                                                  config.y_traps[detector]) /
                                 (2.0 * scale * scale);
    if (log_value > std::log(std::numeric_limits<double>::max())) {
      return std::numeric_limits<double>::infinity();
    }
    return std::exp(log_value);
  };

  auto rebuild_lambda = [&](double scale, double intercept,
                            double habitat_effect) {
    std::vector<double> values(detectors, 0.0);
    for (int detector = 0; detector < detectors; ++detector) {
      for (int individual = 0; individual < abundance; ++individual) {
        values[detector] += contribution(
            centres_x[individual], centres_y[individual], detector, scale,
            intercept, habitat_effect);
      }
    }
    return values;
  };

  auto log_likelihood = [&](const std::vector<double>& lambda) {
    double value = 0.0;
    for (int detector = 0; detector < detectors; ++detector) {
      double intensity = lambda[detector];
      if (intensity < 0.0 && intensity > -1e-12) intensity = 0.0;
      if (!std::isfinite(intensity) || intensity < 0.0) {
        return -std::numeric_limits<double>::infinity();
      }
      if (intensity == 0.0) {
        if (config.y_obs[detector] > 0) {
          return -std::numeric_limits<double>::infinity();
        }
        continue;
      }
      value += config.y_obs[detector] * std::log(intensity) - intensity -
               log_factorial(config.y_obs[detector]);
    }
    return value;
  };

  std::vector<double> lambda = rebuild_lambda(sigma, beta0, beta1);
  double current_log_likelihood = log_likelihood(lambda);
  int save_index = 0;

  for (int iteration = 0; iteration < config.n_iter; ++iteration) {
    ++result.proposed.sigma;
    const double log_sigma = std::log(sigma);
    const double proposed_log_sigma =
        log_sigma + config.log_sigma_step * standard_normal(generator);
    const double proposed_sigma = std::exp(proposed_log_sigma);
    if (std::isfinite(proposed_sigma) && proposed_sigma > 0.0) {
      std::vector<double> proposed_lambda =
          rebuild_lambda(proposed_sigma, beta0, beta1);
      const double proposed_log_likelihood = log_likelihood(proposed_lambda);
      const double log_ratio =
          proposed_log_likelihood - current_log_likelihood +
          normal_log_density(proposed_log_sigma,
                             config.prior_log_sigma_mean,
                             config.prior_log_sigma_sd) -
          normal_log_density(log_sigma, config.prior_log_sigma_mean,
                             config.prior_log_sigma_sd);
      if (accept_log_ratio(log_ratio)) {
        sigma = proposed_sigma;
        lambda.swap(proposed_lambda);
        current_log_likelihood = proposed_log_likelihood;
        ++result.accepted.sigma;
      }
    }

    ++result.proposed.beta0;
    const double proposed_beta0 =
        beta0 + config.beta0_step * standard_normal(generator);
    {
      const double factor = std::exp(proposed_beta0 - beta0);
      std::vector<double> proposed_lambda(lambda.size());
      for (int detector = 0; detector < detectors; ++detector) {
        proposed_lambda[detector] = lambda[detector] * factor;
      }
      const double proposed_log_likelihood = log_likelihood(proposed_lambda);
      const double log_ratio =
          proposed_log_likelihood - current_log_likelihood +
          normal_log_density(proposed_beta0, config.prior_beta0_mean,
                             config.prior_beta0_sd) -
          normal_log_density(beta0, config.prior_beta0_mean,
                             config.prior_beta0_sd);
      if (accept_log_ratio(log_ratio)) {
        beta0 = proposed_beta0;
        lambda.swap(proposed_lambda);
        current_log_likelihood = proposed_log_likelihood;
        ++result.accepted.beta0;
      }
    }

    ++result.proposed.beta1;
    const double proposed_beta1 =
        beta1 + config.beta1_step * standard_normal(generator);
    {
      const double factor = std::exp(proposed_beta1 - beta1);
      std::vector<double> proposed_lambda = lambda;
      for (int detector = 0; detector < detectors; ++detector) {
        if (config.habitat[detector] == 1) {
          proposed_lambda[detector] *= factor;
        }
      }
      const double proposed_log_likelihood = log_likelihood(proposed_lambda);
      const double log_ratio =
          proposed_log_likelihood - current_log_likelihood +
          normal_log_density(proposed_beta1, config.prior_beta1_mean,
                             config.prior_beta1_sd) -
          normal_log_density(beta1, config.prior_beta1_mean,
                             config.prior_beta1_sd);
      if (accept_log_ratio(log_ratio)) {
        beta1 = proposed_beta1;
        lambda.swap(proposed_lambda);
        current_log_likelihood = proposed_log_likelihood;
        ++result.accepted.beta1;
      }
    }

    for (int move = 0;
         move < config.n_location_updates && abundance > 0; ++move) {
      ++result.proposed.location;
      std::uniform_int_distribution<int> individual_distribution(0,
                                                                  abundance - 1);
      const int individual = individual_distribution(generator);
      const double proposed_x =
          centres_x[individual] + config.location_step * standard_normal(generator);
      const double proposed_y =
          centres_y[individual] + config.location_step * standard_normal(generator);
      if (proposed_x <= xmin || proposed_x >= xmax || proposed_y <= ymin ||
          proposed_y >= ymax) {
        continue;
      }
      std::vector<double> proposed_lambda = lambda;
      for (int detector = 0; detector < detectors; ++detector) {
        proposed_lambda[detector] +=
            contribution(proposed_x, proposed_y, detector, sigma, beta0, beta1) -
            contribution(centres_x[individual], centres_y[individual], detector,
                         sigma, beta0, beta1);
      }
      const double proposed_log_likelihood = log_likelihood(proposed_lambda);
      if (accept_log_ratio(proposed_log_likelihood - current_log_likelihood)) {
        centres_x[individual] = proposed_x;
        centres_y[individual] = proposed_y;
        lambda.swap(proposed_lambda);
        current_log_likelihood = proposed_log_likelihood;
        ++result.accepted.location;
      }
    }

    const bool birth = abundance == 0 || uniform01(generator) < 0.5;
    if (birth) {
      ++result.proposed.birth;
      const double new_x = uniform_x(generator);
      const double new_y = uniform_y(generator);
      std::vector<double> proposed_lambda = lambda;
      for (int detector = 0; detector < detectors; ++detector) {
        proposed_lambda[detector] +=
            contribution(new_x, new_y, detector, sigma, beta0, beta1);
      }
      const double proposed_log_likelihood = log_likelihood(proposed_lambda);
      const double log_ratio =
          proposed_log_likelihood - current_log_likelihood +
          std::log(config.prior_mean_N) - std::log(abundance + 1.0) +
          std::log(death_probability(abundance + 1)) -
          std::log(birth_probability(abundance));
      if (accept_log_ratio(log_ratio)) {
        centres_x.push_back(new_x);
        centres_y.push_back(new_y);
        ++abundance;
        lambda.swap(proposed_lambda);
        current_log_likelihood = proposed_log_likelihood;
        ++result.accepted.birth;
      }
    } else {
      ++result.proposed.death;
      std::uniform_int_distribution<int> individual_distribution(0,
                                                                  abundance - 1);
      const int individual = individual_distribution(generator);
      std::vector<double> proposed_lambda = lambda;
      for (int detector = 0; detector < detectors; ++detector) {
        proposed_lambda[detector] -= contribution(
            centres_x[individual], centres_y[individual], detector, sigma, beta0,
            beta1);
      }
      const double proposed_log_likelihood = log_likelihood(proposed_lambda);
      const double log_ratio =
          proposed_log_likelihood - current_log_likelihood +
          std::log(static_cast<double>(abundance)) -
          std::log(config.prior_mean_N) +
          std::log(birth_probability(abundance - 1)) -
          std::log(death_probability(abundance));
      if (accept_log_ratio(log_ratio)) {
        centres_x[individual] = centres_x.back();
        centres_y[individual] = centres_y.back();
        centres_x.pop_back();
        centres_y.pop_back();
        --abundance;
        lambda.swap(proposed_lambda);
        current_log_likelihood = proposed_log_likelihood;
        ++result.accepted.death;
      }
    }

    if (config.rebuild_every > 0 &&
        (iteration + 1) % config.rebuild_every == 0) {
      std::vector<double> rebuilt = rebuild_lambda(sigma, beta0, beta1);
      for (int detector = 0; detector < detectors; ++detector) {
        result.cache_max_abs_error =
            std::max(result.cache_max_abs_error,
                     std::abs(rebuilt[detector] - lambda[detector]));
      }
      lambda.swap(rebuilt);
      current_log_likelihood = log_likelihood(lambda);
    }

    if (iteration >= config.burnin &&
        (iteration - config.burnin) % config.thin == 0) {
      result.N[save_index] = abundance;
      result.sigma[save_index] = sigma;
      result.beta0[save_index] = beta0;
      result.beta1[save_index] = beta1;
      result.log_likelihood[save_index] = current_log_likelihood;
      if (config.store_density) {
        for (int mask = 0; mask < mask_size; ++mask) {
          double value = 0.0;
          for (int individual = 0; individual < abundance; ++individual) {
            value +=
                std::exp(-squared_distance(centres_x[individual],
                                           centres_y[individual],
                                           config.mask_x[mask],
                                           config.mask_y[mask]) /
                         (2.0 * sigma * sigma)) /
                (2.0 * M_PI * sigma * sigma);
          }
          result.density[static_cast<std::size_t>(save_index) * mask_size + mask] =
              value;
        }
      }
      if (config.store_ppc && save_index % 60 == 0) {
        double total=0, zeros=0, observed=0, replicated=0;
        for (int j=0;j<detectors;++j) {
          const double mu=std::max(0.0, lambda[j]);
          std::poisson_distribution<int> poisson(mu);
          const int yrep=poisson(predictive_generator);
          total+=yrep; zeros+=(yrep==0);
          observed+=sqr(config.y_obs[j]-mu)/(mu+0.5);
          replicated+=sqr(yrep-mu)/(mu+0.5);
        }
        result.ppc_total.push_back(total); result.ppc_zeros.push_back(zeros);
        result.ppc_pearson_observed.push_back(observed);
        result.ppc_pearson_replicate.push_back(replicated);
      }
      ++save_index;
    }
  }

  result.final_N = abundance;
  result.final_sigma = sigma;
  result.final_beta0 = beta0;
  result.final_beta1 = beta1;
  result.final_log_likelihood = current_log_likelihood;
  return result;
}

inline double acceptance_rate(std::int64_t accepted, std::int64_t proposed) {
  return proposed > 0 ? static_cast<double>(accepted) / proposed : 0.0;
}

}  // namespace bdsecr

#endif
