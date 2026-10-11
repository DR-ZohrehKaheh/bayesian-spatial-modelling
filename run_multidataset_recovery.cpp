// Project supervisor: Dr. Zohreh Kaheh.
// Undergraduate project collaborator: Jingxuan Hu.
// Initial student coding was developed with supervisory teaching and improvement guidance.
#include "bd_secr_core.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Truth {
  int N = 20;
  double sigma = 0.10;
  double beta0 = -0.35;
  double beta1 = 0.45;
  double res_primary = 1.0;
  double res_secondary = 1.25;
};

struct Dataset {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<int> habitat;
  std::vector<int> counts;
  std::vector<double> effort;
  std::vector<double> centre_x;
  std::vector<double> centre_y;
};

double individual_contribution(double centre_x, double centre_y,
                               double trap_x, double trap_y, int habitat,
                               const Truth& truth) {
  const double multiplier = habitat == 1 ? truth.res_primary
                                          : truth.res_secondary;
  const double distance2 = bdsecr::squared_distance(
      centre_x, centre_y, trap_x, trap_y);
  return std::exp(truth.beta0 + truth.beta1 * habitat -
                  multiplier * multiplier * distance2 /
                      (2.0 * truth.sigma * truth.sigma));
}

Dataset simulate_dataset(std::uint64_t seed, const Truth& truth,
                         int grid_side, int occasions) {
  std::mt19937_64 generator(seed);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  Dataset data;
  data.centre_x.resize(truth.N);
  data.centre_y.resize(truth.N);
  for (int individual = 0; individual < truth.N; ++individual) {
    data.centre_x[individual] = uniform(generator);
    data.centre_y[individual] = uniform(generator);
  }

  for (int occasion = 0; occasion < occasions; ++occasion) {
    for (int row = 0; row < grid_side; ++row) {
      for (int column = 0; column < grid_side; ++column) {
        const double trap_x = 0.10 + 0.80 * column / (grid_side - 1.0);
        const double trap_y = 0.10 + 0.80 * row / (grid_side - 1.0);
        const int habitat = (row + column) % 2;
        double intensity = 0.0;
        for (int individual = 0; individual < truth.N; ++individual) {
          intensity += individual_contribution(
              data.centre_x[individual], data.centre_y[individual], trap_x,
              trap_y, habitat, truth);
        }
        std::poisson_distribution<int> poisson(intensity);
        data.x.push_back(trap_x);
        data.y.push_back(trap_y);
        data.habitat.push_back(habitat);
        data.counts.push_back(poisson(generator));
        data.effort.push_back(1.0);
      }
    }
  }
  return data;
}

void write_header(std::ofstream& stream, const std::string& header) {
  if (!stream) throw std::runtime_error("could not open output file");
  stream << header << '\n';
  stream << std::setprecision(12);
}

}  // namespace

int main(int argc, char** argv) {
  const fs::path output_directory =
      argc > 1 ? fs::path(argv[1]) : fs::path("results/full_recovery");
  fs::create_directories(output_directory);

  const Truth truth;
  const int datasets = 6;
  const int chains = 2;
  const int grid_side = 7;
  const int occasions = 5;
  const int iterations = 96000;
  const int burnin = 24000;
  const int thin = 12;

  std::ofstream traces(output_directory / "recovery_traces.csv");
  std::ofstream acceptance(output_directory / "recovery_acceptance.csv");
  std::ofstream counts(output_directory / "synthetic_counts.csv");
  std::ofstream centres(output_directory / "synthetic_centres.csv");
  std::ofstream design(output_directory / "simulation_design.csv");
  write_header(traces,
               "dataset,chain,draw,N,sigma,beta0,beta1,log_likelihood");
  write_header(acceptance,
               "dataset,chain,sigma,beta0,beta1,location,birth,death,"
               "birth_death,cache_max_abs_error");
  write_header(counts,
               "dataset,occasion,trap,x,y,habitat,count,effort");
  write_header(centres, "dataset,individual,x,y");
  write_header(design, "quantity,value");

  design << "datasets," << datasets << '\n'
         << "chains_per_dataset," << chains << '\n'
         << "grid_side," << grid_side << '\n'
         << "detectors_per_occasion," << grid_side * grid_side << '\n'
         << "occasions," << occasions << '\n'
         << "rows_per_dataset," << grid_side * grid_side * occasions << '\n'
         << "iterations," << iterations << '\n'
         << "burnin," << burnin << '\n'
         << "thin," << thin << '\n'
         << "retained_draws_per_chain,"
         << 1 + (iterations - 1 - burnin) / thin << '\n'
         << "true_N," << truth.N << '\n'
         << "true_sigma," << truth.sigma << '\n'
         << "true_beta0," << truth.beta0 << '\n'
         << "true_beta1," << truth.beta1 << '\n'
         << "res_primary," << truth.res_primary << '\n'
         << "res_secondary," << truth.res_secondary << '\n'
         << "prior_mean_N,22\n"
         << "prior_log_sigma_mean," << std::log(0.12) << '\n'
         << "prior_log_sigma_sd,0.45\n"
         << "prior_beta_sd,1.5\n";

  for (int dataset = 1; dataset <= datasets; ++dataset) {
    const std::uint64_t dataset_seed = 2026080900ULL + dataset;
    const Dataset simulated =
        simulate_dataset(dataset_seed, truth, grid_side, occasions);
    for (int individual = 0; individual < truth.N; ++individual) {
      centres << dataset << ',' << individual + 1 << ','
              << simulated.centre_x[individual] << ','
              << simulated.centre_y[individual] << '\n';
    }
    for (std::size_t row = 0; row < simulated.counts.size(); ++row) {
      const int occasion = static_cast<int>(row) / (grid_side * grid_side) + 1;
      const int trap = static_cast<int>(row) % (grid_side * grid_side) + 1;
      counts << dataset << ',' << occasion << ',' << trap << ','
             << simulated.x[row] << ',' << simulated.y[row] << ','
             << simulated.habitat[row] << ',' << simulated.counts[row] << ','
             << simulated.effort[row] << '\n';
    }

    for (int chain = 1; chain <= chains; ++chain) {
      bdsecr::SamplerConfig config;
      config.n_iter = iterations;
      config.burnin = burnin;
      config.thin = thin;
      config.x_traps = simulated.x;
      config.y_traps = simulated.y;
      config.y_obs = simulated.counts;
      config.habitat = simulated.habitat;
      config.effort = simulated.effort;
      config.mask_x = {0.0, 0.0, 1.0, 1.0};
      config.mask_y = {0.0, 1.0, 0.0, 1.0};
      config.res_primary = truth.res_primary;
      config.res_secondary = truth.res_secondary;
      config.init_N = chain == 1 ? 12 : 30;
      config.init_sigma = chain == 1 ? 0.07 : 0.16;
      config.init_beta0 = chain == 1 ? -0.8 : 0.2;
      config.init_beta1 = chain == 1 ? 0.0 : 0.8;
      config.prior_mean_N = 22.0;
      config.prior_log_sigma_mean = std::log(0.12);
      config.prior_log_sigma_sd = 0.45;
      config.prior_beta0_mean = 0.0;
      config.prior_beta1_mean = 0.0;
      config.prior_beta0_sd = 1.5;
      config.prior_beta1_sd = 1.5;
      config.seed = dataset_seed * 10 + chain;
      config.n_location_updates = 12;
      config.log_sigma_step = 0.070;
      config.beta0_step = 0.100;
      config.beta1_step = 0.110;
      config.location_step = 0.070;
      config.rebuild_every = 200;

      const bdsecr::SamplerResult fit = bdsecr::fit(config);
      for (std::size_t draw = 0; draw < fit.N.size(); ++draw) {
        traces << dataset << ',' << chain << ',' << draw + 1 << ','
               << fit.N[draw] << ',' << fit.sigma[draw] << ','
               << fit.beta0[draw] << ',' << fit.beta1[draw] << ','
               << fit.log_likelihood[draw] << '\n';
      }
      const auto& proposed = fit.proposed;
      const auto& accepted = fit.accepted;
      acceptance
          << dataset << ',' << chain << ','
          << bdsecr::acceptance_rate(accepted.sigma, proposed.sigma) << ','
          << bdsecr::acceptance_rate(accepted.beta0, proposed.beta0) << ','
          << bdsecr::acceptance_rate(accepted.beta1, proposed.beta1) << ','
          << bdsecr::acceptance_rate(accepted.location, proposed.location) << ','
          << bdsecr::acceptance_rate(accepted.birth, proposed.birth) << ','
          << bdsecr::acceptance_rate(accepted.death, proposed.death) << ','
          << bdsecr::acceptance_rate(accepted.birth + accepted.death,
                                    proposed.birth + proposed.death)
          << ',' << fit.cache_max_abs_error << '\n';
      std::cout << "dataset=" << dataset << " chain=" << chain
                << " retained=" << fit.N.size() << '\n';
    }
  }

  std::cout << "multidataset_recovery=COMPLETE\n";
  std::cout << "output_directory=" << output_directory.string() << '\n';
  return 0;
}
