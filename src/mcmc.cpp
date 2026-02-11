/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "tools/types.h"

#include "mcmc.h" /// UPDATE: should remove I/O and Rcpp dependencies from MCMC...
#include "mcmc_impl.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** MCMC : class ***  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

MCMC::MCMC(
    const std::vector<std::string> &data_wide_keys,
    const std::vector<std::string> &hh_size_keys,
    const std::vector<std::string> &hh_size_vec_keys,
    const std::vector<std::string> &hh_double_vec_keys,
    const std::vector<std::string> &hh_double_mtx_keys,
    const std::vector<std::string> &data_long_keys,
    const std::vector<std::string> &param_keys,
    const std::vector<std::string> &group_rate_keys,
    const std::vector<std::string> &augm_keys)
    : pImpl_(std::make_unique<MCMCImpl>(
          data_wide_keys,
          hh_size_keys,
          hh_size_vec_keys,
          hh_double_vec_keys,
          hh_double_mtx_keys,
          data_long_keys,
          param_keys,
          group_rate_keys,
          augm_keys,
          [this]()
          { this->init_model_for_evaluation(); },
          [this](Household &h)
          { this->init_house_for_evaluation(h); },
          [this](Household &h) -> double
          { return this->compute_house_obs_log_lik(h); },
          [this](Household &h, Individual &i)
          { this->update_cache_dependecies(h, i); },
          [this](Household &h, const Individual &k, const Individual &i) -> std::string
          { return this->init_contact_rates(h, k, i); })) {}
MCMC::~MCMC() = default;

double MCMC::uniform() { return pImpl_->uniform(); }
double MCMC::par(size_t p) const { return pImpl_->par(p); }

size_t MCMC::get_max_hh_cases() const noexcept { return pImpl_->get_max_hh_cases(); }
void MCMC::assert_no_missing_status() const { pImpl_->assert_no_missing_status(); }

void MCMC::build_contact_rates() { pImpl_->build_contact_rates(); }

void MCMC::init_infec_profile_template(hsp::dict<std::string> x) { pImpl_->init_infec_profile_template(x); }
double MCMC::infec_profile_inst(const AugmContainer &i, double t) const { return pImpl_->infec_profile_inst(i, t); }
double MCMC::infec_profile_cumul(const AugmContainer &i, double t) const { return pImpl_->infec_profile_cumul(i, t); }

void MCMC::init_missing_augm(Individual &member) { pImpl_->init_missing_augm(member); }

void MCMC::request_discrete_param(size_t p) { pImpl_->request_discrete_param(p); }
void MCMC::request_discrete_augm(size_t a) { pImpl_->request_discrete_augm(a); }
void MCMC::request_discrete_infec_profile() { pImpl_->request_discrete_infec_profile(); }

double MCMC::evaluate_log_lik(std::string_view house_id) { return pImpl_->evaluate_log_lik(house_id); }
void MCMC::run_inference() { pImpl_->run_inference(); }
