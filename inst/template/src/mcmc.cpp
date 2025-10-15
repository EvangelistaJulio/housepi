/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * mcmc.cpp
 *
 */

#include "tools/error.h" // UPDATE: required by contact rates
#include "tools/rates.h" // UPDATE: required by contact rates
#include "mcmc_impl.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : constructor  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

MCMC::MCMC()
    : pImpl_(std::make_unique<Impl>(
          [this]()
          { this->init_model_for_inference(); },
          [this]()
          { this->init_model_for_simulation(); },
          [this](Household &h) -> double
          { return this->compute_household_obs_loglike(h); },
          [this](Household &h)
          { return this->simulate_single_household(h); },
          [this](Household &h, Individual &i)
          { this->update_cache_dependecies(h, i); })) {}
MCMC::~MCMC() = default;

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////  MCMC : data  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

const ReportIndividual &MCMC::get_report_indiv() const { return pImpl_->dataset.report_indiv; }
const ReportHousehold &MCMC::get_report_house() const { return pImpl_->dataset.report_house; }
const ReportHouseSizes &MCMC::get_report_hh_size() const { return pImpl_->dataset.report_hh_size; }
void MCMC::assert_no_missing_status() const
{
    if (pImpl_->dataset.report_indiv.missing > 0)
        hsp::msg::error(
            "Individual with missing status are incompatible with current model.");
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  MCMC : parameters  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double MCMC::param(const std::string &p) const { return pImpl_->param.get(p); }

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : augmented variables  /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::init_missing_augm(Individual &member)
{
    pImpl_->augm.build_missing(member);
    pImpl_->augm.user_reset(member);
}
void MCMC::clear_all_augm(Individual &member) { pImpl_->augm.clear_all(member); }
void MCMC::init_user_list_discrete(const std::string &x) { pImpl_->augm.init_user_list_discrete(x); }

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  MCMC : cache attributes  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::cache_house_integer_vector(const std::string &e) { pImpl_->dataset.init_cache_integer_vector(e); }
void MCMC::cache_house_size_t_vector(const std::string &e) { pImpl_->dataset.init_cache_size_t_vector(e); }
void MCMC::cache_house_double_vector(const std::string &e) { pImpl_->dataset.init_cache_double_vector(e); }
void MCMC::cache_house_double_matrix(const std::string &e) { pImpl_->dataset.init_cache_double_matrix(e); }

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : infectivity profile  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::init_infec_profile_template(u_dict<std::string> x) { pImpl_->infec_profile.init_template(x); }
double MCMC::infec_profile_inst(const Individual &i, double t) const { return pImpl_->infec_profile.inst(&i, t); }
double MCMC::infec_profile_cumul(const Individual &i, double t) const { return pImpl_->infec_profile.cumul(&i, t); }

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : general methods  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double MCMC::uniform() { return pImpl_->setup.unif(); }
double MCMC::evaluate() { return pImpl_->evaluate(); }
void MCMC::infer() { pImpl_->infer(); }
void MCMC::clear_memory() { pImpl_->clear_memory(); }

#ifdef _HOUSEPI_RCPP
Rcpp::DataFrame MCMC::data_wide_for_simulation() const
{
    return pImpl_->dataset.rcpp_export(pImpl_->dataset.headers() - pImpl_->augm.list_export, false);
}
Rcpp::List MCMC::simulate(size_t n_rep) { return pImpl_->simulate(n_rep); }
Rcpp::List MCMC::logevidence(const std::string &scheme, unsigned long n_rungs) { return pImpl_->logevidence(scheme, n_rungs); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : imports from Rcpp ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::import_setup(SEXP input) { pImpl_->setup.rcpp_import(input); }
void MCMC::import_data_wide(SEXP input) { pImpl_->dataset.rcpp_import(input); }
void MCMC::import_multiplicity(SEXP input) { pImpl_->dataset.rcpp_import_multiplicity(input); }
void MCMC::import_data_long(SEXP input) { pImpl_->datalong.rcpp_import(input); }
void MCMC::import_group_rates(SEXP input) { pImpl_->group_rates.rcpp_import(input); }
void MCMC::import_param(SEXP input) { pImpl_->param.rcpp_import(input); }
void MCMC::import_infec_profile(SEXP input) { pImpl_->infec_profile.rcpp_import(input); }
void MCMC::import_augm(SEXP input) { pImpl_->augm.rcpp_import(input); }

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : exports to Rcpp ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Rcpp::DataFrame MCMC::export_data_wide() const { return pImpl_->dataset.rcpp_export(pImpl_->dataset.headers(), true); }
Rcpp::DataFrame MCMC::export_data_long() const { return pImpl_->datalong.rcpp_export(); }
Rcpp::DataFrame MCMC::summarise_cases() const { return pImpl_->dataset.rcpp_summarise_cases(); }
Rcpp::DataFrame MCMC::summarise_house_sizes() const { return pImpl_->dataset.rcpp_summarise_house_sizes(); }
Rcpp::List MCMC::export_setup() const { return pImpl_->setup.rcpp_export(); }
Rcpp::CharacterVector MCMC::export_param_list_const() const { return Rcpp::wrap(pImpl_->param.list_const); }
Rcpp::CharacterVector MCMC::export_param_list_sample_standard() const { return Rcpp::wrap(pImpl_->param.list_sample_standard); }
Rcpp::CharacterVector MCMC::export_param_list_sample_cache() const { return Rcpp::wrap(pImpl_->param.list_sample_cache); }
Rcpp::List MCMC::export_infec_profile() const { return pImpl_->infec_profile.rcpp_export(); }
Rcpp::CharacterVector MCMC::export_augm_list_draws() const { return Rcpp::wrap(pImpl_->augm.list_draws); }
Rcpp::CharacterVector MCMC::export_augm_list_full() const { return Rcpp::wrap(pImpl_->augm.list_full); }
Rcpp::CharacterVector MCMC::export_augm_list_export() const { return Rcpp::wrap(pImpl_->augm.list_export); }
Rcpp::DataFrame MCMC::export_augmented_data_wide() const { return pImpl_->augm.rcpp_export_data_wide(); }
Rcpp::NumericMatrix MCMC::export_draws() const { return pImpl_->rcpp_export_draws(); }
Rcpp::NumericMatrix MCMC::export_house_loglike_draws() const { return pImpl_->dataset.rcpp_export_draws(); }
Rcpp::List MCMC::export_house_snapshot(const std::string &house_id) const { return pImpl_->rcpp_export_house_snapshot(house_id); }

#else
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  MCMC : C++ prints  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::print_simulated_data_wide(const std::string &f) { pImpl_->print_simulated_data_wide(f); }
double MCMC::logevidence(const std::string &scheme, unsigned long n_rungs) { return pImpl_->logevidence(scheme, n_rungs); }

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  MCMC : read  //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::read_setup(const std::string &f) { pImpl_->setup.read(f); }
void MCMC::read_data_wide(const std::string &f) { pImpl_->dataset.read(f); }
void MCMC::read_multiplicity(const std::string &f) { pImpl_->dataset.read_multiplicity(f); }
void MCMC::read_data_long(const std::string &f) { pImpl_->datalong.read(f); }
void MCMC::read_group_rates(const std::string &f) { pImpl_->group_rates.read(f); }
void MCMC::read_param(const std::string &f) { pImpl_->param.read(f); }
void MCMC::read_augm(const std::string &f) { pImpl_->augm.read(f); }
void MCMC::read_infec_profile(const std::string &f) { pImpl_->infec_profile.read(f); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : show, C++ print  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::print_infer_state(const std::string &f) const { pImpl_->print_infer_state(f); }
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : contact rates ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Build all group rates and contact rates and set rate per group. // UPDATE
void MCMC::build_all_rates()
{
    try
    {
        for (Household &house : pImpl_->dataset)
            init_contact_rates_(house);
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While initializing contact rates:\n" +
                   std::string(e.what()));
    }
}

/// Set contact rate according to user assignement. // UPDATE
void MCMC::init_contact_rates_(Household &house)
{
    // loop over all couples
    for (const Individual &member_k : house)     // from : member_k
        for (const Individual &member_i : house) // -> to : member_i
            if (member_k != member_i)
            {
                _pair_contact_rate = hsp::na::STRING;
                init_household_contact_rates(house, member_k, member_i);

                if (!hsp::na::is(_pair_contact_rate))
                {
                    observ_dbl_ptr y = &hsp::rates::VALUE_REFERENCE; // default reference rate

                    if (_pair_contact_rate == hsp::rates::MISSING)
                        throw std::invalid_argument(
                            "Missing contact rate for individuals \"" +
                            member_k.id() + "\" and \"" + member_i.id() +
                            "\" (" + house.err_name() + ").");

                    else if (_pair_contact_rate != hsp::rates::REFERENCE)
                        y = pImpl_->param.get_ptr(_pair_contact_rate); // parameter rate

                    house.init_contact_rates(static_cast<size_t>(member_k), static_cast<size_t>(member_i), y);
                }
            }
}
