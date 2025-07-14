/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * mcmc.cpp
 *
 */

#include "mcmc_impl.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : implementations  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//===== CONSTRUCTOR =====//
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

//===== DATA =====//
ulong MCMC::get_report_indiv(const std::string &s) const { return pImpl_->dataset.get_report_indiv(s); }
ulong MCMC::get_report_house(const std::string &s) const { return pImpl_->dataset.get_report_house(s); }
std::pair<ushort, ulong> MCMC::get_report_hh_size(const std::string &s) const { return pImpl_->dataset.get_report_hh_size(s); }

//===== RAND NUMBER GENERATOR =====//
double MCMC::uniform() { return pImpl_->setup.unif(); }

//===== PARAMETERS =====//
double MCMC::param(const par_key &p) const { return pImpl_->param.get(p); }
const double *MCMC::param_pointer(const par_key &p) const { return pImpl_->param.get_ptr(p); }

//===== AUGMENTED VARIABLES =====//
void MCMC::init_all_augm(Individual &member) { pImpl_->augm.init_all(member); }
void MCMC::clear_all_augm(Individual &member) { pImpl_->augm.clear_all(member); }

//===== CACHE ATTRIBUTES =====//
void MCMC::cache_house_integer_vector(const cache_key &e) { pImpl_->dataset.init_cache_integer_vector(e); }
void MCMC::cache_house_double_vector(const cache_key &e) { pImpl_->dataset.init_cache_double_vector(e); }
void MCMC::cache_house_double_matrix(const cache_key &e) { pImpl_->dataset.init_cache_double_matrix(e); }

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Build all group rates and contact rates and set rate per group */
void MCMC::build_all_rates()
{
    try
    {
        for (Household &house : pImpl_->dataset)
            init_contact_rates_(house);
    }
    catch (const std::exception &e)
    {
        tools::error("While initializing contact rates:\n" +
                     std::string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : set rate per contact  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Set contact rate according to user assignement. */
void MCMC::init_contact_rates_(Household &house)
{
    try
    {
        // loop over all couples
        for (const Individual &member_k : house)     // from : member_k
            for (const Individual &member_i : house) // -> to : member_i
                if (member_k != member_i)
                {
                    _pair_contact_rate = HSP_NA_STRING;
                    init_household_contact_rates(house, member_k, member_i);

                    if (!tools::is_na(_pair_contact_rate))
                    {
                        const double *y = &REFERENCE_VALUE; // default reference rate

                        if (_pair_contact_rate == MISSING_RATE)
                            throw std::invalid_argument(
                                "Missing contact rate for individuals " + member_k.id() +
                                " and " + member_i.id() + " (household " + house.house_id() + ").");

                        else if (_pair_contact_rate != REFERENCE)
                            // parameter rate
                            y = pImpl_->param.get_ptr(_pair_contact_rate);

                        house.init_contact_rates(static_cast<size_t>(member_k), static_cast<size_t>(member_i), y);
                    }
                }
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : infectivity profile  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::init_infec_profile_template(dict_str_str x) { pImpl_->infec_profile.init_template(x); }
double MCMC::infec_profile_inst(const Individual &i, double t) const { return pImpl_->infec_profile.inst(&i, t); }
double MCMC::infec_profile_cumul(const Individual &i, double t) const { return pImpl_->infec_profile.cumul(&i, t); }

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : general methods  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double MCMC::evaluate() { return pImpl_->evaluate(); }
void MCMC::infer() { pImpl_->infer(); }

#ifdef _HOUSEPI_RCPP
Rcpp::DataFrame MCMC::simulate() { return pImpl_->simulate(); }
Rcpp::List MCMC::evidence(const std::string &scheme, ulong n_rungs) { return pImpl_->evidence(scheme, n_rungs); }

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : exports to Rcpp ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Rcpp::List MCMC::export_setup() const { return pImpl_->setup.rcpp_export(); }
Rcpp::List MCMC::export_infec_profile() const { return pImpl_->infec_profile.rcpp_export(); }
Rcpp::DataFrame MCMC::export_param() const { return pImpl_->param.rcpp_export(); }
Rcpp::DataFrame MCMC::export_dataset() const { return pImpl_->dataset.rcpp_export(); }
Rcpp::DataFrame MCMC::export_dataset_report() const { return pImpl_->dataset.rcpp_export_report(); }
Rcpp::DataFrame MCMC::export_house_size_report() const { return pImpl_->dataset.rcpp_export_size_report(); }
Rcpp::NumericMatrix MCMC::export_house_loglike_traces() const { return pImpl_->dataset.rcpp_export_traces(); }
Rcpp::DataFrame MCMC::export_datalist() const { return pImpl_->datalist.rcpp_export(); }
Rcpp::NumericMatrix MCMC::export_traces() const { return pImpl_->rcpp_export_traces(); }
Rcpp::List MCMC::export_house_snapshot(const std::string &house_id) const { return pImpl_->rcpp_export_house_snapshot(house_id); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : imports from Rcpp ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::import_setup(SEXP input) { pImpl_->setup.rcpp_import(input); }
void MCMC::import_dataset(SEXP input) { pImpl_->dataset.rcpp_import(input); }
void MCMC::import_multiplicity(SEXP input) { pImpl_->dataset.rcpp_import_multiplicity(input); }
void MCMC::import_datalist(SEXP input) { pImpl_->datalist.rcpp_import(input); }
void MCMC::import_group_rates(SEXP input) { pImpl_->group_rates.rcpp_import(input); }
void MCMC::import_param(SEXP input) { pImpl_->param.rcpp_import(input); }
void MCMC::import_infec_profile(SEXP input) { pImpl_->infec_profile.rcpp_import(input); }
void MCMC::import_augm(SEXP input) { pImpl_->augm.rcpp_import(input); }

#else
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  MCMC : C++ prints  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::print_simulated_dataset(const std::string &f) { pImpl_->print_simulated_dataset(f); }
void MCMC::print_evidence(const std::string &scheme, ulong n_rungs) { pImpl_->print_evidence(scheme, n_rungs); }

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  MCMC : read  //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::read_setup(const std::string &f) { pImpl_->setup.read(f); }
void MCMC::read_dataset(const std::string &f) { pImpl_->dataset.read(f); }
void MCMC::read_multiplicity(const std::string &f) { pImpl_->dataset.read_multiplicity(f); }
void MCMC::read_datalist(const std::string &f) { pImpl_->datalist.read(f); }
void MCMC::read_group_rates(const std::string &f) { pImpl_->group_rates.read(f); }
void MCMC::read_param(const std::string &f) { pImpl_->param.read(f); }
void MCMC::read_augm(const std::string &f) { pImpl_->augm.read(f); }
void MCMC::read_infec_profile(const std::string &f) { pImpl_->infec_profile.read(f); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : show, C++ print  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::show_setup() const { pImpl_->setup.show(); }
void MCMC::show_all_params() const { pImpl_->param.show(); }
void MCMC::show_results() const { pImpl_->param.show_results(); }
void MCMC::show_summary() const { pImpl_->dataset.show_summary(); }
void MCMC::print_infer_state(const std::string &f) const { pImpl_->print_infer_state(f); }
#endif
