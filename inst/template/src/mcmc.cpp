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

#include <atomic>

#include "mcmc_impl.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : implementations  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//===== CONSTRUCTOR =====//
MCMC::MCMC()
    : pImpl(std::make_unique<Impl>(
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
ulong MCMC::get_report_indiv(const std::string &key) const { return pImpl->get_report_indiv(key); }
ulong MCMC::get_report_house(const std::string &key) const { return pImpl->get_report_house(key); }
std::pair<ushort, ulong> MCMC::get_report_hh_size(const std::string &key) const { return pImpl->get_report_hh_size(key); }

//===== RAND NUMBER GENERATOR =====//
double MCMC::uniform() { return pImpl->setup.unif(); }

//===== PARAMETERS =====//
double MCMC::param(const par_key &p) const { return pImpl->get_param(p); }
const double *MCMC::param_pointer(const par_key &p) const { return pImpl->get_param_pointer(p); }

//===== AUGMENTED VARIABLES =====//
void MCMC::assign_all_augm(Individual &member) { pImpl->assign_all_augm(member); }
void MCMC::clear_all_augm(Individual &member) { pImpl->clear_all_augm(member); }

//===== EXTRA ATTRIBUTES =====//

/* Extra integer vectors for all households. */
void MCMC::assign_cache_house_integer_vector(const cache_key &e)
{
    for (size_t H : pImpl->restrict_hh_indexes)
        pImpl->house_vector[H].assign_cache_integer_vector(e);
}
/* Extra double vectors for all households. */
void MCMC::assign_cache_house_double_vector(const cache_key &e)
{
    for (size_t H : pImpl->restrict_hh_indexes)
        pImpl->house_vector[H].assign_cache_double_vector(e);
}
/* Extra double square matrix for all households. */
void MCMC::assign_cache_house_double_matrix(const cache_key &e)
{
    for (size_t H : pImpl->restrict_hh_indexes)
        pImpl->house_vector[H].assign_cache_double_matrix(e);
}

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
        for (Household &house : pImpl->house_vector)
            _assign_contact_rate(house);
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
void MCMC::_assign_contact_rate(Household &house)
{
    try
    {
        // loop over all couples
        for (const Individual &member_k : house)     // from : member_k
            for (const Individual &member_i : house) // -> to : member_i
                if (member_k != member_i)
                {
                    _pair_contact_rate = HSP_NA_STRING;
                    assign_contact_rate_to_household(house, member_k, member_i);

                    if (!tools::is_na(_pair_contact_rate))
                    {
                        const double *y = &REFERENCE_VALUE; // default reference rate

                        if (_pair_contact_rate == MISSING_RATE)
                            throw std::invalid_argument(
                                "Missing contact rate for individuals " + member_k.id() +
                                " and " + member_i.id() + " (household " + house.house_id() + ").");

                        else if (_pair_contact_rate != REFERENCE)
                            // parameter rate
                            y = pImpl->get_param_pointer(_pair_contact_rate);

                        house.assign_contact_rate(static_cast<size_t>(member_k), static_cast<size_t>(member_i), y);
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

/* Assign template infectivity profile to use if no user defined. */
void MCMC::assign_template_infec_profile(const dict_str_str x)
{
    pImpl->infec_profile.deflist_template = x;
}

/* Value of the instantaneous infectivity profile. */
double MCMC::infec_profile_inst(const Individual &infector, double t) const
{
    return pImpl->infec_profile.inst(infector, t);
}

/* Value of the cumulative infectivity profile. */
double MCMC::infec_profile_cumul(const Individual &infector, double t) const
{
    return pImpl->infec_profile.cumul(infector, t);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : general methods  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Default method to only compute the loglikelihood and display it. */
double MCMC::evaluate() { return pImpl->evaluate(); }

/* Performs the inference using Metropolis-Hastings. */
void MCMC::infer() { pImpl->infer(); }

#ifdef _HOUSEPI_RCPP
/* Perform simulation and return simulated dataset. */
Rcpp::DataFrame MCMC::simulate() { return pImpl->simulate(); }

Rcpp::List MCMC::estimate_evidence(const std::string &scheme, ulong n_rungs)
{
    return pImpl->estimate_evidence(scheme, n_rungs);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : exports to Rcpp ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Rcpp::List MCMC::export_setup() const { return pImpl->export_setup(); }
Rcpp::DataFrame MCMC::export_param() const { return pImpl->export_param(); }
Rcpp::DataFrame MCMC::export_dataset() const { return pImpl->export_dataset(); }
Rcpp::DataFrame MCMC::export_datalist() const { return pImpl->export_datalist(); }
Rcpp::DataFrame MCMC::export_dataset_report() const { return pImpl->export_dataset_report(); }
Rcpp::DataFrame MCMC::export_house_size_report() const { return pImpl->export_house_size_report(); }
Rcpp::NumericMatrix MCMC::export_traces() const { return pImpl->export_traces(); }
Rcpp::NumericMatrix MCMC::export_house_loglike_mtx() const { return pImpl->export_house_loglike_mtx(); }
Rcpp::List MCMC::showcase_infec_profile(int n_bins) const { return pImpl->showcase_infec_profile(n_bins); }
Rcpp::List MCMC::export_for_debug(const std::string &house_id) const { return pImpl->export_for_debug(house_id); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : imports from Rcpp ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::import_setup(const Rcpp::List &x) { pImpl->import_setup(x); }
void MCMC::import_dataset(const Rcpp::DataFrame &x) { pImpl->import_dataset(x); }
void MCMC::import_datalist(const Rcpp::DataFrame &x) { pImpl->import_datalist(x); }
void MCMC::import_multiplicity(const Rcpp::DataFrame &x) { pImpl->import_multiplicity(x); }
void MCMC::import_group_rates(const Rcpp::DataFrame &x) { pImpl->import_group_rates(x); }
void MCMC::import_param(const Rcpp::DataFrame &x) { pImpl->import_param(x); }
void MCMC::import_infec_profile(const Rcpp::List &x) { pImpl->import_infec_profile(x); }
void MCMC::import_augm(const Rcpp::DataFrame &x) { pImpl->import_augm(x); }

#else
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  MCMC : C++ prints  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Print the simulated dataset to the provided file. */
void MCMC::print_simulated_dataset(const std::string &f) { pImpl->print_simulated_dataset(f); }

/* Trapezoidal rule to estimate logarithm of marginal likelihood (evidence) with power posteriors.
    -----
    Ref: Friel, N., Hurn, M. & Wyse, J.
    Improving power posterior estimation of statistical evidence.
    Stat Comput 24, 709–723 (2014).
*/
void MCMC::print_evidence(const std::string &scheme, ulong n_rungs) { pImpl->print_evidence(scheme, n_rungs); }

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  MCMC : read  //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::read_setup(const std::string &f) { pImpl->read_setup(f); }
void MCMC::read_dataset(const std::string &f) { pImpl->read_dataset(f); }
void MCMC::read_datalist(const std::string &f) { pImpl->read_datalist(f); }
void MCMC::read_multiplicity(const std::string &f) { pImpl->read_multiplicity(f); }
void MCMC::read_group_rates(const std::string &f) { pImpl->read_group_rates(f); }
void MCMC::read_param(const std::string &f) { pImpl->read_param(f); }
void MCMC::read_augm(const std::string &f) { pImpl->read_augm(f); }
void MCMC::read_infec_profile(const std::string &f) { pImpl->read_infec_profile(f); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : show, C++ print  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::show_setup() const { pImpl->show_setup(); }
void MCMC::show_all_params() const { pImpl->show_all_params(); }
void MCMC::show_results() const { pImpl->show_results(); }
void MCMC::show_summary() const { pImpl->show_summary(); }
void MCMC::print_infer_state(const std::string &f) const { pImpl->print_infer_state(f); }

#endif
