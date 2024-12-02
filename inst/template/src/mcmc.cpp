/*
 * mcmc.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "mcmc.h"
#include "build.h"

using std::string;
using std::vector;

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : implementations  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

MCMC::MCMC() : pBuild(std::make_unique<Monte_Carlo_Build>()) {}
MCMC::~MCMC() = default;

void MCMC::read_setup(const string &pathtofile) { pBuild->read_setup(pathtofile); }
void MCMC::read_params(const string &pathtofile) { pBuild->read_params(pathtofile); }
void MCMC::read_augm(const string &pathtofile) { pBuild->read_augm(pathtofile); }

#ifdef _HOUSEPI_RCPP
//===== Rcpp IMPORTS =====//
void MCMC::import_setup(const Rcpp::List &x) { pBuild->import_setup(x); }
void MCMC::import_params(const Rcpp::DataFrame &x) { pBuild->import_params(x); }
void MCMC::import_augm(const Rcpp::DataFrame &x) { pBuild->import_augm(x); }

//===== Rcpp EXPORTS =====//
Rcpp::List MCMC::export_setup() const { return pBuild->export_setup(); }
Rcpp::DataFrame MCMC::export_params() const { return pBuild->export_params(); }
Rcpp::DataFrame MCMC::export_data() const { return pBuild->export_data(); }
Rcpp::DataFrame MCMC::export_traces() const { return pBuild->export_traces(); }
Rcpp::NumericMatrix MCMC::export_loglike_matrix() const { return pBuild->export_loglike_matrix(); }
#else

//===== Cpp PRINT =====//
void MCMC::print_infer_state(const string &pathtofile) const { pBuild->print_infer_state(pathtofile); }
void MCMC::print_augm_variables(const string &pathtofile) const { pBuild->print_augm_variables(pathtofile); }
#endif

//===== SHOW =====//
void MCMC::show_full_build() const { pBuild->show_full_build(); }
void MCMC::show_results() const { pBuild->show_results(); }
void MCMC::show_summary() const { pBuild->show_summary(); }
void MCMC::show_group_classes() const { pBuild->show_group_classes(); }

//===== DATA =====//
size_t MCMC::n_houses() const { return pBuild->n_households; }
size_t MCMC::report_data(const string &key) const { return pBuild->report_data(key); }
double MCMC::uniform() { return pBuild->uniform(); }

//===== AUGMENTED =====//
double MCMC::param(const string &parID) const { return pBuild->get_param(parID); }
const double *MCMC::param_pointer(const string &parID) const { return pBuild->get_param_pointer(parID); }
void MCMC::assign_all_augm(Household &house, indiv i) { pBuild->assign_all_augm(house, i); }
void MCMC::clear_all_augm(Household &house, indiv i) { pBuild->clear_all_augm(house, i); }

//===== EXTRA ATTRIBUTES =====//

void MCMC::assign_extra_integer_to_all_houses(const vector<string> &names, int init_value)
{
    for (size_t H = 0; H < n_houses(); ++H)
    {
        pBuild->house_vector[H].allocate_extra_int(names, init_value);
    }
}
void MCMC::assign_extra_double_to_all_houses(const vector<string> &names, double init_value)
{
    for (size_t H = 0; H < n_houses(); ++H)
    {
        pBuild->house_vector[H].allocate_extra_double(names, init_value);
    }
}
void MCMC::assign_extra_vector_to_all_houses(const vector<string> &names, double init_value)
{
    for (size_t H = 0; H < n_houses(); ++H)
    {
        pBuild->house_vector[H].allocate_extra_vector(names, init_value);
    }
}
void MCMC::assign_extra_matrix_to_all_houses(const vector<string> &names, double init_value)
{
    for (size_t H = 0; H < n_houses(); ++H)
    {
        pBuild->house_vector[H].allocate_extra_matrix(names, init_value);
    }
}

//===== INFECTIVITY PROFILE =====//
void MCMC::assign_infectivity_profile(const string &name)
{
    _infec_profile_ = infectivity_profile_switchboard(name);
}
double MCMC::infec_profile_inst(double serial, const vector<double> &x)
{
    return _infec_profile_->pdf(serial, x);
}

double MCMC::infec_profile_cumul(double serial, const vector<double> &x)
{
    return _infec_profile_->cdf(serial, x);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : methods (user defined)  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Update cumulative and instantaneous person-to-person infectivity profiles according to an infector. */
void MCMC::update_augm_dependecies(Household &house, indiv i)
{
    UNUSED(house);
    UNUSED(i);
}

/* Compute the loglikelihood for a single household. */
double MCMC::compute_house_loglike(Household &house) const
{
    UNUSED(house);
    return 0.0;
};

/* Assign group rate to an individual  */
void MCMC::assign_group_rate_to_individual(Household &house, indiv i)
{
    // this_indiv.assign_group_rate(EMPTY_RATE, EMPTY_RATE);
    UNUSED(house);
    UNUSED(i);
}

/* Assign the contact patterns between a couple of individuals. */
void MCMC::assign_contact_rate_to_household(Household &house, indiv from, indiv to)
{
    house.assign_contact_rate(from, to, EMPTY_RATE);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Monte_Carlo_Build : build data  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read data file and then run assignment. */
void MCMC::read_data(const string &pathtofile)
{
    std::map<string, vector<string>> datatable = {};
    tools::read_csv_table(pathtofile, datatable);
    build_data(datatable);

    pBuild->_report_pathtofile_["data"] = pathtofile;
}

#ifdef _HOUSEPI_RCPP
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  Monte_Carlo_Build : imports from Rcpp ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Import data from R DataFrame and run build. */
void MCMC::import_data(const Rcpp::DataFrame &data)
{
    std::map<string, vector<string>> datatable = {};
    tools::rdataframe_to_map(data, datatable);
    build_data(datatable);
    if (pBuild->_report_pathtofile_.find("data") != pBuild->_report_pathtofile_.end())
    {
        pBuild->_report_pathtofile_.erase("data");
    }
}
#endif

/* Assign data after pre-processing. */
void MCMC::build_data(const std::map<string, vector<string>> &datatable)
{
    try
    {
        if (datatable.empty())
        {
            // If nothing was previously defined
            if (pBuild->house_vector.empty())
            {
                throw std::length_error("No database was found.");
            }
        }
        else
        {
            //=====  Drop all Previous Data and Prepare for New  =====//
            pBuild->data_headers.clear();
            for (const auto &x : datatable)
            {
                pBuild->data_headers.push_back(x.first);
            }
            pBuild->n_individuals = datatable.at(pBuild->data_headers[0]).size();
            pBuild->n_households = 0;
            pBuild->house_vector.clear();
            pBuild->_report_group_classes_.clear();

            //===== Assess the presence of house IDs =====//
            if (std::find(pBuild->data_headers.begin(), pBuild->data_headers.end(), "house_id") == pBuild->data_headers.end())
            {
                tools::warning(
                    "Unknown household structure: missing data column 'house_id'.");
            }

            //===== Safe Load of Households =====//
            std::unordered_map<string, size_t> hh_ID_to_index = {};
            for (size_t i = 0; i < pBuild->n_individuals; ++i)
            {
                //----- Read individual data -----//
                string this_indiv_ID = std::to_string(i);     // default value
                string house_ID = std::to_string(n_houses()); // default value
                int this_status = MISSING_STATUS;
                dictargs args = {};

                for (const string &key : pBuild->data_headers)
                {
                    if (key == "id")
                    {
                        if (tools::is_na(datatable.at("id")[i]))
                        {
                            throw std::invalid_argument(
                                "Missing id for individual in row n." + std::to_string(i + 2) + ".");
                        }
                        else
                        {
                            this_indiv_ID = datatable.at("id")[i];
                        }
                    }
                    else if (key == "house_id")
                    {
                        if (tools::is_na(datatable.at("house_id")[i]))
                        {
                            throw std::invalid_argument(
                                "Missing house_id for individual in row n." + std::to_string(i + 2) + ".");
                        }
                        else
                        {
                            house_ID = datatable.at("house_id")[i];
                        }
                    }
                    else if (key == "status")
                    {
                        // NA status is a missing status
                        if (!tools::is_na(datatable.at("status")[i]))
                        {
                            this_status = tools::stoul(datatable.at("status")[i], "individual status");
                        }
                    }
                    else if (tools::is_na(datatable.at(key)[i]))
                    {
                        args.emplace(key, std::numeric_limits<double>::quiet_NaN());
                        tools::warning("Found NA value(s) in the database.");
                    }
                    else
                    {
                        args.emplace(key, tools::stod(datatable.at(key)[i], "data '" + key + "'"));
                    }
                }

                // if household index doesn't exist, add new household
                // note : this is redundant if no structure
                if (hh_ID_to_index.find(house_ID) == hh_ID_to_index.end())
                {
                    // map data index to household number
                    hh_ID_to_index.emplace(house_ID, n_houses());
                    // create a new empty household
                    pBuild->assign_house(house_ID);
                }

                Household &house = pBuild->house_vector[hh_ID_to_index.at(house_ID)];
                // Assign the individual the respective household
                house.assign_individual(this_indiv_ID, this_status, args);
                // Assign and set rate per group
                set_rate_per_group(house, house.size() - 1);
            }

            //===== Set Contact Rates =====//
            pBuild->_n_contacts_ = 0;
            for (Household &house : pBuild->house_vector)
            {
                set_rate_per_contact(house);
            }

            //===== Update Report for Data =====//
            pBuild->update_report_data();

            //===== Clean up data entries for output files =====//
            // eventually remove standard headers to reorder them at the end
            for (const string &col : {"house_id", "id", "status"})
            {
                pBuild->data_headers.erase(
                    std::remove(pBuild->data_headers.begin(), pBuild->data_headers.end(), col),
                    pBuild->data_headers.end());
            }
        }
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Monte_Carlo_Build : set rate per group  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Update model report and set group rates associate to new individual. */
void MCMC::set_rate_per_group(Household &house, indiv i)
{
    try
    {
        assign_group_rate_to_individual(house, i);

        for (const auto &pair : house.get_group_rate_map(i))
        {
            string rateID = pair.first;
            if (rateID == CONTACT_RATE)
            {
                throw std::invalid_argument(
                    "The group rate name '" + CONTACT_RATE + "' is reserved. Choose a different one.");
            }

            string parID = pair.second;

            //=====initialize model report=====//
            // if first individual, initialize rateID-report
            pBuild->_report_group_classes_.emplace(rateID, std::map<string, size_t>());
            pBuild->_report_group_classes_.at(rateID).emplace("[reference]", 0);

            //=====set rates and update report=====//
            if (tools::is_na(parID))
            {
                // error : undefined rate
                std::invalid_argument("Rate '" + rateID + "' is not provided for individual " + house.ID(i) + ".");
            }
            else if (parID == REFERENCE)
            {
                // rate set to reference
                house.set_rate(i, rateID, &REFERENCE_VALUE);
                pBuild->_report_group_classes_.at(rateID).at("[reference]") += 1;
            }
            else if (parID != REFERENCE)
            {
                // rate associated to parameter
                pBuild->assert_param(parID);
                house.set_rate(i, rateID, pBuild->get_param_pointer(parID));
                if (pBuild->_report_group_classes_.at(rateID).find(parID) == pBuild->_report_group_classes_.at(rateID).end())
                {
                    pBuild->_report_group_classes_.at(rateID).emplace(parID, 1);
                }
                else
                {
                    pBuild->_report_group_classes_.at(rateID).at(parID) += 1;
                }
            }
        }
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Monte_Carlo_Build : set rate per contact  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Set contact rate according to user assignement. */
void MCMC::set_rate_per_contact(Household &house)
{
    try
    {
        // allocate contact patterns
        house.init_contact_matrix();

        if (house.index() == 0)
        {
            // if first household, initialize contact rate report
            pBuild->_report_group_classes_.emplace(CONTACT_RATE, std::map<string, size_t>());
            pBuild->_report_group_classes_.at(CONTACT_RATE).emplace("[reference]", 0);
        }
        else
        {
            // all household must have a contact rate
            if (pBuild->_report_group_classes_.find(CONTACT_RATE) == pBuild->_report_group_classes_.end())
            {
                throw std::invalid_argument(
                    "Contact rate has not been provided for all households.");
            }
        }

        // loop over all couples
        for (indiv k = 0; k < house.size(); ++k) // from
        {
            for (indiv i = 0; i < house.size(); ++i) // to
            {
                // note : from 'k' --> to 'i'
                if (k != i)
                {
                    assign_contact_rate_to_household(house, k, i);
                    string parID = house.get_contact_rate_map(k, i);

                    if (tools::is_na(parID))
                    {
                        string sms = "Contact rate not provided for individuals " + house.ID(k);
                        sms += " and " + house.ID(i) + " within household " + house.house_ID() + ".";
                        throw std::invalid_argument(sms);
                    }
                    else if (parID == REFERENCE)
                    {
                        house.set_contact_rate(k, i, &REFERENCE_VALUE); // reference rate
                        pBuild->_report_group_classes_.at(CONTACT_RATE).at("[reference]") += 1;
                    }
                    else if (parID != EMPTY_RATE)
                    {
                        pBuild->assert_param(parID);
                        house.set_contact_rate(k, i, pBuild->get_param_pointer(parID)); // parameter rate
                        if (pBuild->_report_group_classes_.at(CONTACT_RATE).find(parID) == pBuild->_report_group_classes_.at(CONTACT_RATE).end())
                        {
                            pBuild->_report_group_classes_.at(CONTACT_RATE).emplace(parID, 1);
                        }
                        else
                        {
                            pBuild->_report_group_classes_.at(CONTACT_RATE).at(parID) += 1;
                        }
                    }

                    // count the total number of pairs
                    ++pBuild->_n_contacts_;
                }
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

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC : methods (general)  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute new loglikelihood, store old value in buffer. Same for households. */
void MCMC::update_global_loglike()
{
    buffer = global_loglike;
    global_loglike = 0.0;

#pragma omp parallel for reduction(+ : global_loglike)
    for (size_t H = 0; H < n_houses(); ++H)
    {
        if (_error_occurred_)
        {
            continue;
        }
        try
        {
            Household &house = pBuild->house_vector[H];
            // update house loglikelihood
            house.set_loglike(compute_house_loglike(house));
            global_loglike += house.get_loglike() + pBuild->compute_augm_model_loglike(house);
        }
        catch (const std::exception &e)
        {
#pragma omp critical
            {
                if (!_error_occurred_)
                {
                    _error_occurred_ = true;
                    _error_message_ = e.what();
                }
            }
        }
    }
    if (_error_occurred_)
    {
        tools::error(_error_message_);
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : parameter move  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Attempt to move a single parameter for the given inverse temperature. */
void MCMC::move_single_param(const string &parID)
{
    bool is_accepted = false;

    // proceed only if the new value is consistent, otherwise automatically execute refuse_move(parID)
    if (pBuild->propose_valid_param_move(parID))
    {
        if ((this->*_func_param_metropolis_)(parID))
        {
            is_accepted = true;
        }
        else
        {
            // restore buffered value of parameter
            pBuild->restore_param_value(parID);

            // restore buffered value of global loglike
            restore_global_loglike();

            // restore buffered loglikelihood in households
            for (size_t H = 0; H < n_houses(); ++H)
            {
                pBuild->house_vector[H].restore_loglike();
            }
        }
    }

    // update counts for acceptance rate
    pBuild->update_param_accept_rate(parID, is_accepted);
}

bool MCMC::_param_metropolis_1_(const string &parID)
{
    update_global_loglike();

    if (uniform() < exp(get_delta_global_loglike() + pBuild->param_move_logratio(parID)))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool MCMC::_param_metropolis_0_(const string &parID)
{
    if (uniform() < exp(pBuild->param_move_logratio(parID)))
    {
        update_global_loglike();
        return true;
    }
    else
    {
        return false;
    }
}

bool MCMC::_param_metropolis_T_(const string &parID)
{
    update_global_loglike();

    if (uniform() < exp(pBuild->get_inv_T() * get_delta_global_loglike() + pBuild->param_move_logratio(parID)))
    {
        return true;
    }
    else
    {
        return false;
    }
}

/* Move all model parameters. */
void MCMC::move_all_parameters()
{
    for (const string &parID : pBuild->get_param_to_be_inferred())
    {
        move_single_param(parID);
    }
}

/* Move all model parameters and adapt their rule arguments. */
void MCMC::move_and_adapt_all_parameters(ulong step)
{
    double adapt_ratio = static_cast<double>(step) / pBuild->get_n_adapt();

    for (const string &parID : pBuild->get_param_to_be_inferred())
    {
        move_single_param(parID);
        pBuild->adapt_param_rule(parID, adapt_ratio);
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : augmented move  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize augmented dependecies for all households. */
void MCMC::init_augmented_variables()
{
    if (!pBuild->augm_list_full.empty())
    {
        if (pBuild->is_augm_infer)
        {
            _func_move_all_augm_ = &MCMC::_move_all_augm_;

            if (!pBuild->augm_list_infer.empty())
            {
                _func_move_augm_indiv_ = &MCMC::_move_augm_indiv_generic_;
            }
        }

        for (size_t H = 0; H < n_houses(); ++H)
        {
            Household &house = pBuild->house_vector[H];
            for (const indiv &i : house.infectious_set())
            {
                update_augm_dependecies(house, i);
            }
        }
    }
}

/* Move augmented variables from all households. */
void MCMC::_move_all_augm_()
{
    double delta_ll = 0.;

#pragma omp parallel for reduction(+ : delta_ll)
    for (size_t H = 0; H < n_houses(); ++H)
    {
        if (_error_occurred_)
        {
            continue;
        }
        try
        {
            delta_ll += move_augm_single_house(pBuild->house_vector[H]);
        }
        catch (const std::exception &e)
        {
#pragma omp critical
            {
                if (!_error_occurred_)
                {
                    _error_occurred_ = true;
                    _error_message_ = e.what();
                }
            }
        }
    }
    if (_error_occurred_)
    {
        tools::error(_error_message_);
    }
    else
    {
        increment_global_loglike(delta_ll);
    }
}

/* Move augmented variables for infectious members of a single household.
The move is performed one member at a time, one augmented variable at a time.
*/
double MCMC::move_augm_single_house(Household &house)
{
    double delta_global_loglike = 0.0;

    // Loop over infectious residents of the househould
    for (const indiv &i : house.infectious_set())
    {
        (this->*_func_move_augm_indiv_)(house, i, delta_global_loglike);
    }

    return delta_global_loglike;
}

/* Generic move augmented variables one by one. */
void MCMC::_move_augm_indiv_generic_(Household &house, indiv i, double &delta)
{
    // Loop over augmented variables of infectious
    for (const string &augmID : pBuild->augm_list_infer)
    {
        metropolis_augm_single_move(house, i, augmID, delta);
    }
}

/* Move augmented variables by status one variable at a time. */
void MCMC::_move_augm_indiv_by_status_(Household &house, indiv i, double &delta)
{
    // Loop over augmented variables of infectious
    for (const string &augmID : pBuild->augm_status_map_infer.at(house.status(i)))
    {
        metropolis_augm_single_move(house, i, augmID, delta);
    }
}

/* Take care of a a single augmented Metropolis move. */
void MCMC::metropolis_augm_single_move(Household &house, indiv i, const string &augmID, double &delta)
{
    // proceed only if the new value is consistent
    if (house.propose_valid_augm_move(i, augmID, pBuild->get_prng()))
    {
        // update dependencies to perform computations
        house._buffer_all_extra();
        update_augm_dependecies(house, i);

        // update house loglikelihood and store old value in buffer
        house.set_loglike(compute_house_loglike(house));

        // Metropolis-Hasting algorithm
        if (!(this->*_func_augm_metropolis_)(house, i, augmID, delta))
        {
            // restore value and its dependencies
            house.restore_augm_value(i, augmID);
            house._restore_all_extra();
            house.restore_loglike();
        }
    }
}

/* Metropolis algorithm for single augmented variable move (default, inverse temperature 1). */
bool MCMC::_augm_metropolis_1_(const Household &house, indiv i, const string &augmID, double &delta)
{
    double D_hh_ll = house.get_delta_loglike();
    if (uniform() < exp(D_hh_ll + house.augm_move_logratio(i, augmID)))
    {
        delta += D_hh_ll + house.augm_model_logratio(i, augmID);
        return true;
    }
    else
    {
        return false;
    }
}

/* Metropolis algorithm for single augmented variable move (inverse temperature 0). */
bool MCMC::_augm_metropolis_0_(const Household &house, indiv i, const string &augmID, double &delta)
{
    if (uniform() < exp(-house.augm_propose_logratio(i, augmID)))
    {
        delta += house.get_delta_loglike() + house.augm_model_logratio(i, augmID);
        return true;
    }
    else
    {
        return false;
    }
}

/* Metropolis algorithm for single augmented variable move (all inverse temperatures). */
bool MCMC::_augm_metropolis_T_(const Household &house, indiv i, const string &augmID, double &delta)
{
    double D_ll = house.get_delta_loglike() + house.augm_model_logratio(i, augmID);
    if (uniform() < exp(pBuild->get_inv_T() * D_ll - house.augm_propose_logratio(i, augmID)))
    {
        delta += D_ll;
        return true;
    }
    else
    {
        return false;
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : display likelihood  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Default method to only compute the loglikelihood and display it. */
double MCMC::evaluate()
{
    try
    {
        //----- Check no undefined parameters -----//
        pBuild->set_rnd_init_params(false);
        pBuild->init_rand_params();

        //----- Initialize model -----//
        init_model_for_inference();

        //----- Initialize augmented dependencies -----//
        init_augmented_variables();

        //----- Initialize global loglikelihood -----//
        update_global_loglike();
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }

    return get_global_loglike();
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : run markov chain  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* The text to print in the Markov chian progress bar. */
const string MCMC::fmt_progtxt(const string &txt, bool is_power_post)
{
    if (is_power_post)
    {
        return "power post " + std::to_string(_rung_) + "/" + std::to_string(_n_rungs_);
    }
    else if (pBuild->get_n_chains() > 1)
    {
        return txt + " " + std::to_string(_chain_) + "/" + std::to_string(pBuild->get_n_chains());
    }
    else
    {
        return txt;
    }
}

/* Markov chain core, used for inference and power posteriors.
    note: sweep is the number of previously completed thin sweeps. */
void MCMC::run_monte_carlo(CustomProgressBar &progressbar, ulong sweep, bool is_power_post)
{
    //==========================//
    //===== Initialization =====//
    //==========================//

    //----- Initialize model -----//
    init_model_for_inference();

    //----- Initialize diagnostics -----//
    pBuild->init_diagnostics(is_power_post);

    //===== Run Multiple MCMC Chains =====//
    for (_chain_ = 1; _chain_ <= pBuild->get_n_chains(); ++_chain_)
    {
        // update progressbar text
        progressbar.set_text(fmt_progtxt("infer", is_power_post));

        // initialize parameter inference
        pBuild->init_rand_params();

        // initialize augmented dependencies
        init_augmented_variables();

        // initialize global loglikelihood
        update_global_loglike();

        // avoid asserting user definitions after first round
        pBuild->disable_assert();

        //=========================================//
        //======== (Finite) Adaptative MCMC =======//
        //=========================================//

        // initialize the param move counts for adaptive
        pBuild->init_move_counts();

        for (ulong i = 0; i < pBuild->get_n_adapt(); i += pBuild->get_thin())
        {
            for (ulong t = 1; t <= pBuild->get_thin(); ++t)
            {
                (this->*_func_move_all_augm_)();
                move_and_adapt_all_parameters(i + t);
            }
            ++sweep;

            progressbar.update(static_cast<double>(sweep) / _tot_thin_sweeps_);
        }

        //============================================//
        //========= MCMC Burn-in Iterations ==========//
        //============================================//

        for (ulong i = 0; i < pBuild->get_n_burnin(); i += pBuild->get_thin())
        {
            for (ulong t = 1; t <= pBuild->get_thin(); ++t)
            {
                (this->*_func_move_all_augm_)();
                move_all_parameters();
            }
            ++sweep;

            progressbar.update(static_cast<double>(sweep) / _tot_thin_sweeps_);
        }

        //============================================//
        //========= MCMC Posterior Sampling ==========//
        //============================================//

        for (ulong i = pBuild->get_n_burnin(); i < pBuild->get_n_steps(); i += pBuild->get_thin())
        {
            // thinning mcmc moves
            for (ulong t = 1; t <= pBuild->get_thin(); ++t)
            {
                (this->*_func_move_all_augm_)();
                move_all_parameters();
            }
            ++sweep;

            // compute diagnostics for the previous iterations
            pBuild->update_diagnostics(get_global_loglike(), is_power_post);
            progressbar.update(static_cast<double>(sweep) / _tot_thin_sweeps_);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  MCMC : infer  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Performs the inference using Metropolis-Hastings. */
void MCMC::infer(bool is_progbar)
{
    // UPDATE : Check if data is meaningful
    /*
    if (report_data("n.miss.stat") == static_cast<int>(pBuild->n_individuals))
    {
        tools::error("No infectious status has been provided in the data file.");
    }
    */

    //----- total number of mcmc sweeps -----//
    _chain_ = 1;
    _tot_thin_sweeps_ = pBuild->get_c_n_sweeps() * pBuild->get_n_chains();

    //----- Start progress bar and chrono -----//
    CustomProgressBar progressbar(is_progbar);
    progressbar.initialize(fmt_progtxt("infer"));

    //----- Run Markov Chain -----//
    run_monte_carlo(progressbar);

    //----- Close progress bar -----//
    progressbar.finalize();
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Set inverse temperature and its relative functional methods. */
void MCMC::update_inver_temper(double inv_T)
{
    if ((inv_T < 0.) || (inv_T > 1.))
    {
        throw std::range_error(
            "Inverse temperature " + std::to_string(inv_T) + " out of range [0, 1].");
    }
    else if (inv_T == 1.0)
    {
        pBuild->set_inv_T(1.0);
        _func_param_metropolis_ = &MCMC::_param_metropolis_1_;
        _func_augm_metropolis_ = &MCMC::_augm_metropolis_1_;
    }
    else if (inv_T == 0.0)
    {
        pBuild->set_inv_T(0.0);
        _func_param_metropolis_ = &MCMC::_param_metropolis_0_;
        _func_augm_metropolis_ = &MCMC::_augm_metropolis_0_;
    }
    else
    {
        pBuild->set_inv_T(inv_T);
        _func_param_metropolis_ = &MCMC::_param_metropolis_T_;
        _func_augm_metropolis_ = &MCMC::_augm_metropolis_T_;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : power posteriors  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Evaluate the logarithm of marginal likelihood with power posteriors. */
double MCMC::evaluate_marginal(int n_rungs, bool adaptive, bool is_progbar)
{
    /* Ref: Friel, N., Hurn, M. & Wyse, J.
    Improving power posterior estimation of statistical evidence.
    Stat Comput 24, 709–723 (2014).*/

    if (n_rungs < 2)
    {
        tools::error("Number of rungs must be greater than 1.");
    }
    else
    {
        // fix internal number of rungs (for progress bar)
        _n_rungs_ = n_rungs;
    }

    //----- Backup user setup to restore them later -----//
    bool bkp_rnd_init_params = pBuild->get_rnd_init_params();
    double bkp_inv_T = pBuild->get_inv_T();

    //----- Thermodynamic ladder estimates -----//
    bool is_right = false;
    double avg, var;
    vector<double> inv_t = {};   // inverse temperature
    vector<double> avg_vec = {}; // mean loglikelihood
    vector<double> var_vec = {}; // variance loglikelihood

    //----- Total number of mcmc sweeps -----//
    ulong rung_sweeps = pBuild->get_c_n_sweeps() * pBuild->get_n_chains();
    // internal state variables
    _tot_thin_sweeps_ = rung_sweeps * (_n_rungs_ + 1);
    _rung_ = 0;

    //----- Start progress bar and chrono -----//
    CustomProgressBar progressbar(is_progbar);
    progressbar.initialize(fmt_progtxt("", true));

    //==================================//
    //===== INVERSE TEMPEREATURE 0 =====//
    //==================================//

    update_inver_temper(0.0);
    pBuild->set_rnd_init_params(true);
    run_monte_carlo(progressbar, 0, true);
    pBuild->dgn_glob_ll_stats(avg, var);

    inv_t.push_back(0.0);
    avg_vec.push_back(avg);
    var_vec.push_back(var);

    // thermodynamic ladder: use modified power fractions or adaptive scheme
    if (adaptive)
    {
        //===========================//
        //===== ADAPTIVE LADDER =====//
        //===========================//

        //----- Inverse Temperature 1 -----//
        _rung_ = 1;
        update_inver_temper(1.0);
        pBuild->set_rnd_init_params(false);
        run_monte_carlo(progressbar, rung_sweeps, true);
        pBuild->dgn_glob_ll_stats(avg, var);

        inv_t.push_back(1.0);
        avg_vec.push_back(avg);
        var_vec.push_back(var);

        //----- Inverse Temperatures (0, 1) -----//
        size_t k = 0; // index to interval lower bound
        while (_rung_ < _n_rungs_)
        {
            ++_rung_;

            //----- adaptive choice of inverse temperature -----//
            // intersection between tangents
            double new_inv_T = (avg_vec[k + 1] - inv_t[k + 1] * var_vec[k + 1] - avg_vec[k] + inv_t[k] * var_vec[k]) / (var_vec[k] - var_vec[k + 1]);

            // otherwise weighted average
            if ((new_inv_T <= inv_t[k]) || (new_inv_T >= inv_t[k + 1]))
            {
                new_inv_T = inv_t[k] + var_vec[k + 1] * (inv_t[k + 1] - inv_t[k]) / (var_vec[k] + var_vec[k + 1]);
            }

            //----- power posterior sampling -----//
            update_inver_temper(new_inv_T);
            run_monte_carlo(progressbar, _rung_ * rung_sweeps, true);
            pBuild->dgn_glob_ll_stats(avg, var);

            //----- find larger rectangular discretization error -----//
            is_right = (inv_t[k + 1] - pBuild->get_inv_T()) * (avg_vec[k + 1] - avg) > (new_inv_T - inv_t[k]) * (avg - avg_vec[k]);

            //----- insert in-between for all lists -----//
            inv_t.insert(inv_t.begin() + k + 1, pBuild->get_inv_T());
            avg_vec.insert(avg_vec.begin() + k + 1, avg);
            var_vec.insert(var_vec.begin() + k + 1, var);

            // if right interval is chosen, index points to inserted as lower bound
            k += (is_right);
        }
    }
    else
    {
        //==================================//
        //===== POWER FRACTIONS LADDER =====//
        //==================================//

        //----- Inverse Temperatures (0, 1] -----//
        pBuild->set_rnd_init_params(false);
        for (_rung_ = 1; _rung_ < _n_rungs_; ++_rung_)
        {
            update_inver_temper(std::pow(static_cast<double>(_rung_) / (_n_rungs_ - 1), 5));
            run_monte_carlo(progressbar, _rung_ * rung_sweeps, true);
            pBuild->dgn_glob_ll_stats(avg, var);

            inv_t.push_back(pBuild->get_inv_T());
            avg_vec.push_back(avg);
            var_vec.push_back(var);
        }
    }

    //----- trapezoidal rule to compute marginal likelihood -----//
    double d_inv_t, output = 0.;
    for (int k = 0; k < n_rungs; ++k)
    {
        d_inv_t = inv_t[k + 1] - inv_t[k];
        output += (avg_vec[k + 1] + avg_vec[k]) * d_inv_t / 2. - (var_vec[k + 1] - var_vec[k]) * std::pow(d_inv_t, 2) / 12.;
    }

    //----- Close progress bar -----//
    progressbar.finalize();

    //----- Restore user setup -----//
    update_inver_temper(bkp_inv_T);
    pBuild->set_rnd_init_params(bkp_rnd_init_params);

    return output;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : methods (simulation)  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Return the time of infection of a patient. */
void MCMC::simulate_single_house(Household &house)
{
    UNUSED(house);
}

/* Perform simulation as provided by the user */
void MCMC::simulate(bool is_progbar)
{
    //----- Check no undefined parameters -----//
    pBuild->set_rnd_init_params(false);
    pBuild->init_rand_params();

    //----- Initialize model -----//
    init_model_for_simulation();

    //----- Initialize augmented dependencies -----//
    init_augmented_variables();
    CustomProgressBar progressbar(is_progbar);
    progressbar.initialize("simulate");

#ifdef _OPENMP
    int thin_n_house = n_houses() / pBuild->get_n_threads();
#pragma omp parallel for default(shared)
    for (size_t H = 0; H < n_houses(); ++H)
    {
        if (_error_occurred_)
        {
            continue;
        }
        try
        {
            simulate_single_house(pBuild->house_vector[H]);
        }
        catch (const std::exception &e)
        {
#pragma omp critical
            {
                if (!_error_occurred_)
                {
                    _error_occurred_ = true;
                    _error_message_ = e.what();
                }
            }
        }
        if (omp_get_thread_num() == 0)
        {
            progressbar.update(static_cast<double>(H) / thin_n_house);
        }
    }
#else
    for (size_t H = 0; H < n_houses(); ++H)
    {
        if (_error_occurred_)
        {
            continue;
        }
        try
        {
            simulate_single_house(pBuild->house_vector[H]);
        }
        catch (const std::exception &e)
        {
            _error_occurred_ = true;
            _error_message_ = e.what();
        }
        progressbar.update(static_cast<double>(H) / n_houses());
    }
#endif

    if (_error_occurred_)
    {
        tools::error(_error_message_);
    }
    else
    {
        progressbar.finalize();
        pBuild->update_report_data();
    }
}
