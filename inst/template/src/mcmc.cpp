/*
 * mcmc.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "mcmc.h"
#include "build.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : implementations  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

MCMC::MCMC() : pImpl(std::make_unique<Monte_Carlo_Build>()) {}
MCMC::~MCMC() = default;

void MCMC::read_setup(const std::string &pathtofile) { pImpl->read_setup(pathtofile); }
void MCMC::read_params(const std::string &pathtofile) { pImpl->read_params(pathtofile); }
void MCMC::read_augm(const std::string &pathtofile) { pImpl->read_augm(pathtofile); }

#ifdef _HOUSEPI_RCPP
//===== Rcpp IMPORTS =====//
void MCMC::import_setup(const Rcpp::List &x) { pImpl->import_setup(x); }
void MCMC::import_params(const Rcpp::DataFrame &x) { pImpl->import_params(x); }
void MCMC::import_augm(const Rcpp::DataFrame &x) { pImpl->import_augm(x); }

//===== Rcpp EXPORTS =====//
Rcpp::List MCMC::export_setup() const { return pImpl->export_setup(); }
Rcpp::DataFrame MCMC::export_params() const { return pImpl->export_params(); }
Rcpp::DataFrame MCMC::export_data() const { return pImpl->export_data(); }
Rcpp::DataFrame MCMC::export_traces() const { return pImpl->export_traces(); }
Rcpp::NumericMatrix MCMC::export_loglike_matrix() const { return pImpl->export_loglike_matrix(); }
#else

//===== Cpp PRINT =====//
void MCMC::print_infer_state(const std::string &pathtofile) const { pImpl->print_infer_state(pathtofile); }
void MCMC::print_augm_variables(const std::string &pathtofile) const { pImpl->print_augm_variables(pathtofile); }
#endif

//===== SHOW =====//
void MCMC::show_full_build() const { pImpl->show_full_build(); }
void MCMC::show_results() const { pImpl->show_results(); }
void MCMC::show_summary() const { pImpl->show_summary(); }
void MCMC::show_group_classes() const { pImpl->show_group_classes(); }

//===== DATA =====//
int MCMC::house_size(int H) const { return pImpl->house_vec[H].size(); }
size_t MCMC::n_houses() const { return pImpl->n_house; }
int MCMC::report_data(const std::string &key) const { return pImpl->report_data(key); }
double MCMC::uniform() { return pImpl->uniform(); }

//===== AUGMENTED =====//
double MCMC::param(const std::string &parID) const { return pImpl->get_param(parID); }
const double *MCMC::param_pointer(const std::string &parID) const { return pImpl->get_param_pointer(parID); }
void MCMC::assign_all_augm(Household &this_house, int i) { pImpl->assign_all_augm(this_house, i); }
void MCMC::clear_all_augm(Household &this_house, int i) { pImpl->clear_all_augm(this_house, i); }

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
void MCMC::update_augm_dependecies(Household &this_house, int i)
{
    UNUSED(this_house);
    UNUSED(i);
}

/* Compute the loglikelihood for a single household. */
double MCMC::compute_house_loglike(const Household &this_house) const
{
    UNUSED(this_house);
    return 0.0;
};

/* Assign group rate to an individual  */
void MCMC::assign_group_rate_to_individual(Household &this_house, int i)
{
    // this_indiv.assign_group_rate(EMPTY_RATE, EMPTY_RATE);
    UNUSED(this_house);
    UNUSED(i);
}

/* Assign the contact patterns between a couple of individuals. */
void MCMC::assign_contact_rate_to_household(Household &this_house, int from, int to)
{
    this_house.assign_contact_rate(from, to, EMPTY_RATE);
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
void MCMC::read_data(const std::string &pathtofile)
{
    std::map<std::string, std::vector<std::string>> datatable = {};
    tools::read_csv_table(pathtofile, datatable);
    build_data(datatable);

    pImpl->_report_pathtofile_["data"] = pathtofile;
}

#ifdef _HOUSEPI_RCPP
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  Monte_Carlo_Build : imports from Rcpp ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Import data from R DataFrame and run build. */
void MCMC::import_data(const Rcpp::DataFrame &data)
{
    std::map<std::string, std::vector<std::string>> datatable = {};
    tools::rdataframe_to_map(data, datatable);
    build_data(datatable);
    if (pImpl->_report_pathtofile_.find("data") != pImpl->_report_pathtofile_.end())
    {
        pImpl->_report_pathtofile_.erase("data");
    }
}
#endif

/* Assign data after pre-processing. */
void MCMC::build_data(const std::map<std::string, std::vector<std::string>> &datatable)
{
    try
    {
        if (datatable.empty())
        {
            // If nothing was previously defined
            if (pImpl->house_vec.empty())
            {
                throw std::length_error("No database was found.");
            }
        }
        else
        {
            //=====  Drop all Previous Data and Prepare for New  =====//
            pImpl->data_headers.clear();
            for (const auto &x : datatable)
            {
                pImpl->data_headers.push_back(x.first);
            }
            pImpl->n_indiv = datatable.at(pImpl->data_headers[0]).size();
            pImpl->n_house = 0;
            pImpl->house_vec.clear();
            pImpl->_report_group_classes_.clear();

            //===== Assess the presence of house IDs =====//
            if (std::find(pImpl->data_headers.begin(), pImpl->data_headers.end(), "house_id") == pImpl->data_headers.end())
            {
                tools::warning(
                    "Unknown household structure: missing data column 'house_id'. If intentional, ignore this warning.");
            }

            //===== Safe Load of Households =====//
            std::unordered_map<std::string, size_t> hh_ID_to_index = {};
            for (size_t i = 0; i < pImpl->n_indiv; ++i)
            {
                //----- Read individual data -----//
                std::string this_indiv_ID = std::to_string(i);            // default value
                std::string this_house_ID = std::to_string(n_houses()); // default value
                int this_status = MISSING_STATUS;
                dictargs args = {};

                for (const std::string &key : pImpl->data_headers)
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
                            this_house_ID = datatable.at("house_id")[i];
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
                        tools::warning("Found NA value in database. If intentional, ignore this warning.");
                    }
                    else
                    {
                        args.emplace(key, tools::stod(datatable.at(key)[i], "data '" + key + "'"));
                    }
                }

                // if household index doesn't exist, add new household
                // note : this is redundant if no structure
                if (hh_ID_to_index.find(this_house_ID) == hh_ID_to_index.end())
                {
                    // map data index to household number
                    hh_ID_to_index.emplace(this_house_ID, n_houses());
                    // create a new empty household
                    pImpl->assign_house(this_house_ID);
                }

                Household &this_house = pImpl->house_vec[hh_ID_to_index.at(this_house_ID)];
                // Assign the individual the respective household
                this_house.assign_individual(this_indiv_ID, this_status, args);
                // Assign and set rate per group
                set_rate_per_group(this_house, this_house.size() - 1);
            }

            //===== Set Contact Rates =====//
            pImpl->_n_contacts_ = 0;
            for (Household &this_house : pImpl->house_vec)
            {
                set_rate_per_contact(this_house);
            }

            //===== Update Report for Data =====//
            pImpl->update_report_data();

            //===== Clean up data entries for output files =====//
            // eventually remove standard headers to reorder them at the end
            for (const std::string &col : {"house_id", "id", "status"})
            {
                pImpl->data_headers.erase(
                    std::remove(pImpl->data_headers.begin(), pImpl->data_headers.end(), col),
                    pImpl->data_headers.end());
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
void MCMC::set_rate_per_group(Household &this_house, int i)
{
    try
    {
        assign_group_rate_to_individual(this_house, i);

        for (const auto &pair : this_house.get_group_rate_map(i))
        {
            std::string rateID = pair.first;
            if (rateID == CONTACT_RATE)
            {
                throw std::invalid_argument(
                    "The group rate name '" + CONTACT_RATE + "' is reserved. Choose a different one.");
            }

            std::string parID = pair.second;

            //=====initialize model report=====//
            // if first individual, initialize rateID-report
            pImpl->_report_group_classes_.emplace(rateID, std::map<std::string, int>());
            pImpl->_report_group_classes_.at(rateID).emplace("[reference]", 0);

            //=====set rates and update report=====//
            if (tools::is_na(parID))
            {
                // error : undefined rate
                std::invalid_argument("Rate '" + rateID + "' is not provided for individual " + this_house.ID(i) + ".");
            }
            else if (parID == REFERENCE)
            {
                // rate set to reference
                this_house.set_rate(i, rateID, &REFERENCE_VALUE);
                pImpl->_report_group_classes_.at(rateID).at("[reference]") += 1;
            }
            else if (parID != REFERENCE)
            {
                // rate associated to parameter
                pImpl->assert_param(parID);
                this_house.set_rate(i, rateID, pImpl->get_param_pointer(parID));
                if (pImpl->_report_group_classes_.at(rateID).find(parID) == pImpl->_report_group_classes_.at(rateID).end())
                {
                    pImpl->_report_group_classes_.at(rateID).emplace(parID, 1);
                }
                else
                {
                    pImpl->_report_group_classes_.at(rateID).at(parID) += 1;
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
void MCMC::set_rate_per_contact(Household &this_house)
{
    try
    {
        // allocate contact patterns
        int hh_size = this_house.size();
        this_house.init_contact_matrix();

        if (this_house.index() == 0)
        {
            // if first household, initialize contact rate report
            pImpl->_report_group_classes_.emplace(CONTACT_RATE, std::map<std::string, int>());
            pImpl->_report_group_classes_.at(CONTACT_RATE).emplace("[reference]", 0);
        }
        else
        {
            // all household must have a contact rate
            if (pImpl->_report_group_classes_.find(CONTACT_RATE) == pImpl->_report_group_classes_.end())
            {
                throw std::invalid_argument(
                    "Contact rate has not been provided for all households.");
            }
        }

        // loop over all couples
        for (int k = 0; k < hh_size; ++k) // from
        {
            for (int i = 0; i < hh_size; ++i) // to
            {
                // note : from 'k' --> to 'i'
                if (k != i)
                {
                    assign_contact_rate_to_household(this_house, k, i);
                    std::string parID = this_house.get_contact_rate_map(k, i);

                    if (tools::is_na(parID))
                    {
                        std::string sms = "Contact rate not provided for individuals " + this_house.ID(k);
                        sms += " and " + this_house.ID(i) + " within household " + this_house.house_ID() + ".";
                        throw std::invalid_argument(sms);
                    }
                    else if (parID == REFERENCE)
                    {
                        this_house.set_contact_rate(k, i, &REFERENCE_VALUE); // reference rate
                        pImpl->_report_group_classes_.at(CONTACT_RATE).at("[reference]") += 1;
                    }
                    else if (parID != EMPTY_RATE)
                    {
                        pImpl->assert_param(parID);
                        this_house.set_contact_rate(k, i, pImpl->get_param_pointer(parID)); // parameter rate
                        if (pImpl->_report_group_classes_.at(CONTACT_RATE).find(parID) == pImpl->_report_group_classes_.at(CONTACT_RATE).end())
                        {
                            pImpl->_report_group_classes_.at(CONTACT_RATE).emplace(parID, 1);
                        }
                        else
                        {
                            pImpl->_report_group_classes_.at(CONTACT_RATE).at(parID) += 1;
                        }
                    }

                    // count the total number of pairs
                    ++pImpl->_n_contacts_;
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
            Household &this_house = pImpl->house_vec[H];
            // update house loglikelihood
            this_house.set_loglike(compute_house_loglike(this_house));
            global_loglike += this_house.get_loglike() + pImpl->compute_augm_model_loglike(this_house);
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
void MCMC::move_single_param(const std::string &parID)
{
    bool is_accepted = false;

    // proceed only if the new value is consistent, otherwise automatically execute refuse_move(parID)
    if (pImpl->propose_valid_param_move(parID))
    {
        if ((this->*_func_param_metropolis_)(parID))
        {
            is_accepted = true;
        }
        else
        {
            // restore buffered value of parameter
            pImpl->restore_param_value(parID);

            // restore buffered value of global loglike
            restore_global_loglike();

            // restore buffered loglikelihood in households
            for (size_t H = 0; H < n_houses(); ++H)
            {
                pImpl->house_vec[H].restore_loglike();
            }
        }
    }

    // update counts for acceptance rate
    pImpl->update_param_accept_rate(parID, is_accepted);
}

bool MCMC::_param_metropolis_1_(const std::string &parID)
{
    update_global_loglike();

    if (uniform() < exp(get_delta_global_loglike() + pImpl->param_move_logratio(parID)))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool MCMC::_param_metropolis_0_(const std::string &parID)
{
    if (uniform() < exp(pImpl->param_move_logratio(parID)))
    {
        update_global_loglike();
        return true;
    }
    else
    {
        return false;
    }
}

bool MCMC::_param_metropolis_T_(const std::string &parID)
{
    update_global_loglike();

    if (uniform() < exp(pImpl->get_inv_T() * get_delta_global_loglike() + pImpl->param_move_logratio(parID)))
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
    for (const std::string &parID : pImpl->get_param_to_be_inferred())
    {
        move_single_param(parID);
    }
}

/* Move all model parameters and adapt their rule arguments. */
void MCMC::move_and_adapt_all_parameters(ulong step)
{
    double adapt_ratio = static_cast<double>(step) / pImpl->get_n_adapt();

    for (const std::string &parID : pImpl->get_param_to_be_inferred())
    {
        move_single_param(parID);
        pImpl->adapt_param_rule(parID, adapt_ratio);
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
    if (!pImpl->augm_list_full.empty())
    {
        if (pImpl->is_augm_infer)
        {
            _func_move_all_augm_ = &MCMC::_move_all_augm_;

            if (!pImpl->augm_list_infer.empty())
            {
                _func_move_augm_indiv_ = &MCMC::_move_augm_indiv_generic_;
            }
        }

        for (size_t H = 0; H < n_houses(); ++H)
        {
            Household &this_house = pImpl->house_vec[H];
            for (int i : this_house.infectious_set())
            {
                update_augm_dependecies(this_house, i);
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
            delta_ll += move_augm_single_house(pImpl->house_vec[H]);
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
double MCMC::move_augm_single_house(Household &this_house)
{
    double delta_global_loglike = 0.0;

    // Loop over infectious residents of the househould
    for (const int &i : this_house.infectious_set())
    {
        (this->*_func_move_augm_indiv_)(this_house, i, delta_global_loglike);
    }

    return delta_global_loglike;
}

/* Generic move augmented variables one by one. */
void MCMC::_move_augm_indiv_generic_(Household &this_house, int i, double &delta)
{
    // Loop over augmented variables of infectious
    for (const std::string &augmID : pImpl->augm_list_infer)
    {
        metropolis_augm_single_move(this_house, i, augmID, delta);
    }
}

/* Move augmented variables by status one variable at a time. */
void MCMC::_move_augm_indiv_by_status_(Household &this_house, int i, double &delta)
{
    // Loop over augmented variables of infectious
    for (const std::string &augmID : pImpl->augm_status_map_infer.at(this_house.status(i)))
    {
        metropolis_augm_single_move(this_house, i, augmID, delta);
    }
}

/* Take care of a a single augmented Metropolis move. */
void MCMC::metropolis_augm_single_move(Household &this_house, int i, const std::string &augmID, double &delta)
{
    // proceed only if the new value is consistent
    if (this_house.propose_valid_augm_move(i, augmID, pImpl->get_prng()))
    {
        // update dependencies to perform computations
        update_augm_dependecies(this_house, i);

        // update house loglikelihood and store old value in buffer
        this_house.set_loglike(compute_house_loglike(this_house));

        // Metropolis-Hasting algorithm
        if (!(this->*_func_augm_metropolis_)(this_house, i, augmID, delta))
        {
            // restore value and its dependencies
            this_house.augm_restore_value(i, augmID);
            update_augm_dependecies(this_house, i);
            this_house.restore_loglike();
        }
    }
}

/* Metropolis algorithm for single augmented variable move (default, inverse temperature 1). */
bool MCMC::_augm_metropolis_1_(const Household &this_house, int i, const std::string &augmID, double &delta)
{
    double D_hh_ll = this_house.get_delta_loglike();
    if (uniform() < exp(D_hh_ll + this_house.augm_move_logratio(i, augmID)))
    {
        delta += D_hh_ll + this_house.augm_model_logratio(i, augmID);
        return true;
    }
    else
    {
        return false;
    }
}

/* Metropolis algorithm for single augmented variable move (inverse temperature 0). */
bool MCMC::_augm_metropolis_0_(const Household &this_house, int i, const std::string &augmID, double &delta)
{
    if (uniform() < exp(-this_house.augm_propose_logratio(i, augmID)))
    {
        delta += this_house.get_delta_loglike() + this_house.augm_model_logratio(i, augmID);
        return true;
    }
    else
    {
        return false;
    }
}

/* Metropolis algorithm for single augmented variable move (all inverse temperatures). */
bool MCMC::_augm_metropolis_T_(const Household &this_house, int i, const std::string &augmID, double &delta)
{
    double D_ll = this_house.get_delta_loglike() + this_house.augm_model_logratio(i, augmID);
    if (uniform() < exp(pImpl->get_inv_T() * D_ll - this_house.augm_propose_logratio(i, augmID)))
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
    //----- Check no undefined parameters -----//
    pImpl->set_rnd_init_params(false);
    pImpl->init_rand_params();

    //----- Initialize model -----//
    init_model_for_inference();

    //----- Initialize augmented dependencies -----//
    init_augmented_variables();

    //----- Initialize global loglikelihood -----//
    update_global_loglike();

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
// UPDATE : this could be functional too
const std::string MCMC::fmt_progtxt(const std::string &txt, bool is_power_post)
{
    if (is_power_post)
    {
        return "power post " + std::to_string(_rung_) + "/" + std::to_string(_n_rungs_);
    }
    else if (pImpl->get_n_chains() > 1)
    {
        return txt + " " + std::to_string(_chain_) + "/" + std::to_string(pImpl->get_n_chains());
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
    pImpl->init_diagnostics(is_power_post);

    //===== Run Multiple MCMC Chains =====//
    for (_chain_ = 1; _chain_ <= pImpl->get_n_chains(); ++_chain_)
    {
        //-----Initialize parameter inference-----//
        pImpl->init_rand_params();

        //----- Initialize augmented dependencies -----//
        init_augmented_variables();

        //----- Initialize global loglikelihood -----//
        update_global_loglike();

        // avoid asserting user definitions after first round
        // TODO : this can be improved, e. g. with a debug mode
        pImpl->disable_assert();

        //=========================================//
        //======== (Finite) Adaptative MCMC =======//
        //=========================================//

        progressbar.set_text(fmt_progtxt("adapt", is_power_post));

        // initialize the param move counts for adaptive
        pImpl->init_move_counts();

        for (ulong i = 0; i < pImpl->get_n_adapt(); i += pImpl->get_thin())
        {
            for (ulong t = 1; t <= pImpl->get_thin(); ++t)
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

        progressbar.set_text(fmt_progtxt("burn", is_power_post));

        for (ulong i = 0; i < pImpl->get_n_burnin(); i += pImpl->get_thin())
        {
            for (ulong t = 1; t <= pImpl->get_thin(); ++t)
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

        progressbar.set_text(fmt_progtxt("sample", is_power_post));

        for (ulong i = pImpl->get_n_burnin(); i < pImpl->get_n_steps(); i += pImpl->get_thin())
        {
            // thinning mcmc moves
            for (ulong t = 1; t <= pImpl->get_thin(); ++t)
            {
                (this->*_func_move_all_augm_)();
                move_all_parameters();
            }
            ++sweep;

            // compute diagnostics for the previous iterations
            pImpl->update_diagnostics(get_global_loglike(), is_power_post);
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
    if (report_data("n.miss.stat") == static_cast<int>(pImpl->n_indiv))
    {
        tools::error("No infectious status has been provided in the data file.");
    }
    */

    //----- total number of mcmc sweeps -----//
    _tot_thin_sweeps_ = pImpl->get_c_n_sweeps() * pImpl->get_n_chains();

    //----- Start progress bar and chrono -----//
    CustomProgressBar progressbar(is_progbar);
    progressbar.initialize();

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
        pImpl->set_inv_T(1.0);
        _func_param_metropolis_ = &MCMC::_param_metropolis_1_;
        _func_augm_metropolis_ = &MCMC::_augm_metropolis_1_;
    }
    else if (inv_T == 0.0)
    {
        pImpl->set_inv_T(0.0);
        _func_param_metropolis_ = &MCMC::_param_metropolis_0_;
        _func_augm_metropolis_ = &MCMC::_augm_metropolis_0_;
    }
    else
    {
        pImpl->set_inv_T(inv_T);
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
    bool bkp_rnd_init_params = pImpl->get_rnd_init_params();
    double bkp_inv_T = pImpl->get_inv_T();

    //----- Thermodynamic ladder estimates -----//
    bool is_right = false;
    double avg, var;
    std::vector<double> inv_t = {};   // inverse temperature
    std::vector<double> avg_vec = {}; // mean loglikelihood
    std::vector<double> var_vec = {}; // variance loglikelihood

    //----- Total number of mcmc sweeps -----//
    ulong rung_sweeps = pImpl->get_c_n_sweeps() * pImpl->get_n_chains();
    _tot_thin_sweeps_ = rung_sweeps * (_n_rungs_ + 1);

    //----- Start progress bar and chrono -----//
    CustomProgressBar progressbar(is_progbar);

    //==================================//
    //===== INVERSE TEMPEREATURE 0 =====//
    //==================================//

    _rung_ = 0;
    update_inver_temper(0.0);
    pImpl->set_rnd_init_params(true);
    run_monte_carlo(progressbar, 0, true);
    pImpl->dgn_glob_ll_stats(avg, var);

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
        pImpl->set_rnd_init_params(false);
        run_monte_carlo(progressbar, rung_sweeps, true);
        pImpl->dgn_glob_ll_stats(avg, var);

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
            pImpl->dgn_glob_ll_stats(avg, var);

            //----- find larger rectangular discretization error -----//
            is_right = (inv_t[k + 1] - pImpl->get_inv_T()) * (avg_vec[k + 1] - avg) > (new_inv_T - inv_t[k]) * (avg - avg_vec[k]);

            //----- insert in-between for all lists -----//
            inv_t.insert(inv_t.begin() + k + 1, pImpl->get_inv_T());
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
        pImpl->set_rnd_init_params(false);
        for (_rung_ = 1; _rung_ < _n_rungs_; ++_rung_)
        {
            update_inver_temper(std::pow(static_cast<double>(_rung_) / (_n_rungs_ - 1), 5));
            run_monte_carlo(progressbar, _rung_ * rung_sweeps, true);
            pImpl->dgn_glob_ll_stats(avg, var);

            inv_t.push_back(pImpl->get_inv_T());
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
    pImpl->set_rnd_init_params(bkp_rnd_init_params);

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
void MCMC::simulate_single_house(Household &this_house)
{
    UNUSED(this_house);
}

/* Perform simulation as provided by the user */
void MCMC::simulate(bool is_progbar)
{
    //----- Check no undefined parameters -----//
    pImpl->set_rnd_init_params(false);
    pImpl->init_rand_params();

    //----- Initialize model -----//
    init_model_for_simulation();

    //----- Initialize augmented dependencies -----//
    init_augmented_variables();
    CustomProgressBar progressbar(is_progbar);
    progressbar.set_text("simul");

#ifdef _OPENMP
    int thin_n_house = n_houses() / pImpl->get_n_threads();
#pragma omp parallel for default(shared)
    for (size_t H = 0; H < n_houses(); ++H)
    {
        if (_error_occurred_)
        {
            continue;
        }
        try
        {
            simulate_single_house(pImpl->house_vec[H]);
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
            simulate_single_house(pImpl->house_vec[H]);
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
        pImpl->update_report_data();
    }
}
