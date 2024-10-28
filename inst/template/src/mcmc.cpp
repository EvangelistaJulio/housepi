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

#ifdef USE_RCPP
//===== Rcpp IMPORTS =====//
void MCMC::import_setup(Rcpp::List x) { pImpl->import_setup(x); }
void MCMC::import_params(Rcpp::DataFrame x) { pImpl->import_params(x); }
void MCMC::import_augm(Rcpp::DataFrame x) { pImpl->import_augm(x); }

//===== Rcpp EXPORTS =====//
Rcpp::List MCMC::export_setup() const { return pImpl->export_setup(); }
Rcpp::DataFrame MCMC::export_params() const { return pImpl->export_params(); }
Rcpp::DataFrame MCMC::export_data() const { return pImpl->export_data(); }
Rcpp::DataFrame MCMC::export_traces() const { return pImpl->export_traces(); }
Rcpp::NumericMatrix MCMC::export_loglike_matrix() const { return pImpl->export_loglike_matrix(); }
#else
//===== Cpp SHOW =====//
void MCMC::show_full_build() const { pImpl->show_full_build(); }
void MCMC::show_results() const { pImpl->show_results(); }

//===== Cpp PRINT =====//
void MCMC::print_infer_state(const std::string &pathtofile) const { pImpl->print_infer_state(pathtofile); }
void MCMC::print_augm_variables(const std::string &pathtofile) const { pImpl->print_augm_variables(pathtofile); }
#endif

//===== SHOW =====//
void MCMC::show_summary() const { pImpl->show_summary(); }
void MCMC::show_group_classes() const { pImpl->show_group_classes(); }

//===== DATA =====//
size_t MCMC::get_num_house() const { return pImpl->num_house; }
Household &MCMC::house(size_t H) { return pImpl->house[H]; }
int MCMC::report_data(const std::string &key) const { return pImpl->report_data(key); }

double MCMC::uniform() { return pImpl->uniform(); }

//===== AUGMENTED =====//
double MCMC::get_param(const std::string &parID) const { return pImpl->get_param(parID); }
const double *MCMC::get_param_pointer(const std::string &parID) const { return pImpl->get_param_pointer(parID); }
void MCMC::assign_augm(Individual &this_indiv) { pImpl->assign_augm(this_indiv); }
void MCMC::clear_augm(Individual &this_indiv) { pImpl->clear_augm(this_indiv); }

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
void MCMC::assign_group_rate_to_individual(Individual &this_indiv)
{    
    // this_indiv.assign_group_rate(EMPTY_RATE, EMPTY_RATE);
    UNUSED(this_indiv);
}

/* Assign the contact patterns between a couple of individuals. */
void MCMC::assign_contact_rate_to_household(Household &this_house, int from, int to)
{
    this_house.assign_rate_per_contact(from, to, EMPTY_RATE);
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

#ifdef USE_RCPP
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  Monte_Carlo_Build : imports from Rcpp ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Import data from R DataFrame and run build. */
void MCMC::import_data(Rcpp::DataFrame data)
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
            if (pImpl->house.empty())
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
            pImpl->num_indiv = datatable.at(pImpl->data_headers[0]).size();
            pImpl->num_house = 0;
            pImpl->house.clear();
            pImpl->_report_group_classes_.clear();

            //===== Assess the presence of house IDs =====//
            if (std::find(pImpl->data_headers.begin(), pImpl->data_headers.end(), "house_id") == pImpl->data_headers.end())
            {
                std::ostringstream oss;
                oss << "Impossible to associate a household structure: missing data column 'house_id'."
                    << " If intentional, ignore this message.";
                tools::warning(oss.str());
            }

            //===== Safe Load of Households =====//
            std::unordered_map<std::string, int> aux_hh_idx = {};
            for (size_t i = 0; i < pImpl->num_indiv; ++i)
            {
                //----- Read individual data -----//
                std::string this_indiv_id = std::to_string(i);
                std::string this_house_id = std::to_string(pImpl->num_house);
                int this_status = MISSING_STATUS;
                dictargs args = {};
                for (const std::string &key : pImpl->data_headers)
                {
                    if (key == "id")
                    {
                        if (tools::is_na(datatable.at("id")[i]))
                        {
                            throw std::invalid_argument("Missing id for individual in row n." + std::to_string(i + 2) + ".");
                        }
                        else
                        {
                            this_indiv_id = datatable.at("id")[i];
                        }
                    }
                    else if (key == "house_id")
                    {
                        if (tools::is_na(datatable.at("house_id")[i]))
                        {
                            throw std::invalid_argument("Missing house_id for individual in row n." + std::to_string(i + 2) + ".");
                        }
                        else
                        {
                            this_house_id = datatable.at("house_id")[i];
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
                // create individual instance
                Individual this_indiv(this_indiv_id, this_status, args);

                //----- Assign and set rate per group -----//
                set_rate_per_group(this_indiv);

                // if household index doesn't exist, add new household
                // note : this is redundant if no structure
                if (aux_hh_idx.find(this_house_id) == aux_hh_idx.end())
                {
                    // map data index to household number
                    aux_hh_idx.emplace(this_house_id, pImpl->num_house);
                    // create a new empty household
                    pImpl->house.emplace_back();
                    // assign the user data id
                    pImpl->house[pImpl->num_house].set_ID(this_house_id);
                    // set household number as its index
                    pImpl->house[pImpl->num_house].set_index(pImpl->num_house);
                    // increment household number
                    ++pImpl->num_house;
                }
                // assign the individual the respective household
                pImpl->house[aux_hh_idx.at(this_house_id)].assign_individual(this_indiv);
            }

            //===== Set Contact Rates =====//
            pImpl->_num_contacts_ = 0;
            for (Household &this_house : pImpl->house)
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
void MCMC::set_rate_per_group(Individual &this_indiv)
{
    try
    {
        assign_group_rate_to_individual(this_indiv);

        for (const auto &which : this_indiv.which_rate)
        {
            std::string rateID = which.first;
            if (rateID == CONTACT_RATE)
            {
                throw std::invalid_argument("The group rate name '" + CONTACT_RATE + "' is reserved. Choose a different one.");
            }

            std::string parID = which.second;

            //=====initialize model report=====//
            // if first individual, initialize rateID-report
            pImpl->_report_group_classes_.emplace(rateID, std::map<std::string, int>());
            pImpl->_report_group_classes_.at(rateID).emplace("[reference]", 0);

            //=====set rates and update report=====//
            if (tools::is_na(parID))
            {
                // error : undefined rate
                std::invalid_argument("Rate '" + rateID + "' is not provided for individual " + this_indiv.get_ID() + ".");
            }
            else if (parID == REFERENCE)
            {
                // rate set to reference
                this_indiv.set_rate(rateID, &REFERENCE_VALUE);
                pImpl->_report_group_classes_.at(rateID).at("[reference]") += 1;
            }
            else if (parID != REFERENCE)
            {
                // rate associated to parameter
                pImpl->assert_param(parID);
                this_indiv.set_rate(rateID, pImpl->get_param_pointer(parID));
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
        this_house.which_rate.resize(hh_size, std::vector<std::string>(hh_size, UNDEFINED));
        this_house.contact_rates.resize(hh_size, std::vector<const double *>(hh_size, &REFERENCE_VALUE));

        if (this_house.get_index() == 0)
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
                throw std::invalid_argument("Contact rate has not been provided for all households.");
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
                    std::string parID = this_house.which_rate[k][i];

                    if (tools::is_na(parID))
                    {
                        std::ostringstream oss;
                        oss << "Contact rate not provided for individuals " << this_house.indiv[k].get_ID()
                            << " and " << this_house.indiv[i].get_ID() << " within household " << this_house.get_ID() << ".";
                        throw std::invalid_argument(oss.str());
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
                    ++pImpl->_num_contacts_;
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
    for (size_t H = 0; H < get_num_house(); ++H)
    {
        if (_error_occurred_)
        {
            continue;
        }
        try
        {
            Household &this_house = house(H);
            update_house_loglike(this_house);
            global_loglike += this_house.get_loglike() + pImpl->compute_augm_loglike(this_house);
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
        if (metropolis_param_move(parID))
        {
            is_accepted = true;
        }
        else
        {
            // restore buffered value of global loglike
            restore_global_loglike();

            // restore buffered value of parameter
            pImpl->restore_param_value(parID);

            // restore buffered loglikelihood in households
            for (size_t H = 0; H < get_num_house(); ++H)
            {
                house(H).restore_loglike();
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
    double adapt_ratio = static_cast<double>(step) / pImpl->get_num_adapt();

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
    for (size_t H = 0; H < get_num_house(); ++H)
    {
        Household &this_house = house(H);
        for (int i : this_house.infectious_set)
        {
            update_augm_dependecies(this_house, i);
        }
    }
}

/* Move augmented variables from all households. */
void MCMC::move_all_augm_variables()
{
    if (pImpl->is_move_augm)
    {
        double delta_ll = 0.;

#pragma omp parallel for reduction(+ : delta_ll)
        for (size_t H = 0; H < get_num_house(); ++H)
        {
            if (_error_occurred_)
            {
                continue;
            }
            try
            {
                delta_ll += move_single_house_augm(house(H));
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
        std::ostringstream oss;
        oss << "Inverse temperature " << inv_T << " out of range [0, 1].";
        throw std::range_error(oss.str());
    }
    else if (inv_T == 1.0)
    {
        pImpl->set_inv_T(1.0);
        _param_metropolis_method_ = &MCMC::_param_metropolis_1_;
        _augm_metropolis_method_ = &MCMC::_augm_metropolis_1_;
    }
    else if (inv_T == 0.0)
    {
        pImpl->set_inv_T(0.0);
        _param_metropolis_method_ = &MCMC::_param_metropolis_0_;
        _augm_metropolis_method_ = &MCMC::_augm_metropolis_0_;
    }
    else
    {
        pImpl->set_inv_T(inv_T);
        _param_metropolis_method_ = &MCMC::_param_metropolis_T_;
        _augm_metropolis_method_ = &MCMC::_augm_metropolis_T_;
    }
}

/* Move augmented variables for infectious members of a single household.
The move is performed one member at a time, one augmented variable at a time.
*/
double MCMC::move_single_house_augm(Household &this_house)
{
    double delta_global_loglike = 0.0;

    // Loop over infectious residents of the househould
    for (const int &i : this_house.infectious_set)
    {
        Individual &this_indiv = this_house.indiv[i];

        // Loop over augmented variables of infectious
        for (const std::string &augmID : pImpl->augm_to_be_inferred)
        {
            bool is_accepted = false;

            // proceed only if the new value is consistent
            if (this_indiv.propose_valid_augm_move(augmID, pImpl->get_prng()))
            {
                // update dependencies to perform computations
                update_augm_dependecies(this_house, i);
                update_house_loglike(this_house);

                // Metropolis-Hasting algorithm
                if (metropolis_augm_move(this_house, i, augmID, delta_global_loglike))
                {
                    is_accepted = true;
                }
                else
                {
                    // restore value and its dependencies
                    this_indiv.augm_restore_value(augmID);
                    update_augm_dependecies(this_house, i);
                    this_house.restore_loglike();
                }
            }

            this_house.update_augm_accept_rate(augmID, is_accepted);
        }
    }

    return delta_global_loglike;
}

/* Metropolis algorithm for single augmented variable move (default, inverse temperature 1). */
bool MCMC::_augm_metropolis_1_(const Household &this_house, int i, const std::string &augmID, double &delta)
{
    double D_hh_ll = this_house.get_delta_loglike();
    if (uniform() < exp(D_hh_ll + this_house.indiv[i].augm_move_logratio(augmID)))
    {
        delta += D_hh_ll + this_house.indiv[i].augm_model_logratio(augmID);
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
    if (uniform() < exp(-this_house.indiv[i].augm_propose_logratio(augmID)))
    {
        delta += this_house.get_delta_loglike() + this_house.indiv[i].augm_model_logratio(augmID);
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
    double D_ll = this_house.get_delta_loglike() + this_house.indiv[i].augm_model_logratio(augmID);
    if (uniform() < exp(pImpl->get_inv_T() * D_ll - this_house.indiv[i].augm_propose_logratio(augmID)))
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
std::string MCMC::get_progtxt(const std::string &txt, bool is_power_post)
{
    std::ostringstream oss;
    if (is_power_post)
    {
        oss << "power post " << _rung_ << "/" << _num_rungs_;
    }
    else if (pImpl->get_num_chains() > 1)
    {
        oss << txt << " " << _chain_ << "/" << pImpl->get_num_chains();
    }
    else
    {
        oss << txt;
    }
    return oss.str();
}

/* Markov chain core, used for inference and power posteriors.
    note: sweep is the number of previously completed thin sweeps. */
void MCMC::run_monte_carlo(bool progbar, ulong sweep, bool is_power_post)
{
    using namespace tools;

    //==========================//
    //===== Initialization =====//
    //==========================//

    //----- Initialize model -----//
    init_model_for_inference();

    //----- Initialize diagnostics -----//
    pImpl->init_diagnostics(is_power_post);

    //===== Run Multiple MCMC Chains =====//
    for (_chain_ = 1; _chain_ <= pImpl->get_num_chains(); ++_chain_)
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

        // initialize the acceptance counts
        pImpl->init_move_counts();

        //=========================================//
        //======== (Finite) Adaptative MCMC =======//
        //=========================================//

        for (ulong i = 0; i < pImpl->get_num_adapt(); i += pImpl->get_thinning())
        {
            for (ulong t = 1; t <= pImpl->get_thinning(); ++t)
            {
                move_all_augm_variables();
                move_and_adapt_all_parameters(i + t);
            }
            ++sweep;
            progress_bar_run(get_progtxt("adapt", is_power_post),
                             sweep, _tot_thin_sweeps_, progbar);
        }

        // reset the acceptance counts
        // UPDATE : move counts become useless since they are not even printed : improve.
        if (pImpl->get_num_adapt() > 0)
        {
            pImpl->init_move_counts();
        }

        //============================================//
        //========= MCMC Burn-in Iterations ==========//
        //============================================//

        for (ulong i = 0; i < pImpl->get_burn_in(); i += pImpl->get_thinning())
        {
            for (ulong t = 1; t <= pImpl->get_thinning(); ++t)
            {
                move_all_augm_variables();
                move_all_parameters();
            }
            ++sweep;
            progress_bar_run(get_progtxt("burn", is_power_post),
                             sweep, _tot_thin_sweeps_, progbar);
        }

        //============================================//
        //========= MCMC Posterior Sampling ==========//
        //============================================//

        for (ulong i = pImpl->get_burn_in(); i < pImpl->get_num_iter(); i += pImpl->get_thinning())
        {
            // thinning mcmc moves
            for (ulong t = 1; t <= pImpl->get_thinning(); ++t)
            {
                move_all_augm_variables();
                move_all_parameters();
            }
            ++sweep;

            // compute diagnostics for the previous iterations
            pImpl->update_diagnostics(get_global_loglike(), is_power_post);
            progress_bar_run(get_progtxt("sample", is_power_post),
                             sweep, _tot_thin_sweeps_, progbar);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  MCMC : infer  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Performs the inference using Metropolis-Hastings. */
void MCMC::infer(bool progbar)
{
    //----- Check if data is meaningful -----//
    if (report_data("n.miss.stat") == static_cast<int>(pImpl->num_indiv))
    {
        tools::error("No infectious status has been provided in the data file.");
    }

    //----- total number of mcmc sweeps -----//
    _tot_thin_sweeps_ = pImpl->get_c_num_sweeps() * pImpl->get_num_chains();

    //----- Start progress bar and chrono -----//
    tools::progress_bar_init("init", progbar);

    //----- Run Markov Chain -----//
    run_monte_carlo(progbar);

    //----- Close progress bar -----//
    tools::progress_bar_end("end", progbar);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : power posteriors  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Evaluate the logarithm of marginal likelihood with power posteriors. */
double MCMC::evaluate_marginal(int num_rungs, bool adaptive, bool progbar)
{
    /* Ref: Friel, N., Hurn, M. & Wyse, J.
    Improving power posterior estimation of statistical evidence.
    Stat Comput 24, 709–723 (2014).*/

    if (num_rungs < 2)
    {
        tools::error("Number of rungs must be greater than 1.");
    }
    else
    {
        // fix internal number of rungs (for progress bar)
        _num_rungs_ = num_rungs;
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
    ulong rung_sweeps = pImpl->get_c_num_sweeps() * pImpl->get_num_chains();
    _tot_thin_sweeps_ = rung_sweeps * (_num_rungs_ + 1);

    //----- Start progress bar and chrono -----//
    tools::progress_bar_init("init", progbar);

    //==================================//
    //===== INVERSE TEMPEREATURE 0 =====//
    //==================================//

    _rung_ = 0;
    update_inver_temper(0.0);
    pImpl->set_rnd_init_params(true);
    run_monte_carlo(progbar, 0, true);
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
        run_monte_carlo(progbar, rung_sweeps, true);
        pImpl->dgn_glob_ll_stats(avg, var);

        inv_t.push_back(1.0);
        avg_vec.push_back(avg);
        var_vec.push_back(var);

        //----- Inverse Temperatures (0, 1) -----//
        size_t k = 0; // index to interval lower bound
        while (_rung_ < _num_rungs_)
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
            run_monte_carlo(progbar, _rung_ * rung_sweeps, true);
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
        for (_rung_ = 1; _rung_ < _num_rungs_; ++_rung_)
        {
            update_inver_temper(std::pow(static_cast<double>(_rung_) / (_num_rungs_ - 1), 5));
            run_monte_carlo(progbar, _rung_ * rung_sweeps, true);
            pImpl->dgn_glob_ll_stats(avg, var);

            inv_t.push_back(pImpl->get_inv_T());
            avg_vec.push_back(avg);
            var_vec.push_back(var);
        }
    }

    //----- trapezoidal rule to compute marginal likelihood -----//
    double d_inv_t, output = 0.;
    for (int k = 0; k < num_rungs; ++k)
    {
        d_inv_t = inv_t[k + 1] - inv_t[k];
        output += (avg_vec[k + 1] + avg_vec[k]) * d_inv_t / 2. - (var_vec[k + 1] - var_vec[k]) * std::pow(d_inv_t, 2) / 12.;
    }

    //----- Close progress bar -----//
    tools::progress_bar_end("end", progbar);

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
void MCMC::simulate(bool progbar)
{
    //----- Check no undefined parameters -----//
    pImpl->set_rnd_init_params(false);
    pImpl->init_rand_params();

    //----- Initialize model -----//
    init_model_for_simulation();

    //----- Initialize augmented dependencies -----//
    init_augmented_variables();

    tools::progress_bar_init("simul", progbar);

#ifdef _OPENMP
    int thin_num_house = get_num_house() / pImpl->get_num_threads();
#pragma omp parallel for default(shared)
    for (size_t H = 0; H < get_num_house(); ++H)
    {
        if (_error_occurred_)
        {
            continue;
        }
        try
        {
            simulate_single_house(house(H));
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
            tools::progress_bar_run("simul", H, thin_num_house, progbar);
        }
    }
#else
    for (size_t H = 0; H < get_num_house(); ++H)
    {
        if (_error_occurred_)
        {
            continue;
        }
        try
        {
            simulate_single_house(house(H));
        }
        catch (const std::exception &e)
        {
            _error_occurred_ = true;
            _error_message_ = e.what();
        }
        tools::progress_bar_run("simul", H, get_num_house(), progbar);
    }
#endif

    if (_error_occurred_)
    {
        tools::error(_error_message_);
    }
    else
    {
        tools::progress_bar_end("simul", progbar);
        pImpl->update_report_data();
    }
}