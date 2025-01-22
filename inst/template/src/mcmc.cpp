/*
 * mcmc.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#include "mcmc_build.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : implementations  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//===== CONSTRUCTOR =====//
MCMC::MCMC() : pImpl(std::make_unique<Impl>()) {}
MCMC::~MCMC() = default;

//===== READ =====//
void MCMC::read_setup(const std::string &pathtofile) { pImpl->read_setup(pathtofile); }
void MCMC::read_params(const std::string &pathtofile) { pImpl->read_params(pathtofile); }
void MCMC::read_augm(const std::string &pathtofile) { pImpl->read_augm(pathtofile); }

//===== C++ PRINT =====//
void MCMC::print_infer_state(const std::string &pathtofile) const { pImpl->print_infer_state(pathtofile); }
void MCMC::print_augm_variables(const std::string &pathtofile) const { pImpl->print_augm_variables(pathtofile); }

//===== SHOW =====//
void MCMC::show_full_build() const { pImpl->show_full_build(); }
void MCMC::show_results() const { pImpl->show_results(); }
void MCMC::show_summary() const { pImpl->show_summary(); }

//===== DATA =====//
size_t MCMC::n_houses() const { return pImpl->n_households; }
size_t MCMC::report_data(const std::string &key) const { return pImpl->report_data(key); }
double MCMC::uniform() { return pImpl->uniform(); }

//===== AUGMENTED =====//
double MCMC::param(const parID &p) const { return pImpl->get_param(p); }
const double *MCMC::param_pointer(const parID &p) const { return pImpl->get_param_pointer(p); }
void MCMC::assign_all_augm(Household &house, indiv i) { pImpl->assign_all_augm(house, i); }
void MCMC::clear_all_augm(Household &house, indiv i) { pImpl->clear_all_augm(house, i); }

//===== EXTRA ATTRIBUTES =====//

/* Extra single integers for all households. */
void MCMC::assign_extra_integer_to_all_houses(const std::vector<extraID> &names, int xinit)
{
    for (size_t H = 0; H < n_houses(); ++H)
        pImpl->house_vector[H].allocate_extra_integer(names, xinit);
}
/* Extra single doubles for all households. */
void MCMC::assign_extra_double_to_all_houses(const std::vector<extraID> &names, double xinit)
{
    for (size_t H = 0; H < n_houses(); ++H)
        pImpl->house_vector[H].allocate_extra_double(names, xinit);
}
/* Extra double vectors for all households. */
void MCMC::assign_extra_vector_to_all_houses(const std::vector<extraID> &names, double xinit)
{
    for (size_t H = 0; H < n_houses(); ++H)
        pImpl->house_vector[H].allocate_extra_vector(names, xinit);
}
/* Extra double square matrix for all households. */
void MCMC::assign_extra_matrix_to_all_houses(const std::vector<extraID> &names, double xinit)
{
    for (size_t H = 0; H < n_houses(); ++H)
        pImpl->house_vector[H].allocate_extra_matrix(names, xinit);
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
    UNUSED(house);
    UNUSED(i);
}

/* Assign the contact patterns between a couple of individuals. */
void MCMC::assign_contact_rate_to_household(Household &house, indiv from, indiv to)
{
    UNUSED(house);
    UNUSED(from);
    UNUSED(to);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  MCMC : build data  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read data file and then run assignment. */
void MCMC::read_data(const std::string &pathtofile)
{
    dict_str_vec datatable = {};
    tools::read_csv_table(pathtofile, datatable);
    build_data(datatable);

    pImpl->_report_pathtofile_["data"] = pathtofile;
}

/* Assign data after pre-processing. */
void MCMC::build_data(const dict_str_vec &datatable)
{
    try
    {
        if (datatable.empty())
        {
            // If nothing was previously defined
            if (pImpl->house_vector.empty())
                throw std::length_error(
                    "No database was found.");
        }
        else
        {
            //=====  Drop all Previous Data and Prepare for New  =====//
            pImpl->data_headers.clear();
            for (const auto &x : datatable)
                pImpl->data_headers.push_back(x.first);
            pImpl->n_individuals = datatable.at(pImpl->data_headers[0]).size();
            pImpl->n_households = 0;
            pImpl->house_vector.clear();

            //===== Assess the presence of house IDs =====//
            if (std::find(pImpl->data_headers.begin(), pImpl->data_headers.end(), "house_id") == pImpl->data_headers.end())
                tools::warning(
                    "Unknown household structure: missing data column \"house_id\".");

            //===== Safe Load of Households =====//
            std::unordered_map<std::string, size_t> hh_ID_to_index = {};
            for (size_t i = 0; i < pImpl->n_individuals; ++i)
            {
                //----- Read individual data -----//
                std::string ID = std::to_string(i);                // default value
                std::string house_ID = std::to_string(n_houses()); // default value
                int status = MISSING_STATUS;
                dataIDmap<double> data = {};

                for (const std::string &key : pImpl->data_headers)
                {
                    if (key == "id")
                    {
                        if (tools::is_na(datatable.at("id")[i]))
                            throw std::invalid_argument(
                                "Missing id for individual in row n." + std::to_string(i + 2) + ".");
                        else
                            ID = datatable.at("id")[i];
                    }
                    else if (key == "house_id")
                    {
                        if (tools::is_na(datatable.at("house_id")[i]))
                            throw std::invalid_argument(
                                "Missing house_id for individual in row n." + std::to_string(i + 2) + ".");
                        else
                            house_ID = datatable.at("house_id")[i];
                    }
                    else if (key == "status")
                    {
                        // NA status is a missing status
                        if (!tools::is_na(datatable.at("status")[i]))
                            status = tools::stoul(datatable.at("status")[i], "individual status");
                    }
                    else if (tools::is_na(datatable.at(key)[i]))
                    {
                        data.emplace(key, std::numeric_limits<double>::quiet_NaN());
                        tools::warning("Found NA value(s) in the database.");
                    }
                    else
                    {
                        data.emplace(key, tools::stod(datatable.at(key)[i], "data \"" + key + "\""));
                    }
                }

                // if household index doesn't exist, add new household
                // note : this is redundant if no structure
                if (hh_ID_to_index.find(house_ID) == hh_ID_to_index.end())
                {
                    // map data index to household number
                    hh_ID_to_index.emplace(house_ID, n_houses());
                    // create a new empty household
                    pImpl->house_vector.emplace_back(house_ID, pImpl->n_households++);
                }

                // extract the household reference
                Household &house = pImpl->house_vector[hh_ID_to_index.at(house_ID)];

                // Assign the individual the respective household
                house.assign_individual(ID, status, data);

                // Assign and set rate per group
                set_rate_per_group(house, house.size() - 1);
            }

            //===== Set Contact Rates =====//
            pImpl->_n_contacts_ = 0;
            for (Household &house : pImpl->house_vector)
                set_rate_per_contact(house);

            //===== Update Report for Data =====//
            pImpl->update_report_data();

            //===== Clean up data entries for output files =====//
            // eventually remove standard headers to reorder them at the end
            for (const std::string &col : {"house_id", "id", "status"})
                pImpl->data_headers.erase(
                    std::remove(pImpl->data_headers.begin(), pImpl->data_headers.end(), col),
                    pImpl->data_headers.end());
        }
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : set rate per group  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Auxiliary wrapper to emplace_back */
void MCMC::emplace_group_rate(rateID r, parID p)
{
    _individual_group_rates_[r].emplace_back(p);
}

/* Update model report and set group rates associate to new individual. */
void MCMC::set_rate_per_group(Household &house, indiv i)
{
    try
    {
        rateIDmapvec<const double *> group_rates_ptrs = {};

        // call user defined group rate assignmen
        _individual_group_rates_.clear();
        assign_group_rate_to_individual(house, i);

        // convert to dict_rates_ptr
        for (const auto &pair : _individual_group_rates_)
        {
            rateID r = pair.first;

            if (r == CONTACT_RATE)
                throw std::invalid_argument(
                    "The group rate name \"" + CONTACT_RATE + "\" is reserved. Choose a different one.");
            else if (!tools::is_na(r))
            {
                // return parameter pointers
                bool is_all_reference = true; // all reference state
                std::vector<const double *> vec_par_ptr = {};

                for (const parID &p : pair.second)
                {
                    if ((tools::is_na(p)) || (p == MISSING_RATE))
                        std::invalid_argument(
                            "Missing rate \"" + r + "\" for individual " + house.ID(i) + ".");
                    else if (p == REFERENCE)
                        vec_par_ptr.emplace_back(&REFERENCE_VALUE);
                    else
                    {
                        is_all_reference = false;
                        vec_par_ptr.emplace_back(pImpl->get_param_pointer(p));
                    }
                }

                // simplify for the combined reference values
                if (is_all_reference)
                {
                    vec_par_ptr.clear();
                    vec_par_ptr = {&REFERENCE_VALUE};
                }

                // emplace
                group_rates_ptrs[r] = vec_par_ptr;
            }
        }

        // pass the pointers to the individual
        house.assign_all_group_rates(i, group_rates_ptrs);
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : set rate per contact  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Auxiliary wrapper to emplace_back */
void MCMC::emplace_contact_rate(parID p) { _pair_contact_rate_ = p; }

/* Set contact rate according to user assignement. */
void MCMC::set_rate_per_contact(Household &house)
{
    try
    {
        // loop over all couples
        for (indiv k = 0; k < house.size(); ++k) // from
        {
            for (indiv i = 0; i < house.size(); ++i) // to
            {
                // note : from 'k' --> to 'i'
                if (k != i)
                {
                    _pair_contact_rate_ = NULL_STRING;
                    assign_contact_rate_to_household(house, k, i);

                    if (!tools::is_na(_pair_contact_rate_))
                    {
                        if (_pair_contact_rate_ == MISSING_RATE)
                            throw std::invalid_argument(
                                "Missing contact rate for individuals " + house.ID(k) + " and " + house.ID(i) + " (household " + house.house_ID() + ").");
                        else if (_pair_contact_rate_ == REFERENCE)
                            house.assign_contact_rate(k, i, &REFERENCE_VALUE); // reference rate
                        else
                        {
                            house.assign_contact_rate(
                                k, i,
                                pImpl->get_param_pointer(_pair_contact_rate_) // parameter rate
                            );
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
    double new_ll = 0.0;

#pragma omp parallel for reduction(+ : new_ll)
    for (size_t H = 0; H < n_houses(); ++H)
    {

        if (!_error_message_.empty())
            continue;

        try
        {
            Household &house = pImpl->house_vector[H];

            // update house loglikelihood
            house.set_loglike(compute_house_loglike(house));
            new_ll += house.get_loglike();
            if (pImpl->is_augm_model())
                new_ll += pImpl->compute_augm_model_loglike(house);
        }
        catch (const std::exception &e)
        {
#pragma omp critical
            {
                if (_error_message_.empty())
                    _error_message_ = e.what();
            }
        }
    }
    if (!_error_message_.empty())
        tools::error(_error_message_);

    pImpl->set_global_loglike(new_ll);
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
void MCMC::move_single_param(const parID &p)
{
    bool is_accepted = false;
    // proceed only if the new value is consistent, otherwise automatically execute refuse_move(parID)
    if (pImpl->propose_valid_param_move(p))
    {
        update_global_loglike();

        if (uniform() < exp(pImpl->setup_inv_T * pImpl->get_delta_global_loglike() + pImpl->param_move_logratio(p)))
            is_accepted = true;
        else
        {
            // restore buffered value of parameter
            pImpl->restore_param_value(p);

            // restore buffered value of global loglike
            pImpl->restore_global_loglike();

            // restore buffered loglikelihood in households
            for (size_t H = 0; H < n_houses(); ++H)
                pImpl->house_vector[H].restore_loglike();
        }
    }

    // update counts for acceptance rate
    pImpl->update_param_accept_rate(p, is_accepted);
}

/* Move all model parameters. */
void MCMC::move_all()
{
    // move all augmented variables
    (this->*_func_move_all_augm_)();

    // move all model parameters
    for (const parID &p : pImpl->param_to_be_inferred)
        move_single_param(p);
}

/* Move all augmented variables and model parameters and adapt their rule arguments. */
void MCMC::move_all_and_adapt_parameters(ulong step)
{
    // move all augmented variables
    (this->*_func_move_all_augm_)();

    // move all model parameters and adapt their rule arguments
    double adapt_ratio = static_cast<double>(step) / pImpl->setup_n_adapt;

    for (const parID &p : pImpl->param_to_be_inferred)
    {
        move_single_param(p);
        pImpl->adapt_param_rule(p, adapt_ratio);
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : init augmented variables  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize augmented dependecies for all households. */
void MCMC::init_augmented_variables()
{
    // wheter to perform data augmentation or not
    if (pImpl->is_augm_infer())
    {
        _func_move_all_augm_ = &MCMC::_move_all_augm_;

        // wheter augmented variable are generic or status specific
        if (pImpl->is_augm_generic())
            _func_move_augm_indiv_ = &MCMC::_move_augm_indiv_generic_;
        else
            _func_move_augm_indiv_ = &MCMC::_move_augm_indiv_by_status_;

        // assign augmented variables
        for (size_t H = 0; H < n_houses(); ++H)
        {
            // assign augmented to preallocated variables
            Household &house = pImpl->house_vector[H];
            for (indiv i = 0; i < house.size(); ++i)
                assign_all_augm(house, i);

            for (const indiv &i : house.infectious_set())
                update_augm_dependecies(house, i);
        }
    }
    else
    {
        _func_move_all_augm_ = &MCMC::_move_no_augm_;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : augmented move  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Don't move anything. */
void MCMC::_move_no_augm_() {}

/* Move augmented variables from all households. */
void MCMC::_move_all_augm_()
{
    double delta_ll = 0.;

#pragma omp parallel for reduction(+ : delta_ll)
    for (size_t H = 0; H < n_houses(); ++H)
    {
        if (!_error_message_.empty())
            continue;

        try
        {
            delta_ll += move_augm_single_house(pImpl->house_vector[H]);
        }
        catch (const std::exception &e)
        {
#pragma omp critical
            {
                if (_error_message_.empty())
                    _error_message_ = e.what();
            }
        }
    }

    if (!_error_message_.empty())
        tools::error(_error_message_);

    pImpl->increment_global_loglike(delta_ll);
}

/* Move augmented variables for infectious members of a single household.
The move is performed one member at a time, one augmented variable at a time.
*/
double MCMC::move_augm_single_house(Household &house)
{
    double delta_global_loglike = 0.0;

    // Loop over infectious residents of the househould
    for (const indiv &i : house.infectious_set())
        (this->*_func_move_augm_indiv_)(house, i, delta_global_loglike);

    return delta_global_loglike;
}

/* Generic move augmented variables one by one. */
void MCMC::_move_augm_indiv_generic_(Household &house, indiv i, double &delta)
{
    // Loop over augmented variables of infectious
    for (const augmID &a : pImpl->augm_list_to_be_inferred)
        metropolis_augm_single_move(house, i, a, delta);
}

/* Move augmented variables by status one variable at a time. */
void MCMC::_move_augm_indiv_by_status_(Household &house, indiv i, double &delta)
{
    // Loop over augmented variables of infectious
    for (const augmID &a : pImpl->augm_status_map_infer.at(house.status(i)))
        metropolis_augm_single_move(house, i, a, delta);
}

/* Take care of a a single augmented Metropolis move. */
void MCMC::metropolis_augm_single_move(Household &house, indiv i, const augmID &a, double &delta)
{
    // proceed only if the new value is consistent
    if (house.propose_valid_augm_move(i, a, pImpl->get_prng()))
    {
        // update dependencies to perform computations
        house.buffer_all_extra();
        update_augm_dependecies(house, i);

        // update house loglikelihood and store old value in buffer
        house.set_loglike(compute_house_loglike(house));

        // Metropolis-Hasting algorithm
        double D_ll = house.get_delta_loglike() + house.augm_model_logratio(i, a);
        if (uniform() < exp(pImpl->setup_inv_T * D_ll - house.augm_propose_logratio(i, a)))
            delta += D_ll;
        else
            // restore all buffer values
            house.restore_augm_value(i, a);
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
        init_model_for_inference(); // initialize model for inference
        pImpl->init_parameters();  // check no undefined parameters
        init_augmented_variables(); // initialize augmented dependencies
        update_global_loglike();    // initialize global loglikelihood
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }

    return pImpl->get_global_loglike();
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

/* Markov chain core, used for inference and power posteriors.
    note: sweep is the number of previously completed thin sweeps. */
void MCMC::run_monte_carlo(CustomProgressBar &progressbar, ulong sweep)
{
    //==========================//
    //===== Initialization =====//
    //==========================//

    // check if there is something to infer
    if ((pImpl->param_to_be_inferred.empty()) && (!pImpl->is_augm_infer()))
        throw std::length_error(
            "Nothing to infer. Check your parameter and augmented variable definitions.");

    // initialize diagnostics
    pImpl->init_diagnostics();

    // run Multiple MCMC Chains
    for (pImpl->_chain_ = 1; pImpl->_chain_ <= pImpl->setup_n_chains; ++pImpl->_chain_)
    {
        // update progressbar text
        progressbar.set_text(pImpl->fmt_progtxt("infer"));

        // initialize parameter inference
        pImpl->init_parameters();

        // initialize augmented dependencies
        init_augmented_variables();

        // initialize global loglikelihood
        update_global_loglike();

        //=========================================//
        //======== (Finite) Adaptative MCMC =======//
        //=========================================//

        // initialize the param move counts for adaptive
        pImpl->init_move_counts();

        for (ulong i = 0; i < pImpl->setup_n_adapt; i += pImpl->setup_thin)
        {
            for (ulong t = 1; t <= pImpl->setup_thin; ++t)
                move_all_and_adapt_parameters(i + t);

            progressbar.update(static_cast<double>(++sweep) / pImpl->_tot_thin_sweeps_);
        }

        //============================================//
        //========= MCMC Burn-in Iterations ==========//
        //============================================//

        for (ulong i = 0; i < pImpl->setup_n_burnin; i += pImpl->setup_thin)
        {
            for (ulong t = 1; t <= pImpl->setup_thin; ++t)
                move_all();

            progressbar.update(static_cast<double>(++sweep) / pImpl->_tot_thin_sweeps_);
        }

        //============================================//
        //========= MCMC Posterior Sampling ==========//
        //============================================//

        for (ulong i = pImpl->setup_n_burnin; i < pImpl->setup_n_steps; i += pImpl->setup_thin)
        {
            for (ulong t = 1; t <= pImpl->setup_thin; ++t)
                move_all();

            // compute diagnostics for the previous iterations
            pImpl->update_diagnostics();

            progressbar.update(static_cast<double>(++sweep) / pImpl->_tot_thin_sweeps_);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  MCMC : infer  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Model specific initialization (e. g. extra variables, lookup tables, etc). */
void MCMC::init_model_for_inference() {}

/* Performs the inference using Metropolis-Hastings. */
void MCMC::infer(bool is_progbar)
{
    try
    {
        // total number of mcmc sweeps
        pImpl->_chain_ = 1;
        pImpl->_tot_thin_sweeps_ = pImpl->setup_c_n_sweeps * pImpl->setup_n_chains;

        //----- Initialize model for inference -----//
        init_model_for_inference();

        // start progress bar and chrono
        CustomProgressBar progressbar(pImpl->fmt_progtxt("infer"), is_progbar);

        // run Markov Chain
        run_monte_carlo(progressbar);

        // close progress bar
        progressbar.finalize();
    }
    catch (std::exception &e)
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
/////////////////////////////////////  MCMC : power posteriors  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Evaluate the logarithm of marginal likelihood with power posteriors.
    -----
    Ref: Friel, N., Hurn, M. & Wyse, J.
    Improving power posterior estimation of statistical evidence.
    Stat Comput 24, 709–723 (2014).
*/
double MCMC::evaluate_marginal(int n_rungs, std::string scheme, bool is_progbar)
{
    const str_vec options = {"adaptive", "power-fractions"};

    //----- fix internal number of rungs (for progress bar) -----//
    if (n_rungs >= 2)
        pImpl->_n_rungs_ = n_rungs;
    else
        tools::error("Number of rungs must be greater than 1.");

    //----- Backup user setup to restore them later -----//
    bool bkp_rand_init_params = pImpl->setup_rand_init_params;
    double bkp_inv_T = pImpl->setup_inv_T;

    //----- Thermodynamic ladder estimates -----//
    bool is_right = false;
    double avg, var;
    dbl_vec inv_t = {};   // inverse temperature
    dbl_vec avg_vec = {}; // mean loglikelihood
    dbl_vec var_vec = {}; // variance loglikelihood

    //----- Total number of mcmc sweeps -----//
    ulong rung_sweeps = pImpl->setup_c_n_sweeps * pImpl->setup_n_chains;
    // internal state variables
    pImpl->_tot_thin_sweeps_ = rung_sweeps * (pImpl->_n_rungs_ + 1);
    pImpl->_rung_ = 0;

    //----- Initialize model for inference -----//
    init_model_for_inference();

    //----- Start progress bar and chrono -----//
    CustomProgressBar progressbar(pImpl->fmt_progtxt(""), is_progbar);

    //==================================//
    //===== INVERSE TEMPEREATURE 0 =====//
    //==================================//

    pImpl->set_inver_temper(0.0);
    pImpl->setup_rand_init_params = true;
    run_monte_carlo(progressbar, 0);
    pImpl->dgn_glob_ll_stats(avg, var);

    inv_t.push_back(0.0);
    avg_vec.push_back(avg);
    var_vec.push_back(var);

    if (scheme == "adaptive")
    {
        //===========================//
        //===== ADAPTIVE LADDER =====//
        //===========================//

        //----- Inverse Temperature 1 -----//
        pImpl->_rung_ = 1;
        pImpl->set_inver_temper(1.0);
        pImpl->setup_rand_init_params = false;
        run_monte_carlo(progressbar, rung_sweeps);
        pImpl->dgn_glob_ll_stats(avg, var);

        inv_t.push_back(1.0);
        avg_vec.push_back(avg);
        var_vec.push_back(var);

        //----- Inverse Temperatures (0, 1) -----//
        size_t k = 0; // index to interval lower bound
        while (pImpl->_rung_ < pImpl->_n_rungs_)
        {
            ++pImpl->_rung_;

            //----- adaptive choice of inverse temperature -----//
            // intersection between tangents
            double new_inv_T = (avg_vec[k + 1] - inv_t[k + 1] * var_vec[k + 1] - avg_vec[k] + inv_t[k] * var_vec[k]) / (var_vec[k] - var_vec[k + 1]);

            // otherwise weighted average
            if ((new_inv_T <= inv_t[k]) || (new_inv_T >= inv_t[k + 1]))
                new_inv_T = inv_t[k] + var_vec[k + 1] * (inv_t[k + 1] - inv_t[k]) / (var_vec[k] + var_vec[k + 1]);

            //----- power posterior sampling -----//
            pImpl->set_inver_temper(new_inv_T);
            run_monte_carlo(progressbar, pImpl->_rung_ * rung_sweeps);
            pImpl->dgn_glob_ll_stats(avg, var);

            //----- find larger rectangular discretization error -----//
            is_right = (inv_t[k + 1] - pImpl->setup_inv_T) * (avg_vec[k + 1] - avg) > (new_inv_T - inv_t[k]) * (avg - avg_vec[k]);

            //----- insert in-between for all lists -----//
            inv_t.insert(inv_t.begin() + k + 1, pImpl->setup_inv_T);
            avg_vec.insert(avg_vec.begin() + k + 1, avg);
            var_vec.insert(var_vec.begin() + k + 1, var);

            // if right interval is chosen, index points to inserted as lower bound
            k += (is_right);
        }
    }
    else if (scheme == "power-fractions")
    {
        //==================================//
        //===== POWER FRACTIONS LADDER =====//
        //==================================//

        //----- Inverse Temperatures (0, 1] -----//
        pImpl->setup_rand_init_params = false;
        for (pImpl->_rung_ = 1; pImpl->_rung_ < pImpl->_n_rungs_; ++pImpl->_rung_)
        {
            pImpl->set_inver_temper(std::pow(static_cast<double>(pImpl->_rung_) / (pImpl->_n_rungs_ - 1), 5));
            run_monte_carlo(progressbar, pImpl->_rung_ * rung_sweeps);
            pImpl->dgn_glob_ll_stats(avg, var);

            inv_t.push_back(pImpl->setup_inv_T);
            avg_vec.push_back(avg);
            var_vec.push_back(var);
        }
    }
    else
        tools::error(
            "Unknown power posterior scheme \"" + scheme + "\"." + tools::find_closest(scheme, options, 3));

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
    pImpl->set_inver_temper(bkp_inv_T);
    pImpl->setup_rand_init_params = bkp_rand_init_params;

    return output;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : simulation  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Model specific initialization (e. g. extra variables, lookup tables, etc). */
void MCMC::init_model_for_simulation() {}

/* Model specific simulation for a single household. */
void MCMC::simulate_single_house(Household &house) { UNUSED(house); }

/* Perform simulation. */
void MCMC::simulate(bool is_progbar)
{
    //----- Initialize model for simulation -----//
    init_model_for_simulation();

    //----- Check no undefined parameters -----//
    pImpl->setup_rand_init_params = false;
    pImpl->init_parameters();

    //----- Initialize augmented dependencies -----//
    init_augmented_variables();
    CustomProgressBar progressbar("simulate", is_progbar);

#ifdef _OPENMP
    int thin_n_house = n_houses() / pBuild->setup_n_threads;
#pragma omp parallel for default(shared)
    for (size_t H = 0; H < n_houses(); ++H)
    {
        if (!_error_message_.empty())
            continue;

        try
        {
            simulate_single_house(pBuild->house_vector[H]);
        }
        catch (const std::exception &e)
        {
#pragma omp critical
            {
                if (_error_message_.empty())
                    _error_message_ = e.what();
            }
        }

        if (omp_get_thread_num() == 0)
            progressbar.update(static_cast<double>(H) / thin_n_house);
    }
#else
    for (size_t H = 0; H < n_houses(); ++H)
    {
        if (!_error_message_.empty())
            continue;

        try
        {
            simulate_single_house(pImpl->house_vector[H]);
        }
        catch (const std::exception &e)
        {
            _error_message_ = e.what();
        }
        progressbar.update(static_cast<double>(H) / n_houses());
    }
#endif

    if (!_error_message_.empty())
        tools::error(_error_message_);

    progressbar.finalize();
    pImpl->update_report_data();
}
