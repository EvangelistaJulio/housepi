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
#include "mcmc_build.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : implementations  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//===== CONSTRUCTOR =====//
MCMC::MCMC() : pImpl(std::make_unique<Impl>()) {}
MCMC::~MCMC() = default;

//===== DATA =====//
Household &MCMC::get_household(ulong H) const { return pImpl->house_vector.at(H); }
ulong MCMC::n_houses() const { return static_cast<ulong>(pImpl->n_households); }

ulong MCMC::get_report_indiv(const std::string &key) const { return pImpl->get_report_indiv(key); }
ulong MCMC::get_report_house(const std::string &key) const { return pImpl->get_report_house(key); }
std::pair<ushort, ulong> MCMC::get_report_hh_size(const std::string &key) const { return pImpl->get_report_hh_size(key); }

//===== RAND NUMBER GENERATOR =====//
double MCMC::uniform() { return pImpl->uniform(); }

//===== PARAMETERS =====//
double MCMC::param(const parID &p) const { return pImpl->get_param(p); }
const double *MCMC::param_pointer(const parID &p) const { return pImpl->get_param_pointer(p); }

//===== AUGMENTED VARIABLES =====//
bool MCMC::is_augm_any() const { return pImpl->is_augm_any(); }
bool MCMC::is_augm_sample() const { return pImpl->is_augm_sample(); }
bool MCMC::is_augm_model() const { return pImpl->is_augm_model(); }
bool MCMC::is_augm_generic() const { return pImpl->is_augm_generic(); }
void MCMC::assign_all_augm(Household &house, indiv i) { pImpl->assign_all_augm(house, i); }
void MCMC::clear_all_augm(Household &house, indiv i) { pImpl->clear_all_augm(house, i); }
double MCMC::compute_augm_model_loglike(const Household &house) const
{
    return pImpl->compute_augm_model_loglike(house);
}

//===== EXTRA ATTRIBUTES =====//

/* Extra integer vectors for all households. */
void MCMC::allocate_extra_house_integer_vector(const extraID &e)
{
    for (ulong H = 0; H < n_houses(); ++H)
        pImpl->house_vector[H].allocate_extra_integer_vector(e);
}
/* Extra double vectors for all households. */
void MCMC::allocate_extra_house_double_vector(const extraID &e)
{
    for (ulong H = 0; H < n_houses(); ++H)
        pImpl->house_vector[H].allocate_extra_double_vector(e);
}
/* Extra double square matrix for all households. */
void MCMC::allocate_extra_house_double_matrix(const extraID &e)
{
    for (ulong H = 0; H < n_houses(); ++H)
        pImpl->house_vector[H].allocate_extra_double_matrix(e);
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
void MCMC::update_extra_dependecies(Household &house, indiv i)
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

/* Build all group rates and contact rates and set rate per group */
void MCMC::build_all_rates()
{
    try
    {
        for (Household &house : pImpl->house_vector)
        {
            for (indiv i = 0; i < house.size(); ++i)
                set_rate_per_group(house, i);
            set_rate_per_contact(house);
        }
    }
    catch (const std::exception &e)
    {
        tools::error({"While building group and contact rates:",
                      std::string(e.what())});
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
                            "Missing rate \"" + r + "\" for individual " + house.id(i) + ".");

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
        for (indiv k = 0; k < house.size(); ++k)     // from
            for (indiv i = 0; i < house.size(); ++i) // to
                // note : from 'k' --> to 'i'
                if (k != i)
                {
                    _pair_contact_rate_ = NULL_STRING;
                    assign_contact_rate_to_household(house, k, i);

                    if (!tools::is_na(_pair_contact_rate_))
                    {
                        if (_pair_contact_rate_ == MISSING_RATE)
                            throw std::invalid_argument(
                                "Missing contact rate for individuals " + house.id(k) + " and " + house.id(i) + " (household " + house.house_id() + ").");

                        else if (_pair_contact_rate_ == REFERENCE)
                            // reference rate
                            house.assign_contact_rate(k, i, &REFERENCE_VALUE);

                        else
                            // parameter rate
                            house.assign_contact_rate(k, i,
                                                      pImpl->get_param_pointer(_pair_contact_rate_));
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
//////////////////////////////  MCMC : update global loglikelihood  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute new loglikelihood, store old value in buffer. Same for households. */
void MCMC::update_global_loglike(const bool extra_depen)
{
    const bool is_augm_ll = is_augm_model();
    double new_ll = 0.0;
    std::atomic<bool> has_error(false);

    // buffer infectivity profile if lookup tables are active
    if (extra_depen)
        pImpl->update_infec_profile();

#pragma omp parallel for reduction(+ : new_ll)
    for (ulong H = 0; H < n_houses(); ++H)
    {
        if (has_error.load(std::memory_order_relaxed))
            continue;

        try
        {
            Household &house = get_household(H);

            // buffer all household extra dependencies
            if (extra_depen)
            {
                house.buffer_all_extra();
                for (const indiv &i : house.infected_status_list())
                    update_extra_dependecies(house, i);
            }

            // update house loglikelihood
            house.set_loglike(compute_house_loglike(house));
            new_ll += house.get_loglike();
            if (is_augm_ll)
                new_ll += compute_augm_model_loglike(house);
        }
        catch (const std::exception &e)
        {
#pragma omp critical
            {
                if (!has_error.exchange(true))
                    _error_message_ = {"While computing loglikelihood for house " + get_household(H).house_id() + ":",
                                       std::string(e.what())};
            }
        }
    }
    if (!_error_message_.empty())
        tools::error(_error_message_);

    pImpl->set_global_loglike(new_ll);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  MCMC : move all augmented variables  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Move augmented variables from all households. */
void MCMC::move_all_augm()
{
    double delta_ll = 0.;
    std::atomic<bool> has_error(false);

#pragma omp parallel for reduction(+ : delta_ll)
    for (ulong H = 0; H < n_houses(); ++H)
    {
        if (has_error.load(std::memory_order_relaxed))
            continue;

        try
        {
            delta_ll += move_augm_single_house(get_household(H));
        }
        catch (const std::exception &e)
        {
#pragma omp critical
            {
                if (!has_error.exchange(true))
                    _error_message_ = {"While moving augmented variables in house " + get_household(H).house_id() + ":",
                                       std::string(e.what())};
            }
        }
    }

    if (!_error_message_.empty())
        tools::error(_error_message_);

    pImpl->increment_global_loglike(delta_ll);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : parameter move  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Attempt to move a single parameter for the given inverse temperature. */
bool MCMC::move_single_parameter(const parID &p, bool extra_depen)
{
    bool is_accepted = false;
    // proceed only if the new value is consistent
    // otherwise automatically execute refuse_move(parID)
    if (pImpl->propose_valid_param_move(p))
    {
        // update global loglikelihood
        // if extra depending on parameter, update and buffer extra
        update_global_loglike(extra_depen);

        // perform Metropolis-Hastings to accept or not move
        if (pImpl->accept_valid_param_move(p))
            is_accepted = true;
        else
            pImpl->restore_parameter(p, extra_depen);
    }

    return is_accepted;
}

/* Move all model parameters. */
void MCMC::move_all()
{
    // move augmented variables
    if (pImpl->is_augm_sample())
        move_all_augm();

    // move all simple model parameters
    for (const parID &p : pImpl->param_standard_infer_list)
        move_single_parameter(p, false);

    // move parameters with extra dependecies
    for (const parID &p : pImpl->param_extra_infer_list)
        move_single_parameter(p, true);
}

/* Move all augmented variables and model parameters and adapt their rule arguments. */
void MCMC::move_all_and_adapt_parameters(ulong step)
{
    // move augmented variables
    if (pImpl->is_augm_sample())
        move_all_augm();

    // move all model parameters and adapt their rule arguments
    double adapt_ratio = static_cast<double>(step) / pImpl->setup_n_adapt;

    // move all simple model parameters
    for (const parID &p : pImpl->param_standard_infer_list)
    {
        pImpl->update_param_accept_rate(p, move_single_parameter(p, false));
        pImpl->adapt_param_rule(p, adapt_ratio);
    }

    // move parameters with extra dependecies
    for (const parID &p : pImpl->param_extra_infer_list)
    {
        pImpl->update_param_accept_rate(p, move_single_parameter(p, true));
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
    if (is_augm_any())
    {
        if (is_augm_sample())
        {
            // wheter augmented variable are generic or status specific
            if (is_augm_generic())
                _func_move_augm_indiv_ = &MCMC::_move_augm_indiv_generic_;
            else
                _func_move_augm_indiv_ = &MCMC::_move_augm_indiv_by_status_;
        }

        // assign augmented variables
        for (Household &house : pImpl->house_vector)
            for (indiv i = 0; i < house.size(); ++i)
                assign_all_augm(house, i);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : augmented move  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Move augmented variables of infectious members inside household.
The move is performed one member at a time, one augmented variable at a time.
*/
double MCMC::move_augm_single_house(Household &house)
{
    double delta_global_loglike = 0.0;

    // Loop over infectious residents of the househould
    for (const indiv &i : house.infected_status_list())
        (this->*_func_move_augm_indiv_)(house, i, delta_global_loglike);

    return delta_global_loglike;
}

/* Generic move augmented variables one by one. */
void MCMC::_move_augm_indiv_generic_(Household &house, indiv i, double &delta)
{
    // Loop over augmented variables of infectious
    for (const augmID &a : pImpl->augm_infer_list)
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
        update_extra_dependecies(house, i);

        // update house loglikelihood and store old value in buffer
        house.set_loglike(compute_house_loglike(house));

        // Metropolis-Hasting algorithm
        double D_ll = house.get_delta_loglike() + house.augm_model_logratio(i, a);
        if (uniform() < exp(pImpl->setup_inv_T * D_ll - house.augm_propose_logratio(i, a)))
            delta += D_ll;
        else
            // restore all buffered values
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
        pImpl->assert_proper_dataset();        // check if dataset is empty
        pImpl->setup_n_chains = 0;             // use 0 chains as a state variable
        pImpl->setup_rand_init_params = false; // disable random initialization // UPDATE
        pImpl->init_parameters();              // check no undefined parameters
        init_infec_profile();                  // initialize user infectivity profile
        init_model_for_inference();            // initialize model for inference
        init_augmented_variables();            // initialize augmented dependencies
        update_global_loglike(true);           // initialize extra + global loglikelihood

        return pImpl->get_global_loglike();
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }

    return std::numeric_limits<double>::quiet_NaN();
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
void MCMC::run_monte_carlo(progress::Bar &progressbar)
{
    //==========================//
    //===== Initialization =====//
    //==========================//

    pImpl->init_diagnostics(); // initialize diagnostics

    // implementation setup aliases
    ulong &c = pImpl->_chain_;
    const ulong &n_chains = pImpl->setup_n_chains;
    const ulong &n_adapt = pImpl->setup_n_adapt;
    const ulong &thin = pImpl->setup_thin;
    const ulong &n_burnin = pImpl->setup_n_burnin;
    const ulong &n_steps = pImpl->setup_n_steps;

    // run multiple MCMC chains&
    for (c = 1; c <= n_chains; ++c)
    {
        // update progressbar text
        progressbar.set_text(pImpl->fmt_progtxt());

        pImpl->init_parameters();    // initialize parameter inference
        init_augmented_variables();  // initialize augmented dependencies
        update_global_loglike(true); // initialize extra + global loglikelihood

        //----- (Finite) Adaptative MCMC -----//
        // initialize the param move counts for adaptive
        pImpl->init_move_counts();
        for (ulong i = 0; i < n_adapt; i += thin, progressbar.update())
            for (ulong t = 1; t <= thin; ++t)
                move_all_and_adapt_parameters(i + t);

        //----- MCMC Burn-in Iterations -----//
        for (ulong i = 0; i < n_burnin; i += thin, progressbar.update())
            for (ulong t = 1; t <= thin; ++t)
                move_all();

        //----- MCMC Posterior Sampling -----//
        for (ulong i = n_burnin; i < n_steps; i += thin, progressbar.update())
        {
            for (ulong t = 1; t <= thin; ++t)
                move_all();

            pImpl->update_diagnostics(); // store sample
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  MCMC : infer  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Model specific initialization (e. g. extra variables, lookup tables, etc). */
void MCMC::init_model_for_inference() {}

/* Performs the inference using Metropolis-Hastings. */
void MCMC::infer()
{
    try
    {
        pImpl->assert_proper_dataset();    // check if dataset is empty
        pImpl->assert_proper_mcmc_setup(); // check setup for mcmc sampling
        init_infec_profile();              // initialize user infectivity profile
        init_model_for_inference();        // user initialization model

        // start progress bar and chrono
        progress::Bar progressbar(pImpl->setup_c_n_sweeps * pImpl->setup_n_chains,
                                  pImpl->fmt_progtxt(), pImpl->setup_progress);

        run_monte_carlo(progressbar); // run Monte Carlo
        progressbar.finalize();       // close progress bar
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

/* Power posteriors numerical integration for evidence.
    -----
    Ref: Friel, N., Hurn, M. & Wyse, J.
    Improving power posterior estimation of statistical evidence.
    Stat Comput 24, 709–723 (2014).
*/
const std::vector<dbl_vec> MCMC::power_posteriors(const std::string &scheme, ulong N)
{
    const str_vec options = {"adaptive", "power-fractions"};

    pImpl->assert_proper_dataset();    // check if dataset is empty
    pImpl->assert_proper_mcmc_setup(); // check setup for mcmc sampling

    // fix internal number of rungs (for progress bar)
    if (N >= 2)
        pImpl->setup_n_rungs = N;
    else
        throw std::invalid_argument("N. of rungs must be greater than 1.");

    // buffer user setup to restore them later
    const bool buffer_rand_init_params = pImpl->setup_rand_init_params;
    const double buffer_inv_T = pImpl->setup_inv_T;

    // implementation setup aliases
    ulong &rung = pImpl->_rung_;    // current theromdynamic ladder rung
    double &t = pImpl->setup_inv_T; // current inverse temperature
    const ulong &n_rungs = pImpl->setup_n_rungs;
    const ulong rung_sweeps = pImpl->setup_c_n_sweeps * pImpl->setup_n_chains;

    // thermodynamic ladder estimates
    std::pair<double, double> stats = {}; // results of stats computation
    dbl_vec itm = {};                     // inverse temperature
    dbl_vec avg = {};                     // mean loglikelihood
    dbl_vec var = {};                     // variance loglikelihood

    // initialize model for inference
    init_model_for_inference();

    // start progress bar and chrono
    progress::Bar progressbar(rung_sweeps * (n_rungs + 1),
                              pImpl->fmt_progtxt(), pImpl->setup_progress);

    //=================================//
    //===== INVERSE TEMPERATURE 0 =====//
    //=================================//

    // UPDATE: should I use previous inference for itm = 1 if available?
    // UPDATE: always under-estimate time if starting from itm = 0

    pImpl->setup_rand_init_params = true;

    rung = 0;
    t = 0.0;
    run_monte_carlo(progressbar);
    stats = pImpl->dgn_glob_ll_stats();

    itm.push_back(t);
    avg.push_back(stats.first);
    var.push_back(stats.second);

    pImpl->setup_rand_init_params = false;

    if (scheme == "adaptive")
    {
        //===========================//
        //===== ADAPTIVE LADDER =====//
        //===========================//

        //----- Inverse Temperature 1 -----//
        rung = 1;
        t = 1.0;
        run_monte_carlo(progressbar);
        stats = pImpl->dgn_glob_ll_stats();

        itm.push_back(t);
        avg.push_back(stats.first);
        var.push_back(stats.second);

        //----- Inverse Temperatures (0, 1) -----//
        size_t k = 0; // index to interval lower bound
        bool is_right = false;
        while (rung < n_rungs)
        {
            ++rung;

            //----- adaptive choice of inverse temperature -----//
            // intersection between tangents
            double new_t = (avg[k + 1] - itm[k + 1] * var[k + 1] - avg[k] + itm[k] * var[k]) / (var[k] - var[k + 1]);
            // otherwise weighted average
            if ((new_t <= itm[k]) || (new_t >= itm[k + 1]))
                new_t = itm[k] + var[k + 1] * (itm[k + 1] - itm[k]) / (var[k] + var[k + 1]);

            //----- power posterior sampling -----//
            t = new_t;
            run_monte_carlo(progressbar);
            stats = pImpl->dgn_glob_ll_stats();

            //----- find larger rectangular discretization error -----//
            is_right = (itm[k + 1] - t) * (avg[k + 1] - stats.first) > (new_t - itm[k]) * (stats.first - avg[k]);

            //----- insert in-between for all lists -----//
            itm.insert(itm.begin() + k + 1, t);
            avg.insert(avg.begin() + k + 1, stats.first);
            var.insert(var.begin() + k + 1, stats.second);

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
        for (rung = 1; rung <= n_rungs; ++rung)
        {
            // inverse tempearture = (rung / #rungs)^5
            t = std::pow(static_cast<double>(rung) / n_rungs, 5);
            run_monte_carlo(progressbar);
            stats = pImpl->dgn_glob_ll_stats();

            itm.push_back(t);
            avg.push_back(stats.first);
            var.push_back(stats.second);
        }
    }
    else
        throw std::invalid_argument(
            "Unknown power posterior scheme \"" + scheme + "\"." + tools::find_closest(scheme, options, 3));

    //----- Close progress bar -----//
    progressbar.finalize();

    //----- Restore user setup -----//
    t = buffer_inv_T;
    pImpl->setup_rand_init_params = buffer_rand_init_params;

    return {itm, avg, var};
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

/* Auxiliary method to initialize simulations. */
void MCMC::_aux_init_simulate_()
{
    pImpl->assert_proper_dataset(); // check if dataset is empty
    pImpl->setup_n_chains = 0;      // use 0 chains as a state variable

    pImpl->setup_n_threads = 1; // impose no parallelization
    pImpl->set_n_threads();

    pImpl->setup_rand_init_params = false; // disable random parameter initialization
    pImpl->init_parameters();              // check no undefined parameters
    init_infec_profile();                  // initialize user infectivity profile
    init_model_for_simulation();           // initialize model for simulation
}

#ifdef _HOUSEPI_RCPP
/* Perform simulation and return simulated dataset. */
Rcpp::DataFrame MCMC::simulate()
{
    using namespace Rcpp;

    try
    {
        size_t n_cols = pImpl->n_individuals;
        size_t row = 0;
        str_vec _aux_data_headers = pImpl->data_headers;
        Rcpp::CharacterVector df_index = CharacterVector(n_cols);
        Rcpp::CharacterVector df_house_index = CharacterVector(n_cols);
        Rcpp::IntegerVector df_status = IntegerVector(n_cols);
        std::vector<Rcpp::NumericVector> df_dbl = {};

        // UPDATE: why augm_list_full is still a set?
        // remove data columns if they have same name than augmented
        str_vec &v = _aux_data_headers;
        const std::unordered_set<augmID> &y = pImpl->augm_list_full;
        auto new_end = std::remove_if(v.begin(), v.end(),
                                      [&y](const dataID &elem)
                                      {
                                          return y.find(elem) != y.end();
                                      });
        v.erase(new_end, v.end());

        // other data and augmented variable columns
        for (size_t s = 0; s < _aux_data_headers.size() + pImpl->augm_list_full.size(); ++s)
            df_dbl.emplace_back(n_cols);

        //----- initialization -----//
        _aux_init_simulate_();
        progress::Bar progressbar(n_houses(), "simulate",
                                  pImpl->setup_progress);

        //----- run simulation -----//
        for (ulong H = 0; H < n_houses(); ++H, progressbar.update())
            try
            {
                // deep copy of the household
                Household house = get_household(H);
                // inplace simlation
                simulate_single_house(house);

                // store simulation
                for (indiv i = 0; i < house.size(); ++i)
                {
                    df_index[row] = house.id(i);
                    df_house_index[row] = house.house_id();
                    df_status[row] = house.status(i);

                    size_t col = 0;
                    for (const dataID &d : _aux_data_headers)
                        df_dbl[col++][row] = house.data(i, d);

                    for (const augmID &a : pImpl->augm_list_full)
                        df_dbl[col++][row] = house.augm(i, a);

                    ++row;
                }
            }
            catch (const std::exception &e)
            {
                tools::error({"While simulating house " + get_household(H).house_id() + " :",
                              e.what()});
            }

        // close progress bar
        progressbar.finalize();

        // assign vectors to dataframe
        DataFrame output = DataFrame::create(
            _["id"] = std::move(df_index),
            _["house_id"] = std::move(df_house_index),
            _["status"] = std::move(df_status));

        size_t col = 0;
        for (const dataID &d : _aux_data_headers)
            output.push_back(std::move(df_dbl[col++]), d);

        for (const augmID &a : pImpl->augm_list_full)
            output.push_back(std::move(df_dbl[col++]), a);

        return output;
    }
    catch (std::exception &e)
    {
        tools::error(e.what());
    }

    return DataFrame::create();
}

Rcpp::List MCMC::estimate_evidence(const std::string &scheme,
                                   ulong n_rungs)
{
    using namespace Rcpp;

    try
    {
        double output = 0;
        const std::vector<dbl_vec> pp = power_posteriors(scheme, n_rungs);
        const dbl_vec &inv_tem = pp[0];
        const dbl_vec &avg_vec = pp[1];
        const dbl_vec &var_vec = pp[2];

        for (int k = 0; k < n_rungs; ++k)
        {
            double d_inv_tem = inv_tem[k + 1] - inv_tem[k];
            output += (avg_vec[k + 1] + avg_vec[k]) * d_inv_tem / 2.;
            output += (var_vec[k] - var_vec[k + 1]) * d_inv_tem * d_inv_tem / 12.;
        }

        return List::create(
            _["estimate"] = output,
            _["scheme"] = scheme,
            _["n_rungs"] = n_rungs,
            _["inv_temper"] = wrap(inv_tem),
            _["mean_loglike"] = wrap(avg_vec),
            _["var_loglike"] = wrap(var_vec));
    }
    catch (const std::exception &e)
    {
        tools::error({"While estimating evidence with power posteriors:",
                      std::string(e.what())});
    }
    return List();
}

#else

/* Print the simulated dataset to the provided file. */
void MCMC::print_simulated_dataset(const std::string &pathtofile)
{
    pImpl->_aux_print_dataset_init_(pathtofile);

    try
    {
        //----- initialization -----//
        _aux_init_simulate_();
        progress::Bar progressbar(n_houses(), "simulate",
                                  pImpl->setup_progress);

        //----- run simulation -----//
        for (ulong H = 0; H < n_houses(); ++H, progressbar.update())
            try
            {
                Household house = get_household(H);       // deep copy of the household
                simulate_single_house(house);             // inplace simulation
                pImpl->_aux_print_dataset_update_(house); // print the simulated dataset
            }
            catch (const std::exception &e)
            {
                tools::error({"While simulating house " + get_household(H).house_id() + " :",
                              e.what()});
            }

        // close progress bar
        progressbar.finalize();

        // finalize output
        pImpl->_aux_print_dataset_finalize_();
    }
    catch (std::exception &e)
    {
        tools::error(e.what());
    }
}

/* Trapezoidal rule to estimate logarithm of marginal likelihood (evidence) with power posteriors.
    -----
    Ref: Friel, N., Hurn, M. & Wyse, J.
    Improving power posterior estimation of statistical evidence.
    Stat Comput 24, 709–723 (2014).
*/
void MCMC::print_evidence(const std::string &scheme, ulong n_rungs)
{
    try
    {
        double output = 0;
        const std::vector<dbl_vec> pp = power_posteriors(scheme, n_rungs);
        const dbl_vec &inv_tem = pp[0];
        const dbl_vec &avg_vec = pp[1];
        const dbl_vec &var_vec = pp[2];

        for (ulong k = 0; k < n_rungs; ++k)
        {
            double d_inv_tem = inv_tem[k + 1] - inv_tem[k];
            output += (avg_vec[k + 1] + avg_vec[k]) * d_inv_tem / 2.;
            output += (var_vec[k] - var_vec[k + 1]) * d_inv_tem * d_inv_tem / 12.;
        }

        // print result to stdout
        tools::print_stdout(std::to_string(output) + "\n");
    }
    catch (const std::exception &e)
    {
        tools::error({"While estimating evidence with power posteriors:",
                      std::string(e.what())});
    }
}

#endif
