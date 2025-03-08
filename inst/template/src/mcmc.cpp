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

//===== DATA =====//
Household &MCMC::get_household(size_t H) const { return pImpl->house_vector.at(H); }
size_t MCMC::n_houses() const { return pImpl->n_households; }

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

/* Initialize extra dependencies. */
void MCMC::init_extra_dependecies()
{
    if (pImpl->house_vector[0].is_any_extra())
        for (Household &house : pImpl->house_vector)
            for (const indiv &i : house.infected_status_list())
                update_extra_dependecies(house, i);
}

/* Extra integer vectors for all households. */
void MCMC::allocate_extra_house_integer_vector(const extraID &e)
{
    for (size_t H = 0; H < n_houses(); ++H)
        pImpl->house_vector[H].allocate_extra_integer_vector(e);
}
/* Extra double vectors for all households. */
void MCMC::allocate_extra_house_double_vector(const extraID &e)
{
    for (size_t H = 0; H < n_houses(); ++H)
        pImpl->house_vector[H].allocate_extra_double_vector(e);
}
/* Extra double square matrix for all households. */
void MCMC::allocate_extra_house_double_matrix(const extraID &e)
{
    for (size_t H = 0; H < n_houses(); ++H)
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
        tools::error(
            "While building group and contact rates: " + std::string(e.what()));
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
///////////////////////////////////  MCMC : methods (general)  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

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
bool MCMC::move_independent_param(const parID &p)
{
    bool is_accepted = false;
    // proceed only if the new value is consistent, otherwise automatically execute refuse_move(parID)
    if (pImpl->propose_valid_param_move(p))
    {
        update_simple_global_loglike();

        if (pImpl->accept_valid_param_move(p))
            is_accepted = true;
        else
            pImpl->restore_standard_param_value(p);
    }

    return is_accepted;
}

/* Attempt to move parameter with extra dependencies. */
bool MCMC::move_extra_depen_param(const parID &p)
{
    bool is_accepted = false;

    // proceed only if the new value is consistent, otherwise automatically execute refuse_move(parID)
    if (pImpl->propose_valid_param_move(p))
    {
        update_extra_global_loglike();

        if (pImpl->accept_valid_param_move(p))
            is_accepted = true;
        else
            pImpl->restore_extra_param_value(p);
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
        move_independent_param(p);

    // move parameters with extra dependecies
    for (const parID &p : pImpl->param_extra_infer_list)
        move_extra_depen_param(p);
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
        pImpl->update_param_accept_rate(p, move_independent_param(p));
        pImpl->adapt_param_rule(p, adapt_ratio);
    }

    // move parameters with extra dependecies
    for (const parID &p : pImpl->param_extra_infer_list)
    {
        pImpl->update_param_accept_rate(p, move_extra_depen_param(p));
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
        pImpl->setup_rand_init_params = false; // disable random initialization
        pImpl->init_parameters();              // check no undefined parameters
        init_infec_profile();                  // initialize user infectivity profile
        init_model_for_inference();            // initialize model for inference
        init_augmented_variables();            // initialize augmented dependencies
        init_extra_dependecies();              // initialize extra dependencies
        update_simple_global_loglike();        // initialize global loglikelihood
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
void MCMC::run_monte_carlo(progress::Bar &progressbar, ulong sweep)
{
    //==========================//
    //===== Initialization =====//
    //==========================//

    pImpl->init_diagnostics(); // initialize diagnostics

    // run multiple MCMC chains
    for (pImpl->_chain_ = 1; pImpl->_chain_ <= pImpl->setup_n_chains; ++pImpl->_chain_)
    {
        // update progressbar text
        progressbar.set_text(pImpl->fmt_progtxt());

        pImpl->init_parameters();       // initialize parameter inference
        init_augmented_variables();     // initialize augmented dependencies
        init_extra_dependecies();       // initialize extra dependencies
        update_simple_global_loglike(); // initialize global loglikelihood

        //=========================================//
        //======== (Finite) Adaptative MCMC =======//
        //=========================================//

        // initialize the param move counts for adaptive
        pImpl->init_move_counts();

        for (ulong i = 0; i < pImpl->setup_n_adapt; i += pImpl->setup_thin)
        {
            for (ulong t = 1; t <= pImpl->setup_thin; ++t)
                move_all_and_adapt_parameters(i + t);

            progressbar.update(++sweep);
        }

        //============================================//
        //========= MCMC Burn-in Iterations ==========//
        //============================================//

        for (ulong i = 0; i < pImpl->setup_n_burnin; i += pImpl->setup_thin)
        {
            for (ulong t = 1; t <= pImpl->setup_thin; ++t)
                move_all();

            progressbar.update(++sweep);
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

            progressbar.update(++sweep);
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
                                  pImpl->fmt_progtxt(), pImpl->setup_is_progbar);

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

/* Evaluate the logarithm of marginal likelihood with power posteriors.
    -----
    Ref: Friel, N., Hurn, M. & Wyse, J.
    Improving power posterior estimation of statistical evidence.
    Stat Comput 24, 709–723 (2014).
*/
double MCMC::estimate_evidence(const std::string &scheme, int n_rungs)
{
    const str_vec options = {"adaptive", "power-fractions"};
    double output = 0.;

    try
    {
        pImpl->assert_proper_dataset();    // check if dataset is empty
        pImpl->assert_proper_mcmc_setup(); // check setup for mcmc sampling

        // fix internal number of rungs (for progress bar)
        if (n_rungs >= 2)
            pImpl->_n_rungs_ = n_rungs;
        else
            throw std::invalid_argument("N. of rungs must be greater than 1.");

        // backup user setup to restore them later
        bool bkp_rand_init_params = pImpl->setup_rand_init_params;
        double bkp_inv_T = pImpl->setup_inv_T;

        // thermodynamic ladder estimates
        bool is_right = false;
        double avg, var;
        dbl_vec inv_t = {};   // inverse temperature
        dbl_vec avg_vec = {}; // mean loglikelihood
        dbl_vec var_vec = {}; // variance loglikelihood

        // total number of mcmc sweeps
        ulong rung_sweeps = pImpl->setup_c_n_sweeps * pImpl->setup_n_chains;
        // internal state variables
        pImpl->_rung_ = 0;

        // initialize model for inference
        init_model_for_inference();

        // start progress bar and chrono
        progress::Bar progressbar(rung_sweeps * (pImpl->_n_rungs_ + 1),
                                  pImpl->fmt_progtxt(), pImpl->setup_is_progbar);

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
            throw std::invalid_argument(
                "Unknown power posterior scheme \"" + scheme + "\"." + tools::find_closest(scheme, options, 3));

        //----- trapezoidal rule to compute marginal likelihood -----//
        double d_inv_t;
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
    }
    catch (std::exception &e)
    {
        tools::error(e.what());
    }

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
