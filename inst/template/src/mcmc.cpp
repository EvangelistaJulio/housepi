/*
 * mcmc.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "mcmc.h"

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

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC : methods (general)  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* House loglikelihood including augmented modl contribution. */
double MCMC::compute_augm_loglike(const Household &this_house) const
{
    double hh_loglike = 0.;

    // Contribution of model associated to augmented variables
    for (const std::string &augmID : augm_with_model_set)
    {
        for (int i : this_house.infectious_set)
        {
            hh_loglike += this_house.indiv[i].augm_model_logpdf(augmID);
        }
    }

    return hh_loglike;
}

/* Compute new loglikelihood, store old value in buffer. Same for households. */
void MCMC::update_global_loglike()
{
    std::atomic<bool> error_occurred = false;
    std::string error_message;

    buffer = global_loglike;
    global_loglike = 0.0;

#pragma omp parallel for reduction(+ : global_loglike)
    for (int h = 0; h < num_house; ++h)
    {
        if (_error_occurred_)
        {
            continue;
        }
        try
        {
            Household &this_house = house[h];
            update_house_loglike(this_house);
            global_loglike += this_house.get_loglike() + compute_augm_loglike(this_house);
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

/* Set accept rates to zero. */
void MCMC::init_move_counts()
{
    // ... for parameters
    for (const std::string &parID : param_to_be_inferred)
    {
        reset_accept_rate(parID);
    }

    // ... for augmented variables
    for (Household &this_house : house)
    {
        for (const std::string &augmID : augm_to_be_inferred)
        {
            this_house.reset_augm_accept_rate(augmID);
        }
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
    if (propose_valid_param_move(parID))
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
            restore_param_value(parID);

            // restore buffered loglikelihood in households
            for (int h = 0; h < num_house; ++h)
            {
                house[h].restore_loglike();
            }
        }
    }

    // update counts for acceptance rate
    update_param_accept_rate(parID, is_accepted);
}

bool MCMC::_param_metropolis_1_(const std::string &parID)
{
    update_global_loglike();

    if (uniform() < exp(get_delta_global_loglike() + param_move_logratio(parID)))
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
    if (uniform() < exp(param_move_logratio(parID)))
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

    if (uniform() < exp(setup_inv_T * get_delta_global_loglike() + param_move_logratio(parID)))
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
    for (const std::string &parID : param_to_be_inferred)
    {
        move_single_param(parID);
    }
}

/* Move all model parameters and adapt their rule arguments. */
void MCMC::move_and_adapt_all_parameters(ulong step)
{
    double adapt_ratio = static_cast<double>(step) / setup_num_adapt;

    for (const std::string &parID : param_to_be_inferred)
    {
        move_single_param(parID);
        adapt(parID, adapt_ratio);
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
    for (Household &this_house : house)
    {
        for (int i : this_house.infectious_set)
        {
            update_augm_dependecies(this_house, i);
        }
    }
}

/* Move augmented variables from all households. */
void MCMC::move_all_augm_variables()
{
    if (_is_move_augm_)
    {
        double delta_ll = 0.;

#pragma omp parallel for reduction(+ : delta_ll)
        for (int h = 0; h < num_house; ++h)
        {
            if (_error_occurred_)
            {
                continue;
            }
            try
            {
                delta_ll += move_single_house_augm(house[h]);
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
        for (const std::string &augmID : augm_to_be_inferred)
        {
            bool is_accepted = false;

            // proceed only if the new value is consistent
            // FIXME : Gibbs sampling should be always valid but it's not now
            if (this_indiv.propose_valid_augm_move(augmID, get_prng()))
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

// FIXME developing
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

// FIXME developing
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

// FIXME developing
bool MCMC::_augm_metropolis_T_(const Household &this_house, int i, const std::string &augmID, double &delta)
{
    double D_ll = this_house.get_delta_loglike() + this_house.indiv[i].augm_model_logratio(augmID);
    if (uniform() < exp(setup_inv_T * D_ll - this_house.indiv[i].augm_propose_logratio(augmID)))
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
    //===== Initialization =====//

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
/////////////////////////////////  MCMC : methods (diagnostics) ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Init Bayesian diagnostics. */
void MCMC::init_diagnostics()
{
    ulong sample_size = (setup_num_iter - setup_burn_in) / setup_thinning;
    sample_size *= setup_num_chains;

    dgn_count = 0;
    dgn_glob_loglike = std::vector<double>(sample_size);

    if (!_is_power_post_)
    {
        // diagnostics for all parameters
        init_param_diagnostics(sample_size);

        // diagnostics over households
        for (int h = 0; h < num_house; ++h)
        {
            house[h].init_diagnostics(sample_size);
        }
    }
}

/* Update Bayesian diagnostics */
void MCMC::update_diagnostics()
{
    // diagnostics for global loglike
    dgn_glob_loglike[dgn_count++] = get_global_loglike();

    if (!_is_power_post_)
    {
        // diagnostics for all parameters
        update_param_diagnostics();

        // diagnostics over households
        for (int h = 0; h < num_house; ++h)
        {
            house[h].update_diagnostics();
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : run markov chain  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* The text to print in the Markov chian progress bar. */
std::string MCMC::get_progtxt(const std::string &txt)
{
    std::ostringstream oss;
    if (_is_power_post_)
    {
        oss << "power post " << _rung_ << "/" << _num_rungs_;
    }
    else if (setup_num_chains > 1)
    {
        oss << txt << " " << _chain_ << "/" << setup_num_chains;
    }
    else
    {
        oss << txt;
    }
    return oss.str();
}

/* Disable all assert methods for user definitions. */
void MCMC::disable_assert()
{
    disable_assert_param();

    // disable augmented variable assertion for all individuals
    for (Household &this_house : house)
    {
        for (Individual &this_indiv : this_house.indiv)
        {
            this_indiv.disable_assert();
        }
    }
}

/* Markov chain core, used for inference and power posteriors.
    note: sweep is the number of previously completed thin sweeps. */
void MCMC::run_monte_carlo(bool progbar, ulong sweep)
{
    using namespace tools;

    //==========================//
    //===== Initialization =====//
    //==========================//

    //----- Initialize model -----//
    init_model_for_inference();

    //----- Initialize diagnostics -----//
    init_diagnostics();

    //===== Run Multiple MCMC Chains =====//
    for (_chain_ = 1; _chain_ <= setup_num_chains; ++_chain_)
    {
        //-----Initialize parameter inference-----//
        init_rand_params(setup_rnd_init_params);

        //----- Initialize augmented dependencies -----//
        init_augmented_variables();

        //----- Initialize global loglikelihood -----//
        update_global_loglike();

        // avoid asserting user definitions after first round
        // UPDATE : this can be improved, e. g. with a debug mode
        disable_assert();

        // initialize the acceptance counts
        init_move_counts();

        //=========================================//
        //======== (Finite) Adaptative MCMC =======//
        //=========================================//

        for (ulong i = 0; i < setup_num_adapt; i += setup_thinning)
        {
            for (ulong t = 1; t <= setup_thinning; ++t)
            {
                move_all_augm_variables();
                move_and_adapt_all_parameters(i + t);
            }
            ++sweep;
            progress_bar_run(get_progtxt("adapt"), sweep, _tot_thin_sweeps_, progbar);
        }

        // reset the acceptance counts
        // UPDATE : move counts become useless since they are not even printed : improve.
        if (setup_num_adapt > 0)
        {
            init_move_counts();
        }

        //============================================//
        //========= MCMC Burn-in Iterations ==========//
        //============================================//

        for (ulong i = 0; i < setup_burn_in; i += setup_thinning)
        {
            for (ulong t = 1; t <= setup_thinning; ++t)
            {
                move_all_augm_variables();
                move_all_parameters();
            }
            ++sweep;
            progress_bar_run(get_progtxt("burn"), sweep, _tot_thin_sweeps_, progbar);
        }

        //============================================//
        //========= MCMC Posterior Sampling ==========//
        //============================================//

        for (ulong i = setup_burn_in; i < setup_num_iter; i += setup_thinning)
        {
            // thinning mcmc moves
            for (ulong t = 1; t <= setup_thinning; ++t)
            {
                move_all_augm_variables();
                move_all_parameters();
            }
            ++sweep;

            // compute diagnostics for the previous iterations
            update_diagnostics();
            progress_bar_run(get_progtxt("sample"), sweep, _tot_thin_sweeps_, progbar);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  MCMC : infer  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Performs the inference using Metropolis-Hastings. */
// FIXME : Gelman Rubin could be integrated
void MCMC::infer(bool progbar)
{
    //----- Check if data is meaningful -----//
    if (report_data["n.miss.stat"] == num_indiv)
    {
        tools::error("No infectious status has been provided in the data file.");
    }

    //----- total number of mcmc sweeps -----//
    _tot_thin_sweeps_ = _thin_sweeps_;

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
    bool bkp_rnd_init_params = setup_rnd_init_params;
    double bkp_inv_T = setup_inv_T;

    //----- Thermodynamic ladder estimates -----//
    bool is_right = false;
    double avg, var;
    std::vector<double> inv_t = {};   // inverse temperature
    std::vector<double> avg_vec = {}; // mean loglikelihood
    std::vector<double> var_vec = {}; // variance loglikelihood

    //----- Fix internal states -----//
    _is_power_post_ = true;

    //----- Total number of mcmc sweeps -----//
    _tot_thin_sweeps_ = _thin_sweeps_ * (_num_rungs_ + 1);

    //----- Start progress bar and chrono -----//
    tools::progress_bar_init("init", progbar);

    //==================================//
    //===== INVERSE TEMPEREATURE 0 =====//
    //==================================//

    setup_rnd_init_params = true;
    _rung_ = 0;
    set_inver_temper(0.0);
    run_monte_carlo(progbar, 0);
    tools::statistics(dgn_glob_loglike, avg, var);

    inv_t.push_back(setup_inv_T);
    avg_vec.push_back(avg);
    var_vec.push_back(var);

    // thermodynamic ladder: use modified power fractions or adaptive scheme
    if (adaptive)
    {
        //===========================//
        //===== ADAPTIVE LADDER =====//
        //===========================//

        //----- Inverse Temperature 1 -----//
        setup_rnd_init_params = false;
        _rung_ = 1;
        set_inver_temper(1.0);
        run_monte_carlo(progbar, _thin_sweeps_);
        tools::statistics(dgn_glob_loglike, avg, var);

        inv_t.push_back(setup_inv_T);
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
            set_inver_temper(new_inv_T);
            run_monte_carlo(progbar, _rung_ * _thin_sweeps_);
            tools::statistics(dgn_glob_loglike, avg, var);

            //----- find larger rectangular discretization error -----//
            is_right = (inv_t[k + 1] - setup_inv_T) * (avg_vec[k + 1] - avg) > (new_inv_T - inv_t[k]) * (avg - avg_vec[k]);

            //----- insert in-between for all lists -----//
            inv_t.insert(inv_t.begin() + k + 1, setup_inv_T);
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
        setup_rnd_init_params = false;
        for (_rung_ = 1; _rung_ < _num_rungs_; ++_rung_)
        {
            set_inver_temper(std::pow(static_cast<double>(_rung_) / (_num_rungs_ - 1), 5));
            run_monte_carlo(progbar, _rung_ * _thin_sweeps_);
            tools::statistics(dgn_glob_loglike, avg, var);

            inv_t.push_back(setup_inv_T);
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
    set_inver_temper(bkp_inv_T);
    setup_rnd_init_params = bkp_rnd_init_params;

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

    //===== Initialization =====//

    //----- Initialize model -----//
    init_model_for_simulation();

    //----- Initialize augmented dependencies -----//
    init_augmented_variables();

    tools::progress_bar_init("simul", progbar);

#ifdef _OPENMP
    int thin_num_house = num_house / get_num_threads();
#pragma omp parallel for default(shared)
    for (int h = 0; h < num_house; ++h)
    {
        if (_error_occurred_)
        {
            continue;
        }
        try
        {
            simulate_single_house(house[h]);
            if (omp_get_thread_num() == 0)
            {
                tools::progress_bar_run("simul", h, thin_num_house, progbar);
            }
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
#else
    for (int h = 0; h < num_house; ++h)
    {
        if (_error_occurred_)
        {
            continue;
        }
        try
        {
            simulate_single_house(house[h]);
            tools::progress_bar_run("simul", h, num_house, progbar);
        }
        catch (const std::exception &e)
        {
            _error_occurred_ = true;
            _error_message_ = e.what();
        }
    }
#endif

    if (_error_occurred_)
    {
        tools::error(_error_message_);
    }
    else
    {
        tools::progress_bar_end("simul", progbar);
        update_report_data();
    }
}