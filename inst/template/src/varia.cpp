/*
 * varia.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "varia.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base_Proposer : constructor  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Base_Proposer::Base_Proposer(Base_Distribution *distr,
                             double rule_start,
                             std::vector<int> args_to_adapt) : distr_ptr(distr)
{
    _args_to_adapt_ = args_to_adapt;

    if (rule_start == INF_DOUBLE)
    {
        _is_start_from_zero_ = false;
    }
    else
    {
        _is_start_from_zero_ = true;
        set_rule_start_value(rule_start);
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Base_Proposer : methods  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Increase spread if acceptation rate is high and viceversa. */
void Base_Proposer::adapt(double ratio)
{
    double factor = 1. + (1. - ratio) * (get_accept_rate() - TARGET_ACCEPT_RATE);

    for (const int &argID : _args_to_adapt_)
    {
        (distr_ptr->set_arg)(argID, ((distr_ptr->get_arg)(argID)) * factor);
    }
}

/* Reset acceptation rate. */
void Base_Proposer::reset_accept_rate()
{
    accept_rate = 0.0;
    num_proposed_moves = 0;
};

/* Update the parameter accept rate. */
void Base_Proposer::update_accept_rate(bool is_accepted)
{
    accept_rate = (double)(accept_rate * num_proposed_moves + is_accepted) / (num_proposed_moves + 1);
    num_proposed_moves += 1;
}

/* Switchboard to select the custom metropolis proposer rule. */
Base_Proposer *proposer_switchboard(std::string rule_name, dictargs rule_args, double rule_from)
{
    Base_Proposer *output = new Base_Proposer(new Base_Distribution());

    if (rule_name == "unif")
    {
        output = new Uniform_Proposer(rule_args, rule_from);
    }
    else if (rule_name == "norm")
    {
        output = new Normal_Proposer(rule_args, rule_from);
    }
    else if (rule_name == "lnorm")
    {
        output = new Lognorm_Proposer(rule_args, rule_from);
    }
    else if (!tools::is_na(rule_name))
    {
        throw std::invalid_argument("Unknown proposer rule '" + rule_name + "' for Metropolis-Hastings.");
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
/////////////////////////////////////  Parameter : methods  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Init the values for the Bayesian diagnostics. */
void Parameter::init_diagnostics(ulong sample_size)
{
    dgn_count = 0;
    dgn_value.clear();
    dgn_value = std::vector<double>(sample_size, std::numeric_limits<double>::quiet_NaN());
}

/* Update the values for the Bayesian diagnostics. */
void Parameter::update_diagnostics()
{
    dgn_value[dgn_count++] = get_value();
}

/* Message to be used at the moment of display diagnostics. */
std::string Parameter::show_diagnostics() const
{
    using namespace tools;
    std::ostringstream oss;

    double avg = 0, var = 0;
    statistics(dgn_value, avg, var);

    oss << avg << C_DELIM << std::sqrt(var) << C_DELIM << truth(infer)
        << C_DELIM << show_distr_name() << C_DELIM << show_distr_arg_names() << C_DELIM << show_distr_arg_values()
        << C_DELIM << show_prop_name() << C_DELIM << show_prop_arg_names() << C_DELIM << show_prop_arg_values()
        << C_DELIM << show_range_min() << C_DELIM << show_range_max();
    return oss.str();
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Monte_Carlo_Core : methods (getter)  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Functional method to get a parameter value. */
double Monte_Carlo_Core::get_param(const std::string &parID) const
{
    return (this->*_get_param_method_)(parID);
}

/* Default method to get a value. */
double Monte_Carlo_Core::_get_param_default_(const std::string &parID) const
{
    return params.at(parID).get_value();
}

/* Safe method to get a value with previous assertion. */
double Monte_Carlo_Core::_get_param_assert_(const std::string &parID) const
{
    assert_param(parID);
    return _get_param_default_(parID);
}

/* Functional method to get a pointer to a parameter. */
const double *Monte_Carlo_Core::get_param_pointer(const std::string &parID) const
{
    return (this->*_get_param_pointer_method_)(parID);
}

/* Default method to get a pointer to a parameter. */
const double *Monte_Carlo_Core::_get_param_pointer_default_(const std::string &parID) const
{
    return params.at(parID).get_pointer();
}

/* Safe method to get a pointer to a parameter with previous assertion. */
const double *Monte_Carlo_Core::_get_param_pointer_assert_(const std::string &parID) const
{
    assert_param(parID);
    return _get_param_pointer_default_(parID);
}

/* Disable parameter assertion. */
void Monte_Carlo_Core::disable_assert_param()
{
    _get_param_method_ = &Monte_Carlo_Core::_get_param_default_;
    _get_param_pointer_method_ = &Monte_Carlo_Core::_get_param_pointer_default_;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Monte_Carlo_Core : methods (setter)  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Drop all params dictionary. */
void Monte_Carlo_Core::drop_all_params()
{
    params.clear();
    param_to_be_static.clear();
    param_to_be_inferred.clear();
}

/* Assign a new parameter to the dictionary. */
void Monte_Carlo_Core::assign_param(const std::string &parID, const Parameter &new_param)
{
    // create the new parameter if it's not already defined
    if (params.find(parID) != params.end())
    {
        throw std::invalid_argument("Parameter '" + parID + "' has multiple definitions.");
    }
    else
    {
        // emplace the new parameter
        params.emplace(parID, new_param);

        // update the lists of parameters to be inferred, and these to be static
        if (new_param.get_infer())
        {
            param_to_be_inferred.push_back(parID);
        }
        else
        {
            param_to_be_static.push_back(parID);
        }
    }
}

/* Set the number of threads for parrell computations (default is 1). */
void Monte_Carlo_Core::set_num_threads(unsigned n)
{
    if (n == 0)
    {
        throw std::invalid_argument("Number of threads must be equal or greater than 1.");
    }

#ifdef _OPENMP
    if (n > _max_num_available_threads)
    {
        std::ostringstream oss;
        oss << "Number of threads (" << n << ") changed to number of available threads (" << _max_num_available_threads << ").";
        tools::warning(oss.str());
        n = _max_num_available_threads;
    }

    // declare the global number of omp threads
    omp_set_num_threads(n);
    setup_num_threads = n;

    // create different generators for different parallel threads
    // UPDATE: they are used only to parallelize data augmentantion
    if (n > 1)
    {
        prng_parallel.reserve(n);
        _prng_method_ = &Monte_Carlo_Core::_prng_parallel_threads_;
    }
    else
    {
        _prng_method_ = &Monte_Carlo_Core::_prng_unique_thread_;
    }
#else
    if (n != 1)
    {
        tools::warning("Parallelization with OpenMP not supported.");
    }
    setup_num_threads = 1;
    _prng_method_ = &Monte_Carlo_Core::_prng_unique_thread_;
#endif
}

/* Assign seed for random generators. */
void Monte_Carlo_Core::set_prng_seed(int seed)
{
    // seed for the unique thread rand num generator
    prng.seed(seed);

    // seeds for the case of parallelization
    if (setup_num_threads > 1)
    {
        std::uniform_int_distribution<unsigned long int> seed_generator(0, 0XFFFFFFFF);
        for (unsigned n = 0; n < setup_num_threads; ++n)
        {
            prng_parallel[n].seed(seed_generator(prng));
        }
    }
}

/* Add a new scalar for augmented typedict FIXME : temporary. */
const double *Monte_Carlo_Core::assign_scalar_pointer(double new_val)
{
    if (_scalar_type_params_.find(new_val) == _scalar_type_params_.end())
    {
        Parameter new_param(UNDEFINED, new_val, UNDEFINED, {}, UNDEFINED, {}, {});
        _scalar_type_params_.emplace(new_val, new_param);
    }
    return _scalar_type_params_.at(new_val).get_pointer();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Monte_Carlo_Core : methods (setter)  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Check if the given parameter is defined, otherwise raise error. */
void Monte_Carlo_Core::assert_param(const std::string &parID) const
{
    if (params.find(parID) == params.end())
    {
        throw std::invalid_argument("Parameter '" + parID + "' not found. Check your model and parameter definitions.");
    }
}

/* randomly initialize values of the parameters to be inferred */
void Monte_Carlo_Core::init_rand_params(bool is_rand_params)
{
    if (is_rand_params)
    {
        for (const std::string &parID : param_to_be_inferred)
        {
            params.at(parID).rand(prng);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////  Monte_Carlo_Core : methods (parameter moves)  //////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* propose a valid move (retrun true), otherwise update refusal (retrun false) */
// FIXME : Gibbs sampling is missing here!
// 1) always valid 2) move_logratio == 0
bool Monte_Carlo_Core::propose_valid_param_move(const std::string &parID)
{
    // set new value of the parameter and store the old to buffer
    params.at(parID).propose(prng);

    // proceed only if valid move, otherwise restore buffer
    if (params.at(parID).valid_move())
    {
        return true;
    }
    else
    {
        // reassign buffered value to param
        restore_param_value(parID);
        return false;
    }
}

/* Init diagnostics for all parameters */
void Monte_Carlo_Core::init_param_diagnostics(ulong sample_size)
{
    for (const std::string &parID : param_to_be_inferred)
    {
        params.at(parID).init_diagnostics(sample_size);
    }
}

/* Update diagnostics for all parameters */
void Monte_Carlo_Core::update_param_diagnostics()
{
    for (const std::string &parID : param_to_be_inferred)
    {
        params.at(parID).update_diagnostics();
    }
}

/* Format parameter to show. */
std::string Monte_Carlo_Core::show_param(const std::string &parID) const
{
    std::ostringstream oss;
    oss << parID << C_DELIM << params.at(parID).show();
    return oss.str();
}

/* Format parameter diagnostics to show. */
std::string Monte_Carlo_Core::show_param_diagnostics(const std::string &parID) const
{
    std::ostringstream oss;
    oss << parID << C_DELIM << params.at(parID).show_diagnostics();
    return oss.str();
}