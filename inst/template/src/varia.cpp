// varia.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#include "varia.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  WARNINGS  /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

bool Invalid_Proposer_Range_Warning = false; // warning for invalid proposer range declaration
bool Disabled_OpenMP_Warning = false;        // warning for attempt to use openmp while disabled
bool Above_Max_Threads_Warning = false;      // warning for over the threshold n. of threads

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base_Proposer : constructor  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Base_Proposer::Base_Proposer(
    Distribution_1D *distr,
    double strategy_from_zero,
    sdargs value_range,
    std::vector<int> list_of_args_to_adapt) : distr_ptr(distr)
{
    args_to_adapt = list_of_args_to_adapt;
    set_range(value_range);

    if (strategy_from_zero == INF_DOUBLE)
    {
        zero_start = false;
    }
    else
    {
        if (!valid(strategy_from_zero))
        {
            std::ostringstream oss;
            oss << "Initial value " << strategy_from_zero << " and value range"
                << " {min:" << _range_min_ << ", max:" << _range_max_ << "} are inconsistent."
                << "\nIf it was intentional, ignore this message.";
            screen::raise_warning(oss.str(), Invalid_Proposer_Range_Warning);
        }
        zero_start = true;
        set_strat_zero(strategy_from_zero);
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Base_Proposer : methods  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Increase spread if acceptation rate is high and viceversa. */
void Base_Proposer::adapt(double ratio)
{
    double factor = 1. + (1. - ratio) * (get_accept_rate() - TARGET_ACCEPT_RATE);

    for (int argID : args_to_adapt)
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

/* Set the range of validity of the value, and the corresponding validation function. */
void Base_Proposer::set_range(sdargs new_range)
{
    int condition = 0;

    for (auto const &x : new_range)
    {
        if (x.first == "min")
        {
            _range_min_ = x.second;
            condition += 2;
        }
        else if (x.first == "max")
        {
            _range_max_ = x.second;
            condition += 3;
        }
    }

    if (condition == 0)
    {
        _valid_func_ = &Base_Proposer::_valid_dummy_;
    }
    else if (condition == 2)
    {
        _valid_func_ = &Base_Proposer::_valid_range_min_;
    }
    else if (condition == 3)
    {
        _valid_func_ = &Base_Proposer::_valid_range_max_;
    }
    else if (condition == 5)
    {
        if (_range_min_ >= _range_max_)
        {
            std::ostringstream oss;
            oss << "Provided an inconsistent range: {min:" << _range_min_ << ", max:" << _range_max_ << "}.";
            screen::raise_error(oss.str());
        }
        else
        {
            _valid_func_ = &Base_Proposer::_valid_range_both_;
        }
    }
    else
    {
        screen::raise_error("Provided multiple definition for the range.");
    }
}

/* Switchboard to select the custom metropolis strategy. */
Base_Proposer *proposer_switchboard(std::string name, sdargs args, double strat_zero, sdargs value_range)
{
    Base_Proposer *output = new Undefined_Proposer();

    if (name == "uniform")
    {
        output = new Uniform_Proposer(args, strat_zero, value_range);
    }
    else if (name == "normal")
    {
        output = new Normal_Proposer(args, strat_zero, value_range);
    }
    else if (name == "lognorm")
    {
        output = new Lognorm_Proposer(args, strat_zero, value_range);
    }
    else if (name != UNDEFINED)
    {
        screen::raise_error("Unknown metropolis strategy `" + name + "`.");
    }

    return output;
}

/* Display method for proposer (building on the distribution method). */
std::ostringstream Base_Proposer::message_to_display()
{
    std::ostringstream oss;
    oss << distr_ptr->message_to_display().str() << screen::C_DELIM;
    if (_range_min_ > -INF_DOUBLE)
    {
        oss << "min" << screen::C_DICT_DELIM << _range_min_;
        if (_range_max_ < INF_DOUBLE)
        {
            oss << screen::C_DICT_DELIM;
        }
    }
    if (_range_max_ < INF_DOUBLE)
    {
        oss << "max" << screen::C_DICT_DELIM << _range_max_;
    }

    return oss;
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Parameter : methods  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Init the values for the Bayesian diagnostics. */
void Parameter::init_diagnostics(int sample_size)
{
    dgn_sample_size = sample_size;
    dgn_mean = 0.;
    dgn_mean_sqr = 0;
}

/* Update the values for the Bayesian diagnostics. */
void Parameter::update_diagnostics()
{
    dgn_mean += get_value() / dgn_sample_size;
    dgn_mean_sqr += std::pow(get_value(), 2) / dgn_sample_size;
}

/* Message to be used at the moment of display diagnostics. */
std::string Parameter::message_to_diagnostics()
{
    std::ostringstream oss;
    oss << dgn_mean << screen::C_DELIM
        << std::sqrt((dgn_mean_sqr - std::pow(dgn_mean, 2)) * (dgn_sample_size / (dgn_sample_size - 1.)))
        << screen::C_DELIM << ((infer) ? "TRUE" : "FALSE");
    if (infer)
    {
        oss << screen::C_DELIM << distr_ptr->message_to_display().str()
            << screen::C_DELIM << prop_ptr->message_to_display().str();
    }

    return oss.str();
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Monte_Carlo_Core : methods (getter)  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Monte_Carlo_Core::disable_assert_param()
{
    _get_value_method_ = &Monte_Carlo_Core::_get_value_default_;
}

double Monte_Carlo_Core::_get_value_default_(const std::string k) const
{
    return params.at(k).get_value();
}

double Monte_Carlo_Core::_get_value_assert_(const std::string k) const
{
    assert_param(k);
    return _get_value_default_(k);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Monte_Carlo_Core : methods (setter)  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Drop all params dictionary. */
void Monte_Carlo_Core::drop_all_params()
{
    params.clear();
    param_static.clear();
    param_to_be_inferred.clear();
}

/* Assign a new parameter to the dictionary. */
void Monte_Carlo_Core::assign_parameter(
    const std::string parID, double new_val,
    std::string prior_name, sdargs prior_args,
    std::string strat_name, sdargs strat_args,
    sdargs value_range)
{
    Parameter new_param(new_val,
                        prior_name, prior_args,
                        strat_name, strat_args,
                        value_range);

    if (params.find(parID) == params.end())
    {
        params.emplace(parID, new_param);
    }
    else
    {
        screen::raise_error("Parameter `" + parID + "` has multiple definitions.");
    }

    // update the lists of parameters to be inferred, and these to be static
    if (new_param.get_infer())
    {
        param_to_be_inferred.push_back(parID);
    }
    else
    {
        param_static.push_back(parID);
    }
}

/* Set the number of threads for parrell computations (default is 1). */
void Monte_Carlo_Core::set_num_threads(unsigned n)
{
    if (n == 0)
    {
        screen::raise_error("Number of threads must be equal or greater than 1.");
    }

#ifdef _OPENMP
    unsigned max_detected_threads = omp_get_max_threads();
    if (n > max_detected_threads)
    {
        std::ostringstream oss;
        oss << "Number of threads (" << n << ") changed to max amount of detected threads (" << max_detected_threads << ").";
        screen::raise_warning(oss.str(), Above_Max_Threads_Warning);
        n = max_detected_threads;
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
        screen::raise_warning("Parallelization with OpenMP not supported.", Disabled_OpenMP_Warning);
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
        for (unsigned n = 0; n < setup_num_threads; n++)
        {
            prng_parallel[n].seed(seed_generator(prng));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Monte_Carlo_Core : methods (setter)  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Check if the given parameter is defined, otherwise raise error. */
void Monte_Carlo_Core::assert_param(const std::string k) const
{
    if (params.find(k) == params.end())
    {
        screen::raise_error("Parameter `" + k + "` not found. Check your model and/or parameter definition.");
    }
}

/* randomly initialize values of the parameters to be inferred */
void Monte_Carlo_Core::init_rand_params(bool is_rand_params)
{
    if (is_rand_params)
    {
        int num_attempts = 0;
        double new_value;

        for (std::string parID : param_to_be_inferred)
        {
            new_value = params.at(parID).prior_rand(prng);
            while (params.at(parID).valid(new_value) == false)
            {
                if (num_attempts >= MAX_NUMBER_ATTEMPTS)
                {
                    screen::raise_error("Number of attempts exceeded in random initialization of `" + parID + "`.");
                }
                num_attempts++;
                new_value = params.at(parID).prior_rand(prng);
            }
            set_value(parID, new_value);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////  Monte_Carlo_Core : methods (parameter moves)  //////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* propose a valid move (retrun true), otherwise update refusal (retrun false) */
bool Monte_Carlo_Core::propose_valid_param_move(const std::string parID)
{
    double new_value = params.at(parID).propose(prng); // propose new value

    // proceed only if valid
    if (params.at(parID).valid(new_value))
    {
        // store old value into buffer
        buffer_value = get_value(parID);
        // add the prior logratio to buffer
        buffer_logratio = params.at(parID).prior_logratio(buffer_value, new_value);
        // add the proposal logratio to buffer
        buffer_logratio -= params.at(parID).propose_logratio(buffer_value, new_value);
        // set the new value of the parameter for computations
        set_value(parID, new_value);
        return true;
    }
    else
    {
        // update accept rate with the refusal
        params.at(parID).update_accept_rate(false);
        return false;
    }
}

/* move is accepted */
void Monte_Carlo_Core::accept_param_move(const std::string parID)
{
    params.at(parID).update_accept_rate(true);
}

/* move is refused: restore old value in the buffer */
void Monte_Carlo_Core::refuse_param_move(const std::string parID)
{
    set_value(parID, buffer_value);
    params.at(parID).update_accept_rate(false);
}

/* Init diagnostics for all parameters */
void Monte_Carlo_Core::init_param_diagnostics(int sample_size)
{
    for (std::string parID : param_to_be_inferred)
    {
        params.at(parID).init_diagnostics(sample_size);
    }
}

/* Update diagnostics for all parameters */
void Monte_Carlo_Core::update_param_diagnostics()
{
    for (std::string parID : param_to_be_inferred)
    {
        params.at(parID).update_diagnostics();
    }
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//