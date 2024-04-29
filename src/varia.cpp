// varia.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#include "varia.h"

// Warning for status
bool Invalid_Proposer_Warning = true;

///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////  Base_Proposer : constructor ////////////////////////////////// //
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
                << " {low:" << _range_min_ << ", up:" << _range_max_ << "} are inconsistent."
                << "\nIf it was intentional, ignore this message.";
            raise_warning_message(oss.str(), &Invalid_Proposer_Warning);
        }
        zero_start = true;
        set_strat_zero(strategy_from_zero);
    }
};

/* increase spread if acceptation rate is high and viceversa */
void Base_Proposer::adapt(double ratio, double target)
{
    // linear
    double factor = 1. + (1. - ratio) * (get_accept_rate() - target);
    
    // quadratic 
    //double factor = 1. + sqrt(1. - ratio) * (get_accept_rate() - target);

    for (int argID : args_to_adapt)
    {
        (distr_ptr->set_arg)(argID, ((distr_ptr->get_arg)(argID)) * factor);
    }
}

/* reset acceptation rate */
void Base_Proposer::reset_accept_rate()
{
    accept_rate = 0.0;
    num_proposed_moves = 0;
};

/* update accept rate */
void Base_Proposer::update_accept_rate(bool is_accepted)
{
    accept_rate = (double)(accept_rate * num_proposed_moves + is_accepted) / (num_proposed_moves + 1);
    num_proposed_moves += 1;
}

/* set range validation functions */
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
            std::string err_text = "Provided an inconsistent range: {low:" + std::to_string(_range_min_) + ", up:" + std::to_string(_range_max_) + "}.";
            raise_error_message(err_text);
            std::exit(EXIT_FAILURE);
        }
        else
        {
            _valid_func_ = &Base_Proposer::_valid_range_both_;
        }
    }
    else
    {
        std::string err_text = "Provided multiple definition for the range.";
        raise_error_message(err_text);
        std::exit(EXIT_FAILURE);
    }
}

/* Switchboard to select the custom metropolis strategy */
Base_Proposer *proposer_switchboard(std::string name, sdargs args, double strat_zero, sdargs value_range)
{
    if (name == UNDEFINED)
    {
        return new Undefined_Proposer();
    }
    else if (name == "uniform")
    {
        return new Uniform_Proposer(args, strat_zero, value_range);
    }
    else if (name == "normal")
    {
        return new Normal_Proposer(args, strat_zero, value_range);
    }
    else if (name == "lognorm")
    {
        return new Lognorm_Proposer(args, strat_zero, value_range);
    }
    else
    {
        std::string err_text = "Unknown metropolis strategy `" + name + "`.";
        raise_error_message(err_text);
        std::exit(EXIT_FAILURE);
    }
}

/* Display method for proposer (building on the distribution method) */
std::string Base_Proposer::message_to_display()
{
    std::string output = distr_ptr->message_to_display();
    output += DELIM;
    if (_range_min_ > -INF_DOUBLE)
    {
        std::string aux = "min";
        output += aux + DICT_DELIM + std::to_string(_range_min_);
        if (_range_max_ < INF_DOUBLE)
        {
            output += DICT_DELIM;
        }
    }
    if (_range_max_ < INF_DOUBLE)
    {
        std::string aux = "max";
        output += aux + DICT_DELIM + std::to_string(_range_max_);
    }

    return output;
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////
//////// Parameter : methods ///////////////////
///////////////////////////////////////////////////////////

/* Update the values for the Bayesian diagnostics */
void Parameter::update_diagnostics()
{
    dgn_sum += get_value();
    dgn_sqr_sum += std::pow(get_value(), 2);
    dgn_sample_size++;
}

/* Message to be used at the moment of display diagnostics */
std::string Parameter::message_to_display_diagnostics()
{
    std::string output = std::to_string(dgn_sum / dgn_sample_size) + DELIM;
    output += std::to_string(std::sqrt(est_sample_variance(dgn_sum, dgn_sqr_sum, dgn_sample_size)));
    output += DELIM + (std::string)((infer) ? "true" : "false");
    if (infer)
    {
        output += DELIM + distr_ptr->message_to_display() + DELIM + prop_ptr->message_to_display();
    }

    return output;
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////
//////// Monte_Carlo_Core : methods ///////////////////
///////////////////////////////////////////////////////////

/* Check if the given parameter is defined, otherwise exits. */
void Monte_Carlo_Core::assert_param(const std::string k)
{
    if (params.find(k) == params.end())
    {
        std::string err_text = "Parameter `" + k + "` is not defined.\n Check your model and parameter file.";
        raise_error_message(err_text);
        std::exit(EXIT_FAILURE);
    }
}

/* add_parameter */
void Monte_Carlo_Core::add_parameter(
    const std::string parID,
    double new_val,
    std::string prior_name,
    sdargs prior_args,
    std::string strat_name,
    sdargs strat_args,
    sdargs value_range)
{
    Parameter new_param(new_val,
                        prior_name,
                        prior_args,
                        strat_name,
                        strat_args,
                        value_range);

    if (params.find(parID) == params.end())
    {
        params.emplace(parID, new_param);
    }
    else
    {
        std::string err_text = "Parameter `" + parID + "` has multiple definitions.";
        raise_error_message(err_text);
        std::exit(EXIT_FAILURE);
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

/* randomly initialize values of the parameters to be inferred */
void Monte_Carlo_Core::init_rand_params()
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
                std::string err_text = "Number of attempts exceeded in random initialization of `" + parID + "`.";
                raise_error_message(err_text);
                std::exit(EXIT_FAILURE);
            }
            num_attempts++;
            new_value = params.at(parID).prior_rand(prng);
        }
        set_value(parID, new_value);
    }
}

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

/* Update the diagnostics of all parameters */
void Monte_Carlo_Core::update_diagnostics()
{
    for (std::string parID : param_to_be_inferred)
    {
        params.at(parID).update_diagnostics();
    }
}

///////////////////////////////////////////////////////////
//////// Monte_Carlo_Core : display ///////////////////
///////////////////////////////////////////////////////////

/* display the user model parameters */
void Monte_Carlo_Core::display_parameters()
{
    if (Verbose)
    {
        std::ostringstream oss;
        oss << param_file_required_headers[0];
        for (size_t k = 1; k < param_file_required_headers.size(); k++)
        {
            oss << DELIM << param_file_required_headers[k];
        }
        for (std::string parID : param_to_be_inferred)
        {
            oss << "\n" << parID << DELIM << params.at(parID).message_to_display();
        }
        for (std::string parID : param_static)
        {
            oss << "\n" << parID << DELIM << params.at(parID).message_to_display();
        }
        display_stdout(oss.str());
    }
}

/* display the diagnostics */
void Monte_Carlo_Core::display_param_diagnostics()
{
    std::ostringstream oss;
    oss << "name" << DELIM << "value" << DELIM << "value_sd";
    for (size_t k = 2; k < param_file_required_headers.size(); k++)
    {
        oss << DELIM << param_file_required_headers[k];
    }
    for (std::string parID : param_to_be_inferred)
    {
        oss << "\n" << parID << DELIM << params.at(parID).message_to_display_diagnostics();
    }
    display_stdout(oss.str());
}