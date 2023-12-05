// varia.cpp //

#include "varia.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////  Base_Proposer : constructor ////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

Base_Proposer::Base_Proposer(
    Distribution_1D *distr,
    double strategy_from_zero,
    kwargs value_range,
    std::vector<size_t> list_of_args_to_adapt) : distr_ptr(distr)
{
    args_to_adapt = list_of_args_to_adapt;
    set_range(value_range);

    if (strategy_from_zero == INFINITY)
    {
        zero_start = false;
    }
    else
    {
        zero_start = true;
        set_strat_zero(strategy_from_zero);
        /* 
        NOTE : I sometimes initialized to not valid values
        // UPDATE : place a warning maybe
        if (valid(strategy_from_zero))
        {
            zero_start = true;
            set_strat_zero(strategy_from_zero);
        }
        else
        {
            std::string err_text = "Initial value " + std::to_string(strategy_from_zero) + " and value range ";
            err_text += " {low:" + std::to_string(_range_min_) + ", up:" + std::to_string(_range_max_) + "} are inconsistent.";
            exit_with_failure_error(err_text);
        }
        */
    }
};

/* increase spread if acceptation rate is high and viceversa */
void Base_Proposer::adapt(double ratio)
{
    double factor = 1. + (1. - ratio) * (get_accept_rate() - DEFAULT_ACCEPT_RATE);

    for (size_t argID : args_to_adapt)
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
void Base_Proposer::set_range(kwargs new_range)
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
            exit_with_failure_error(err_text);
        }
        else
        {
            _valid_func_ = &Base_Proposer::_valid_range_both_;
        }
    }
    else
    {
        std::string err_text = "Provided multiple definition for the range.";
        exit_with_failure_error(err_text);
    }
}

/* Switchboard to select the custom metropolis strategy */
Base_Proposer *proposer_switchboard(std::string name, kwargs args, double strat_zero, kwargs value_range)
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
        exit_with_failure_error(err_text);
        return nullptr;
    }
}

std::string Base_Proposer::display(char delim, char dict_delim)
{
    std::string output = distr_ptr->display(delim, dict_delim);
    output += delim;
    if (_range_min_ > -INFINITY)
    {
        std::string aux = "min";
        output += aux + dict_delim + std::to_string(_range_min_);
        if (_range_max_ < INFINITY)
        {
            output += dict_delim;
        }
    }
    if (_range_max_ < INFINITY)
    {
        std::string aux = "max";
        output += aux + dict_delim + std::to_string(_range_max_);
    }

    return output;
}

///////////////////////////////////////////////////////////
//////// Parameter_Dictionary : methods ///////////////////
///////////////////////////////////////////////////////////

/* Check if the given parameter is defined, otherwise exits. */
void Parameter_Dictionary::check(const std::string k)
{
    if (params.find(k) == params.end())
    {
        std::string err_text = "Parameter `" + k + "` is not defined.\n Check your model and parameter file.";
        exit_with_failure_error(err_text);
    }
}

/* add_param */
void Parameter_Dictionary::add_param(
    const std::string k,
    double new_val,
    std::string prior_name,
    kwargs prior_args,
    std::string strat_name,
    kwargs strat_args,
    kwargs value_range)
{
    Parameter new_param(new_val, prior_name, prior_args, strat_name, strat_args, value_range);

    if (params.find(k) == params.end())
    {
        params.emplace(k, new_param);
    }
    else
    {
        std::string err_text = "Parameter `" + k + "` has multiple definitions.";
        exit_with_failure_error(err_text);
    }

    // update the lists of parameters to be inferred or not
    if (new_param.get_infer())
    {
        param_to_be_inferred.push_back(k);
        num_param_to_infer = param_to_be_inferred.size();
    }
    else
    {
        param_static.push_back(k);
    }
}

/* display the user model parameters */
void Parameter_Dictionary::display_parameters(char delim, char dict_delim)
{
    for (std::string k : param_to_be_inferred)
    {
        std::cout << k << delim << params.at(k).display(delim, dict_delim) << endl;
    }
    for (std::string k : param_static)
    {
        std::cout << k << delim << params.at(k).display(delim, dict_delim) << endl;
    }
}

/* randomly initialize values of the parameters to be inferred */
void Parameter_Dictionary::init_rand_params(std::mt19937_64 &rng)
{
    double new_val;
    for (std::string k : param_to_be_inferred)
    {
        new_val = params.at(k).prior_rand(rng);
        while (params.at(k).valid(new_val) == false)
        {
            new_val = params.at(k).prior_rand(rng);
        }
        set_value(k, new_val);
    }
}

/* propose a valid move */
bool Parameter_Dictionary::propose_valid_move(const std::string k, std::mt19937_64 &rng)
{
    double new_val = params.at(k).propose(rng); // propose new value
    bool is_valid = false;

    // proceed only if valid
    if (params.at(k).valid(new_val))
    {
        is_valid = true;
        buffer_value = get_value(k); // store old value into buffer
        buffer_log_ratio = params.at(k).prior_log_ratio(buffer_value, new_val);
        buffer_log_ratio -= params.at(k).propose_log_ratio(buffer_value, new_val);
        set_value(k, new_val);
    }
    else
    {
        params.at(k).update_accept_rate(false);
    }

    return is_valid;
}

/* move is accepted */
void Parameter_Dictionary::accept_param_move(const std::string k)
{
    params.at(k).update_accept_rate(true);
}

/* move is refused: restore old value in the buffer */
void Parameter_Dictionary::refuse_param_move(const std::string k)
{
    set_value(k, buffer_value);
    params.at(k).update_accept_rate(false);
}