// varia.cpp //

#include "varia.h"

using namespace std;

//----------------------------------------------------------------------
// Strategy : functional methods
//----------------------------------------------------------------------

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
        if (valid(strategy_from_zero))
        {
            zero_start = true;
            set_strat_zero(strategy_from_zero);
        }
        else
        {
            std::cerr << "ERROR: initial value " << strategy_from_zero << " and value range {low:" << _range_low_ << ", up:" << _range_low_ << "} are inconsistent." << std::endl;
            std::cerr << "Exit...\n"
                      << std::endl;
            std::exit(EXIT_FAILURE);
        }
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
        if (x.first == "low")
        {
            _range_low_ = x.second;
            condition += 2;
        }
        else if (x.first == "up")
        {
            _range_up_ = x.second;
            condition += 3;
        }
    }

    if (condition == 0)
    {
        _valid_func_ = &Base_Proposer::_valid_dummy_;
    }
    else if (condition == 2)
    {
        _valid_func_ = &Base_Proposer::_valid_range_low_;
    }
    else if (condition == 3)
    {
        _valid_func_ = &Base_Proposer::_valid_range_up_;
    }
    else if (condition == 5)
    {
        if (_range_low_ >= _range_up_)
        {
            std::cerr << "ERROR: Provided an inconsistent range: {" << _range_low_ << ", " << _range_up_ << std::endl;
            std::cerr << "Exit...\n"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
        else
        {
            _valid_func_ = &Base_Proposer::_valid_range_both_;
        }
    }
    else
    {
        std::cerr << "ERROR: Provided multiple definition for the range." << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        exit(EXIT_FAILURE);
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
        std::cerr << "ERROR: Unknown metropolis strategy " << name << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

std::string Base_Proposer::display(char delim, char dict_delim)
{
    std::string output = distr_ptr->display(delim, dict_delim);
    output += delim;
    if (_range_low_ > -INFINITY)
    {
        std::string aux = "low";
        output += aux + dict_delim + std::to_string(_range_low_);
        if (_range_up_ < INFINITY)
        {
            output += dict_delim;
        }
    }
    if (_range_up_ < INFINITY)
    {
        std::string aux = "up";
        output += aux + dict_delim + std::to_string(_range_up_);
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
        std::cerr << "ERROR: Parameter `" << k << "` is not defined." << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        std::exit(EXIT_FAILURE);
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
        std::cerr << "ERROR: Parameter `" << k << "` has multiple definitions." << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // update the lists of parameters to be inferred or not
    if (new_param.get_infer())
    {
        param_to_be_inferred.push_back(k);
    }
    else
    {
        param_not_to_be_inferred.push_back(k);
    }
}

/* display the user model parameters */
void Parameter_Dictionary::display_parameters(char delim, char dict_delim)
{
    for (std::string k : param_to_be_inferred)
    {
        std::cout << k << delim << params.at(k).display(delim, dict_delim) << endl;
    }
    for (std::string k : param_not_to_be_inferred)
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

/* FIXME */
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
void Parameter_Dictionary::accept_move(const std::string k)
{
    params.at(k).update_accept_rate(true);
}

/* move is refused: restore old value in the buffer */
void Parameter_Dictionary::refuse_move(const std::string k)
{
    set_value(k, buffer_value);
    params.at(k).update_accept_rate(false);
}