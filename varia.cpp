// varia.cpp //

#include "varia.h"

using namespace std;

//----------------------------------------------------------------------
// Strategy : functional methods
//----------------------------------------------------------------------

Base_Proposer::Base_Proposer(
    Distribution_1D *distr,
    double strat_zero,
    kwargs value_range,
    std::vector<std::string> args_to_adapt) : distr_ptr(distr)
{
    this->args_to_adapt = args_to_adapt;
    set_range(value_range);
    if (strat_zero == INFINITY)
    {
        this->zero_start = false;
    }
    else
    {
        if (valid(strat_zero))
        {
            this->zero_start = true;
            set_strat_zero(strat_zero);
        }
        else
        {
            std::cerr << "ERROR: initial value " << strat_zero << " and value range {";
            if (value_range.find("low") != value_range.end())
            {
                std::cerr << "low:" << value_range.at("low") << ", ";
            }
            if (value_range.find("up") != value_range.end())
            {
                std::cerr << "up:" << value_range.at("up");
            }
            std::cerr << "} are inconsistent." << std::endl;
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

    for (std::string argID : args_to_adapt)
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

    if (new_range.find("low") != new_range.end())
    {
        _range_low_ = new_range.at("low");
        condition += 1;
    }

    if (new_range.find("up") != new_range.end())
    {
        _range_up_ = new_range.at("up");
        condition += 2;
    }

    if (condition == 0)
    {
        _valid_func_ = &Base_Proposer::_valid_dummy_;
    }
    else if (condition == 1)
    {
        _valid_func_ = &Base_Proposer::_valid_range_low_;
    }
    else if (condition == 2)
    {
        _valid_func_ = &Base_Proposer::_valid_range_up_;
    }
    else
    {
        if (_range_low_ >= _range_up_)
        {
            std::cerr << "ERROR: Provided an inconsistent range: {" << new_range.at("low") << ", " << new_range.at("up") << std::endl;
            std::cerr << "Exit...\n"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
        else
        {
            _valid_func_ = &Base_Proposer::_valid_range_both_;
        }
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

///////////////////////////////////////////////////////////
//////// Parameter_Dictionary : methods ///////////////////
///////////////////////////////////////////////////////////

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
void Parameter_Dictionary::display_parameters(char delim)
{
    std::cout << "-NAME-" << delim << "-INFER-" << delim << "-VALUE-" << endl;
    for (std::string k : param_to_be_inferred)
    {
        std::cout << k << delim << "true" << delim << get_value(k) << endl;
    }
    for (std::string k : param_not_to_be_inferred)
    {
        std::cout << k << delim << "false" << delim << get_value(k) << endl;
    }
}

/* randomly initialize values of the parameters to be inferred */
void Parameter_Dictionary::init_rand_params(std::mt19937_64 &rng)
{
    for (std::string k : param_to_be_inferred)
    {
        set_value(k, params.at(k).prior_rand(rng));
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