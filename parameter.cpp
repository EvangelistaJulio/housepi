// parameter.cpp //

#include "parameter.hpp"

using namespace std;

//----------------------------------------------------------------------
// Strategy : functional methods
//----------------------------------------------------------------------

// --- prior --- //
// FIXME : missing

/* adapt : increase spread argument accept_rate is small and viceversa */
double Metropolis_Hastings::_adapter_(int step, int num_adapt_iter)
{
    return 1. + (1. - step / num_adapt_iter) * (DEFAULT_ACCEPT_RATE - get_accept_rate());
}

/* update accept rate */
void Metropolis_Hastings::update_accept_rate(bool is_accepted, int step)
{
    accept_rate = ((is_accepted + accept_rate * (step - 1)) / step);
}

/* reset moves */
void Metropolis_Hastings::reset_moves()
{
    num_propos = 0;
    num_accept = 0;
    accept_rate = 0.0;
}

/* set range */
void Metropolis_Hastings::set_range(kwargs new_range)
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
        _valid_func_ = &Metropolis_Hastings::_valid_dummy_;
    }
    else if (condition == 1)
    {
        _valid_func_ = &Metropolis_Hastings::_valid_range_low_;
    }
    else if (condition == 2)
    {
        _valid_func_ = &Metropolis_Hastings::_valid_range_up_;
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
            _valid_func_ = &Metropolis_Hastings::_valid_range_both_;
        }
    }
}

/////////////////////////////////////////////////////
//  Parameter : class definition
/////////////////////////////////////////////////////

Parameter::Parameter(
    double * new_value_ptr,
    std::string strat_name,
    kwargs strat_args,
    kwargs value_range,
    std::string prior_name,
    kwargs prior_args) : _distribution_(), _proposer_()
{
    value_ptr = new_value_ptr;
    value_0 = get_value();

    // --- metropolis strategy --- //

    if (strat_name == "_none_")
    {
        infer = false;
    }
    else if (strat_name == "normal")
    {
        infer = true;
        _proposer_ = Normal_Proposal(strat_args);
    }
    else if (strat_name == "lognorm")
    {
        infer = true;
        _proposer_ = Lognorm_Proposal(strat_args);
    }
    else
    {
        std::cerr << "ERROR: Unknown metropolis strategy " << strat_name << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    _proposer_.set_range(value_range);

    // --- prior function --- //
    if (infer == true && prior_name != "_none_")
    {
        if (prior_name == "unif")
        {
            _distribution_ = Uniform_Distr(prior_args);
        }
        else if (prior_name == "normal")
        {
            _distribution_ = Normal_Distr(prior_args);
        }
        else if (prior_name == "lognorm")
        {
            _distribution_ = Lognorm_Distr(prior_args);
        }
        else if (prior_name == "gamma")
        {
            _distribution_ = Gamma_Distr(prior_args);
        }
        else if (prior_name == "exp")
        {
            _distribution_ = Exponential_Distr(prior_args);
        }
        else
        {
            std::cerr << "\nERROR: Unknown distribution " << prior_name << std::endl;
            std::cerr << "Exit...\n"
                      << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        prior_name = "_none_";
        prior_args = {};
    }
}

///////////////////////////////////////////////////////////
//////// Parameter_Dictionary : methods ///////////////////
///////////////////////////////////////////////////////////

/* add_param */
void Parameter_Dictionary::add_param(const std::string k, Parameter new_param)
{
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

    if (new_param.get_infer() == true)
    {
        to_be_inferred.push_back(k);
    }
}

/* print */
void Parameter_Dictionary::print()
{
    std::cout << "_NAME_" << TAB << "_VALUE_" << TAB << "_INFER_" << endl;
    for (const auto &entry : params)
    {
        Parameter param = entry.second;
        std::cout << entry.first << TAB << param.get_value();
        std::cout << TAB << (param.get_infer() > 0 ? "true" : "false") << endl;
    }
}