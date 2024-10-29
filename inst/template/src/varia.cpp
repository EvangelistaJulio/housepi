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
                             double rule_start_value,
                             std::vector<int> args_to_adapt) : distr_ptr(distr)
{
    _args_to_adapt_ = args_to_adapt;

    if (rule_start_value == INF_DOUBLE)
    {
        _is_start_from_zero_ = false;
    }
    else
    {
        _is_start_from_zero_ = true;
        set_rule_start_value(rule_start_value);
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

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Uniform_Proposer : class  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Uniform_Proposer : public Base_Proposer
{
public:
	//=====Constructors=====//
	Uniform_Proposer(
		dictargs rule_args = {{"min", (double)-0.5}, {"max", (double)0.5}},
		double rule_from = INF_DOUBLE) : Base_Proposer(distribution_ptr("unif", rule_args), rule_from, {"min", "max"})
	{
		_propose_method_ = &Uniform_Proposer::_propose_sum_to_previous_;
		_logratio_method_ = &Uniform_Proposer::_logratio_symm_;
	}
	~Uniform_Proposer() {}
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Normal_Proposer : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Normal_Proposer : public Base_Proposer
{
public:
	//=====Constructor=====//
	Normal_Proposer(
		dictargs rule_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
		double rule_from = INF_DOUBLE) : Base_Proposer(distribution_ptr("norm", rule_args), rule_from, {1})
	{
		_propose_method_ = &Normal_Proposer::_propose_sum_to_previous_;
		_logratio_method_ = is_rule_start_from_zero() ? &Normal_Proposer::_logratio_sum_to_zero_ : &Normal_Proposer::_logratio_symm_;
	}
	~Normal_Proposer() {}
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Lognorm_Proposer : class  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Lognorm_Proposer : public Base_Proposer
{
public:
	//=====Constructor=====//
	Lognorm_Proposer(
		dictargs rule_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
		double rule_from = INF_DOUBLE) : Base_Proposer(distribution_ptr("lnorm", rule_args), rule_from, {1})
	{
		_propose_method_ = &Lognorm_Proposer::_propose_times_previous_;
		_logratio_method_ = is_rule_start_from_zero() ? &Lognorm_Proposer::_logratio_times_zero_ : &Lognorm_Proposer::_logratio_lognorm_symm_;
	}
	~Lognorm_Proposer() {}
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Gibbs_Proposer : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Gibbs_Proposer : public Base_Proposer
{
public:
	//=====Constructor=====//
	Gibbs_Proposer(Base_Distribution *distribution) : Base_Proposer(distribution)
	{
		_propose_method_ = &Gibbs_Proposer::_propose_draw_rand_;
		_logratio_method_ = &Gibbs_Proposer::_logratio_pdf_;
	}
	~Gibbs_Proposer() {}
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

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
//////////////////////////////////  Base_Variable : methods  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Base_Variable::assign_proposer(const std::string &rule_name,
                                    dictargs rule_args,
                                    double rule_from)
{
    // check that rule args have no mismatch
    if (rule_args.find(IO_MISMATCH) != rule_args.end())
    {
        tools::error("Mismatching definitions in rule arguments of '" + name + "'.");
    }

    if (tools::is_na(rule_name))
    {
        // immutable variable if no rule is defined
        infer = false;
        _propose_method_ = &Base_Variable::_propose_error_;
    }
    else
    {
        infer = true;

        // define proposer rule with functional programming
        if (rule_name == "gibbs")
        {
            prop_ptr = new Gibbs_Proposer(distr_ptr);
            _propose_method_ = &Base_Variable::_propose_gibbs_sample_;
            _move_logratio_method_ = &Base_Variable::_move_logratio_gibbs_;
        }
        else
        {
            if ((rule_from != INF_DOUBLE) && (!distr_ptr->valid(rule_from)))
            {
                std::string min_val = distr_ptr->show_range_min();
                std::string max_val = distr_ptr->show_range_max();
                std::ostringstream oss;
                oss << "Initial value " << rule_from << " outside of value range"
                    << " {min:" << (tools::is_na(min_val) ? std::numeric_limits<double>::quiet_NaN() : tools::stod(min_val))
                    << ", max:" << (tools::is_na(max_val) ? std::numeric_limits<double>::quiet_NaN() : tools::stod(max_val))
                    << "}."
                    << " If intentional, ignore this message.";
                tools::warning(oss.str());
            }

            // assigne the pointer to the proposer
            prop_ptr = proposer_switchboard(rule_name, rule_args, rule_from);
            _propose_method_ = is_from_zero() ? &Base_Variable::_propose_from_zero_ : &Base_Variable::_propose_from_previous_;
        }
    }
}
