/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "tools/atomic.h"
#include "tools/loglike.h"
#include "tools/nan.h"

#include "mcmc_base.h"
#include "setup.h"
#include "data_wide.h"
#include "data_long.h"
#include "parameter.h"
#include "group_rate.h"
#include "contact_rate.h"
#include "infec_profile.h"
#include "augmented.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  *** MCMC : class ***  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

template <class Model>
class MCMC : public MCMCBase
{
public:
	///////////////////////////
	//////  CONSTRUCTOR  //////
	///////////////////////////

	Model model;
	ContactRate contact_rate; // UPDATE: move to base

	MCMC()
		: model(),
		  contact_rate(
			  [this](Household &h, const Individual &k, const Individual &i) -> std::string
			  { return model.init_contact_rates(h, k, i); },
			  &data_wide, &param)
	{
		init_registry(
			model.get_DataKey(),
			model.get_HHsizeKey(),
			model.get_HHsizevecKey(),
			model.get_HHdoublevecKey(),
			model.get_HHdoublemtxKey(),
			model.get_LongKey(),
			model.get_ParamKey(),
			model.get_GroupRateKey(),
			model.get_AugmKey());
		model.bind_solver(&param, &data_wide, &setup, &infec_profile, &augm);
		model.init_advanced_setup();
	}

	//////////////////////////
	//////  EVALUATION  //////
	//////////////////////////

	/// Initialize model and evaluate log-likelihood.
	double evaluate_log_lik(std::string_view house_id = "")
	{
		setup.n_chains = 0; // use 0 chains as a state variable
		init_evaluate_();	// initialize for evaluation

		// mode: single house evaluation
		if (!hsp::na::is(house_id))
		{
			atomic_error.restore(); // reset atomic error
			infec_profile.update(); // eventually recompute infectivity profile
			double output = atomic_materialize_house_log_lik_(data_wide.at(house_id));
			atomic_error.throw_exception(); // catch and throw if previous failed
			return output;
		}

		// mode: total evaluation
		return update_total_log_lik_();
	}

	///////////////////////////////////
	//////  METROPOLIS-HASTINGS  //////
	///////////////////////////////////

	/// Performs the inference using Metropolis-Hastings.
	void run_inference()
	{
		init_evaluate_(false);			 // initialize for evaluation
		init_monte_carlo_();			 // initialize mcmc sampling
		setup.prog_declare();			 // declare current progress bar
		setup.prog_open_for_inference(); // open progress bar
		run_monte_carlo_();				 // run Monte Carlo
		setup.prog_close();				 // close progress bar
	}

	/*
	/// Trapezoidal rule to estimate logarithm of marginal likelihood (evidence) with power posteriors.
	///
	/// Ref
	/// -----
	/// Friel, N., Hurn, M. & Wyse, J.
	/// Improving power posterior estimation of statistical evidence.
	/// Stat Comput 24, 709–723 (2014).
	std::vector<std::vector<double>> power_posteriors(std::string_view scheme, unsigned long N);
	*/

private:
	//////////////////////////
	//////  EVALUATION  //////
	//////////////////////////

	/// Assert setup definitions, then proceed with initialization for evaluation.
	void init_evaluate_(bool is_reset = true)
	{
		data_wide.assert_non_empty();			   // check if data is empty
		data_wide.assert_not_all_missing_status(); // warn if all missing status
		augm.init();							   // initialize individual augmented relations

		if (is_reset)
		{
			param.user_reset(setup.prng()); // initialize non constant parameters
			augm.user_reset(setup.prng());	// initialize non constant augmented
		}

		model.init_model_for_evaluation(); // user initialization of model for inference

		// init cache dependencies for all households and individuals
		data_wide.allocate_all_cache();
		for (size_t R = 0; R < data_wide.restricted_size(); ++R)
		{
			Household &house = data_wide.restricted(R);
			model.init_house_for_evaluation(house);
		}
	}

	/// Compute new log-likelihood, store old value in buffer. Same for households.
	double update_total_log_lik_(bool is_cache_depen = true)
	{
		double new_ll = 0.0;

		// buffer infectivity profile if lookup tables are active
		if (is_cache_depen)
		{
			data_wide.buffer_all_cache();
			infec_profile.update();
		}

		atomic_error.restore();
#pragma omp parallel for reduction(+ : new_ll)
		for (size_t R = 0; R < data_wide.restricted_size(); ++R)
		{
			if (atomic_error.is_pending())
				continue;

			Household &house = data_wide.restricted(R);
			new_ll += data_wide.multiplicity(R) *
					  atomic_materialize_house_log_lik_(house, is_cache_depen);
		}
		atomic_error.throw_exception();

		return new_ll;
	}

	///////////////////////////////////
	//////  METROPOLIS-HASTINGS  //////
	///////////////////////////////////

	/// Markov chain core, used for inference and power posteriors.
	/// NOTE: sweep is the number of previously completed thin sweeps.
	void run_monte_carlo_()
	{
		init_posterior_sample_(); // initialize posterior sample

		// run multiple MCMC chains&
		for (setup.chain = 1; setup.chain <= setup.n_chains; ++setup.chain)
		{
			try
			{
				init_monte_carlo_chain_();
				tot_log_lik.set(update_total_log_lik_()); // init cached, set total log-likelihood

				// assert that initial log-likelihood value is finite
				if (tot_log_lik.is_inf())
					data_wide.assert_finite_log_lik();
			}
			catch (const std::exception &e)
			{
				throw std::runtime_error(
					"While initializing Markov chain:\n" + std::string(e.what()));
			}

			// Finite adaptative iterations
			for (unsigned long i = 0; i < setup.n_adapt; i += setup.n_thin, setup.prog_update())
				for (unsigned long t = 1; t <= setup.n_thin; ++t)
				{
					move_all_();
					param.adapt_proposal_of_sampled(i + t, setup.n_adapt);
				}

			// Burn-in iterations
			for (unsigned long i = 0; i < setup.n_burnin; i += setup.n_thin, setup.prog_update())
				for (unsigned long t = 1; t <= setup.n_thin; ++t)
					move_all_();

			// Posterior sampling
			for (unsigned long i = setup.n_burnin; i < setup.n_steps; i += setup.n_thin, setup.prog_update())
			{
				for (unsigned long t = 1; t <= setup.n_thin; ++t)
					move_all_();
				store_posterior_sample_();
			}
		}
	}

	/// Update and return contribution of given household to full logliklehood.
	double atomic_materialize_house_log_lik_(Household &house, bool is_cache_depen = true)
	{
		double output = 0;
		try
		{
			// buffer all household cache dependencies
			if (is_cache_depen)
				model.update_house_cache_dependecies(house);

			// compute and set household observation log-likelihood
			output = model.compute_house_obs_log_lik(house);
			house.obs_log_lik_.set(output);

			// compute and set individual augmented log-likelihood
			if (augm.is_any_model())
				for (size_t i : house.infected_hhindexes())
					output += augm.materialize_augm_model_log_lik(house.member(i));
		}
		catch (const std::exception &e)
		{
			atomic_error.set_if_first("While computing log-likelihood for " +
									  house.err_name() + ":\n" + std::string(e.what()));
		}

		return output;
	}

	/// Move all model parameters and augmented variables.
	void move_all_()
	{
		// move augmented variables to sample, if any
		if (augm.is_any_sample())
			move_all_augm_();

		// move all standard model parameters
		for (size_t p : param.list_adaptable_sample_standard)
			param.at(p).update_accept_rate(move_single_parameter_(p, false));
		for (size_t p : param.list_non_adapt_sample_standard)
			move_single_parameter_(p, false);

		// move parameters with cached dependecies
		for (size_t p : param.list_adaptable_sample_cached)
			param.at(p).update_accept_rate(move_single_parameter_(p, true));
		for (size_t p : param.list_non_adapt_sample_cached)
			move_single_parameter_(p, true);
	}

	/// Attempt to move a single parameter for the given inverse temperature.
	bool move_single_parameter_(size_t p, bool is_cache_depen)
	{
		// proceed only if the new value is consistent
		// otherwise automatically refuse move
		if (param.propose_valid_move(p, setup.prng()))
		{
			// if parameter is cache dependent, update cache
			tot_log_lik.set(update_total_log_lik_(is_cache_depen));

			// perform Metropolis-Hastings to accept or not move
			const double D_ll = tot_log_lik.logratio();
			const double move_logratio = param.at(p).get_move_logratio();

			if (setup.exponential() > -(setup.inv_T * D_ll + move_logratio))
				return true; // accepted move
			recover_parameter_(p, is_cache_depen);
		}
		return false; // refused move
	}

	/// Move augmented variables of all infected individuals.
	void move_all_augm_()
	{
		double delta_ll = 0.;
		atomic_error.restore();

#pragma omp parallel for reduction(+ : delta_ll)
		for (size_t H = 0; H < data_wide.report_house.total; ++H)
		{
			if (atomic_error.is_pending())
				continue;

			delta_ll += atomic_move_augm_single_house_(H);
		}
		atomic_error.throw_exception();
		tot_log_lik.increment(delta_ll);
	}

	/// Move augmented variables of infectious members inside household.
	/// The move is performed one member at a time, one augmented variable at a time.
	double atomic_move_augm_single_house_(size_t H)
	{
		Household &house = data_wide.at(H);
		double delta_ll = 0.;

		try
		{
			// loop over augmented variables of infectious
			for (size_t i : house.infected_hhindexes())
			{
				Individual &member_i = house.member(i);
				for (size_t a : augm.list_sample(member_i))
					delta_ll += metropolis_augm_single_move_(house, member_i, a);
			}
		}
		catch (const std::exception &e)
		{
			atomic_error.set_if_first("While sampling augmented variables for " +
									  house.err_name() + ":\n" + std::string(e.what()));
		}

		return delta_ll;
	}

	/// Move a single augmented variable.
	double metropolis_augm_single_move_(Household &house, Individual &member, size_t a)
	{
		Variable &xaugm = member.get_augm_variable(a);
		xaugm.propose(setup.prng()); // propose new value and buffer old

		// proceed only if the new value is consistent
		if (xaugm.is_valid_move())
		{
			size_t i = member.hhindex();
			// update dependencies to perform computations
			data_wide.buffer_member_cache(house, i);
			model.update_member_cache_dependecies(house, member);

			// Metropolis-Hastings
			double delta_ll = augm.materialize_augm_model_logratio(member, a);
			// compute and set household observation log-likelihood
			house.obs_log_lik_.set(model.compute_house_obs_log_lik(house));
			delta_ll += house.obs_log_lik_.logratio();

			if (setup.exponential() > xaugm.get_propose_logratio() - setup.inv_T * delta_ll)
				return delta_ll;

			// else recover cached values
			house.obs_log_lik_.recover();
			member.augm_log_lik_.recover();
			data_wide.recover_member_buffer(house, i);
		}
		// else recover buffered value
		xaugm.recover();

		return 0.;
	}

	/* --FUTURE-LOG-EVIDENCE--
	/// Power posteriors at 0 inverse temperature.
	void power_post_zero_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var);
	/// Power posteriors with power fractions scheme for the thermodynamic ladder.
	void power_post_fractions_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var);
	void power_post_adaptive_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var);
	*/
};
