/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * housepi.cpp
 *
 */

#include "tools/error.h" // UPDATE: required by contact rates
#include "tools/rates.h"
#include "tools/params.h"
#include "variable.h"
#include "household.h"
#include "mcmc_impl.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  AugmentedMap : setter  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Assign definitions for the augmented variable model.
/// NOTE: a Variable must already be allocated with clear_augm.
void AugmentedMap::init_augm(const std::string &a,
							 double value,
							 bool is_const,
							 const DefinitionParser &func_name,
							 bool is_prim,
							 const DefinitionParser &model_name)
{
	Variable &this_augm = safe_at(augm_map_, a, "augmented variable");
	this_augm.init(value, is_const, func_name, is_prim);
	this_augm.init_distr(model_name);
}

/// A string containing a summary of the augmented map.
std::string AugmentedMap::str_augm_summary() const
{
	std::string output;
	for (const auto &pair : augm_map_)
	{
		output += "- " + pair.first + " : ";
		if (pair.second.is_function() && !pair.second.get_function()->is_implicit_callable())
			output += pair.second.get_function()->str_summary();
		else
			output += pair.second.str_summary();
		output += "\n";
	}
	return output;
}

/// Associate functional methods for augmented model pdfs.
void AugmentedMap::init_augm_model_logratios(const u_dict<str_vec> &models_ordered)
{
	// >> clean pre-existing
	for (auto &pair : model_logratios_)
	{
		pair.second.clear();
		pair.second.shrink_to_fit();
	}
	model_logratios_.clear();
	model_logratios_.rehash(0);
	// <<

	// associate models if applicable
	for (const auto &pair : models_ordered)
	{
		const std::string &a = pair.first;
		model_logratios_[a] = {};

		// repeat for all dependent models
		for (const std::string &b : pair.second)
		{
			const Variable &xaugm = augm_map_[b];
			model_logratios_[a].emplace_back([&xaugm]()
											 { return xaugm.get_distr_logratio(); });
		}
	}
}

/// Value of loglikelihood associated to augmented variable models.
/// Store the loglikelihood of individual augmented variables.
double AugmentedMap::compute_and_set_augm_model_loglike()
{
	buffer_ = loglike_; // store old value in buffer
	loglike_ = 0.0;

	for (const auto &pair : model_logratios_)
		loglike_ += augm_map_[pair.first].get_distr_logpdf();

	if (hsp::na::is(loglike_))
		throw std::runtime_error(
			"Augmented loglikelihood returned NaN for " + err_name() + ".");

	return loglike_;
}

/// Change in loglikelihood associated to given augmented.
double AugmentedMap::compute_and_set_augm_model_logratio(const std::string &a)
{
	buffer_ = loglike_; // store old value in buffer

	double delta_ll = 0.;
	for (const std::function<double()> &f : model_logratios_[a])
		delta_ll += f();
	loglike_ += delta_ll;

	if (hsp::na::is(loglike_))
		throw std::runtime_error(
			"Augmented loglikelihood returned NaN for " + err_name() +
			" after change of " + augm_map_[a].err_name() + ".");

	return delta_ll;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Individual : constructor  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Individual constructor.
Individual::Individual(size_t index,
					   size_t hhindex,
					   const std::string &id,
					   int status,
					   const u_dict<double> &data)
{
	index_ = index;
	hhindex_ = hhindex;
	id_ = id;
	status_ = status;
	data_map_ = data;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Individual : setter  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize group rate according to given attribute.
void Individual::init_group_rate(const std::string &g, const Attribute &attr)
{
	group_rate_map_[g].transfer(attr);
	group_rate_map_[g].name(g);
}

/// Initialize group rate to reference.
void Individual::init_group_rate_to_reference(const std::string &g)
{
	group_rate_map_[g] = hsp::rates::VALUE_REFERENCE;
	group_rate_map_[g].name(g);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Household : constructor  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Household constructor.
Household::Household(const std::string &house_id, size_t H)
{
	house_id_ = house_id;
	index_ = H;
	house_multiplicity_ = 1;
	multiplicity_class_ = H;
	house_size_ = 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Household : setter  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Assign a multiplicity to the household.
void Household::set_multiplicity(unsigned long m, size_t class_idx)
{
	house_multiplicity_ = m;
	multiplicity_class_ = class_idx;
}

/// Link an individual to the household.
void Household::init_member(size_t index,
							const std::string &id,
							int status,
							const u_dict<double> &data)
{
	// emplace new individual
	members_.emplace_back(index, house_size_, id, status, data);

	// update status lists within the household
	if (members_.back().is_infected_status())
		infected_list_.push_back(house_size_);
	else if (members_.back().is_healthy_status())
		healthy_list_.push_back(house_size_);
	else
		missing_list_.push_back(house_size_);

	++house_size_;
}

/// Change_status of individual i to new status (without controls).
///	NOTE : this introduces augmented variables -> it requires RJMCMC in case of inference.
void Household::set_status_(Individual &member, int new_status)
{
	if (new_status != member.status())
	{
		const size_t i = static_cast<size_t>(member); // get within house index

		if (member.is_infected_status())
			infected_list_.erase(std::find(infected_list_.begin(), infected_list_.end(), i));
		else if (member.is_healthy_status())
			healthy_list_.erase(std::find(healthy_list_.begin(), healthy_list_.end(), i));
		else
			missing_list_.erase(std::find(missing_list_.begin(), missing_list_.end(), i));

		member.status_ = new_status; // change member status
		if (member.is_infected_status())
			infected_list_.push_back(i);
		else if (member.is_healthy_status())
			healthy_list_.push_back(i);
		else
			missing_list_.push_back(i);
	}
}

/// Associate a pointer to parameter value for the contact from->to.
void Household::init_contact_rates(size_t from, size_t to, observ_dbl_ptr x)
{
	// initialize contact rate matrix if necessary
	if (contact_rates_.empty())
		contact_rates_.resize(size(), std::vector<observ_dbl_ptr>(size(), &hsp::rates::VALUE_REFERENCE));

	contact_rates_[from][to] = x;
}

/// Store the loglikelihood of household observation.
void Household::set_obs_loglike(double x)
{
	if (hsp::na::is(x))
		throw std::runtime_error(
			"Observation loglikelihood returned NaN.");

	obs_buffer_ = obs_loglike_;
	obs_loglike_ = x;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Household : cache attributes  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Recover previous values of the loglikelihood and cached depedencies.
void Household::recover_buffer_values(bool is_cache_depen)
{
	obs_loglike_ = obs_buffer_;
	if (is_cache_depen)
	{
		// swap buffer into  cache attributes
		cache_int_vec_.swap(buffer_cache_int_vec_);
		cache_size_t_vec_.swap(buffer_cache_size_t_vec_);
		cache_dbl_vec_.swap(buffer_cache_dbl_vec_);
		cache_dbl_mtx_.swap(buffer_cache_dbl_mtx_);
	}
}

/// Deep copies of cache attributes in the buffer
void Household::buffer_all_cache()
{
	buffer_cache_int_vec_ = cache_int_vec_;
	buffer_cache_size_t_vec_ = cache_size_t_vec_;
	buffer_cache_dbl_vec_ = cache_dbl_vec_;
	buffer_cache_dbl_mtx_ = cache_dbl_mtx_;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : constructor  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

MCMC::MCMC()
    : pImpl_(std::make_unique<Impl>(
          [this]()
          { this->init_model_for_inference(); },
          [this]()
          { this->init_model_for_simulation(); },
          [this](Household &h) -> double
          { return this->compute_household_obs_loglike(h); },
          [this](Household &h)
          { return this->simulate_single_household(h); },
          [this](Household &h, Individual &i)
          { this->update_cache_dependecies(h, i); })) {}
MCMC::~MCMC() = default;

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////  MCMC : data  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

const ReportIndividual &MCMC::get_report_indiv() const { return pImpl_->dataset.report_indiv; }
const ReportHousehold &MCMC::get_report_house() const { return pImpl_->dataset.report_house; }
const ReportHouseSizes &MCMC::get_report_hh_size() const { return pImpl_->dataset.report_hh_size; }
void MCMC::assert_no_missing_status() const
{
    if (pImpl_->dataset.report_indiv.missing > 0)
        hsp::msg::error(
            "Individual with missing status are incompatible with current model.");
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  MCMC : parameters  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double MCMC::param(const std::string &p) const { return pImpl_->param.get(p); }

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : augmented variables  /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::init_missing_augm(Individual &member)
{
    pImpl_->augm.build_missing(member);
    pImpl_->augm.user_reset(member);
}
void MCMC::clear_all_augm(Individual &member) { pImpl_->augm.clear_all(member); }
void MCMC::request_augm_discrete(const std::string &x) { pImpl_->augm.request_augm_discrete(x); }

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  MCMC : cache attributes  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::cache_house_integer_vector(const std::string &e) { pImpl_->dataset.init_cache_integer_vector(e); }
void MCMC::cache_house_size_t_vector(const std::string &e) { pImpl_->dataset.init_cache_size_t_vector(e); }
void MCMC::cache_house_double_vector(const std::string &e) { pImpl_->dataset.init_cache_double_vector(e); }
void MCMC::cache_house_double_matrix(const std::string &e) { pImpl_->dataset.init_cache_double_matrix(e); }

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : infectivity profile  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::init_infec_profile_template(u_dict<std::string> x) { pImpl_->infec_profile.init_template(x); }
double MCMC::infec_profile_inst(const Individual &i, double t) const { return pImpl_->infec_profile.inst(&i, t); }
double MCMC::infec_profile_cumul(const Individual &i, double t) const { return pImpl_->infec_profile.cumul(&i, t); }

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : general methods  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double MCMC::uniform() { return pImpl_->setup.unif(); }
double MCMC::evaluate() { return pImpl_->evaluate(); }
void MCMC::infer() { pImpl_->infer(); }
void MCMC::clear_memory() { pImpl_->clear_memory(); }

#ifdef _HOUSEPI_RCPP
Rcpp::DataFrame MCMC::data_wide_for_simulation() const
{
    return pImpl_->dataset.rcpp_export(pImpl_->dataset.headers() - pImpl_->augm.list_export, false);
}
Rcpp::List MCMC::simulate(size_t n_rep) { return pImpl_->simulate(n_rep); }
Rcpp::List MCMC::logevidence(const std::string &scheme, unsigned long n_rungs) { return pImpl_->logevidence(scheme, n_rungs); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : imports from Rcpp ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::import_setup(SEXP input) { pImpl_->setup.rcpp_import(input); }
void MCMC::import_data_wide(SEXP input) { pImpl_->dataset.rcpp_import(input); }
void MCMC::import_multiplicity(SEXP input) { pImpl_->dataset.rcpp_import_multiplicity(input); }
void MCMC::import_data_long(SEXP input) { pImpl_->datalong.rcpp_import(input); }
void MCMC::import_group_rates(SEXP input) { pImpl_->group_rates.rcpp_import(input); }
void MCMC::import_param(SEXP input) { pImpl_->param.rcpp_import(input); }
void MCMC::import_infec_profile(SEXP input) { pImpl_->infec_profile.rcpp_import(input); }
void MCMC::import_augm(SEXP input) { pImpl_->augm.rcpp_import(input); }

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : exports to Rcpp ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Rcpp::DataFrame MCMC::export_data_wide() const { return pImpl_->dataset.rcpp_export(pImpl_->dataset.headers(), true); }
Rcpp::DataFrame MCMC::export_data_long() const { return pImpl_->datalong.rcpp_export(); }
Rcpp::DataFrame MCMC::summarise_cases() const { return pImpl_->dataset.rcpp_summarise_cases(); }
Rcpp::DataFrame MCMC::summarise_house_sizes() const { return pImpl_->dataset.rcpp_summarise_house_sizes(); }
Rcpp::List MCMC::export_setup() const { return pImpl_->setup.rcpp_export(); }
Rcpp::CharacterVector MCMC::export_param_list_const() const { return Rcpp::wrap(pImpl_->param.list_const); }
Rcpp::CharacterVector MCMC::export_param_list_sample_standard() const { return Rcpp::wrap(pImpl_->param.list_sample_standard); }
Rcpp::CharacterVector MCMC::export_param_list_sample_cache() const { return Rcpp::wrap(pImpl_->param.list_sample_cached); }
Rcpp::List MCMC::export_infec_profile() const { return pImpl_->infec_profile.rcpp_export(); }
Rcpp::CharacterVector MCMC::export_augm_list_draws() const { return Rcpp::wrap(pImpl_->augm.list_draws); }
Rcpp::CharacterVector MCMC::export_augm_list_full() const { return Rcpp::wrap(pImpl_->augm.list_full); }
Rcpp::CharacterVector MCMC::export_augm_list_export() const { return Rcpp::wrap(pImpl_->augm.list_export); }
Rcpp::CharacterVector MCMC::export_augm_list_discrete() const { return Rcpp::wrap(pImpl_->augm.user_list_discrete); }
Rcpp::DataFrame MCMC::export_augmented_data_wide() const { return pImpl_->augm.rcpp_export_data_wide(); }
Rcpp::NumericMatrix MCMC::export_draws() const { return pImpl_->rcpp_export_draws(); }
Rcpp::NumericMatrix MCMC::export_house_loglike_draws() const { return pImpl_->dataset.rcpp_export_draws(); }
Rcpp::List MCMC::export_house_snapshot(const std::string &house_id) const { return pImpl_->rcpp_export_house_snapshot(house_id); }

#else
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  MCMC : C++ prints  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::print_simulated_data_wide(const std::string &f) { pImpl_->print_simulated_data_wide(f); }
double MCMC::logevidence(const std::string &scheme, unsigned long n_rungs) { return pImpl_->logevidence(scheme, n_rungs); }

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  MCMC : read  //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::read_setup(const std::string &f) { pImpl_->setup.read(f); }
void MCMC::read_data_wide(const std::string &f) { pImpl_->dataset.read(f); }
void MCMC::read_multiplicity(const std::string &f) { pImpl_->dataset.read_multiplicity(f); }
void MCMC::read_data_long(const std::string &f) { pImpl_->datalong.read(f); }
void MCMC::read_group_rates(const std::string &f) { pImpl_->group_rates.read(f); }
void MCMC::read_param(const std::string &f) { pImpl_->param.read(f); }
void MCMC::read_augm(const std::string &f) { pImpl_->augm.read(f); }
void MCMC::read_infec_profile(const std::string &f) { pImpl_->infec_profile.read(f); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : show, C++ print  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::print_infer_state(const std::string &f) const { pImpl_->print_infer_state(f); }
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : contact rates ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Build all group rates and contact rates and set rate per group. // UPDATE
void MCMC::build_all_rates()
{
    try
    {
        for (Household &house : pImpl_->dataset)
            init_contact_rates_(house);
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While initializing contact rates:\n" +
                        std::string(e.what()));
    }
}

/// Set contact rate according to user assignement. // UPDATE
void MCMC::init_contact_rates_(Household &house)
{
    // loop over all couples
    for (const Individual &member_k : house)     // from : member_k
        for (const Individual &member_i : house) // -> to : member_i
            if (member_k != member_i)
            {
                _pair_contact_rate = hsp::na::STRING;
                init_household_contact_rates(house, member_k, member_i);

                if (!hsp::na::is(_pair_contact_rate))
                {
                    observ_dbl_ptr y = &hsp::rates::VALUE_REFERENCE; // default reference rate

                    if (_pair_contact_rate == hsp::rates::MISSING)
                        throw std::invalid_argument(
                            "Missing contact rate for individuals \"" +
                            member_k.id() + "\" and \"" + member_i.id() +
                            "\" (" + house.err_name() + ").");

                    else if (_pair_contact_rate != hsp::rates::REFERENCE)
                        y = pImpl_->param.get_ptr(_pair_contact_rate); // parameter rate

                    house.init_contact_rates(static_cast<size_t>(member_k), static_cast<size_t>(member_i), y);
                }
            }
}
