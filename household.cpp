#include "household.hpp"

using namespace std;

/////////////////////////////////////////////////////
//  Individual : Class
/////////////////////////////////////////////////////

// ----------------------------
//  Constructors
// ----------------------------

Individual::Individual(
	int indid,
	int status,
	int age,
	int sex,
	bool isolation,
	bool recovered,
	int vacc_status,
	double study_period,
	double vacc_delay,
	double onset_time)
{
	// Data variables
	i_indid = indid;
	i_status = status;
	i_age = age;
	i_sex = sex;
	i_isolation = isolation;
	i_recovered = recovered;
	i_vacc_status = vacc_status;
	i_study_period = study_period;
	i_vacc_delay = vacc_delay;
	i_onset_time = onset_time;

	// Initialize augmented variables
	i_augm_infec_time = onset_time;
	i_augm_onset_time = onset_time;
}

// ----------------------------
//  Methods : Computations
// ----------------------------

// FIXME : these are probabilities, but we want rates.

double Individual::infec_inst(double t_infec, double study_period)
{
	double output = 0.0, tost, t_incub;
	std::vector<double> args = {
		glob_params.value("tost_tau"), glob_params.value("tost_mu"),
		glob_params.value("tost_sigma"), glob_params.value("tost_alpha")};

	if (i_status == SYMPTOMATIC)
	{
		if (i_augm_infec_time < study_period && i_augm_infec_time < t_infec)
		{
			tost = t_infec - i_augm_onset_time;
			t_incub = i_augm_onset_time - i_augm_infec_time;
			output = t_onset_sympt_to_transm_pdf(tost, t_incub, args);
		}
	}

	return output;
}

double Individual::infec_cumul(double t_infec, double study_period)
{
	double output = 0.0, tost, t_incub;
	std::vector<double> args = {
		glob_params.value("tost_tau"), glob_params.value("tost_mu"),
		glob_params.value("tost_sigma"), glob_params.value("tost_alpha")};

	if (i_status == SYMPTOMATIC)
	{
		if (i_augm_infec_time < study_period && i_augm_infec_time < t_infec)
		{
			tost = t_infec - i_augm_onset_time;
			t_incub = i_augm_onset_time - i_augm_infec_time;
			output = t_onset_sympt_to_transm_cdf(tost, t_incub, args);
		}
	}

	return output;
}

// ////////////// AUGMENTED ////////////// //
// UPDATE : generalize augmented

// Initialize augmented variables
void Individual::init_augm_data(std::mt19937_64 &rand_gen)
{
	std::vector<double> args = {glob_params.value("incub_shape"), glob_params.value("incub_scale")};

	if (i_status == NOT_A_CASE)
	{
		i_augm_infec_time = INFINITY;
		i_augm_onset_time = INFINITY;
	}
	else if (i_status == SYMPTOMATIC)
	{
		// FIXME : this initialization looks too arbitrary
		i_augm_infec_time = i_onset_time - gamma_rand(rand_gen, args);
		i_augm_onset_time = i_onset_time + unif_rand(rand_gen);
	}
	else
	{
		std::cerr << "ERROR: status not in the model." << endl;
		std::cerr << "Exit...\n"
				  << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

// Propose new value for symptom onset time
double Individual::propose_augm_onset_time(std::mt19937_64 &rand_gen)
{
	double new_onset_time = 0.0;

	if (i_status == NOT_A_CASE)
	{
		new_onset_time = INFINITY;
	}
	else
	{
		new_onset_time = i_onset_time + unif_rand(rand_gen);
	}

	return new_onset_time;
}
// Propose new value for infection time
double Individual::propose_augm_infec_time(std::mt19937_64 &rand_gen)
{
	double new_infec_time = 0.0, incub_period = 0.0;
	std::vector<double> args = {glob_params.value("incub_shape"), glob_params.value("incub_scale")};

	if (i_status == SYMPTOMATIC)
	{
		incub_period = gamma_rand(rand_gen, args);
		new_infec_time = get_augm_onset_time() - incub_period;
	}
	else
	{ // Resident without infection or with unknown outcome
		new_infec_time = INFINITY;
	}

	return new_infec_time;
}

/////////////////////////////////////////////////////
//  Household : Class
/////////////////////////////////////////////////////

// Add individual to household
void Household::add_individual(Individual this_indiv)
{
	h_residents.push_back(this_indiv);
	h_size++;
	update_h_infec_time();
}

// Initialize augmented variables
void Household::model_init_augm_data(std::mt19937_64 &rand_gen)
{

	// loop over residents within the house
	for (size_t i = 0; i < h_size; i++)
	{
		h_residents[i].init_augm_data(rand_gen);
	}
	update_h_infec_time();
	init_infec_profiles();
}

// Change symptom onset
void Household::set_augm_onset_time(size_t idx, double new_onset_time)
{
	h_residents[idx].set_augm_onset_time(new_onset_time);
	update_infec_profiles(idx);
}

// Change infection time
void Household::update_h_infec_time()
{
	double new_t0 = INFINITY;
	size_t new_case0 = h_size; // NOTE: init to impossible value

	for (size_t i = 0; i < h_size; i++)
	{
		// Find first infection time
		new_t0 = std::min(h_residents[i].get_augm_infec_time(), new_t0);
		// Identify first case
		if (h_residents[i].get_augm_infec_time() == new_t0)
		{
			new_case0 = i;
		}
	}

	h_first_case = new_case0;
	h_infec_time = new_t0;
}

// Change individual infection time
void Household::set_augm_infec_time(size_t idx, double new_infec_time)
{
	h_residents[idx].set_augm_infec_time(new_infec_time);
	update_h_infec_time();
	update_infec_profiles(idx);
}

// ----------------------------
//  Methods : Computations
// ----------------------------

// Propose new value for symptom onset time
double Household::propose_augm_onset_time(size_t idx, std::mt19937_64 &rand_gen)
{
	return h_residents[idx].propose_augm_onset_time(rand_gen);
}

// Propose new value for infection time
double Household::propose_augm_infec_time(size_t idx, std::mt19937_64 &rand_gen)
{
	return h_residents[idx].propose_augm_infec_time(rand_gen);
}

// Compute the risk of infection
void Household::init_infec_profiles()
{

	// initialize 2-dim vectors h_cumul_rate and h_inst_rate
	h_cumul_rate.resize(h_size);
	h_inst_rate.resize(h_size);
	for (size_t i = 0; i < h_size; i++)
	{
		h_cumul_rate[i].resize(h_size);
		h_inst_rate[i].resize(h_size);
	}

	// Compute cumulative and instantaneous person-to-person transmission rates
	for (size_t infector = 0; infector < h_size; infector++)
	{
		// FIXME : is this working ?
		update_infec_profiles(infector);
	}
}

// Update h_cumul_rate and h_inst_rate matrices during data augmentation
void Household::update_infec_profiles(size_t idx)
{
	double t1;
	// NOTE : h_cumul_rate[idx][i] : from `idx` to `i`

	// resident `idx` is considered as an infector
	for (size_t i = 0; i < h_size; i++)
	{											   // resident `i` : infectee
		t1 = h_residents[i].get_augm_infec_time(); // time of infection of resident `i`
		if (get_augm_infec_time(idx) < t1)
		{ // NOTE : this ensures idx != i

			// Cumulative transmission rate for all household contacts
			h_cumul_rate[idx][i] = h_residents[idx].infec_cumul(t1, h_residents[i].i_study_period);

			// Instantaneous transmission rate for secondary cases
			if (h_residents[i].i_status > 0)
			{ // resident i got infectedd
				h_inst_rate[idx][i] = h_residents[idx].infec_inst(t1, h_residents[i].i_study_period);
			}
		}
	}

	// resident `idx` is considered as an infectee
	if (h_residents[idx].i_status > 0)
	{
		t1 = get_augm_infec_time(idx);
	}
	else
	{
		t1 = h_residents[idx].i_study_period;
	}

	for (size_t i = 0; i < h_size; i++)
	{ // resident `i` : infector
		if (h_residents[i].get_augm_infec_time() < t1)
		{ // NOTE : this ensures idx != i

			// Cumulative transmission rate for all household contacts
			h_cumul_rate[i][idx] = h_residents[i].infec_cumul(t1, h_residents[idx].i_study_period);

			// Instantaneous transmission rate for secondary cases
			if (h_residents[idx].i_status > 0)
			{ // resident `idx` got infected
				h_inst_rate[i][idx] = h_residents[i].infec_inst(t1, h_residents[idx].i_study_period);
			}
		}
	}
}

double Household::log_incub_proba(size_t idx)
{
	// logarithm of the pdf of the incubation period
	double out = 0.0, incub_period = 0.0;
	std::vector<double> args = {glob_params.value("incub_shape"), glob_params.value("incub_scale")};

	if (h_residents[idx].i_status == SYMPTOMATIC)
	{
		incub_period = get_augm_onset_time(idx) - get_augm_infec_time(idx);
		out = log_gamma_pdf(incub_period, args);
	}
	return out;
}

// Instaneous risk of infection from community
double Household::model_community_risk_inst(double t1)
{
	double risk_inst = glob_params.value("alpha") * Inst_Community_Infec[(int)t1];

	return risk_inst;
}

// Cumulative risk of infection from community
double Household::model_community_risk_cumul(double t0, double t1)
{
	double risk_cumul = glob_params.value("alpha") * (Cumul_Community_Infec[(int)t1] - Cumul_Community_Infec[(int)t0]);

	return risk_cumul;
}

// Relative infectivity of symptomatic child cases versus symptomatic adult cases
double Household::model_contact_modulation(size_t infectee, size_t idx)
{
	double factor = 1.0;
	size_t i = infectee;

	// NOTE : this call is missleading : use only if to be inferred...

	// Child (2) - Child (2)
	if (glob_params.infer("rCC") == true && h_residents[idx].i_age == 0 && h_residents[i].i_age == 0)
	{
		factor *= glob_params.value("rCC");

	} // Child (2) - Mother (0)
	else if (glob_params.infer("rMC") == true && ((h_residents[idx].i_age == 0 && h_residents[i].i_age == 1 && h_residents[i].i_sex == 1) || (h_residents[i].i_age == 0 && h_residents[idx].i_age == 1 && h_residents[idx].i_sex == 1)))
	{
		factor *= glob_params.value("rMC");

	} // Child (2) - Father (1)
	else if (glob_params.infer("rFC") == true && ((h_residents[idx].i_age == 0 && h_residents[i].i_age == 1 && h_residents[i].i_sex == 2) || (h_residents[i].i_age == 0 && h_residents[idx].i_age == 1 && h_residents[idx].i_sex == 2)))
	{
		factor *= glob_params.value("rFC");
	}

	return factor;
}

double Household::model_susceptibility(size_t idx)
{
	double factor = 1.0;

	// Suceptibility of different vaccination classes
	if (glob_params.infer("rU") == true && h_residents[idx].i_vacc_status == 1)
	{
		factor *= glob_params.value("rU");
	}

	if (glob_params.infer("r1V") == true && h_residents[idx].i_vacc_status == 2 && h_residents[idx].i_vacc_delay < VACC_ONSET)
	{
		factor *= glob_params.value("r1V");
	}

	if (glob_params.infer("r1V90") == true && h_residents[idx].i_vacc_status == 2 && h_residents[idx].i_vacc_delay >= VACC_ONSET)
	{
		factor *= glob_params.value("r1V90");
	}

	if (glob_params.infer("r2V") == true && h_residents[idx].i_vacc_status == 3 && h_residents[idx].i_vacc_delay < VACC_ONSET)
	{
		factor *= glob_params.value("r2V");
	}

	if (glob_params.infer("r3V") == true && h_residents[idx].i_vacc_status == 4)
	{
		factor *= glob_params.value("r3V");
	}

	// Suceptibility of children
	if (glob_params.infer("rC") == true && h_residents[idx].i_age == 0)
	{
		factor *= glob_params.value("rC");
	}

	// Suceptibility of isolated people
	if (glob_params.infer("rQ") == true && h_residents[idx].i_isolation == true)
	{
		factor *= glob_params.value("rQ");
	}

	// Suceptibility of recovered people
	if (glob_params.infer("rR") == true && h_residents[idx].i_recovered == true && (h_residents[idx].i_vacc_status == 1 || h_residents[idx].i_vacc_status == 2))
	{
		factor *= glob_params.value("rR");
	}

	return factor;
}

double Household::model_transmissability()
{
	double factor = 1.0;

	// Instantaneous per capita transmission rate
	factor *= glob_params.value("beta");

	// If transmission hazard depends on household size
	if (glob_params.infer("delta") == true)
	{
		factor /= pow(h_size / glob_params.value("typical_size"), glob_params.value("delta"));
	}

	return factor;
}

// UPDATE : to generalize with Household::log_infec_proba, log_surv_proba

/* logarithm of the infection probability */
double Household::log_infec_proba(size_t idx)
{
	// Instantaneous risk of infection from other household members
	double beta_i = 0.0, beta = 0.0, output = 0.0;
	double t1 = get_augm_infec_time(idx);

	// Instantaneous risk of infection from other household members
	for (size_t i = 0; i < h_size; ++i)
	{
		if (idx != i)
		{
			beta_i = h_inst_rate[i][idx];
			beta_i *= model_contact_modulation(i, idx);
			beta += beta_i;
		}
	}

	beta *= model_transmissability();
	output = model_susceptibility(idx) * (model_community_risk_inst(t1) + beta);
	output = std::log(1. - std::exp(-output));

	return output;
}

/* Survival probability from time t0 of resident idx */
double Household::log_surv_proba(size_t idx, double t0)
{
	// Survival from t0 to t1
	double t1, output = 0.0, beta = 0.0, beta_i = 0.0;

	if (h_residents[idx].i_status > 0)
	{ // infected case
		t1 = get_augm_infec_time(idx);
	}
	else
	{ // not a case
		t1 = h_residents[idx].i_study_period;
	}

	// Cumulative risk of infection from other household members
	for (size_t i = 0; i < h_size; ++i)
	{
		if (idx != i)
		{
			beta_i = h_cumul_rate[i][idx];
			beta_i *= model_contact_modulation(i, idx);
			beta += beta_i;
		}
	}

	beta *= model_transmissability();
	output = model_susceptibility(idx) * (model_community_risk_cumul(t0, t1) + beta);
	output *= (-1.);

	return output;
}

/* Compute the log-likelihood within household */
double Household::compute_loglike()
{
	double LL = 0.0; // Log Likelihood

	// All residents
	for (size_t i = 0; i < h_size; i++)
	{

		if (h_residents[i].i_status == NOT_A_CASE)
		{ // Contribution of not infected cases
			LL += log_surv_proba(i, h_infec_time);
		}
		else if (i == h_first_case)
		{ // Contribution of first case
			LL += log_incub_proba(i);
		}
		else
		{ // Contribution of secondary cases
			LL += log_incub_proba(i);
			LL += log_surv_proba(i, h_infec_time);
			LL += log_infec_proba(i);
		}
	}

	return LL;
}