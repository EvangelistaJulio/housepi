// individual.cpp //

#include "individual.hpp"

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
}

// ----------------------------
//  Methods : Infectivity Profiles
// ----------------------------

// FIXME : these are probabilities, but we want rates.

double Individual::infec_inst(double t_transm)
{
	double output = 0.0, tost;
	kwargs args = {
		{"alpha", glob_params.value("tost_alpha")},
		{"mu", glob_params.value("tost_mu")},
		{"sigma", glob_params.value("tost_sigma")},
		{"tau", glob_params.value("tost_tau")}};

	tost = t_transm - get_augm_onset_time();
	output = t_onset_sympt_to_transm_pdf(tost, get_augm_incub_period(), args);

	return output;
}

double Individual::infec_cumul(double t_transm)
{
	double output = 0.0, tost;
	kwargs args = {
		{"alpha", glob_params.value("tost_alpha")},
		{"mu", glob_params.value("tost_mu")},
		{"sigma", glob_params.value("tost_sigma")},
		{"tau", glob_params.value("tost_tau")}};

	tost = t_transm - get_augm_onset_time();
	output = t_onset_sympt_to_transm_cdf(tost, get_augm_incub_period(), args);

	return output;
}

// ----------------------------
//  Methods : Incubation Period
// ----------------------------

// Logorithm of the pdf of the incubation period
// UPDATE: Functional programming + consistency checks :
// to be safer one could just load pointer to function
// now I just call the function upon condition
double Individual::log_incub_period_pdf(double x)
{
	kwargs args = {
		{"shape", glob_params.value("incub_shape")},
		{"scale", glob_params.value("incub_scale")}};

	return log_gamma_pdf(x, args);
}

// ----------------------------
//  Methods : Augmented Variables
// ----------------------------

// UPDATE : generalize augmented

// Initialize augmented variables
// FIXME : this is the only control on i_status being in {HEALTHY,SYMPTOMATIC}
void Individual::init_augm_data(std::mt19937_64 &rand_gen)
{
	if (i_status == HEALTHY)
	{
		i_augm_incub_period = 0.;
		i_augm_onset_time = i_study_period;
	}
	else if (i_status == SYMPTOMATIC)
	{
		i_augm_incub_period = 0.;
		i_augm_onset_time = i_onset_time;
		i_augm_incub_period = propose_augm_incub_period(rand_gen);
		i_augm_onset_time = propose_augm_onset_time(rand_gen);
	}
	else
	{
		std::cerr << "ERROR: status not in the model." << endl;
		std::cerr << "Exit...\n"
				  << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

// Propose new value for incubation period
double Individual::propose_augm_incub_period(std::mt19937_64 &rand_gen)
{
	kwargs args = {
		{"shape", glob_params.value("incub_shape")},
		{"scale", glob_params.value("incub_scale")}};

	return gamma_rand(rand_gen, args);
}

// Propose new value for symptom onset time
double Individual::propose_augm_onset_time(std::mt19937_64 &rand_gen)
{
	return get_augm_onset_time() + unif_rand(rand_gen);
}
