#include "individual.hpp"

using namespace std;

/*
/////////////////////////////////////////////////////
//  Augmented_Variable : Class
/////////////////////////////////////////////////////

double Augmented_Variable::propose(std::mt19937_64 &rand_gen)
{
	return old_value() + rand(rand_gen);
}


void Augmented_Variable::init(std::mt19937_64 &rand_gen)
{

}
*/

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
// UPDATE : functional programming to avoid checking i_status

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

// ----------------------------
//  Methods : Incubation Period
// ----------------------------

// UPDATE: these could be compacted from utils

// Incubation period random generator
double Individual::incub_period_rand(std::mt19937_64 &rand_gen)
{
	std::vector<double> args = {
		glob_params.value("incub_shape"),
		glob_params.value("incub_scale")};

	return gamma_rand(rand_gen, args);
}

// pdf of the ncubation period
double Individual::incub_period_pdf(double x)
{
	std::vector<double> args = {
		glob_params.value("incub_shape"),
		glob_params.value("incub_scale")};

	if (i_status == NOT_A_CASE)
	{
		return 0.0;
	}
	else if (i_status == SYMPTOMATIC)
	{
		return gamma_pdf(x, args);
	}
	else
	{
		std::cerr << "ERROR: status not in the model." << endl;
		std::cerr << "Exit...\n"
				  << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

// Logorithm of the pdf of the incubation period
double Individual::log_incub_period_pdf(double x)
{
	std::vector<double> args = {
		glob_params.value("incub_shape"),
		glob_params.value("incub_scale")};

	if (i_status == NOT_A_CASE)
	{
		return -INFINITY;
	}
	else if (i_status == SYMPTOMATIC)
	{
		return log_gamma_pdf(x, args);
	}
	else
	{
		std::cerr << "ERROR: status not in the model." << endl;
		std::cerr << "Exit...\n"
				  << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

// ----------------------------
//  Methods : Augmented Variables
// ----------------------------

// ////////////// AUGMENTED ////////////// //
// FIXME : generalize augmented
// Functional programming + consistency checks


// Initialize augmented variables
void Individual::init_augm_data(std::mt19937_64 &rand_gen)
{
	i_augm_infec_time = i_onset_time;
	i_augm_onset_time = i_onset_time;

	i_augm_infec_time = propose_augm_infec_time(rand_gen);
	i_augm_onset_time = propose_augm_onset_time(rand_gen);
}

// Propose new value for symptom onset time
double Individual::propose_augm_onset_time(std::mt19937_64 &rand_gen)
{
	if (i_status == NOT_A_CASE)
	{
		return INFINITY;
	}
	else if (i_status == SYMPTOMATIC)
	{
		return get_augm_onset_time() + unif_rand(rand_gen);
	}
	else
	{
		std::cerr << "ERROR: status not in the model." << endl;
		std::cerr << "Exit...\n"
				  << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

// Propose new value for infection time
double Individual::propose_augm_infec_time(std::mt19937_64 &rand_gen)
{
	if (i_status == NOT_A_CASE)
	{
		return INFINITY;
	}
	else if (i_status == SYMPTOMATIC)
	{
		return get_augm_onset_time() - incub_period_rand(rand_gen);
	}
	else
	{
		std::cerr << "ERROR: status not in the model." << endl;
		std::cerr << "Exit...\n"
				  << std::endl;
		std::exit(EXIT_FAILURE);
	}
}