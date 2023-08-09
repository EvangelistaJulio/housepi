// model.cpp //

#include "household.hpp"

using namespace std;

// Shared Parameter Dictionary
extern std::vector<double> Cumul_Community_Infec;
extern std::vector<double> Inst_Community_Infec;

/////////////////////////////////////////////////////
//  Household : Methods (model)
/////////////////////////////////////////////////////

// UPDATE : this could be interpolated instead

// Instaneous risk of infection from community
double Household::model_community_risk_inst(double t1)
{
	return glob_params.value("alpha") * Inst_Community_Infec[(int) t1];
}

// Cumulative risk of infection from community
double Household::model_community_risk_cumul(double t0, double t1)
{
	return glob_params.value("alpha") * (Cumul_Community_Infec[(int) t1] - Cumul_Community_Infec[(int) t0]);
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