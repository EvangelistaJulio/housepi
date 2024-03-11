// saturation.cpp //

#include "saturation.h"

// auxiliary global variable
std::map<int, std::vector<gentype>> Saturation_Lookup_Table = {};

///////////////////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////// Auxiliary Functions //////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Standard definition for the factorial */
unsigned long long factorial(int n)
{
	if (n <= 1)
	{
		return 1;
	}
	else
	{
		return n * factorial(n - 1);
	}
}

/* Count-compression of all possible ranked saturated offpsring generations */
std::vector<gentype> saturation_count_span(int n)
{
	std::vector<gentype> output = {};

	for (int w = 0; w < n; w++)
	{
		gentype counts = {n - w, w};
		output.push_back(counts);
	}

	return output;
}

/* Initialize the Lookup Table */
void Model::initialize_model()
{
	size_t max_num = 0, aux = 0;

	// find maximum possible number of sero-positives
	for (size_t h = 0; h < num_househ; h++)
	{
		aux = househ[h].infectors_list.size() + househ[h].unknown_status_list.size();
		max_num = (aux > max_num) ? aux : max_num;
	}

	// initialize LookupTable to largest possible value
	Saturation_Lookup_Table.emplace(0, std::vector<gentype>());
	for (size_t n = 1; n <= max_num; n++)
	{
		Saturation_Lookup_Table.emplace(n, saturation_count_span(n));
	}
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////// Individual : Methods for Model ////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Time of infection */
// FIXME : this is not user-friendly
double Individual::get_infec_time() const
{
	return 0.;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////// Model : Methods (likelihood)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood within household. */
double Model::compute_single_house_loglike(std::size_t h)
{
	// >>>> FIXME : temporary
	if (househ[h].unknown_status_list.size() > 0)
	{
		std::string err_txt = "Missing status found for house " + std::to_string(h) + ". Cse not modeled yet.";
		exit_with_failure_error(err_txt);
	}
	// <<<< <<<< <<<< <<<<

	//=====community escape log-probability=====//
	double log_commu_esc_proba = -get_value("alpha");

	//=====within-houshold escape probability=====//
	double log_within_esc_proba = -get_value("beta");

	//=====compute the likelihood according to the Reed-Frost model=====//
	int hh_size = househ[h].size;
	int n_positive = househ[h].infectors_list.size();
	int multiplicity;

	double log_commu_infec_proba = log(1. - exp(log_commu_esc_proba));
	double log_within_infec_proba = log(1. - exp(log_within_esc_proba));
	double tmp, log_sero_negative_term = 0., sero_positive_terms = 0.;

	if (n_positive == 0)
	{
		//-----contribution of sero-negatives-----//
		return hh_size * log_commu_esc_proba;
	}
	else
	{
		//-----contribution of sero-negatives-----//
		log_sero_negative_term = (hh_size - n_positive) * (log_commu_esc_proba + log_within_esc_proba);

		//-----contribution of sero-positives-----//
		// sum over all possible saturated generations
		for (gentype sat_counts : Saturation_Lookup_Table.at(n_positive))
		{
			// number of equivalent configurations returning the same kind of generation
			// euclidean division... normally
			multiplicity = factorial(n_positive) / (factorial(sat_counts[0]) * factorial(sat_counts[1]));

			// GENERATION (0) : infected by community
			tmp = sat_counts[0] * log_commu_infec_proba;
			// GENERATION (0) : escaped by community
			tmp += sat_counts[1] * log_commu_esc_proba;

			// GENERATION (1) : infected within the house
			tmp += sat_counts[1] * log_within_infec_proba;

			sero_positive_terms += multiplicity * exp(tmp);
		}

		return log(sero_positive_terms) + log_sero_negative_term;
	}
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////// Model : Methods (simulation)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Poisson point process to simulate data. */
void Model::simulate_infection()
{
	bool is_infec_by_community;
	std::vector<std::size_t> copy_list = {};

	//=====community escape probability=====//
	double commu_esc_proba = exp(-get_value("alpha"));

	//=====within-household escape probability=====//
	double within_esc_proba = exp(-get_value("beta"));

	//=====Simulate the infection for each household=====//
	for (std::size_t h = 0; h < num_househ; h++)
	{
		//-----INITIALIZE-----//
		// set all cases to healthy (change infectors and unknown to susceptibles)
		copy_list = househ[h].infectors_list;
		for (std::size_t i : copy_list)
		{
			househ[h].change_status(i, HEALTHY);
		}
		copy_list = househ[h].unknown_status_list;
		for (std::size_t i : copy_list)
		{
			househ[h].change_status(i, HEALTHY);
		}
		is_infec_by_community = false;
		
		//-----GENERATION 0 : infection by community-----//
		copy_list = househ[h].susceptibles_list;
		for (std::size_t i : copy_list)
		{
			if (uniform(rng) > commu_esc_proba)
			{
				househ[h].change_status(i, SYMPTOMATIC);
				is_infec_by_community = true;
			}
		}

		//-----GENERATION =1 : infection by previous generation only-----//
		if (is_infec_by_community)
		{
			copy_list = househ[h].susceptibles_list;
			for (std::size_t i : copy_list)
			{
				if (uniform(rng) > within_esc_proba)
				{
					househ[h].change_status(i, SYMPTOMATIC);
				}
			}
		}
	}
}