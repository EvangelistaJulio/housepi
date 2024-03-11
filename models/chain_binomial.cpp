// chain_binomial.cpp //

#include "chain_binomial.h"

// auxiliary global variable
std::map<int, std::vector<gentype>> Generations_Lookup_Table = {};

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

/* All possible ranked offspring generations for a given number of infected individuals n > 0.

	e.g.
	----
	n_infec = 3
	return [[0, 0, 0], [0, 0, 1], [0, 1, 1], [0, 1, 2]]
*/
std::vector<gentype> generation_span(int n)
{
	int root = n - 1;
	std::vector<gentype> output(1, gentype(n, 0));

	while (root > 0)
	{
		for (gentype g : generation_span(n - root))
		{
			gentype temp(root, 0);
			for (int i = 0; i < n - root; i++)
			{
				g[i] += 1;
			}
			temp.insert(temp.end(), g.begin(), g.end());
			output.push_back(temp);
		}
		root--;
	}

	return output;
}

/* Count-compression of all possible ranked offpsring generations */
std::vector<gentype> generation_count_span(int n)
{
	std::vector<gentype> output = {};
	for (gentype g : generation_span(n))
	{
		// I take advantage of the order granted by the map
		std::map<int, int> freq_map = {};
		for (int val : g)
		{
			if (freq_map.find(val) != freq_map.end())
			{
				freq_map[val]++;
			}
			else
			{
				freq_map.emplace(val, 1);
			}
		}

		gentype counts = {};
		for (const auto &pair : freq_map)
		{
			counts.push_back(pair.second);
		}
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
	Generations_Lookup_Table.emplace(0, std::vector<gentype>());
	for (size_t n = 1; n <= max_num; n++)
	{
		Generations_Lookup_Table.emplace(n, generation_count_span(n));
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
	double log_commu_esc_proba = - get_value("alpha");
	// community infection log-probability
	double log_commu_infec_proba = log(1. - exp(-get_value("alpha")));

	//=====within house person-to-person escape probability=====//
	double log_p2p_esc_proba = -get_value("beta");

	//=====compute the likelihood according to the Reed-Frost model=====//
	double log_likelihood = 0;
	double tmp, sero_positive_term;
	int multiplicity;

	int hh_size = househ[h].size;
	int n_positive = househ[h].infectors_list.size();

	//-----contribution of sero-negatives-----//
	log_likelihood += (hh_size - n_positive) * (log_commu_esc_proba + n_positive * log_p2p_esc_proba);

	//-----contribution of sero-positives-----//
	if (n_positive > 0)
	{
		sero_positive_term = 0.;

		// sum over all possible generations
		for (gentype gen_counts : Generations_Lookup_Table.at(n_positive))
		{
			// number of equivalent configurations returning the same kind of generation
			multiplicity = factorial(n_positive) / factorial(gen_counts[0]); // euclidean division... normally

			// GENERATION (0) : infected by community
			tmp = gen_counts[0] * log_commu_infec_proba;
			// GENERATION (0) : escaped by community
			tmp += (n_positive - gen_counts[0]) * log_commu_esc_proba;

			// GENERATION (i) : within house, i > 0
			for (size_t i = 1; i < gen_counts.size(); i++)
			{
				// GENERATION (i) : infected by previous generation i - 1
				tmp += gen_counts[i] * log(1. - exp(gen_counts[i - 1] * log_p2p_esc_proba));

				// GENERATION (i) : escaped from all others generation s.t. 0 <= generation < i - 1
				for (size_t j = 0; j < i - 1; j++)
				{
					tmp += gen_counts[i] * gen_counts[j] * log_p2p_esc_proba;
				}

				multiplicity /= factorial(gen_counts[i]);
			}

			sero_positive_term += multiplicity * exp(tmp);
		}

		log_likelihood += log(sero_positive_term);
	}

	return log_likelihood;
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
	size_t num_new_infectors = 0, num_old_infectors = 0;
	std::vector<std::size_t> copy_list = {};

	//=====community escape probability=====//
	double commu_esc_proba = exp(-get_value("alpha"));

	//=====person-to-person escape probability=====//
	double p2p_esc_proba = exp(-get_value("beta")); 

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

		//-----GENERATION 0 : infection by community-----//
		num_new_infectors = 0;
		copy_list = househ[h].susceptibles_list;
		for (std::size_t i : copy_list)
		{
			if (uniform(rng) > commu_esc_proba)
			{
				househ[h].change_status(i, SYMPTOMATIC);
				num_new_infectors++;
			}
		}

		//-----GENERATION >=1 : infection by previous generation only-----//
		while ((househ[h].susceptibles_list.size() > 0) && (num_new_infectors > 0))
		{
			num_old_infectors = num_new_infectors;
			num_new_infectors = 0;

			copy_list = househ[h].susceptibles_list;
			for (std::size_t i : copy_list)
			{
				if (uniform(rng) > pow(p2p_esc_proba, num_old_infectors))
				{
					househ[h].change_status(i, SYMPTOMATIC);
					num_new_infectors++;
				}
			}
		}
	}
}