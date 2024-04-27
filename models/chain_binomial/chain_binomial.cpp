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
	for (int h = 0; h < num_househ; h++)
	{
		// >>>> FIXME : temporary
		if (househ[h].unknown_status_list.size() > 0)
		{
			std::string err_txt = "Missing status found for house " + std::to_string(h) + ". Case not modeled yet.";
			raise_error_message(err_txt);
			std::exit(EXIT_FAILURE);
		}
		// <<<< <<<< <<<< <<<<

		aux = househ[h].infectious_list.size() + househ[h].unknown_status_list.size();
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

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////// Model : Methods (likelihood)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood within household. */
double Model::compute_single_house_loglike(int h)
{
	//=====community escape log-probability=====//
	double log_commu_esc_proba = -get_value("alpha");
	// community infection log-probability
	double log_commu_infec_proba = log(1. - exp(-get_value("alpha")));

	//=====within house person-to-person escape probability=====//
	double log_p2p_esc_proba = -get_value("beta");

	//=====compute the likelihood according to the Reed-Frost model=====//
	double log_likelihood = 0;
	double tmp, sero_positive_term;
	int multiplicity;

	int hh_size = househ[h].house_size;
	int n_positive = househ[h].infectious_list.size();

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

/* Simulate the infection within a household. */
void Model::simulate_single_house(int h)
{
	int num_new_infectious = 0, num_old_infectious = 0;
	std::vector<int> store_old_list = {};

	//-----community escape probability-----//
	double commu_esc_proba = exp(-get_value("alpha"));

	//-----person-to-person escape log-probability-----//
	double log_p2p_esc_proba = -get_value("beta");

	//=====INITIALIZE=====//
	// set all cases to healthy (change infectious and unknown to susceptibles)
	store_old_list = househ[h].infectious_list;
	for (int i : store_old_list)
	{
		househ[h].change_status(i, SUSCEPTIBLE);
	}
	store_old_list = househ[h].unknown_status_list;
	for (int i : store_old_list)
	{
		househ[h].change_status(i, SUSCEPTIBLE);
	}
	num_new_infectious = 0;

	//====GENERATION 0 : infection by community=====//
	store_old_list = househ[h].susceptibles_list;
	for (int i : store_old_list)
	{
		if (uniform() > commu_esc_proba)
		{
			househ[h].change_status(i, INFECTIOUS);
			num_new_infectious++;
		}
	}

	//====GENERATION >=1 : infection by previous generatioon=====//
	while ((househ[h].susceptibles_list.size() > 0) && (num_new_infectious > 0))
	{
		num_old_infectious = num_new_infectious;
		num_new_infectious = 0;

		store_old_list = househ[h].susceptibles_list;
		for (int i : store_old_list)
		{
			if (uniform() > exp(num_old_infectious * log_p2p_esc_proba))
			{
				househ[h].change_status(i, INFECTIOUS);
				num_new_infectious++;
			}
		}
	}
}