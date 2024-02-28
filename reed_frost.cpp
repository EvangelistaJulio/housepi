// reed_frost.cpp //

#include "reed_frost.h"

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

/* All possible ranked offspring generations for a given number of infected individuals.
	e.g.
	----
	n_infec = 3
	return [[0, 0, 0], [0, 0, 1], [0, 1, 1], [0, 1, 2]]
*/
std::vector<gentype> generation_span(int n_infec)
{
	int root = n_infec - 1;
	std::vector<gentype> output(1, gentype(n_infec, 0));

	while (root > 0)
	{
		for (gentype g : generation_span(n_infec - root))
		{
			gentype temp(root, 0);
			for (int i = 0; i < n_infec - root; i++)
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
std::vector<gentype> generation_count_span(int n_infec)
{
	std::vector<gentype> output;
	for (gentype g : generation_span(n_infec))
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
	int hh_size = househ[h].size;
	int num_unknown_status = househ[h].unknown_status_list.size();
	int num_infec = househ[h].infectors_list.size();
	int augm_num_infec = 0;

	double output = 0;
	double degener = 0., tmp = 0., infectious_term = 0.;
	double hh_surface = 0.;
	double log_p2p_esc_proba = 0;

	// community escape probability
	double log_commu_esc_proba = -get_value("alpha");

	// augmented number of possibly infected household members
	for (int k = 0; k <= num_unknown_status; k++)
	{
		augm_num_infec = num_infec + k;

		// the case of zero infected is hard coded
		if (num_infec == 0)
		{
			output += hh_size * log_commu_esc_proba;
		}
		else
		{
			// FIXME : this convoluted way of accessing the surface should be revised
			hh_surface = indiv[househ[h].infectors_list[0]].get_data("HH_surface");

			// person-to-person escape probability
			log_p2p_esc_proba = -get_value("beta") - get_value("aero_eff") * get_value("typical_surface") / hh_surface;

			// sum over all possible generations
			for (gentype gen_counts : generations_Look_Up_Table.at(augm_num_infec))
			{
				degener = (double)factorial(augm_num_infec) / factorial(gen_counts[0]);

				tmp = gen_counts[0] * log(1. - exp(log_commu_esc_proba));
				tmp += (augm_num_infec - gen_counts[0]) * log_commu_esc_proba;

				for (size_t i = 1; i < gen_counts.size(); i++)
				{
					tmp += gen_counts[i] * log(1. - pow(exp(log_p2p_esc_proba), gen_counts[i - 1]));
					for (size_t j = 0; j < i - 1; j++)
					{
						tmp += gen_counts[j] * gen_counts[i] * log_p2p_esc_proba;
					}

					degener /= factorial(gen_counts[i]);
				}
				infectious_term += degener * exp(tmp);
			}

			output += log(infectious_term) + (hh_size - augm_num_infec) * (augm_num_infec * log_p2p_esc_proba + log_commu_esc_proba);
		}
	}

	return output;
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
	double commu_infec_proba = 1. - exp(-get_value("alpha"));
	double log_p2p_esc_proba = 0.;
	double hh_surface = 0.;

	// simulate the infection for each household
	for (std::size_t h = 0; h < num_househ; h++)
	{
		//=====INITIALIZE=====//
		// set all cases to healthy (change infectors to susceptibles)
		copy_list = househ[h].infectors_list;
		for (std::size_t i : copy_list)
		{
			indiv[i].set_status(HEALTHY);
			househ[h].change_status(i, HEALTHY);
		}

		// FIXME : this convoluted way of accessing the surface should be revised
		hh_surface = indiv[househ[h].susceptibles_list[0]].get_data("HH_surface");

		//====GENERATION 0 : infection by community=====//
		copy_list = househ[h].susceptibles_list;
		for (std::size_t i : copy_list)
		{
			if (uniform(rng) < commu_infec_proba)
			{
				househ[h].change_status(i, SYMPTOMATIC);
				indiv[i].set_status(SYMPTOMATIC);
			}
		}

		// compute person-to-person transmission
		log_p2p_esc_proba = -get_value("beta") - get_value("aero_eff") * get_value("typical_surface") / hh_surface;

		//====GENERATION >=1 : infection by previous generatioon=====//
		num_new_infectors = househ[h].infectors_list.size();
		while ((househ[h].susceptibles_list.size() > 0) && (num_new_infectors > 0))
		{
			num_old_infectors = num_new_infectors;
			num_new_infectors = 0;

			copy_list = househ[h].susceptibles_list;
			for (std::size_t i : copy_list)
			{
				if (uniform(rng) > exp(num_old_infectors * log_p2p_esc_proba))
				{
					househ[h].change_status(i, SYMPTOMATIC);
					indiv[i].set_status(SYMPTOMATIC);
					num_new_infectors++;
				}
			}
		}
	}
}