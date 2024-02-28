// heterogeneous_reed_frost.cpp //

#include "heterogeneous_reed_frost.h"

// GLOBAL VARIABLES
// Lookup table to store combinations and complements
usi ****Combinations_LookupTable;
// Lookup table to store the number of combinations
unsigned **Num_of_Combs_LookupTable;

///////////////////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////// Auxiliary Functions //////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Binomial coefficient of N elements in k positions*/
unsigned binomial_coefficient(const usi N, const usi k)
{
	if ((N == 0) || (k == 0))
	{
		return 1;
	}
	else
	{
		unsigned recur[k];
		recur[0] = N - k + 1;

		for (usi i = 1; i < k; ++i)
		{
			recur[i] = recur[i - 1] * (N - k + 1 + i) / (i + 1);
		}
		return recur[k - 1];
	}
}

/*
Algorithm to generate all possible combinations of N elements in k positions, plus the N - k complements.
The elements are indexed from 1 to N.
https://www.cs.utexas.edu/users/djimenez/utsa/cs3343/lecture25.html
*/
void _recursive_combinations_span(usi N, usi k, vec_usi &current, usi start, usi loop, std::vector<vec_usi> &results)
{
	if (loop > k)
	{
		//-----complement of the combination-----//
		// Mark elements present in the list
		std::vector<bool> missing(N + 1, true);
		for (usi i = 0; i < k; i++)
		{
			missing[current[i]] = false;
		}
		// store complement
		for (usi i = 1; i <= N; i++)
		{
			if (missing[i])
			{
				current[loop - 1] = i;
				loop++;
			}
		}
		results.push_back(current);
	}
	else
	{
		for (usi i = start; i <= N; i++)
		{
			current[loop - 1] = i;
			_recursive_combinations_span(N, k, current, i + 1, loop + 1, results);
		}
	}
}

/* Store in the Look Up Table all combinations of N elements in k positions. */
void store_combinations(const usi N, const usi k)
{
	// Generate combinations with the recurisve algorithm
	std::vector<vec_usi> results;
	vec_usi current(N, 0); // initialize to impossible value (not in [1, N])
	_recursive_combinations_span(N, k, current, 1, 1, results);

	unsigned c = 0;
	for (vec_usi current : results)
	{
		// store the generated combination + complement
		Combinations_LookupTable[N][k][c] = new unsigned short[N];
		for (short i = 0; i < N; i++)
		{
			Combinations_LookupTable[N][k][c][i] = current[i];
		}
		c++;
	}
}

/* Initialize the Look Up Table */
void initialize_LookupTable(const short max_N)
{
	unsigned num_combin;

	Combinations_LookupTable = new usi ***[max_N + 1];
	Num_of_Combs_LookupTable = new unsigned *[max_N + 1];

	for (usi N = 0; N <= max_N; N++)
	{
		Combinations_LookupTable[N] = new usi **[N];
		Num_of_Combs_LookupTable[N] = new unsigned[N];

		for (usi k = 0; k < N; k++)
		{
			num_combin = binomial_coefficient(N, k);
			Combinations_LookupTable[N][k] = new usi *[num_combin];
			Num_of_Combs_LookupTable[N][k] = num_combin;

			store_combinations(N, k);
		}
	}
}

/*
Compute the likelihood for a heterogeneous Reed-Frost model.
Arguments :
-----
	N : household size
	k : number of seropositive
	indidx : [0 - k) seropositive indexes, [k - N) seronegative indexes. Indexes are [1, N]
	escape_matrix : probability 'col' escaped from 'row' (including community as index 0)
*/
double Likelihood_Heterogeneous_ReedFrost(usi N, usi k, vec_usi &indidx, mtx &escape_matrix)
{
	double output;
	// submatrix has always at least 1 infector: community
	mtx sub_escape_matrix(k + 1, std::vector<double>(k));

	// find the submatrix of the given infectors if any
	usi row = 0, inf = 0;
	while ((k > 0) && (row <= k))
	{
		for (usi col = 0; col < k; col++)
		{
			sub_escape_matrix[row][col] = escape_matrix[inf][indidx[col] - 1];
		}
		// next infector
		inf = indidx[row];
		row++;
	}

	// product over all infectors of escape matrix for the escaped individuals
	double esc_product = 1.;
	for (usi i = k; i < N; i++)
	{
		// from community
		esc_product *= escape_matrix[0][indidx[i] - 1];

		// from infected
		for (usi j = 0; j < k; j++)
		{
			esc_product *= escape_matrix[indidx[j]][indidx[i] - 1];
		}
	}
	output = Recursive_Function(k, sub_escape_matrix) * esc_product;

	return output;
}

/*
Recursive fonction to break likelihood into shorter escape patterns.
 */
double Recursive_Function(usi N, mtx &escape_matrix)
{
	double output = 1.;
	for (usi k = 0; k < N; k++)
	{
		for (unsigned c = 0; c < Num_of_Combs_LookupTable[N][k]; c++)
		{
			// FIXME : temporary
			vec_usi sub_indidx(N, 0);
			for(usi j = 0; j < N; j ++)
			{
				sub_indidx[j] = Combinations_LookupTable[N][k][c][j];
			}
			//

			output -= Likelihood_Heterogeneous_ReedFrost(N, k, sub_indidx, escape_matrix);
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

/* Initialize the Lookup Table. */
void Model::initialize_model()
{
	size_t max_size = 0;

	for (size_t h = 0; h < num_househ; h++)
	{
		max_size = (househ[h].size > max_size) ? househ[h].size : max_size;
	}
	initialize_LookupTable(max_size);
}

/* Compute the log-likelihood within household. */
double Model::compute_single_house_loglike(std::size_t h)
{
	double output = 0;

	int hh_size = househ[h].size;
	int num_infec = househ[h].infectors_list.size();
	int num_unknown_status = househ[h].unknown_status_list.size();

	// allocate memory for escape matrix
	mtx escape_matrix(hh_size + 1, std::vector<double>(hh_size));

	// ---------->> //
	// Assign the values of the escape probabilities here :
	// ---------->> //

	bool is_child;

	// community escape probability
	int c = 0;
	for (usi i : househ[h].house_members)
	{
		is_child = indiv[i].get_data("age") < get_value("child_age");
		escape_matrix[0][c] = exp(-((is_child) ? get_value("beta_child") : get_value("beta_adult")));
		c++;
	}

	// FIXME : temporary, need to distinguish between infectors, unknown, susceptibles...
	// person-to-person within house escape probability
	double p2p_esc_proba = exp(-get_value("lambda") * pow(hh_size, -get_value("gamma")));
	for (usi k = 1; k <= hh_size; k++)
	{
		for (usi i = 0; i < hh_size; i++)
		{
			escape_matrix[k][i] = p2p_esc_proba;
		}
		escape_matrix[k][k - 1] = 1.; // note : this is just a place-holder
	}

	// <<---------- //
	// <<---------- //

	// FIXME : if heterogeneous within house this part must change

	// initial array of individuals
	vec_usi indidx(hh_size);
	for (int i = 0; i < hh_size; i++)
	{
		indidx[i] = i + 1;
	}

	// compute the loglikelihood for all the scenario compatible with the missing data
	for (int i = 0; i <= num_unknown_status; i++)
	{
		output += log(Likelihood_Heterogeneous_ReedFrost(hh_size, num_infec + i, indidx, escape_matrix));
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
	bool is_child;
	size_t num_new_infectors = 0, num_old_infectors = 0;
	std::vector<std::size_t> copy_list = {};
	double commu_esc_proba, log_p2p_esc_proba;

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

		//====GENERATION 0 : infection by community=====//
		copy_list = househ[h].susceptibles_list;
		for (std::size_t i : copy_list)
		{
			is_child = indiv[i].get_data("age") < get_value("child_age");
			commu_esc_proba = exp(-((is_child) ? get_value("beta_child") : get_value("beta_adult")));

			if (uniform(rng) > commu_esc_proba)
			{
				househ[h].change_status(i, SYMPTOMATIC);
				indiv[i].set_status(SYMPTOMATIC);
			}
		}

		// compute person-to-person transmission // FIXME : this is the constant transmission proba case
		log_p2p_esc_proba = -get_value("lambda") * pow(househ[h].size, -get_value("gamma"));

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