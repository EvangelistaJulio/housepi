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
			// FIXME : temporary (LookupTable to be converted to array)
			vec_usi sub_indidx(N, 0);
			for (usi j = 0; j < N; j++)
			{
				sub_indidx[j] = Combinations_LookupTable[N][k][c][j];
			}

			output -= Likelihood_Heterogeneous_ReedFrost(N, k, sub_indidx, escape_matrix);
			//
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
	size_t max_num = 0, aux = 0;

	// find maximum possible number of sero-positives
	for (size_t h = 0; h < num_househ; h++)
	{
		aux = househ[h].infectors_list.size() + househ[h].unknown_status_list.size();
		max_num = (aux > max_num) ? aux : max_num;
	}

	initialize_LookupTable(max_num);
}

/* Compute the log-likelihood within household. */
double Model::compute_single_house_loglike(std::size_t h)
{
	// >>>> FIXME : temporary
	if (househ[h].unknown_status_list.size() > 0)
	{
		std::string err_txt = "Missing status found for house " + std::to_string(h) + ". Case not modeled yet.";
		exit_with_failure_error(err_txt);
	}
	// <<<< <<<< <<<< <<<<

	double log_likelihood = 0;

	int hh_size = househ[h].size;
	int num_infec = househ[h].infectors_list.size();

	// allocate memory for escape matrix
	mtx escape_matrix(hh_size + 1, std::vector<double>(hh_size));

	// ---------->> //
	// Assign the values of the escape probabilities here :
	// ---------->> //

	//=====community escape log-probability=====//
	double commu_esc_proba = get_value("commu_esc_proba"); // exp(-get_value("alpha"));
	// bool is_child;
	// double commu_esc_proba_child = exp(-get_value("alpha_child"));

	for (usi i = 0; i < hh_size; i++)
	{
		// is_child = indiv[househ[h].house_members[i]].get_data("age") < get_value("child_age");
		// escape_matrix[0][i] = (is_child) ? commu_esc_proba : commu_esc_proba_child;
		escape_matrix[0][i] = commu_esc_proba;
	}

	//=====within house person-to-person escape probability=====//
	// FIXME : temporary, need to distinguish between infectors, unknown, susceptibles...
	// double p2p_esc_proba = exp(-get_value("lambda") * pow(hh_size, -get_value("gamma")));
	double p2p_esc_proba = get_value("p2p_esc_proba"); // exp(-get_value("beta"));

	for (usi k = 1; k <= hh_size; k++)
	{
		for (usi i = 0; i < hh_size; i++)
		{
			escape_matrix[k][i] = p2p_esc_proba;
		}
		escape_matrix[k][k - 1] = 1.;
	}

	// <<---------- //
	// <<---------- //

	//=====compute the loglikelihood=====//
	// initialize array of indexes to {1, ..., hh_size}
	vec_usi indidx(hh_size);
	std::iota(indidx.begin(), indidx.end(), 1);

	log_likelihood = log(Likelihood_Heterogeneous_ReedFrost(hh_size, num_infec, indidx, escape_matrix));

	// FLAG
	// std::cout << "log_likelihood : " << log_likelihood << std::endl;
	// std::exit(1);
	// END-of-FLAG

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

	// bool is_child = indiv[i].get_data("age") < get_value("child_age");
	double commu_esc_proba = get_value("commu_esc_proba"); // exp(-((is_child) ? get_value("beta_child") : get_value("beta_adult")));

	// compute person-to-person transmission // FIXME : this is the constant transmission proba case
	double log_p2p_esc_proba = log(get_value("p2p_esc_proba")); // -get_value("lambda") * pow(househ[h].size, -get_value("gamma"));

	// simulate the infection for each household
	for (std::size_t h = 0; h < num_househ; h++)
	{
		//=====INITIALIZE=====//
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
		num_new_infectors = 0;

		//====GENERATION 0 : infection by community=====//
		copy_list = househ[h].susceptibles_list;
		for (std::size_t i : copy_list)
		{
			if (uniform(rng) > commu_esc_proba)
			{
				househ[h].change_status(i, SYMPTOMATIC);
				num_new_infectors++;
			}
		}

		//====GENERATION >=1 : infection by previous generatioon=====//
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
					num_new_infectors++;
				}
			}
		}
	}
}