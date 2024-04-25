// chain_binomial_heterogeneous.cpp //

#include "chain_binomial_heterogeneous.h"

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

/* Initialize the Look Up Table of combinations and their number */
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
double ll_chain_binomial_heterogeneous(usi N, usi k, usi *indidx, std::vector<double> &escape_matrix)
{
	double output;
	// submatrix has always at least 1 infector: community
	std::vector<double> sub_escape_matrix((k + 1) * (k));

	// find the submatrix of the given infectious if any
	usi row = 0, inf = 0; // community (inf=0) is always at position (row=0)
	while ((k > 0) && (row <= k))
	{
		for (usi col = 0; col < k; col++)
		{
			sub_escape_matrix[col + k * row] = escape_matrix[indidx[col] - 1 + N * inf];
		}
		// next infector
		inf = indidx[row];
		row++;
	}

	// product over all infectious of escape matrix for the escaped individuals
	double esc_product = 1.;
	for (usi i = k; i < N; i++)
	{
		// from community
		esc_product *= escape_matrix[indidx[i] - 1];

		// from infected
		for (usi j = 0; j < k; j++)
		{
			esc_product *= escape_matrix[indidx[i] - 1 + N * indidx[j]];
		}
	}

	output = esc_product * _recursive_function(k, sub_escape_matrix);

	return output;
}

/*
Recursive fonction to break likelihood into shorter escape patterns.
 */
double _recursive_function(usi N, std::vector<double> &escape_matrix)
{
	double output = 1.;
	for (usi k = 0; k < N; k++)
	{
		for (unsigned c = 0; c < Num_of_Combs_LookupTable[N][k]; c++)
		{
			output -= ll_chain_binomial_heterogeneous(
				N, k, Combinations_LookupTable[N][k][c], escape_matrix);
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
// ///////////////////////////// Model : Methods (likelihood)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize the Lookup Table. */
void Model::initialize_model()
{
	int max_num = 0, aux = 0;

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

	initialize_LookupTable(max_num);
}

/* Compute the log-likelihood for a single household. */
double Model::compute_single_house_loglike(int h)
{
	usi hh_size = househ[h].house_size;
	usi num_infec = househ[h].infectious_list.size();
	double log_likelihood = 0;

	// allocate memory for escape matrix
	std::vector<double> escape_matrix((hh_size + 1) * hh_size);

	// ---------->>---------->>---------->>---------->> //
	// Assign the values of the escape probabilities here :
	// ---------->>---------->>---------->>---------->> //

	//=====community escape probability=====//
	double commu_esc_proba = exp(-get_value("alpha"));

	for (usi i = 0; i < hh_size; i++)
	{
		escape_matrix[i] = commu_esc_proba;
	}

	//=====person-to-person escape log-probability=====//
	// NOTE : change if to distinguish between infectious, unknown, susceptibles...
	double p2p_esc_proba = exp(-get_value("beta"));

	for (usi k = 1; k <= hh_size; k++)
	{
		for (usi i = 0; i < hh_size; i++)
		{
			escape_matrix[i + k * hh_size] = p2p_esc_proba;
		}
		escape_matrix[k - 1 + k * hh_size] = 1.;
	}

	// <<----------<<----------<<----------<<---------- //
	// <<----------<<----------<<----------<<---------- //

	//=====compute the loglikelihood=====//
	// initialize array of indexes to {1, ..., hh_size}
	usi *indidx = new usi[hh_size];
	for (usi i = 0; i < hh_size; i++)
	{
		indidx[i] = i + 1;
	}

	log_likelihood = log(ll_chain_binomial_heterogeneous(hh_size, num_infec, indidx, escape_matrix));
	delete[] indidx;

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
void Model::simulate_single_house(int h)()
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
		househ[h].change_status(i, HEALTHY);
	}
	store_old_list = househ[h].unknown_status_list;
	for (int i : store_old_list)
	{
		househ[h].change_status(i, HEALTHY);
	}
	num_new_infectious = 0;

	//====GENERATION 0 : infection by community=====//
	store_old_list = househ[h].susceptibles_list;
	for (int i : store_old_list)
	{
		if (uniform(rng) > commu_esc_proba)
		{
			househ[h].change_status(i, SYMPTOMATIC);
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
			if (uniform(rng) > exp(num_old_infectious * log_p2p_esc_proba))
			{
				househ[h].change_status(i, SYMPTOMATIC);
				num_new_infectious++;
			}
		}
	}
	
}