// chain_binomial_heterogeneous.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#include "chain_binomial_heterogeneous.h"

usi ****Combinations_LUT;	 // look-up-table to store combinations and complements
unsigned **Num_of_Combs_LUT; // look-up-table to store the number of combinations

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////// Auxiliary Functions //////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Binomial coefficient of N elements in k positions*/
unsigned binomial_coefficient(const usi N, const usi k)
{
	if ((N == 0) || (k == N))
	{
		return 1;
	}
	else
	{
		unsigned result = N;
		usi stop = (k > N - k) ? k : N - k; // use symmetry
		for (usi i = 1; i < stop; i++)
		{
			result = (result * (N - i)) / (i + 1);
		}
		return result;
	}
}

/*
Algorithm to generate all possible combinations of N elements in k positions, plus the N - k complements.
The elements are indexed from 1 to N.

----
Ref:
----
Daniel Angel Jimenez
https://www.cs.utexas.edu/users/djimenez/utsa/cs3343/lecture25.html
*/
void recursive_combinations_span(usi N, usi k, vec_usi &current, usi start, usi loop, std::vector<vec_usi> &results)
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
			recursive_combinations_span(N, k, current, i + 1, loop + 1, results);
		}
	}
}

/* Store in the Look Up Table all combinations of N elements in k positions. */
void store_combinations(const usi N, const usi k)
{
	// Generate combinations with the recurisve algorithm
	std::vector<vec_usi> results;
	vec_usi current(N, 0); // initialize to impossible value (not in [1, N])
	recursive_combinations_span(N, k, current, 1, 1, results);

	unsigned c = 0;
	for (vec_usi current : results)
	{
		// store the generated combination + complement
		Combinations_LUT[N][k][c] = new unsigned short[N];
		for (short i = 0; i < N; i++)
		{
			Combinations_LUT[N][k][c][i] = current[i];
		}
		c++;
	}
}

/* Initialize the Look Up Table of combinations and their number */
void initialize_LUT(const short max_N)
{
	unsigned num_combin;

	Combinations_LUT = new usi ***[max_N + 1];
	Num_of_Combs_LUT = new unsigned *[max_N + 1];

	for (usi N = 0; N <= max_N; N++)
	{
		Combinations_LUT[N] = new usi **[N];
		Num_of_Combs_LUT[N] = new unsigned[N];

		for (usi k = 0; k < N; k++)
		{
			num_combin = binomial_coefficient(N, k);
			Combinations_LUT[N][k] = new usi *[num_combin];
			Num_of_Combs_LUT[N][k] = num_combin;

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

	return esc_product * recursive_function(k, sub_escape_matrix);
}

/* Recursive function to break likelihood into shorter escape patterns. */
double recursive_function(usi N, std::vector<double> &escape_matrix)
{
	double output = 1.;

	for (usi k = 0; k < N; k++)
	{
		for (unsigned c = 0; c < Num_of_Combs_LUT[N][k]; c++)
		{
			output -= ll_chain_binomial_heterogeneous(
				N, k, Combinations_LUT[N][k][c], escape_matrix);
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
// ///////////////////////////////  Model : initialize model  //////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize the Lookup Table. */
void Model::initialize_model()
{
	// Unknown status are not modeled
	if (report_data["n.unk.stat"] > 0)
	{
		display_error_message("A missing status was found : incompatible with the model.");
		std::exit(EXIT_FAILURE);
	}

	//=====Find maximum possible number of sero-positives=====//
	int max_num = 0, aux = 0;
	for (int h = 0; h < num_househ; h++)
	{
		aux = househ[h].infectious_list.size();
		max_num = (aux > max_num) ? aux : max_num;
	}

	//=====Initialize look-up-table to largest possible value=====//
	initialize_LUT(max_num);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////  Model : compute single house loglike  /////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood for a single household. */
double Model::compute_single_house_loglike(int h)
{
	usi hh_size = househ[h].house_size;
	usi num_infec = househ[h].infectious_list.size();
	double log_likelihood = 0;

	// allocate memory for escape matrix
	std::vector<double> escape_matrix((hh_size + 1) * hh_size, 1.);

	// ---------->>---------->>---------->>---------->>
	// ---------->>---------->>---------->>---------->>

	// Assign the values of the escape probabilities HERE :

	// community escape probability
	double commu_esc_proba = exp(-get_value("alpha"));

	for (usi i = 0; i < hh_size; i++)
	{
		escape_matrix[i] = commu_esc_proba;
	}

	// person-to-person escape log-probability
	// NOTE : change if to distinguish between infectious, unknown, susceptibles...
	double p2p_esc_proba = exp(-get_value("beta"));

	for (usi k = 1; k <= hh_size; k++)
	{
		for (usi i = 0; i < hh_size; i++)
		{
			escape_matrix[i + k * hh_size] = p2p_esc_proba;
		}
	}

	// <<----------<<----------<<----------<<----------
	// <<----------<<----------<<----------<<----------

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
void Model::simulate_single_house(int h)
{
	// FIXME

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

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//