// model.h //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  * *** * CHAIN BINOMIAL HETEROGENEOUS * *** *  /////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"

//------------------------------------------------------------------------------------------------------------//
//---------------------------------------- - - - CHANGE-ME - - - ---------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- \/ -------------------------------------------------------//

/* Assign the parameters controlling infectivity for each individual */
void Model::assign_rate_per_group_to_individual(Individual &this_indiv)
{
	//=====Susceptibility=====//
	if (this_indiv.get_data("age") <= get_value("child_age"))
	{ // child
		this_indiv.assign_rate_per_group("Susceptibility", "sus_child");
	}
	else
	{ // adult
		if (this_indiv.get_data("sex") == 0)
		{ // adult women
			this_indiv.assign_rate_per_group("Susceptibility", "sus_women");
		}
		else
		{ // adult men : Susceptibility [REFERENCE]
			this_indiv.assign_rate_per_group("Susceptibility", REFERENCE);
		}
	}

	//=====Infectivity=====//
	if (this_indiv.get_data("age") > get_value("child_age"))
	{ // adult : Susceptiblity [REFERENCE]
		this_indiv.assign_rate_per_group("Infectivity", REFERENCE);
	}
	else
	{ // child
		this_indiv.assign_rate_per_group("Infectivity", "inf_child");
	}
}

/* The house escape matrix is a flattern vector containing escape probability. */
void Model::fill_house_escape_matrix(int h, std::vector<double> &escape_matrix)
{
	int hh_size = househ[h].house_size;
	double roi;

	// ordered members group rates : 1) infectious, 2) susceptibles
	int count = 0;
	std::vector<double> susc_rates(1., hh_size), infec_rates(1., hh_size);
	for (int i : househ[h].infectious_list)
	{
		susc_rates[count] = househ[h].indiv[i].get_group_rate("Susceptibility");
		infec_rates[count] = househ[h].indiv[i].get_group_rate("Infectivity");
		count++;
	}
	for (int i : househ[h].susceptibles_list)
	{
		susc_rates[count] = househ[h].indiv[i].get_group_rate("Susceptibility");
		infec_rates[count] = househ[h].indiv[i].get_group_rate("Infectivity");
		count++;
	}

	// assign values to matrix
	for (int i = 0; i < hh_size; i++)
	{
		// community terms
		escape_matrix[i] = exp(-susc_rates[i] * get_value("alpha"));

		// person-to-person terms
		for (int k = 0; k < hh_size; k++)
		{
			if (k != i)
			{
				roi = get_value("beta") * susc_rates[i] * infec_rates[k];
				escape_matrix[i + (k + 1) * hh_size] = exp(-roi);
			}
			// else is equal to 1 (by declaration outside)
		}
	}
}

//------------------------------------------------- /\ -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//---------------------------------------- - - - CHANGE-ME - - - ---------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

usi ****Combinations_LUT;	 // look-up-table to store combinations and complements
unsigned **Num_of_Combs_LUT; // look-up-table to store the number of combinations

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
		raise_error_message("A missing status was found : incompatible with the model.");
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

	//=====Escape Matrix=====//
	std::vector<double> escape_matrix((hh_size + 1) * hh_size, 1.);
	fill_house_escape_matrix(h, escape_matrix);

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
	int hh_size = househ[h].house_size;
	int num_new_infectious = 0;
	std::vector<int> store_old_list = {};
	std::vector<int> old_infectors_list = {}, new_infectors_list = {};

	//=====Escape Matrix=====//
	std::vector<double> escape_matrix((hh_size + 1) * hh_size, 1.);
	fill_house_escape_matrix(h, escape_matrix);

	//=====INITIALIZE=====//
	// set all cases to healthy (change infectious to susceptibles)
	store_old_list = househ[h].infectious_list;
	for (int i : store_old_list)
	{
		househ[h].change_status(i, SUSCEPTIBLE);
	}

	//====GENERATION 0 : infection by community=====//
	store_old_list = househ[h].susceptibles_list;
	for (int i : store_old_list)
	{
		if (uniform() > escape_matrix[i])
		{
			househ[h].change_status(i, INFECTIOUS);
			new_infectors_list.push_back(i);
		}
	}

	//====GENERATION >=1 : infection by previous generatioon=====//
	while ((househ[h].susceptibles_list.size() > 0) && (new_infectors_list.size() > 0))
	{
		old_infectors_list = new_infectors_list;
		new_infectors_list.clear();

		store_old_list = househ[h].susceptibles_list;
		for (int i : store_old_list)
		{
			for (int k : old_infectors_list)
			{
				if (househ[h].indiv[i].get_status() == SUSCEPTIBLE)
				{
					if (uniform() > escape_matrix[i + (k + 1) * hh_size])
					{
						househ[h].change_status(i, INFECTIOUS);
						new_infectors_list.push_back(i);
					}
				}
			}
		}
	}
}
