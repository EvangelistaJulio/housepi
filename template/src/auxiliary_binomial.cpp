// auxiliary_binomial.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#include "auxiliary_binomial.h"

// Binomial Auxiliary Namespace
namespace bax
{

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////  GLOBAL VARIABLES  ///////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////

	// lookup-table with generations
	std::vector<std::vector<gentype>> Generations_LookupTable;
	// lookup-table with number of generations
	std::vector<std::vector<unsigned>> Num_Generations_LookupTable;

	// lookup-table to store combinations and complements
	std::vector<std::vector<std::vector<std::vector<ushort>>>> Combinations_LookupTable = {};
	// lookup-table to store the number of combinations
	std::vector<std::vector<unsigned>> Num_Combinations_LookupTable;

	/* Standard definition for the factorial. */
	unsigned factorial(unsigned n)
	{
		unsigned result = 1;
		for (unsigned i = 2; i <= n; i++)
		{
			result *= i;
		}
		return result;
	}

	/* Binomial coefficient of N elements in k positions*/
	unsigned binomial_coefficient(const ushort N, const ushort k)
	{
		if ((N == 0) || (k == N))
		{
			return 1;
		}
		else
		{
			unsigned result = N;
			ushort stop = (k > N - k) ? k : N - k; // use symmetry
			for (ushort i = 1; i < stop; i++)
			{
				result = (result * (N - i)) / (i + 1);
			}
			return result;
		}
	}

	/* All possible ranked offspring generations for a given number of infected individuals n > 0.

		e.g.
		----
		n_infec = 3
		return [[0, 0, 0], [0, 0, 1], [0, 1, 1], [0, 1, 2]]
	*/
	std::vector<gentype> generation_span(ushort n)
	{
		ushort root = n - 1;
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

	/* Initialize the Look Up Table of count-compressed generations */
	void init_generations_lookuptable(const ushort max_N)
	{
		Num_Generations_LookupTable.resize(max_N + 1);
		Generations_LookupTable.resize(max_N + 1);

		// Empty generation
		Num_Generations_LookupTable[0] = {0};
		Generations_LookupTable[0] = std::vector<gentype>();

		// Count-compression of all possible ranked offpsring generations.
		for (int n = 1; n <= max_N; n++)
		{
			for (gentype g : generation_span(n))
			{
				std::map<ushort, ushort> freq_map = {};
				gentype counts = {};
				unsigned multiplicities = factorial(n);

				// I take advantage of the order granted by the map
				for (ushort val : g)
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

				// assign counts
				for (const auto &x : freq_map)
				{
					counts.push_back(x.second);
					multiplicities /= factorial(x.second);
				}

				// add to the list
				Generations_LookupTable[n].push_back(counts);
				Num_Generations_LookupTable[n].push_back(multiplicities);
			}
		}
	}

	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//

	/*
	Algorithm to generate all possible combinations of N elements in k positions, plus the N - k complements.
	The elements are indexed from 1 to N.

	----
	Ref:
	----
	Daniel Angel Jimenez
	https://www.cs.utexas.edu/users/djimenez/utsa/cs3343/lecture25.html
	*/
	void recursive_combinations_span(ushort N, ushort k, vec_ushort &current, ushort start, ushort loop, std::vector<vec_ushort> &results)
	{
		if (loop > k)
		{
			//-----complement of the combination-----//
			// Mark elements present in the list
			std::vector<bool> missing(N + 1, true);
			for (ushort i = 0; i < k; i++)
			{
				missing[current[i]] = false;
			}
			// store complement
			for (ushort i = 1; i <= N; i++)
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
			for (ushort i = start; i <= N; i++)
			{
				current[loop - 1] = i;
				recursive_combinations_span(N, k, current, i + 1, loop + 1, results);
			}
		}
	}

	/* Store all combinations + complement of N elements in k positions. */
	void store_combinations(const ushort N, const ushort k)
	{
		// Generate combinations with the recurisve algorithm
		std::vector<vec_ushort> results;
		vec_ushort current(N, 0); // initialize to impossible value (not in [1, N])
		recursive_combinations_span(N, k, current, 1, 1, results);

		// Store all results inside the Lookup Table
		unsigned c = 0;
		for (vec_ushort current : results)
		{
			Combinations_LookupTable[N][k][c].resize(N);

			// store current combination + complement
			for (ushort i = 0; i < N; i++)
			{
				Combinations_LookupTable[N][k][c][i] = current[i];
			}
			c++;
		}
	}

	/* Initialize the Look Up Table of combinations and their number */
	void init_combinations_lookuptable(const ushort max_N)
	{
		unsigned num_combin;

		Combinations_LookupTable.resize(max_N + 1);
		Num_Combinations_LookupTable.resize(max_N + 1);

		for (ushort N = 0; N <= max_N; N++)
		{
			Combinations_LookupTable[N].resize(N);
			Num_Combinations_LookupTable[N].resize(N);

			for (ushort k = 0; k < N; k++)
			{
				num_combin = binomial_coefficient(N, k);
				Combinations_LookupTable[N][k].resize(num_combin);
				Num_Combinations_LookupTable[N][k] = num_combin;

				store_combinations(N, k);
			}
		}
	}

	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//

	/*
	Compute the likelihood for a heterogeneous Reed-Frost model.
	Arguments :
	-----
		N : household size
		k : number of seropositive
		indidx : [0 - k) seropositive indexes, [k - N) seronegative indexes. Indexes are [1, N]
		escape_matrix : probability 'col' escaped from 'row' (including community as index 0)
	*/
	double ll_chain_binomial_heterogeneous(ushort N, ushort k, vec_ushort &indidx, std::vector<double> &escape_matrix)
	{

		// FLAG
		std::cerr << "INSIDE ll_chain_binomial_heterogeneous "
				  << "N" << N << " k" << k << std::endl;
		// FLAG-END

		// product over all infectious of escape matrix for the escaped individuals
		double esc_product = 1.;
		for (ushort i = k; i < N; i++)
		{
			// from community
			esc_product *= escape_matrix[indidx[i] - 1];

			// from infected
			for (ushort j = 0; j < k; j++)
			{
				esc_product *= escape_matrix[indidx[i] - 1 + N * indidx[j]];
			}
		}

		// submatrix has always at least 1 infector: community
		std::vector<double> sub_escape_matrix((k + 1) * (k));
		// find the submatrix of the given infectious if any
		// community (inf=0) is always at position (row=0)
		ushort row = 0, inf = 0;
		while ((k > 0) && (row <= k))
		{
			for (ushort col = 0; col < k; col++)
			{
				sub_escape_matrix[col + k * row] = escape_matrix[indidx[col] - 1 + N * inf];
			}
			// next infector
			inf = indidx[row];
			row++;
		}

		return esc_product * recursive_function(k, sub_escape_matrix);
	}

	/* Recursive function to break likelihood into shorter escape patterns. */
	double recursive_function(ushort N, std::vector<double> &escape_matrix)
	{
		double output = 1.;

		for (ushort k = 0; k < N; k++)
		{
			for (unsigned c = 0; c < Num_Combinations_LookupTable[N][k]; c++)
			{
				output -= ll_chain_binomial_heterogeneous(
					N, k, Combinations_LookupTable[N][k][c], escape_matrix);
			}
		}

		return output;
	}
}