// auxiliary_binomial.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#include "auxiliary_binomial.h"

// Binomial Auxiliary Namespace
namespace bax
{
	/////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////  GLOBAL VARIABLES  ///////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////

	//===== Transmission Generations : Lookup Tables =====//
	// numbers of generations
	std::vector<unsigned> Num_Gener_LookupTable;
	// n. of multiplicites per generation
	std::vector<std::vector<unsigned>> Multipl_Gener_LookupTable;
	// generation vectors
	std::vector<std::vector<gentype>> Generations_LookupTable;

	//===== Combinations : Lookup Tables =====//
	// lookup-table to store the number of combinations
	std::vector<std::vector<unsigned>> Num_Combin_LookupTable;
	// lookup-table to store combinations and complements
	std::vector<std::vector<std::vector<std::vector<ushort>>>> Combinations_LookupTable;

	//===== EscapePatterns : Lookup Tables =====//
	// lookup-table to store the number of combinations
	std::vector<EscapePattern> ChainBinomial_LookupTable;
	// lookup-table to store combinations and complements
	std::vector<std::vector<EscapePattern>> EscapePatterns_LookupTable;

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////  COMBINATORICS FUNCTIONS  ////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////

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

	/////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////  Escape Pattern : Class   ///////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////

	/* Compute the loglikelihood starting from the current escape pattern.
		Important : it requires initialization of the lookup table. */
	double EscapePattern::likelihood(std::vector<double> &esc_mtx)
	{
		double output = 1.;

		for (int p : pattern_idxs)
		{
			output -= EscapePatterns_LookupTable[N_star][p].likelihood(esc_mtx);
		}

		for (ushort i : entries)
		{
			output *= esc_mtx[i];
		}

		return output;
	}

	/* Print the formula of the loglikelihood. */
	std::string EscapePattern::show()
	{
		std::ostringstream oss;

		for (size_t c = 0; c < entries.size(); c++)
		{
			ushort u = entries[c];
			if (u < N_star)
			{
				oss << "B" << u + 1;
			}
			else
			{
				oss << "Q" << u / N_star << "_" << u % N_star + 1;
			}

			if (c < entries.size() - 1)
			{
				oss << "*";
			}
		}

		if (!pattern_idxs.empty())
		{
			if (!entries.empty())
			{
				oss << "*(";
			}

			oss << "1";
			for (int p : pattern_idxs)
			{
				oss << " - " << EscapePatterns_LookupTable[N_star][p].show();
			}

			if (!entries.empty())
			{
				oss << ")";
			}
		}

		return oss.str();
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////  TRANSMISSION GENERATIONS  ///////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////

	/* All possible ranked offspring generations for a given number of infected individuals n > 0.

		e.g.
		----
		n_infec = 3
		return [[0, 0, 0], [0, 0, 1], [0, 1, 1], [0, 1, 2]]
	*/
	std::vector<gentype> recursive_span_generation(ushort n)
	{
		ushort root = n - 1;
		std::vector<gentype> output(1, gentype(n, 0));

		while (root > 0)
		{
			for (gentype g : recursive_span_generation(n - root))
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
		// Allocate LUT memory
		Num_Gener_LookupTable.resize(max_N + 1);
		Multipl_Gener_LookupTable.resize(max_N + 1);
		Generations_LookupTable.resize(max_N + 1);

		// Empty generation n = 0
		// note : just for access clarity
		Num_Gener_LookupTable[0] = 0;
		Multipl_Gener_LookupTable[0].push_back(0);
		Generations_LookupTable[0].push_back({});

		// Count-compression of all possible ranked offpsring generations.
		size_t tot_num_gen = 0;
		std::vector<unsigned> aux_num_lut;
		std::vector<gentype> aux_gen_lut;

		for (int n = 1; n <= max_N; n++)
		{
			aux_num_lut.clear();
			aux_gen_lut.clear();

			for (gentype g : recursive_span_generation(n))
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
				aux_num_lut.push_back(multiplicities);
				aux_gen_lut.push_back(counts);
			}

			// allocate lookup table
			// note : how to obtain tot_num_gen a priori
			tot_num_gen = aux_gen_lut.size();
			Num_Gener_LookupTable[n] = tot_num_gen;
			Multipl_Gener_LookupTable[n].resize(tot_num_gen);
			Generations_LookupTable[n].resize(tot_num_gen);

			// copy from auxiliary vectors
			for (size_t c = 0; c < tot_num_gen; c++)
			{
				Multipl_Gener_LookupTable[n][c] = aux_num_lut[c];
				Generations_LookupTable[n][c] = aux_gen_lut[c];
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
	Ref: (modified from)
	----
	Daniel Angel Jimenez
	https://www.cs.utexas.edu/users/djimenez/utsa/cs3343/lecture25.html
	*/
	void recursive_span_combinations(ushort N, ushort k, vec_ushort &current, std::vector<vec_ushort> &results, ushort start, ushort loop)
	{
		if (loop > k - 1)
		{
			//----- complement of the combination -----//
			// Mark elements present in the list
			std::vector<bool> missing(N, true);
			for (ushort i = 0; i < k; i++)
			{
				missing[current[i]] = false;
			}
			// store complement
			for (ushort i = 0; i < N; i++)
			{
				if (missing[i])
				{
					current[loop] = i;
					loop++;
				}
			}

			//----- append current combination & complement to results -----//
			results.push_back(current);
		}
		else
		{
			for (ushort i = start; i < N; i++)
			{
				current[loop] = i;
				recursive_span_combinations(N, k, current, results, i + 1, loop + 1);
			}
		}
	}

	/* Store all combinations + complement of N elements in k positions. */
	void store_combinations(const ushort N, const ushort k)
	{
		// Generate combinations with the recurisve algorithm
		std::vector<vec_ushort> results;
		vec_ushort current(N);
		recursive_span_combinations(N, k, current, results, 0, 0);

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

	/* Initialize the Look Up Table of combinations and their counts */
	void init_combinations_lookuptable(const ushort max_N)
	{
		unsigned num_combin;

		Num_Combin_LookupTable.resize(max_N + 1);
		Combinations_LookupTable.resize(max_N + 1);

		for (ushort N = 0; N <= max_N; N++)
		{
			Num_Combin_LookupTable[N].resize(max_N);
			Combinations_LookupTable[N].resize(max_N);

			for (ushort k = 0; k < N; k++)
			{
				num_combin = binomial_coefficient(N, k);
				Num_Combin_LookupTable[N][k] = num_combin;
				Combinations_LookupTable[N][k].resize(num_combin);
				store_combinations(N, k);
			}
		}
	}

	/* Clear the Look Up Table of combinations and their counts */
	void clear_combinations_lookuptable()
	{
		Num_Combin_LookupTable.clear();
		Num_Combin_LookupTable.shrink_to_fit();
		Combinations_LookupTable.clear();
		Combinations_LookupTable.shrink_to_fit();
	}

	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//

	/*  */
	void recursive_escape_patterns_function(ushort N_star, vec_ushort access_matrix, EscapePattern &father, std::vector<EscapePattern> &results)
	{
		ushort N = (std::sqrt(4 * access_matrix.size() + 1) - 1) / 2;

		for (ushort k = 0; k < N; k++)
		{
			for (unsigned c = 0; c < Num_Combin_LookupTable[N][k]; c++)
			{
				EscapePattern new_child(N_star);

				// associate entries for this pattern
				std::vector<ushort> &indidx = Combinations_LookupTable[N][k][c];
				for (ushort i = k; i < N; i++)
				{
					//  community
					new_child.entries.push_back(access_matrix[indidx[i]]);

					// person-to-person
					for (ushort j = 0; j < k; j++)
					{
						new_child.entries.push_back(access_matrix[indidx[i] + N * (indidx[j] + 1)]);
					}
				}

				if (k > 0)
				{
					// find the submatrix of the given infectious
					// note: always at least 1 infector (community)
					vec_ushort sub_access_matrix((k + 1) * k);

					ushort inf = 0; // community (inf=0) always at row=0
					for (ushort row = 0; row <= k; row++)
					{
						for (ushort col = 0; col < k; col++)
						{
							sub_access_matrix[col + k * row] = access_matrix[indidx[col] + N * inf];
						}
						inf = indidx[row] + 1; // next infector
					}

					recursive_escape_patterns_function(N_star, sub_access_matrix, new_child, results);
					sub_access_matrix.clear();
				}

				// attach pattern index
				father.pattern_idxs.push_back(results.size());
				// store pattern
				results.push_back(new_child);
			}
		}
	}

	/* Initialize the Look Up Table of escape patterns and the access object. */
	void init_escape_patterns_lookuptable(const ushort max_N)
	{
		//----- lookuptable with all combinations + complements -----//
		init_combinations_lookuptable(max_N);

		//----- lookuptable with escape patterns linked to chain-binomial -----//
		// case N_star == 0, for consistency
		ChainBinomial_LookupTable.emplace_back(0);
		EscapePatterns_LookupTable.push_back(std::vector<EscapePattern>());

		// case N_star > 0
		for (int N_star = 1; N_star <= max_N; N_star++)
		{
			vec_ushort access_matrix((N_star + 1) * N_star);
			for (int i = 0; i < (N_star + 1) * N_star; i++)
			{
				access_matrix[i] = i;
			}

			// store to lookup table
			ChainBinomial_LookupTable.emplace_back(N_star);
			EscapePatterns_LookupTable.push_back({});
			recursive_escape_patterns_function(N_star, access_matrix,
											   ChainBinomial_LookupTable[N_star],
											   EscapePatterns_LookupTable[N_star]);
		}

		// deallocate memory of combinations lookup table
		clear_combinations_lookuptable();
	}

	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//
	// -----------------------------------------------------------------------------------------------------------//

	/* Heterogeneous chain binomial loglikelihood for ordered probabilities (seropositives, seronegatives).

		Arguments
		---
		N_star : n. of seropositves
		commu_esc_proba : ordered community escape probabilities 
		p2p_esc_proba : ordered person-to-person escape probabilities 
	*/
	double chain_binomial_heter_loglikelihood(ushort N_star,
											  std::vector<double> &commu_esc_proba,
											  std::vector<std::vector<double>> &p2p_esc_proba)
	{
		ushort hh_size = commu_esc_proba.size();
		double loglikelihood = 0.;

		//----- seronegatives contribution -----//
		for (ushort i = N_star; i < hh_size; i++)
		{
			// community : i escaped from community
			loglikelihood += std::log(commu_esc_proba[i]);

			// person-to-person : i escaped from k
			for (ushort k = 0; k < N_star; k++)
			{
				loglikelihood += std::log(p2p_esc_proba[k][i]);
			}
		}

		if (N_star == 1)
		{ //----- single seropositive contribution -----//
			loglikelihood += std::log(1. - commu_esc_proba[0]);
		}
		else if (N_star > 1)
		{ //----- seropositives contribution -----//

			// flatten sub-matrix with only seropositives
			std::vector<double> escape_matrix((N_star + 1) * N_star, 1);
			for (ushort i = 0; i < N_star; i++)
			{
				// community : i escaped from community
				escape_matrix[i] = commu_esc_proba[i];

				// person-to-person : i escaped from k
				for (ushort k = 0; k < N_star; k++)
				{
					escape_matrix[N_star * (k + 1) + i] = p2p_esc_proba[k][i];
				}
			}

			loglikelihood += std::log(ChainBinomial_LookupTable[N_star].likelihood(escape_matrix));
		}

		return loglikelihood;
	}
}

/* Debug mode, to check values of the likelihood and print formulation.
	g++ -Wall -Wextra -std=gnu++11 -g -O2 -o auxiliary_binomial.o auxiliary_binomial.cpp

#include <chrono>
int main(int argc, char *argv[])
{
	int max_num = std::stoi(argv[1]);
	init_escape_patterns_lookuptable(max_num);

	for (int N_star = 0; N_star <= max_num; N_star++)
	{
		std::cout << "\n<<" << N_star << ">>";
		std::cout << " n.terms = " << 1 + (EscapePatterns_LookupTable[N_star].size() - 1) / 2 << std::endl;

		//----- elements for the escape matrix -----//
		std::vector<double> commu_esc_proba(N_star);
		for (ushort i = 0; i < N_star; i++)
		{
			commu_esc_proba[i] = 0.15 + 0.025 * i;
			std::cout << "B" << i + 1 << " = " << commu_esc_proba[i] << std::endl;
		}

		std::vector<std::vector<double>> p2p_esc_proba(N_star);
		for (ushort row = 0; row < N_star; row++)
		{
			p2p_esc_proba[row].resize(N_star);
			for (ushort col = 0; col < N_star; col++)
			{
				p2p_esc_proba[row][col] = 0.25 + 0.0125 * col + 0.0125 * row;
				std::cout << "Q" << row + 1 << "_" << col + 1 << " = " << p2p_esc_proba[row][col] << std::endl;
			}
		}

		if (N_star < 5)
		{
			std::cerr << ChainBinomial_LookupTable[N_star].show() << std::endl;
		}
		auto start = std::chrono::high_resolution_clock::now();
		std::cerr << chain_binomial_heter_loglikelihood(N_star, commu_esc_proba, p2p_esc_proba);
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = end - start;
		std::cout << " | time: " << elapsed.count() << " seconds" << std::endl;
	}
}
*/