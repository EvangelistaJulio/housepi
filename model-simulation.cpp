// model-simulation.cpp //

#include "model.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////// Model : Methods (simulation)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/*
Poisson point process to simulate data.
	Initialize all individuals to healthy
	Except the first case within household
	SImulate the Poisson process
	FIXME : is this correct? first onset maybe...

	Also, what about different in-house study-periods?
*/
void Model::simulate_infection(timetype dt, mt19937_64 &rng)
{
	timetype time_end;
	uniform_real_distribution<double> unif(0.0, 1.0);

	// simulate the infection for each household
	for (size_t h = 0; h < num_househ; h++)
	{
		time_end = 0;

		//=====INITIALIZATION=====//
		// generate new incubation times and first case
		for (size_t i : househ[h].healthy_list)
		{
			indiv[i].init_augm("incub_period", rng);
		}
		for (size_t i : househ[h].symptomatic_list)
		{
			indiv[i].init_augm("incub_period", rng);
			time_end = max(time_end, indiv[i].get("study_period"));
		}
		update_first_case(h);

		// set all secondary cases to healthy
		vector<size_t> copy_symptomatic_list = househ[h].symptomatic_list;
		for (size_t i : copy_symptomatic_list)
		{
			if (i != househ[h].first_case)
			{
				indiv[i].set_status(HEALTHY);
				househ[h].change_status(i, HEALTHY);
			}
		}
		update_infectivity_profiles(househ[h].first_case);

		// FLAG
		cerr << "infec_time = " << househ[h].infec_time << endl;
		// FLAG

		//=====simulate infection within household=====//
		for (timetype t = househ[h].infec_time; t < time_end; t += dt)
		{
			if (househ[h].healthy_list.empty())
			{
				t = time_end;
			}
			else
			{
				vector<size_t> copy_healthy_list = househ[h].healthy_list;
				for (size_t i : copy_healthy_list)
				{
					// NOTE : assign infec_time from onset_time
					indiv[i].set_status(SYMPTOMATIC);
					indiv[i].set_augm("onset_time", t + dt * unif(rng) + indiv[i].get_augm("incub_period"));
					househ[h].change_status(i, SYMPTOMATIC);
					update_infectivity_profiles(i);

					// FIXME : in this way onset times could anticipate...

					// FLAG
					cerr << "t = " << t << " | prob (" << i << ") = " << exp(log_infec_proba(i)) << endl;
					// FLAG

					// refuse
					if (unif(rng) >= exp(log_infec_proba(i)))
					{
						indiv[i].set_status(HEALTHY);
						househ[h].change_status(i, HEALTHY);
						update_infectivity_profiles(i);
					}
				}
			}
			// FLAG
			cerr << endl;
			// FLAG
		}
	}

	// format
	/*
		for (size_t i = 0; i < num_indiv; i++)
	{
		if (indiv[i].get_augm("onset_time") > indiv[i].get("study_period"))
		{
			indiv[i].set_status(HEALTHY);
		}
		else
		{
			indiv[i].set_augm("onset_time", round(indiv[i].get_augm("onset_time")));
		}

		if (indiv[i].get_status() == HEALTHY)
		{
			indiv[i].set_augm("onset_time", 1000);
			indiv[i].set_augm("incub_period", 0);
		}
	}
	*/
}