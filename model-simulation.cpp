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
void Model::simulate_infection(timetype dt, std::mt19937_64 &rng)
{
	timetype time_end;
	std::uniform_real_distribution<double> unif(0.0, 1.0);

	// simulate the infection for each household
	for (size_t h = 0; h < num_househ; h++)
	{
		time_end = 0;
		std::vector<size_t> copy_list = househ[h].symptomatic_list;
		for (size_t i : copy_list)
		{
			// set all secondary cases to healthy
			if (i != househ[h].first_case)
			{
				indiv[i].set_status(HEALTHY);
				indiv[i].set_augm("incub_period", 0);
				indiv[i].set_augm("onset_time", indiv[i].get("study_period"));
				househ[h].change_status(i, HEALTHY);
			}
			time_end = std::max(time_end, indiv[i].get("study_period"));
		}
		update_infectivity_profiles(househ[h].first_case);

		//=====simulate infection within household=====//
		for (timetype t = househ[h].infec_time; t < time_end; t += dt)
		{
			if (househ[h].healthy_list.empty())
			{
				t = time_end;
			}
			else
			{
				std::vector<size_t> copy_list = househ[h].healthy_list;
				for (size_t i : copy_list)
				{
					// NOTE : assign infec_time from onset_time
					indiv[i].set_status(SYMPTOMATIC);
					indiv[i].init_augm("incub_period", rng);
					indiv[i].set_augm("onset_time", t + indiv[i].get_augm("incub_period"));
					househ[h].change_status(i, SYMPTOMATIC);

					if (unif(rng) < std::exp(log_infec_proba(i)))
					{
						update_infectivity_profiles(i);
					}
					else
					{
						indiv[i].set_status(HEALTHY);
						indiv[i].set_augm("incub_period", 0);
						indiv[i].set_augm("onset_time", time_end);
						househ[h].change_status(i, HEALTHY);
					}
				}
			}
		}
	}
}