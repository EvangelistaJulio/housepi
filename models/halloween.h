// halloween.h //

#ifndef MODEL__H
#define MODEL__H

#include "mcmc.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////// Model : class /////////////////////////////////////////// //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

class Model : public MCMC
{
protected:
	double *inst_profile_lookuptab;
	double *cumul_profile_lookuptab;

public:

	void assign_infec_profile_function() override;
	void assign_rate_per_group_to_individual(int i) override;
	void store_augm_var_dependecies(int i) override;
	void initialize_model() override;

	double get_infec_time(int i) const override;
	double compute_single_house_loglike(int h) override;

	void simulate_infection() override;

	//=====Optional Methods=====//
	void assign_rate_per_contact(int i, int k) override;

	//=====Infectivity Profiles=====//
	double compute_infec_profile_inst(int i, double t_transm);
	double compute_infec_profile_cumul(int i, double t_transm);

	//=====LogLikelihood=====//
	double _sum_hazards_(double *profiles, int idx);
};

#endif // MODEL__H