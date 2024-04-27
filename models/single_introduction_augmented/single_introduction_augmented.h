// single_introduction_augmented.h //

#ifndef MODEL__H
#define MODEL__H

#include "mcmc.h"

//#define TEST_EPI_CURVE

//=====Forward declaration of the Infectivity Profile Function=====//
class Infectivity_Profile_Function;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////// Model : class /////////////////////////////////////////// //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

class Model : public MCMC
{
protected:
	double *inst_profile_lookuptab;
	double *cumul_profile_lookuptab;

public:

#ifdef TEST_EPI_CURVE
	void build_extra(std::string extra_pathtofile ) override;
#endif

	//=====Infectivity Profiles=====//
	Infectivity_Profile_Function *infec_profile_func;
	double compute_infec_profile_inst(int i, double t_transm);
	double compute_infec_profile_cumul(int i, double t_transm);

	//=====Rate Choices=====//
	void assign_rate_per_group_to_individual(int i) override;
	std::string assign_rate_per_contact(int i, int k) override;

	void store_augm_var_dependecies(int i) override;
	double get_infec_time(int i) const override;
	double compute_single_house_loglike(int h) override;
	void initialize_model() override;

	void simulate_single_house(int h) override;
	void format_simulation() override;

	//=====LogLikelihood=====//
	double _sum_hazards_(double *profiles, int idx);
};

#endif // MODEL__H