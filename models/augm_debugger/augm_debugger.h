// augm_debugger.h //

#ifndef MODEL__H
#define MODEL__H

#include "mcmc.h"

//=====Forward declaration of the Infectivity Profile Function=====//
class Infectivity_Profile_Function;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////// Model : class /////////////////////////////////////////// //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

class Model : public MCMC
{
protected:
	std::vector<double> inst_profile_LUT = {};
	std::vector<double> cumul_profile_LUT = {};

public:
	//=====Model Specifc Methods =====//
	Infectivity_Profile_Function *infec_profile_func;
	double compute_infec_profile_inst(int i, double t_transm);
	double compute_infec_profile_cumul(int i, double t_transm);
	double hazards_of_infec(std::vector<double> *profiles, int idx);

	void augm_var_store_dependecies(int i) override;
	double get_infec_time(int i) const override;
	double compute_single_house_loglike(int h) override;
	void initialize_model() override;

	void simulate_single_house(int h) override;
	void format_augmented_data() override;

};

#endif // MODEL__H