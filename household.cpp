#include "household.hpp"

using namespace std;

/////////////////////////////////////////////////////
//  Individual : Class
/////////////////////////////////////////////////////

// ----------------------------
//  Constructors
// ----------------------------

Individual::Individual (
	int indid,
	int is_case,
	int age,
	int sex,
	bool isolation,
	bool recovered,
	int vacc_status,
	double study_period,
	double vacc_delay,
	double onset_time
	) {
	// Data variables
	i_indid = indid;
	i_is_case = is_case;
	i_age = age;
	i_sex = sex;
	i_isolation = isolation;
	i_recovered = recovered;
	i_vacc_status = vacc_status;
 	i_study_period = study_period;
	i_vacc_delay = vacc_delay;
	i_onset_time = onset_time;

	// Initialize augmented variables
	i_augm_infec_time = onset_time;
	i_augm_onset_time = onset_time;
}

// ----------------------------
//  Methods : Computations
// ----------------------------

double Individual::infec_inst (double t_infec, double study_period) {
	double out = 0.0, tost, t_incub;
	std::vector<double> args = {
		_param_value("tost_tau"), _param_value("tost_mu"),
		_param_value("tost_sigma"), _param_value("tost_alpha")};

	// NOTE : check for Symptomatic infector inside Household
	if (i_augm_infec_time < study_period && i_augm_infec_time < t_infec) {
		tost = t_infec - i_augm_onset_time;
		t_incub = i_augm_onset_time - i_augm_infec_time;
		out = t_onset_sympt_to_transm_pdf(tost, t_incub, args); 
	} 

	return out;
}

double Individual::infec_cumul (double t_infec, double study_period) {
	double out = 0.0, tost, t_incub;
	std::vector<double> args = {
		_param_value("tost_tau"), _param_value("tost_mu"),
		_param_value("tost_sigma"), _param_value("tost_alpha")};

	// NOTE : check for Symptomatic infector inside Household
	if (i_augm_infec_time < study_period && i_augm_infec_time < t_infec) {
		tost = t_infec - i_augm_onset_time;
		t_incub = i_augm_onset_time - i_augm_infec_time;
		out = t_onset_sympt_to_transm_cdf(tost, t_incub, args); 
	} 

	return out;
}

/////////////////////////////////////////////////////
//  Household : Class
/////////////////////////////////////////////////////

// ----------------------------
//  Methods : Getter
// ----------------------------

size_t Household::get_size () const {
	return h_residents.size();
}

// number of symptomatic
size_t Household::num_symptomatic () const {
	size_t num_symptomatic = 0;

	for (size_t i = 0; i < h_residents.size(); i++) {
		if (h_residents[i].i_onset_time < INFINITY) {
			num_symptomatic++;
		}
	}

	return num_symptomatic;
}

// get augmented onset time for given resident
double Household::get_augm_onset_time(size_t idx) const {
	return h_residents[idx].i_augm_onset_time;
}

// ----------------------------
//  Methods : Setter
// ----------------------------

// Add individual to household
void Household::add_individual (Individual this_indiv) {
	
	h_residents.push_back(this_indiv);
	h_infec_time = std::min(h_infec_time, this_indiv.i_augm_infec_time);
	h_size++; 
}


// Initialize augmented onset time
void Household::init_onset_time (size_t idx, std::mt19937_64& rand_gen) {

	if (h_residents[idx].i_is_case > 0) { 
		h_residents[idx].i_augm_onset_time = h_residents[idx].i_onset_time + unif_rand(rand_gen);
	
	} else {
		// Household residents without infection or with unknown final outcome
		h_residents[idx].i_augm_onset_time = INFINITY; 
	}

	// --- Display --- //	
	if (Verbose == true) {
		std::cout << "Observed onset time: " << h_residents[idx].i_onset_time << std::endl;
		std::cout << "Augmented onset time" << h_residents[idx].i_augm_onset_time << std::endl;
  	} 
	// --- End-of-Display --- //
}

// Initialize augmented infection time
void Household::init_infec_time (size_t idx, std::mt19937_64& rand_gen) {
	double incub_period = 0.0;
	std::vector<double> args = {_param_value("shape_incub"), _param_value("scale_incub")};

	if (h_residents[idx].i_is_case == 1) {   // Symptomatic case 
		incub_period = gamma_rand(rand_gen, args);
		h_residents[idx].i_augm_infec_time = h_residents[idx].i_augm_onset_time - incub_period;
	
	} else {    
		// Household residents without infection or with unknown final outcome
		h_residents[idx].i_augm_infec_time = INFINITY; 
	}

	// --- Display --- //	
	if (Verbose == true) {
		std::cout << "Inf time " << h_residents[idx].i_indid << ": " << h_residents[idx].i_augm_infec_time << std::endl;
	}
	// --- End-of-Display --- //
}

// ----------------------------
//  Methods : Computations
// ----------------------------

// Change symptom onset // FIXME : generalize augmented
void Household::set_augm_onset_time (size_t idx, double onset_time) {
	h_residents[idx].i_augm_onset_time = onset_time;
	update_infec_profiles(idx);
}

// Change infection time
void Household::set_augm_infec_time (size_t idx, double infec_time) {
	h_residents[idx].i_augm_infec_time = infec_time;
	update_infec_profiles(idx);
	h_infec_time = std::min(h_infec_time, infec_time);
}

// Compute the risk of infection
void Household::init_infec_profiles() {

	// initialize 2-dim vectors h_cumul_Lambda and h_inst_Lambda
	h_cumul_Lambda.resize(h_size);
	h_inst_Lambda.resize(h_size);
	for (size_t i = 0; i < h_size; i++) {
		h_cumul_Lambda[i].resize(h_size);
		h_inst_Lambda[i].resize(h_size);
	}

	// Compute cumulative and instantaneous person-to-person transmission rates
	for (size_t infector = 0; infector < h_size; infector++) {
		// FIXME : is this working ?
		update_infec_profiles(infector);
	}

	// --- Display --- //	
	if (Verbose == true) {
		// Cumulative transmission rate
		display_2D_vector(h_cumul_Lambda);

		// Instantaneous transmission rate
		display_2D_vector(h_inst_Lambda);

		// Infection time 
		for (size_t i = 0; i < h_size; i++) {
			std::cout << h_residents[i].i_augm_infec_time << " ";
		}
		std::cout << std::endl;
	}
	// --- End-of-Display --- //	
}

// Update h_cumul_Lambda and h_inst_Lambda matrices during data augmentation
void Household::update_infec_profiles(size_t idx) {
	double t1;
	// FIXME : to be checked

	if (h_residents[idx].i_is_case > 0){				// resident idx is an infector

		for (size_t infectee = 0; infectee < h_size; infectee++) {
			t1 = h_residents[infectee].i_augm_infec_time;

			if (h_residents[idx].i_augm_infec_time < t1){ // NOTE : idx != infectee 
				// Cumulative transmission rate for all household contacts
				h_cumul_Lambda[idx][infectee] = h_residents[idx].infec_cumul(t1, h_residents[infectee].i_study_period);	
				
				// Instantaneous transmission rate for secondary cases
				if (h_residents[infectee].i_is_case > 0) {
					h_inst_Lambda[idx][infectee] = h_residents[idx].infec_inst(t1, h_residents[infectee].i_study_period);
				}
			}			
		}

	} else if (h_residents[idx].i_is_case == 0) {		// resident idx is an infectee
		t1 = h_residents[idx].i_study_period;

	} else { // FIXME: there is not else
		t1 = h_residents[idx].i_augm_infec_time;    
	}

	for (size_t infector = 0; infector < h_size; infector++) {    
		// FIXME : idx != infector

		if (h_residents[infector].i_is_case > 0){
			// Cumulative transmission rate for all household contacts
			if (h_residents[infector].i_augm_infec_time < t1) { 
				h_cumul_Lambda[infector][idx] = h_residents[infector].infec_cumul(t1, h_residents[idx].i_study_period);
			}

			// Instantaneous transmission rate for secondary cases
			if (h_residents[infector].i_augm_infec_time < h_residents[idx].i_augm_infec_time && h_residents[idx].i_is_case > 0) {
				h_inst_Lambda[infector][idx] = h_residents[infector].infec_inst(h_residents[idx].i_augm_infec_time, h_residents[idx].i_study_period);
			}
		}
	}

	// --- Display --- //	
	if (Verbose == true) {
		display_2D_vector(h_cumul_Lambda);

		for (size_t i = 0; i < h_size; i++) {
			std::cout << h_residents[i].i_augm_infec_time << " ";
		}
		std::cout << std::endl;
	}
	// --- End-of-Display --- //	
}


double Household::propose_onset_time (size_t idx, std::mt19937_64& rand_gen) {
	// Propose new value for symptom onset time
	double new_onset_time = 0.0;		

	if (h_residents[idx].i_is_case > 0) {           
		// is_case household member 
		new_onset_time = h_residents[idx].i_onset_time + unif_rand(rand_gen);
		
	} else {    // Resident without infection or with unknown outcome
		new_onset_time = INFINITY;
	}

	return new_onset_time;
}

double Household::propose_infec_time (size_t idx, std::mt19937_64& rand_gen) {
	// Propose new value for infection time
	double new_infec_time = 0.0, incub_period = 0.0;
	std::vector<double> args = {_param_value("shape_incub"), _param_value("scale_incub")};

	if (h_residents[idx].i_is_case == 1) {		// Symptomatic case 
		incub_period = gamma_rand(rand_gen, args);
		new_infec_time = h_residents[idx].i_augm_onset_time - incub_period;

	} else {  	// Resident without infection or with unknown outcome
		new_infec_time = INFINITY;
	}

	return new_infec_time;
}


/* Compute the log-likelihood within household */
double Household::compute_loglike () {
	size_t idx_first_case = 0;
	double LL = 0.0; // Log Likelihood
	double t0 = h_infec_time;

	// Identify idx of first case
	for (size_t i = 0; i < h_size; i++) {
		if (h_residents[i].i_augm_infec_time == t0) {
			idx_first_case = i;
		}
	}

	// All residents
	for (size_t i = 0; i < h_size; i++) { 
		
		if (i == idx_first_case) {     					// Contribution of first case                                      		
			LL += log_incub_proba(i); 

		} else if (h_residents[i].i_is_case > 0) {      // Contribution of secondary cases
			LL += log_incub_proba(i);
	  		LL += log_surv_proba(i, t0);
			LL += log_infec_proba(i);

		} else {										// Contribution of not infected cases
			LL += log_surv_proba(i, t0); 
		}

		// --- Display --- //
		if (Verbose == true) {
			std::cout << i << " " << idx_first_case << " " << h_residents[i].i_augm_onset_time;
			std::cout << " " << h_residents[i].i_augm_infec_time << " " << h_residents[i].i_is_case << std::endl; 
		}
		// --- End-of-Display --- //
	}

	// --- Display --- //
	if (Verbose == true) {
		std::cout << "LL: " << LL << "\n" << std::endl;
	}
	// --- End-of-Display --- //

  return LL;
}


double Household::log_incub_proba (size_t idx) {
	// logarithm of the pdf of the incubation period
	double out = 0.0, incub_period = 0.0;
	std::vector<double> args = {_param_value("shape_incub"), _param_value("scale_incub")};

	if (h_residents[idx].i_is_case == 1) {	// Symptomatic case with known symptom onset 
		incub_period = h_residents[idx].i_augm_onset_time - h_residents[idx].i_augm_infec_time;
		out = log_gamma_pdf(incub_period, args);
	} 

	// --- Display --- //
	if (Verbose == true) {
		std::cout << "log_incub_proba: " << out << std::endl;
	}
	// --- End-of-Display --- //	

	return out;
}

// Instaneous risk of infection from community
double Household::model_community_risk_inst (double t1){
	double risk_inst = _param_value("alpha") * Inst_Community_Infec[(int) t1]; 

	return risk_inst;
}

// Cumulative risk of infection from community
double Household::model_community_risk_cumul (double t0, double t1){
	double risk_cumul = _param_value("alpha") * (Cumul_Community_Infec[(int) t1] - Cumul_Community_Infec[(int) t0]); 

	return risk_cumul;
}

double Household::model_contact_modulation (size_t idx, size_t infectee){
	// Relative infectivity of symptomatic child cases versus symptomatic adult cases
	double factor = 1.0;
	size_t i = infectee;

	// NOTE : this call is missleading : use only if to be inferred...

	// Child (2) - Child (2)
	if (_param_infer("rCC") == true && h_residents[idx].i_age == 0 && h_residents[i].i_age == 0) {
		factor *= _param_value("rCC");

	// Child (2) - Mother (0)
	} else if (_param_infer("rMC") == true && ((h_residents[idx].i_age == 0 && h_residents[i].i_age == 1 && h_residents[i].i_sex == 1) || (h_residents[i].i_age == 0 && h_residents[idx].i_age == 1 && h_residents[idx].i_sex == 1))) {
		factor *= _param_value("rMC");

	// Child (2) - Father (1)
	} else if (_param_infer("rFC") == true && ((h_residents[idx].i_age == 0 && h_residents[i].i_age == 1 && h_residents[i].i_sex == 2) || (h_residents[i].i_age == 0 && h_residents[idx].i_age == 1 && h_residents[idx].i_sex == 2))) {
		factor *= _param_value("rFC");
	}

	return factor;
}


double Household::model_susceptibility (size_t idx){
	double factor = 1.0; 
	
	// Suceptibility of different vaccination classes
	if (_param_infer("rU") == true && h_residents[idx].i_vacc_status == 1) {
		factor *= _param_value("rU");
	}

	if (_param_infer("r1V") == true && h_residents[idx].i_vacc_status == 2 && h_residents[idx].i_vacc_delay < VACC_ONSET) {
		factor *= _param_value("r1V");
	}

	if (_param_infer("r1V90") == true && h_residents[idx].i_vacc_status == 2 && h_residents[idx].i_vacc_delay >= VACC_ONSET) {
		factor *= _param_value("r1V90");
	}

	if (_param_infer("r2V") == true && h_residents[idx].i_vacc_status == 3 && h_residents[idx].i_vacc_delay < VACC_ONSET) {
		factor *= _param_value("r2V");
	}

	if (_param_infer("r3V") == true && h_residents[idx].i_vacc_status == 4) {
		factor *= _param_value("r3V");
	}

	// Suceptibility of children
	if (_param_infer("rC") == true && h_residents[idx].i_age == 0) {
		factor *= _param_value("rC");
	}

	// Suceptibility of isolated people
	if (_param_infer("rQ") == true && h_residents[idx].i_isolation == true) {
		factor *= _param_value("rQ");
	}

	// Suceptibility of recovered people
	if (_param_infer("rR") == true && h_residents[idx].i_recovered == true && (h_residents[idx].i_vacc_status == 1 || h_residents[idx].i_vacc_status == 2)) {
		factor *= _param_value("rR");
	}

	return factor;
}

double Household::model_transmissability (){
	double factor = 1.0; 

	// Instantaneous per capita transmission rate 
	factor *= _param_value("beta");

	// If transmission hazard depends on household size
	if (_param_infer("delta") == true) {
		factor /= pow( h_size / _param_value("typical_size"), _param_value("delta") );
	} 

	return factor;
}
// FIXME : to generalize with Household::log_infec_proba, log_surv_proba

/* logarithm of the infection probability */
double Household::log_infec_proba (size_t idx) {
	// Instantaneous risk of infection from other household members
	double beta_i = 0.0, beta = 0.0, output = 0.0; 
	double t1 = h_residents[idx].i_augm_infec_time;

	// all residents
	for (size_t i = 0; i < h_size; ++i) {
		if (idx != i) {
			beta_i = h_inst_Lambda[i][idx];
			beta_i *= model_contact_modulation(idx, i);
			beta += beta_i;
			
			// --- Display --- //
			if (Verbose == true) {
				std::cout << "pInf: " << idx << " " << i << " " << beta_i << std::endl;
			}
			// --- End-of-Display --- //
	  }
	}

	beta *= model_transmissability();
	output = model_susceptibility(idx) * ( model_community_risk_inst(t1) + beta );
	output = std::log(1. - std::exp(-output));

	// --- Display --- //
	if (Verbose == true) {
		std::cout << "log_infec_proba: " << output << std::endl;
	}
	// --- End-of-Display --- //

	return output;
}

/* Survival probability from time t0 of resident idx */
double Household::log_surv_proba (size_t idx, double t0) {
	// Survival from t0 to t1	
  	double t1, output = 0.0, beta = 0.0, beta_i = 0.0;

	if (h_residents[idx].i_is_case > 0) { 		// infected case
		t1 = h_residents[idx].i_augm_infec_time;

	} else {									// not a case
		t1 = h_residents[idx].i_study_period;
	}

  	// Cumulative risk of infection from other household members
	for (size_t i = 0; i < h_size; ++i) {
		if (idx != i) {
			beta_i = h_cumul_Lambda[i][idx];
			beta_i *= model_contact_modulation(idx, i);
			beta += beta_i;

			// --- Display --- //
			if (Verbose == true){
				std::cout << "pInf: " << idx << " " << i << " " << beta_i << std::endl;
			}
			// --- End-of-Display --- //
		}
	}

	beta *= model_transmissability();
	output = model_susceptibility(idx) * (model_community_risk_cumul(t0, t1) + beta);
	output *= (-1.) ;

	// --- Display --- //
  	if (Verbose == true){
		std::cout << "log_surv_proba: " << output << std::endl;
	}
	// --- End-of-Display --- //

	return output;
}