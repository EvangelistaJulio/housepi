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
	double start_follow_up,
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
	i_start_follow_up = start_follow_up;
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

double Individual::infec_profile (double t_infec, double study_period) {
	double out = 0.0, tost, t_incub;
	std::vector<double> args = {
		_param_value("tost_tau"), _param_value("tost_mu"),
		_param_value("tost_sigma"), _param_value("tost_alpha")};

	// Symptomatic infector
	if (i_is_case > 0) {
		if (i_augm_infec_time < study_period && i_augm_infec_time < t_infec) {
			tost = t_infec - i_augm_onset_time;
			t_incub = i_augm_onset_time - i_augm_infec_time;
			out = t_onset_sympt_to_transm_pdf(tost, t_incub, args); 
		} 
	} 

	return out;
}

double Individual::infec_cumulative (double t_infec, double study_period) {
	double out = 0.0, tost, t_incub;
	std::vector<double> args = {
		_param_value("tost_tau"), _param_value("tost_mu"),
		_param_value("tost_sigma"), _param_value("tost_alpha")};

	// Symptomatic infector
	if (i_is_case > 0) {
		if (i_augm_infec_time < study_period && i_augm_infec_time < t_infec) {
			tost = t_infec - i_augm_onset_time;
			t_incub = i_augm_onset_time - i_augm_infec_time;
			out = t_onset_sympt_to_transm_cdf(tost, t_incub, args); 
		} 
	}

	return out;
}

/////////////////////////////////////////////////////
//  Household : Class
/////////////////////////////////////////////////////

// ----------------------------
//  Constructors
// ----------------------------

// FIXME : check if default constructor is working

Household::Household (bool contact_pattern, std::vector<double> instval, std::vector<double> cumval) {
	h_inst_community_infec = instval;
	h_cumul_community_infec = cumval;
	h_contact_pattern = contact_pattern; 
}

// ----------------------------
//  Methods : Getter
// ----------------------------

size_t Household::get_size () const { return h_residents.size(); }
size_t Household::num_infec_cases () const { return h_confCase.size(); }
std::vector<int> Household::get_infec_list () const { return h_confCase; }

// number of symptomatic
size_t Household::num_symptomatic () const{
	size_t num_symptomatic = 0;

	for (size_t i = 0; i < h_residents.size(); i++) {
		if (h_residents[i].i_onset_time < INFINITY) {
			num_symptomatic += 1;
		}
	}

	return num_symptomatic;
}

// get individual idx
std::string Household::get_resident_index (void) const {
	std::string indid;

	for (size_t i = 0; i < h_residents.size(); i++) {
		indid += std::to_string(h_residents[i].i_indid) + " "; 
	}

	return indid;
}

// get augmented onset time for given resident
double Household::get_augm_onset_time(size_t idx) const{
	return h_residents[idx].i_augm_onset_time;
}

// ----------------------------
//  Methods : Setter
// ----------------------------

// Add individual to household
void Household::add_individual (Individual this_indiv) {
	
	h_residents.push_back(this_indiv);
	h_size += 1; 
	h_infec_time = std::min(h_infec_time, this_indiv.i_augm_infec_time);
	//h_start_follow_up = this_indiv.start_follow_up; // FIXME : what is this?

	// Add confirmed cases
	if (this_indiv.i_is_case == 0) { 
		// FIXME : is this still the case?
		h_confCase.push_back(h_size - 1);
	} 
}


// Initialize augmented onset time
void Household::init_onset_time (size_t idx, std::mt19937_64& rand_gen, bool display) {

	if (h_residents[idx].i_is_case > 0) { 
		h_residents[idx].i_augm_onset_time = h_residents[idx].i_onset_time + unif_rand(rand_gen);
	} else {
		// Household residents without infection or with unknown final outcome
		h_residents[idx].i_augm_onset_time = INFINITY; 
	}

	// --- Display --- //	
	if (display == true) {
		cout << "Observed onset time: " << h_residents[idx].i_onset_time << endl;
		cout << "Augmented onset time" << h_residents[idx].i_augm_onset_time << endl;
  	} 
	// --- End-of-Display --- //
}

// Initialize augmented infection time
void Household::init_infec_time(size_t idx, std::mt19937_64& rand_gen, bool display) {
	double incubPeriod = 0.0;
	std::vector<double> args = {_param_value("shape_incub"), _param_value("scale_incub")};

 	// FIXME : this should be part off the model class
	if (h_residents[idx].i_is_case == 1) {   
		// Symptomatic case 
		incubPeriod = gamma_rand(rand_gen, args);
		h_residents[idx].i_augm_infec_time = h_residents[idx].i_augm_onset_time - incubPeriod;
	} else {    
		// Household residents without infection or with unknown final outcome
		h_residents[idx].i_augm_infec_time = INFINITY; 
	}

	// --- Display --- //	
	if (display == true) {
		cout << "Inf time " << h_residents[idx].i_indid << ": " << h_residents[idx].i_augm_infec_time << endl;
	}
	// --- End-of-Display --- //
}

// ----------------------------
//  Methods : Computations
// ----------------------------

// Change symptom onset // FIXME : augmented
void Household::set_augm_onset_time (size_t idx, double onset_time, bool display) {
	h_residents[idx].i_augm_onset_time = onset_time;
	update_infec_profiles(idx, display);
}

// Change infection time
void Household::set_augm_infec_time (size_t idx, double infec_time, bool display) {
	h_residents[idx].i_augm_infec_time = infec_time;
	h_infec_time = std::min(h_infec_time, infec_time);
	update_infec_profiles(idx, display);
}

// Compute the risk of infection
void Household::init_infec_profiles(bool display) {

	// initialize 2-dim vectors h_cumul_Lambda and h_inst_Lambda
	h_cumul_Lambda.resize(h_size);
	h_inst_Lambda.resize(h_size);
	for (size_t i = 0; i  <h_size; i++) {
		h_cumul_Lambda[i].resize(h_size);
		h_inst_Lambda[i].resize(h_size);
	}

	// Compute cumulative and instantaneous person-to-person transmission rates
	for (size_t infector = 0; infector < h_size; infector++) {
		// FIXME : is this working ?
		update_infec_profiles(infector, false); // 
	}

	// --- Display --- //	
	if (display == true) {
		// Cumulative transmission rate
		for (size_t i = 0; i < h_size;i++) {
			for (size_t j = 0; j < h_size; j++) {
				cout << h_cumul_Lambda[i][j] << " "; 
			}
			cout << endl;
		}
		cout << endl;
		// Instantaneous transmission rate
		for (size_t i = 0; i < h_size;i++) {
			for (size_t j = 0; j < h_size; j++) {
				cout << h_inst_Lambda[i][j] << " "; 
			}
			cout << endl;
		}
		cout << endl;
		// Infection time 
		for (size_t i = 0; i < h_size; i++) {
			cout << h_residents[i].i_augm_infec_time << " ";
		}
		cout << endl;
	}
	// --- End-of-Display --- //	
}

// Update h_cumul_Lambda and h_inst_Lambda matrices during data augmentation
void Household::update_infec_profiles(size_t idx, bool display) {
	double t1;

	// idx is the infector
	for (size_t infectee = 0; infectee < h_size; infectee++) {
		// FIXME : idx != infectee

		if (h_residents[idx].i_is_case > 0){
			// Cumulative transmission rate for all household contacts
			t1 = h_residents[infectee].i_augm_infec_time;
			if (h_residents[idx].i_augm_infec_time < t1 && h_residents[infectee].i_is_case >= 0) { // FIXME : is_case >= 0 always
				h_cumul_Lambda[idx][infectee] = h_residents[idx].infec_cumulative(t1, h_residents[infectee].i_study_period);
			}

			// Instantaneous transmission rate for secondary cases
			if (h_residents[idx].i_augm_infec_time < h_residents[infectee].i_augm_infec_time && h_residents[infectee].i_is_case > 0) {
				h_inst_Lambda[idx][infectee] = h_residents[idx].infec_profile(t1, h_residents[infectee].i_study_period);
			} else {
				t1 = h_residents[infectee].i_study_period;
			}
		}
	}

	// idx is an infectee
	if (h_residents[idx].i_is_case == 0) {
		t1 = h_residents[idx].i_study_period;
	} else {
		t1 = h_residents[idx].i_augm_infec_time;    
	}

	for (size_t infector = 0; infector < h_size; infector++) {    
		// FIXME : idx != infector

		if (h_residents[infector].i_is_case > 0){
			// Cumulative transmission rate for all household contacts
			if (h_residents[infector].i_augm_infec_time < t1 && h_residents[idx].i_is_case >= 0) { // FIXME : is_case >= 0 always
				h_cumul_Lambda[infector][idx] = h_residents[infector].infec_cumulative(t1, h_residents[idx].i_study_period);
			}
			// Instantaneous transmission rate for secondary cases
			if (h_residents[infector].i_augm_infec_time < h_residents[idx].i_augm_infec_time && h_residents[idx].i_is_case > 0) {
				h_inst_Lambda[infector][idx] = h_residents[infector].infec_profile(h_residents[idx].i_augm_infec_time, h_residents[idx].i_study_period);
			}
		}
	}

	// --- Display --- //	
	if (display == true) {
		for (size_t i = 0; i < h_size; i++) {
			for (size_t j = 0; j < h_size; j++) {
				cout << h_cumul_Lambda[i][j] << " "; 
			}
			cout << endl;
		}
		for (size_t i = 0; i < h_size; i++) {
			cout << h_residents[i].i_augm_infec_time << " ";
		}
		cout << endl;
	}
	// --- End-of-Display --- //	
}


double Household::new_onset_time (size_t idx, std::mt19937_64& rand_gen) {
	// Propose new value for symptom onset time
	double new_onset_time = 0.0;		

	if (h_residents[idx].i_is_case > 0) {           
		// is_case household member 
		new_onset_time = h_residents[idx].i_onset_time + unif_rand(rand_gen);
	} else {    
		// Household residents without infection or with unknown final outcome
		new_onset_time = INFINITY;
	}

	return new_onset_time;
}

double Household::new_infec_time (size_t idx, std::mt19937_64& rand_gen) {
	// Propose new value for infection time
	double new_infec_time = 0.0, incubPeriod = 0.0;
	std::vector<double> args = {_param_value("shape_incub"), _param_value("scale_incub")};

	if (h_residents[idx].i_is_case == 1) {           
		// Symptomatic case 
		incubPeriod = gamma_rand(rand_gen, args);
		new_infec_time += h_residents[idx].i_augm_onset_time - incubPeriod;
	} else {    
		// Household residents without infection or with unknown final outcome
		new_infec_time = INFINITY;
	}

	return new_infec_time;
}
// FIXME : change h_infec_time in the Household accordingly!


double Household::compute_loglike (bool display) {
	// FIXME : comment
	size_t firstCase(0);
	double LL = 0.0; // Log Likelihood
	double t0 = h_infec_time;

	// Identify idx case
	for (size_t i = 0; i < h_size; i++) {
		if (h_residents[i].i_augm_infec_time == t0) {
			firstCase = i;
		}
	}

	// All residents
	for (size_t i = 0; i < h_size; i++) { 

		// --- Display --- //
		if (display == true) {
			cout << i << " " << firstCase << " " << h_residents[i].i_augm_onset_time;
			cout << " " << h_residents[i].i_augm_infec_time << " " << h_residents[i].i_is_case << endl; 
		}
		// --- End-of-Display --- //
		
		if (i == firstCase) {     													// Incubation period of 1st is_case                                      		
			LL += log_incub_proba(i, display); 
		} else if (h_residents[i].i_is_case > 0) {             					// Contribution of the secondary cases
			LL += log_incub_proba(i, display);
	  		LL += log_surv_proba(i, t0, display);
			LL += log_infec_proba(i, display);
		} else {																// Non is_case individuals tha contribute to the likelihood
	  		if (h_residents[i].i_is_case == 0) {
				// FIXME : isn't this always true at this point of the logical case
				LL += log_surv_proba(i, t0, display); 
			} 
		}
	}

	// --- Display --- //
	if (display == true) {
		cout << "LL: " << LL << "\n" << endl;
	}
	// --- End-of-Display --- //

  return LL;
}



double Household::log_incub_proba (size_t idx, bool display) {
	// FIXME : comment
	double out = 0.0, incubPeriod = 0.0;
	std::vector<double> args = {_param_value("shape_incub"), _param_value("scale_incub")};

	if (h_residents[idx].i_is_case == 1) {       
		// Symptomatic case with known symptom onset 
		incubPeriod = h_residents[idx].i_augm_onset_time - h_residents[idx].i_augm_infec_time;
		out = log_gamma_pdf(incubPeriod, args);
	} 

	// --- Display --- //
	if (display == true) {
		cout << "log_dIncub: " << out << endl;
	}
	// --- End-of-Display --- //	

	return out;
}

double Household::child_infectivity (size_t curr, size_t infectee){
	// Relative infectivity of symptomatic child cases versus symptomatic adult cases
	double beta_i = 1.0;
	size_t i = infectee;

	// Child (2) - Child (2)
	if (_param_infer("rCC") == true && h_residents[curr].i_age == 0 && h_residents[i].i_age == 0) {
		beta_i *= _param_value("rCC");

	// Child (2) - Mother (0)
	} else if (_param_infer("rMC") == true && ((h_residents[curr].i_age == 0 && h_residents[i].i_age == 1 && h_residents[i].i_sex == 1) || (h_residents[i].i_age == 0 && h_residents[curr].i_age == 1 && h_residents[curr].i_sex == 1))) {
		beta_i *= _param_value("rMC");

	// Child (2) - Father (1)
	} else if (_param_infer("rFC") == true && ((h_residents[curr].i_age == 0 && h_residents[i].i_age == 1 && h_residents[i].i_sex == 2) || (h_residents[i].i_age == 0 && h_residents[curr].i_age == 1 && h_residents[curr].i_sex == 2))) {
		beta_i *= _param_value("rFC");
	}

	return beta_i;
}


double Household::model_susceptibility (size_t curr){
	// FIXME : if you don't want to infer a parameter it just should't be in the model... maybe it's not clean
	double out = 1.0; 

	// FIXME : how to loop this?
	
	// Suceptibility of different vaccination classes
	if (_param_infer("rU") == true && h_residents[curr].i_vacc_status == 1) {
		out *= _param_value("rU");
	}
	if (_param_infer("r1V") == true && h_residents[curr].i_vacc_status == 2 && h_residents[curr].i_vacc_delay < VACC_ONSET) {
		out *= _param_value("r1V");
	}
	if (_param_infer("r1V90") == true && h_residents[curr].i_vacc_status == 2 && h_residents[curr].i_vacc_delay >= VACC_ONSET) {
		out *= _param_value("r1V90");
	}
	if (_param_infer("r2V") == true && h_residents[curr].i_vacc_status == 3 && h_residents[curr].i_vacc_delay < VACC_ONSET) {
		out *= _param_value("r2V");
	}
	if (_param_infer("r3V") == true && h_residents[curr].i_vacc_status == 4) {
		out *= _param_value("r3V");
	}

	// Suceptibility of children
	if (_param_infer("rC") == true && h_residents[curr].i_age == 0) {
		out *= _param_value("rC");
	}
	// suceptibility of isolated people
	if (_param_infer("rQ") == true && h_residents[curr].i_isolation) {
		out *= _param_value("rQ");
	}
	// Suceptibility of recovered people
	if (_param_infer("rR") == true && h_residents[curr].i_recovered && (h_residents[curr].i_vacc_status == 1 || h_residents[curr].i_vacc_status == 2)) {
		out *= _param_value("rR");
	}

	return out;
}

double Household::log_infec_proba (size_t curr, bool display) {
	// FIXME : 
	// Probability of infection
	int t_infec = (int) h_residents[curr].i_augm_infec_time; 
	// Instaneous risk of infection from community
	double alpha = _param_value("alpha") * h_inst_community_infec[t_infec]; 
	// Instantaneous risk of infection from other household members
	double beta_i = 0.0, beta = 0.0; 

	for (size_t i = 0; i < h_size; ++i) {
		if (curr != i) {
			beta_i = h_inst_Lambda[i][curr];

			// --- Display --- //
			if (display == true) {
				cout << "pInf: " << curr << " " << i << " " << beta_i << endl;
			}
			// --- End-of-Display --- //

			if (h_contact_pattern == true) {
				beta_i *= child_infectivity(curr, i);
			}

			beta += beta_i;
	  }
	}

	// Instantaneous per capita transmission rate 
	beta *= _param_value("beta");

	// If transmission hazard depends on household size
	if (_param_infer("delta") == true) {
		beta /= pow( h_size / _param_value("typical_size"), _param_value("delta") );
	} 
	beta = beta+alpha;

	beta *= model_susceptibility(curr);

	// --- Display --- //
	if (display == true) {
		cout << "log_infec_proba: " << log((alpha + beta)) << endl;
	}
	// --- End-of-Display --- //

	return log((beta));
}

double Household::log_surv_proba (size_t curr, double t0, bool display) {
	// Survival from t0 to t1	
  	double t1, alpha = 0.0, beta = 0.0, beta_i = 0.0;
  	int tb = (int) t0;

	// If the individual is not is_case, it's infection
	if ( h_residents[curr].i_is_case == 0 ) {
		t1 = h_residents[curr].i_study_period;
	}else{
		t1 = h_residents[curr].i_augm_infec_time;
	}
  
	// Cumulative risk of infection from community
	alpha = _param_value("alpha") * (h_cumul_community_infec[(int) t1] - h_cumul_community_infec[tb]);

  	// Cumulative risk of infection from other household members
	for (size_t i = 0; i < h_size; ++i) {
		if (curr != i) {
			beta_i = h_cumul_Lambda[i][curr];

			// --- Display --- //
			if (display == true){
				cout << "pInf: " << curr << " " << i << " " << beta_i << endl;
			}
			// --- End-of-Display --- //

			if (h_contact_pattern == true) {
				beta_i *= child_infectivity(curr, i);
			}
			beta += beta_i;
		}
	}

	// FIXME : this is a copy of Household::log_infec_proba !
	// FIXME : to be generalized

	// Instantaneous per capita transmission rate 
	beta *= _param_value("beta");

	// If transmission hazard depends on household size
	if (_param_infer("delta") == true) {
		beta /= pow( h_size / _param_value("typical_size"), _param_value("delta") );
	} 

	beta = beta+alpha;

	beta *= model_susceptibility(curr);

	// --- Display --- //
  	if (display == true){
		cout << beta << " " << alpha << endl;
		cout << "log_surv_proba: " << -(alpha + beta) << endl;
	}
	// --- End-of-Display --- //

	return - beta;
}