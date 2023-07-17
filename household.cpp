#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <algorithm>
#include "household.hpp"
#include "utils.h"

using namespace std;

/////////////////////////////////////////////////////
// Constant parameters
/////////////////////////////////////////////////////

// FIXME : this might not work
const struct const_params{
	double mGamma = 26.1;
	double vGamma = 7.0;
  	double shift = 25.6;
	double shapeInf = pow(mGamma, 2) / vGamma;
  	double scaleInf = vGamma / mGamma;

  	double mIncub_g = 3.9;
  	double sdIncub_g = 2.1;
  	double shape_incub = pow(mIncub_g, 2) / pow(sdIncub_g, 2);
  	double scale_incub = pow(sdIncub_g, 2) / mIncub_g;

  	double shape_gt = 4.0;
  	double scale_gt = 5.0 / 4.0;
	double typical_size = 5.0;
} my_const_params;

/////////////////////////////////////////////////////
//  Individual : Class
/////////////////////////////////////////////////////

// ----------------------------
//  Constructors
// ----------------------------

Individual::Individual (
	int indid,
	double onset_time,
	bool is_case,
	int age,
	int sex,
	int start_follow_up, 
	int study_period,
	int vacc_status,
	int isolation,
	int vacc_delay,
	int recovered
	) {
	// Data Variables
	indid = indid;
	onset_time = onset_time;
	is_case = is_case;
	age = age;
	sex = sex;
	start_follow_up = start_follow_up;
 	study_period = study_period;
	vacc_status = vacc_status;
	isolation = isolation;
	vacc_delay = vacc_delay;
	recovered = recovered;
	// Augmented Variables
	augm_onset_time = onset_time;
}

// ----------------------------
//  Methods : Computations
// ----------------------------

double Individual::infec_profile (double tinf, int study_period) {
	double out = 0.;

	// Symptomatic infector
	if (is_case > 0) {
		if (infec_time < study_period && infec_time < tinf)
		{
			out = dtost(tinf-augm_onset_time, augm_onset_time-infec_time); 
		} 
	}

	return out;
}

double Individual::infec_cumulative (double tinf, int study_period) {
	double out = 0.;

	// Symptomatic infector
	if (is_case > 0) {
		if (infec_time < study_period && infec_time < tinf) {
			out = ptost(tinf-augm_onset_time, augm_onset_time-infec_time); 
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
	h_cum_community_infec = cumval;
	h_contact_pattern = contact_pattern; 
}

// ----------------------------
//  Methods : Getter
// ----------------------------

size_t Household::get_size (void) const {
	return h_size;
}

size_t Household::num_infec_cases (void) const {
	return h_confCase.size();
}

std::vector<int> Household::get_infec_list (void) const {
	return h_confCase;
}

// number of symptomatic // FIXME: should be rather a function of <is_case> for clarity?
size_t Household::num_symptomatic (void) const{
	size_t num_symptomatic = 0;

	for (size_t i = 0; i < h_size; i++) {
		if (h_residents[i].onset_time < INFINITY) {
			num_symptomatic += 1;
		}
	}

	return num_symptomatic;
}

// get individual index
std::string Household::get_resident_index (void) const {
	std::string indid;

	for (size_t i = 0; i < h_size; i++) {
		indid += std::to_string(h_residents[i].indid) + " "; 
	}

	return indid;
}

// get augmented onset time for given resident
double Household::get_augm_onset_time(int index) const{
	return h_residents[index].augm_onset_time;
}

// ----------------------------
//  Methods : Setter
// ----------------------------

// Add individual to household
void Household::add_individual (Individual this_indiv) {
	
	h_residents.push_back(this_indiv);
	h_size += 1; 
	h_infec_time = std::min(h_infec_time, this_indiv.infec_time);
	//h_start_follow_up = this_indiv.start_follow_up; // FIXME : what is this?

	// Add confirmed cases
	if (this_indiv.is_case == 0) { 
		/* 
		1: symptomatic cases
		2: asymptomatic cases
		3: symptomatic cases with unknown symptom onset
		*/
		// FIXME : is this still the case?
		h_confCase.push_back(h_size-1);
	} else {
		h_notis_case += 1;
	}
}

// Change symptom onset // FIXME : augmented
void Household::set_onset_time (int index, double onset_time) {
	h_residents[index].augm_onset_time = onset_time;
}

// Change infection time
void Household::set_infec_time (int index, double infec_time) {
	h_residents[index].infec_time = infec_time;
	h_infec_time = std::min(h_infec_time, infec_time);
}

// Initialize augmented onset time
void Household::init_onset_time (int index, std::mt19937_64& rnd_generator, bool display) {

	if (h_residents[index].is_case > 0) { 
		h_residents[index].augm_onset_time = h_residents[index].onset_time + runif(rnd_generator, 0.0, 1.0);
	} else {
		// Household residents without infection or with unknown final outcome
		h_residents[index].augm_onset_time = INFINITY; 
	}

	// --- Display --- //	
	if (display == true) {
		cout << "Observed onset time: " << h_residents[index].onset_time << endl;
		cout << "Augmented onset time" << h_residents[index].augm_onset_time << endl;
  	} 
	// --- End-of-Display --- //
}

// Initialize augmented infection time
void Household::init_inf_time(int index, std::mt19937_64& rnd_generator, bool display) {
	double incubPeriod = 0.0;

	if (h_residents[index].is_case == 1) {   
		// Symptomatic case 
		incubPeriod = rgamma(rnd_generator, my_const_params.shape_incub, my_const_params.scale_incub);
		h_residents[index].infec_time = h_residents[index].augm_onset_time - incubPeriod;
	} else {    
		// Household residents without infection or with unknown final outcome
		h_residents[index].infec_time = INFINITY; 
	}

	// --- Display --- //	
	if (display == true) {
		cout << "Inf time " << h_residents[index].indid << ": " << h_residents[index].infec_time << endl;
	}
	// --- End-of-Display --- //
}

// ----------------------------
//  Methods : Computations
// ----------------------------

// Compute the risk of infection
void Household::compute_infec_profiles(bool display) {
	double t1;

	// initialize 2dim vectors h_cum_Lambda and h_inst_Lambda
	h_cum_Lambda.resize(h_size);
	h_inst_Lambda.resize(h_size);
	for (size_t i=0; i<h_size; i++) {
		h_cum_Lambda[i].resize(h_size);
		h_inst_Lambda[i].resize(h_size);
	}

	// Compute cumulative and instantaneous person-to-person transmission rates
	for (size_t infector = 0; infector < h_size; infector++) {
		for (size_t infectee = 0; infectee < h_size; infectee++) {

			if ( h_residents[infectee].is_case == 0 ) {
				t1 = h_residents[infectee].study_period;
	  		} else {
				t1 = h_residents[infectee].infec_time;
			}

			if (display == true) {
				cout << "Infector: " << infector << " " << h_residents[infector].infec_time  << " " << h_residents[infector].augm_onset_time << endl;
				cout << "Infectee: " << infectee << " " << t1 << " " << h_residents[infectee].is_case << endl;
			}

			// Cumulative transmission rate for household contacts
			if (h_residents[infector].is_case > 0 && h_residents[infector].infec_time < t1 && h_residents[infectee].is_case >=0) {
				h_cum_Lambda[infector][infectee] += h_residents[infector].infec_cumulative(t1, h_residents[infectee].study_period);
			}

			// Instantaneous transmission rate for secondary cases
			if (h_residents[infector].is_case > 0 && h_residents[infector].infec_time < h_residents[infectee].infec_time && h_residents[infectee].is_case > 0) {
				h_inst_Lambda[infector][infectee] += h_residents[infector].infec_profile(h_residents[infectee].infec_time, h_residents[infectee].study_period);
			}
		}
	}

	// --- Display --- //	
	if (display == true) {
		// Cumulative transmission rate
		for (size_t i = 0; i < h_size;i++) {
			for (size_t j = 0; j < h_size; j++) {
				cout << h_cum_Lambda[i][j] << " "; 
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
			cout << h_residents[i].infec_time << " ";
		}
		cout << endl;
	}
	// --- End-of-Display --- //	
}

void Household::update_infec_profiles(int ind, bool display) {
	// Update h_cum_Lambda and h_inst_Lambda matrices during data augmentation
	double t1;

	// Ind is the infector
	for (size_t infectee = 0; infectee < h_size; infectee++) {
		// Reinitialize
		h_cum_Lambda[ind][infectee] = 0;
		h_inst_Lambda[ind][infectee] = 0;

		// Cumulative transmission rate for all household contacts
		if ( h_residents[infectee].is_case == 0 ) {
			t1 = h_residents[infectee].study_period;
		}else{
			t1 = h_residents[infectee].infec_time;
		}

		if (h_residents[ind].is_case > 0 && h_residents[ind].infec_time < t1 && h_residents[infectee].is_case >=0) {

			h_cum_Lambda[ind][infectee] = h_residents[ind].infec_cumulative(t1, h_residents[infectee].study_period);
		}

		// Instantaneous transmission rate for secondary cases
		if (h_residents[ind].is_case > 0 && h_residents[ind].infec_time < h_residents[infectee].infec_time && h_residents[infectee].is_case > 0) {
			h_inst_Lambda[ind][infectee] = h_residents[ind].infec_profile(t1, h_residents[infectee].study_period);
		}
	}

	// Ind is an infectee
	if (h_residents[ind].is_case == 0) {
		t1 = h_residents[ind].study_period;
	} else {
		t1 = h_residents[ind].infec_time;    
	}

	for (size_t infector = 0; infector < h_size; infector++) {    
		// Reinitialize
		h_cum_Lambda[infector][ind] = 0;
		h_inst_Lambda[infector][ind] = 0;

		// Cumulative transmission rate for all household contacts
		if (h_residents[infector].is_case > 0 && h_residents[infector].infec_time < t1 && h_residents[ind].is_case >=0) {
			h_cum_Lambda[infector][ind] = h_residents[infector].infec_cumulative(t1, h_residents[ind].study_period);
		}
		// Instantaneous transmission rate for secondary cases
		if (h_residents[infector].is_case > 0 && h_residents[infector].infec_time < h_residents[ind].infec_time && h_residents[ind].is_case > 0) {
			h_inst_Lambda[infector][ind] = h_residents[infector].infec_profile(h_residents[ind].infec_time, h_residents[ind].study_period);
		}
	}

	// --- Display --- //	
	if (display==true) {
		for (size_t i = 0; i < h_size;i++) {
			for (size_t j = 0; j < h_size; j++) {
				cout << h_cum_Lambda[i][j] << " "; 
			}
			cout << endl;
		}
		for (size_t i = 0; i < h_size; i++) {
			cout << h_residents[i].infec_time << " ";
		}
		cout << endl;
	}
	// --- End-of-Display --- //	
}


double Household::new_onset_time (int index, std::mt19937_64& rnd_generator) {
	// Propose new value for symptom onset time
	double new_onset_time = 0.0;		

	if (h_residents[index].is_case > 0) {           
		// is_case household member 
		new_onset_time = h_residents[index].onset_time + runif(rnd_generator, 0.0, 1.0);
	} else {    
		// Household residents without infection or with unknown final outcome
		new_onset_time = INFINITY;
	}

	return new_onset_time;
}

double Household::new_infec_time (int index, std::mt19937_64& rnd_generator) {
	// Propose new value for infection time
	double new_infec_time = 0.0, incubPeriod = 0.0;

	if (h_residents[index].is_case == 1) {           
		// Symptomatic case 
		incubPeriod = rgamma(rnd_generator, my_const_params.shape_incub, my_const_params.scale_incub);
		new_infec_time += h_residents[index].augm_onset_time - incubPeriod;
	} else {    
		// Household residents without infection or with unknown final outcome
		new_infec_time = INFINITY;
	}

	return new_infec_time;
}
// FIXME : change h_infec_time in the Household accordingly!

double Household::log_dIncub (int index, bool display) {
	// FIXME : comment
	double out = 0.0, incubPeriod = 0.0;

	if (h_residents[index].is_case == 1) {       
		// Symptomatic case with known symptom onset 
		incubPeriod = h_residents[index].augm_onset_time - h_residents[index].infec_time;
		out = logdlnorm(incubPeriod, mIncub, sdIncub);
	} 

	// --- Display --- //
	if (display == true) {
		cout << "log_dIncub: " << out << endl;
	} // --- End-of-Display --- //	

	return out;
}

double Household::log_dIncub_g (int index, bool display) {
	// FIXME : comment
	double out = 0.0, incubPeriod = 0.0;

	if (h_residents[index].is_case == 1) {       
		// Symptomatic case with known symptom onset 
		incubPeriod = h_residents[index].augm_onset_time - h_residents[index].infec_time;
		out = log(dgamma(incubPeriod, my_const_params.shape_incub, my_const_params.scale_incub));
	} 

	// --- Display --- //
	if (display == true) {
		cout << "log_dIncub: " << out << endl;
	}
	// --- End-of-Display --- //	

	return out;
}


double Household::compute_log_lik (std::vector<double> parameter, std::vector<int> selectedParam, double maxPCRDetectability, bool display) {
	// FIXME : comment
	size_t firstCase(0);
	double LL = 0.0; // Log Likelihood
	double t0 = h_infec_time;

	// Identify index case
	for (size_t i = 0; i < h_size; i++) {
		if (h_residents[i].infec_time == t0) {
			firstCase = i;
		}
	}

	// All residents
	for (size_t i = 0; i < h_size; i++) { 

		// --- Display --- //
		if (display == true) {
			cout << i << " " << firstCase << " " << h_residents[i].augm_onset_time;
			cout << " " << h_residents[i].infec_time << " " << h_residents[i].is_case << endl; 
		}
		// --- End-of-Display --- //
		
		if (i == firstCase) {     												// Incubation period of 1st is_case                                      		
			LL += log_dIncub_g(i, display); 
		} else if (h_residents[i].is_case > 0) {             					// Contribution of the secondary cases
	  		LL += log_S(i, t0, parameter, selectedParam, display);
			LL += log_pInf(i, parameter, selectedParam, display);
			LL += log_dIncub_g(i, display);
		} else {																// Non is_case individuals tha contribute to the likelihood
	  		if (h_residents[i].is_case == 0) {
				// FIXME : isn't this always true at this point of the logical case
				LL += log_S(i, t0, parameter, selectedParam, display); 
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

double Household::log_pInf (size_t curr, std::vector<double> parameter, std::vector<int> selectedParam, bool display) {
	// FIXME : 
	int tinf = (int) h_residents[curr].infec_time; // Probability of infection at t1
	double alpha = parameter[0] * h_inst_community_infec[tinf]; // Instaneous risk of infection from community
	double beta_i = 0.0, beta = 0.0; // Instantaneous risk of infection from other household members

	for (size_t i = 0; i < h_size; ++i) {
		if (curr != i) {
			beta_i = h_inst_Lambda[i][curr];

			// --- Display --- //
			if (display == true) {
				cout << "pInf: " << curr << " " << i << " " << beta_i << endl;
			}
			// --- End-of-Display --- //

			// Relative infectivity of symptomatic child cases versus symptomatic adult cases
			// FIXME : this is saying daughter-mother, father - son 
			// FIXME : there could be else if instead
			if (h_contact_pattern == true) {
				// Child (2) - Mother (0)
				if (((h_residents[curr].age == 0 && h_residents[i].age == 1 && h_residents[i].sex == 1) || (h_residents[curr].age == 1 && h_residents[i].age == 0 && h_residents[curr].sex == 1)) && selectedParam[11]) {
					beta_i *= parameter[11];
				}
				// Child (2) - Father (1)
				if (((h_residents[curr].age == 0 && h_residents[i].age == 1 && h_residents[i].sex == 2) || (h_residents[curr].age == 1 && h_residents[i].age == 0 && h_residents[curr].sex == 2))&& selectedParam[12] ) {
					beta_i *= parameter[12];
				}
				// Child (2) - Child (2)
				if ( h_residents[curr].age == 0 && h_residents[i].age == 0 && selectedParam[10]) {
					beta_i *= parameter[10];
				} 
	   		}
			beta += beta_i;
	  }
	}

	// Instantaneous per capita transmission rate 
	beta *= parameter[1];

	// If transmission hazard depends on household size
	if (selectedParam[2] == 1) {
		beta /= pow( h_size / my_const_params.typical_size, parameter[2] );
	} 
	beta = beta+alpha;

	// Suceptibility of children
	if (h_residents[curr].age==0 && selectedParam[8]) {
		beta *= parameter[8];
	}
	// suceptibility of isolated people
	if (h_residents[curr].isolation && selectedParam[9]) {
		beta *= parameter[9];
	}
	// Suceptibility of recovered people
	if (h_residents[curr].recovered && selectedParam[13] && (h_residents[curr].vacc_status==1 || h_residents[curr].vacc_status==2)) {
		beta *= parameter[13];
	}
	
	// Suceptibility of different vaccination classes
	if (h_residents[curr].vacc_status==1 && selectedParam[3]) {
		beta *= parameter[3];
	}
	if (h_residents[curr].vacc_status==2 && h_residents[curr].vacc_delay<90 && selectedParam[4]) {
		beta *= parameter[4];
	}
	if (h_residents[curr].vacc_status==2 && h_residents[curr].vacc_delay>=90 && selectedParam[5]) {
		beta *= parameter[5];
	}
	if (h_residents[curr].vacc_status==3 && h_residents[curr].vacc_delay<90 && selectedParam[6]) {
		beta *= parameter[6];
	}
	if (h_residents[curr].vacc_status==4 && selectedParam[7]) {
		beta *= parameter[7];
	}

	// --- Display --- //
	if (display == true) {
		cout << "log_pInf: " << log((alpha+beta)) << endl;
	}
	// --- End-of-Display --- //

	return log((beta));
}

double Household::log_S (size_t curr, double t0, std::vector<double> parameter, std::vector<int> selectedParam, bool display) {
	// Survival from t0 to t1	
  	double t1, alpha = 0.0, beta = 0.0, beta_i = 0.0;
  	int tb = (int) t0;
  	int te = (int) t1;

	// If the individual is not is_case, it's infection
	if ( h_residents[curr].is_case == 0 ) {
		t1 = h_residents[curr].study_period;
	}else{
		t1 = h_residents[curr].infec_time;
	}
  
	// Cumulative risk of infection from community
	alpha = parameter[0] * (h_cum_community_infec[te] - h_cum_community_infec[tb]);

  	// Cumulative risk of infection from other household members
	for (size_t i = 0; i < h_size; ++i) {
		if (curr != i) {
			beta_i = h_cum_Lambda[i][curr];

			// --- Display --- //
			if (display == true){
				cout << "pInf: " << curr << " " << i << " " << beta_i << endl;
			}
			// --- End-of-Display --- //

			// Relative infectivity of symptomatic child cases versus symptomatic adult cases
			// FIXME : this is saying daughter-mother, father - son 
			// FIXME : there should be <if child> and <else if>
			if (h_contact_pattern == true) {
				// Child (2) - Mother (0)
				if (((h_residents[curr].age == 0 && h_residents[i].age == 1 && h_residents[i].sex == 1) || (h_residents[curr].age == 1 && h_residents[i].age == 0 && h_residents[curr].sex == 1)) && selectedParam[11]) {
					beta_i *= parameter[11];
				}
				// Child (2) - Father (1)
				if (((h_residents[curr].age == 0 && h_residents[i].age == 1 && h_residents[i].sex == 2) || (h_residents[curr].age == 1 && h_residents[i].age == 0 && h_residents[curr].sex == 2))&& selectedParam[12] ) {
					beta_i *= parameter[12];
				}
				// Child (2) - Child (2)
				if ( h_residents[curr].age == 0 && h_residents[i].age == 0 && selectedParam[10]) {
					beta_i *= parameter[10];
				} 
			}
			beta += beta_i;
		}
	}

	// FIXME : this is a copy of Household::log_pInf !
	// FIXME : to be generalized

	// Instantaneous per capita transmission rate 
	beta *= parameter[1];

	// If transmission hazard depends on household size
	if (selectedParam[2] == 1) {
		beta /= pow( h_size / my_const_params.typical_size, parameter[2] );
	} 

	beta = beta+alpha;

	// Suceptibility of children
	if (h_residents[curr].age==0 && selectedParam[8]) {
		beta *= parameter[8];
	}
	// suceptibility of isolated people
	if (h_residents[curr].isolation && selectedParam[9]) {
		beta *= parameter[9];
	}
	// Suceptibility of recovered people
	if (h_residents[curr].recovered && selectedParam[13] && (h_residents[curr].vacc_status==1 || h_residents[curr].vacc_status==2)) {
		beta *= parameter[13];
	}
	
	// Suceptibility of different vaccination classes
	if (h_residents[curr].vacc_status==1 && selectedParam[3]) {
		beta *= parameter[3];
	}
	if (h_residents[curr].vacc_status==2 && h_residents[curr].vacc_delay<90 && selectedParam[4]) {
		beta *= parameter[4];
	}
	if (h_residents[curr].vacc_status==2 && h_residents[curr].vacc_delay>=90 && selectedParam[5]) {
		beta *= parameter[5];
	}
	if (h_residents[curr].vacc_status==3 && h_residents[curr].vacc_delay<90 && selectedParam[6]) {
		beta *= parameter[6];
	}
	if (h_residents[curr].vacc_status==4 && selectedParam[7]) {
		beta *= parameter[7];
	}

	// --- Display --- //
  	if (display == true){
		cout << beta << " " << alpha << endl;
		cout << "log_S: " << -(alpha+beta) << endl;
	}
	// --- End-of-Display --- //

	return - beta;
}