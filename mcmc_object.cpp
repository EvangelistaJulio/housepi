#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "mcmc_object.hpp"
#include "household.hpp"
#include "utils.h"

using namespace std;

/////////////////////////////////////////////////////
//  McmcObject : Class
/////////////////////////////////////////////////////

// ----------------------------
//  Constructors
// ----------------------------

// default 
McmcObject::McmcObject() {
    m_rnd_generator.seed(123);
    m_globalLogLik = 0.0;
    m_iterations = 0;
    m_iter_time_infec = 0;
    m_accepted_move_data = 0;
	m_proposed_move_data = 0;
    m_rate_for_random_walk.resize(0);
    m_num_accepted_move.resize(0);
    m_num_proposed_move.resize(0);
    m_accept.resize(0);
    m_parameter.resize(0);
    m_selected_param.resize(0);
    m_hhLogLik.resize(0);
    std::vector<Household> m_data;
    m_data.resize(0);
    m_nHH = 0;
    m_lnormPrior.resize(0);
    m_sdLNormPrior.resize(0);
    m_max_PCR_detectability = 10.0;
    m_adaptive_steps = 0.0;
}

// user 
McmcObject::McmcObject(
    unsigned seed,
    int num_iterations,
    int num_adaptive_steps,
    std::vector<Household> data,
    std::vector<double> parameter,
    std::vector<int> selected_parameter,
    std::vector<double> rate_for_random_walk, 
    std::vector<int> lnormPrior,
    std::vector<double> sdLNormPrior,
    int n_iter_time_infec,
    double max_PCR_detectability
    ) {

    // Default values
    m_globalLogLik = 0.0;
    m_num_accepted_move.resize(0);
    m_num_proposed_move.resize(0);
    m_accept.resize(parameter.size());
    m_accepted_move_data = 0;
	m_proposed_move_data = 0;

    // Copy input vectors
    m_rnd_generator.seed(seed);
    m_hhLogLik.resize(data.size());
    m_iterations = num_iterations;
    m_adaptive_steps = num_adaptive_steps;
    m_iter_time_infec = n_iter_time_infec;
    m_rate_for_random_walk = rate_for_random_walk;
    m_lnormPrior = lnormPrior;
    m_sdLNormPrior = sdLNormPrior;
    m_max_PCR_detectability = max_PCR_detectability;
    m_parameter = parameter;
    m_selected_param = selected_parameter;
    m_data = data;
    m_nHH = data.size();
}

// ----------------------------
//  Methods : Setter
// ----------------------------

// Initialize parameter values
void McmcObject::init_param_values () {
    
    for (size_t parID = 0; parID < m_parameter.size(); parID++) {
        // FIXME : what is this if?
        if (m_selected_param[parID]) {
            
            switch (parID) {
            case 0: // alpha - Uniform(0,1)
                m_parameter[parID] = runif(m_rnd_generator); 
                //m_parameter[parID] = rexp(m_rnd_generator, 500); 
                break;

            case 1: // beta - Uniform(0,5)
                m_parameter[parID] = runif(m_rnd_generator, 0.0, 5.0); 
                break;

            case 2: // delta - Uniform(-3,3)
                m_parameter[parID] = runif(m_rnd_generator, -3.0, 3.0); 
                break;

            case 3:
                if ( m_lnormPrior[parID] ) { // rUV - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_rnd_generator, 0.0, m_sdLNormPrior[parID]);   
                } else {                // rUV - Uniform(0,5)
                    m_parameter[parID] = runif(m_rnd_generator, 0.0, 5.0);
                } 
                break;

            case 4:
                if ( m_lnormPrior[parID] ) { // rVC - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_rnd_generator, 0.0, m_sdLNormPrior[parID]);   
                } else {                // rVC - Uniform(0,5)
                    m_parameter[parID] = runif(m_rnd_generator, 0.0, 5.0);
                } 
                break;

            case 5:
                if ( m_lnormPrior[parID] ) { // r1VA - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_rnd_generator, 0.0, m_sdLNormPrior[parID]);   
                } else {                // r1VA - Uniform(0,5)
                    m_parameter[parID] = runif(m_rnd_generator, 0.0, 0.1);
                } 
                break;

            case 6:
                if ( m_lnormPrior[parID] ) { // r2VA - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_rnd_generator, 0.0, m_sdLNormPrior[parID]);   
                } else {                // r2VA - Uniform(0,5)
                    m_parameter[parID] = runif(m_rnd_generator, 0.0, 0.1);
                } 
                break;

            case 7:
                if ( m_lnormPrior[parID] ) { // rUA - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_rnd_generator, 0.0, m_sdLNormPrior[parID]);   
                } else {                // rUA - Uniform(0,5)
                    m_parameter[parID] = runif(m_rnd_generator, 0.0, 5.0);
                } 
                break;

            case 8:
                if ( m_lnormPrior[parID] ) { // rBA - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_rnd_generator, 0.0, m_sdLNormPrior[parID]);   
                } else {                // rBA - Uniform(0,5)
                    m_parameter[parID] = runif(m_rnd_generator, 0.0, 5.0);
                } 
                break;

            case 9:
                if ( m_lnormPrior[parID] ) { // rQ - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_rnd_generator, 0.0, m_sdLNormPrior[parID]);   
                } else {                // rQ - Uniform(0,5)
                    m_parameter[parID] = runif(m_rnd_generator, 0.0, 5.0);
                } 
                break;

            case 10:
                if ( m_lnormPrior[parID] ) { // rCC - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_rnd_generator, 0.0, m_sdLNormPrior[parID]);   
                } else {                // rCC - Uniform(0,5)
                    m_parameter[parID] = runif(m_rnd_generator, 0.0, 5.0);
                } 
                break;

            case 11:
                if ( m_lnormPrior[parID] ) { // rMC - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_rnd_generator, 0.0, m_sdLNormPrior[parID]);   
                } else {                // rPC - Uniform(0,5)
                    m_parameter[parID] = runif(m_rnd_generator, 0.0, 5.0);
                } 
                break;
            
            case 12:
                if ( m_lnormPrior[parID] ) { // rFC - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_rnd_generator, 0.0, m_sdLNormPrior[parID]);   
                } else {                // rPC - Uniform(0,5)
                    m_parameter[parID] = runif(m_rnd_generator, 0.0, 5.0);
                } 
                break;

            case 13:
                if ( m_lnormPrior[parID] ) { // rFC - LogNormal(0,sd)
                    m_parameter[parID] = rlnorm(m_rnd_generator, 0.0, m_sdLNormPrior[parID]);   
                } else {                // rPC - Uniform(0,5)
                    m_parameter[parID] = runif(m_rnd_generator, 0.0, 5.0);
                } 
                break;
            
            }
        }
    }
}

// Modification of attributes
void McmcObject::reset_moves () {
    m_num_proposed_move.clear();
    m_num_proposed_move.resize(m_parameter.size());
    m_num_accepted_move.clear();
    m_num_accepted_move.resize(m_parameter.size());
    m_accepted_move_data = 0;
	m_proposed_move_data = 0;
}

// Initialize infection time
void McmcObject::init_augm_data (bool display) {
    
    for (size_t house = 0; house < m_nHH; house++) {
        for (size_t ind = 0; ind < m_data[house].get_size(); ind++) {
            m_data[house].init_onset_time(ind, m_rnd_generator, display);
            m_data[house].init_inf_time(ind, m_rnd_generator, display);
        }
        // Compute person-to-person transmission rate within household
        m_data[house].compute_infec_profiles(display);
    }
}

// Initialize infection time
void McmcObject::init_log_lik (bool display) {
    int nCases = 0;
    
    // --- Display -- //
    if (display == true) {
        cout << m_parameter[1] << endl;
    }
    // --- End-of-Display -- //

    for (size_t house=0; house < m_nHH; house++) {
        m_hhLogLik[house] = m_data[house].compute_log_lik(m_parameter, m_selected_param, m_max_PCR_detectability, display);
        m_globalLogLik += m_hhLogLik[house];
        nCases += m_data[house].num_symptomatic();
    }
}

// ----------------------------
//  Methods : Computations
// ----------------------------

// Update parameter value in the mcmc chain
void McmcObject::update_parameter (int parID, double step, bool display) {

    // Random walk
	double oldValue = m_parameter[parID];
	double newValue = 0.0;
    double randomWalk = 0.0;
    double logRatioProposal = 0.0;
    double logRatioPrior = 0.0;

    if (parID != 2) {
        randomWalk = m_rate_for_random_walk[parID] * rnorm(m_rnd_generator, 0.0, 1.0);
        newValue = oldValue * exp(randomWalk);
        logRatioProposal = randomWalk;
    } else {
        newValue = oldValue + rnorm(m_rnd_generator, 0.0, m_rate_for_random_walk[parID]);
    }

    // Log ratio of priors 
    // FIXME : if log(0) then skip the move without any further computation
    // FIXME : those are clippers (logspaced), it should be a parameter
    switch (parID) {
    	case 0: // alpha - Uniform(0,1)
    		if (newValue > 1) logRatioPrior = log(0); // FIXME : log of 0
            //logRatioPrior = logdexp(newValue, 500) - logdexp(oldValue, 500); 
    		break;

        case 1: // beta - Uniform(0,10)
        	if (newValue > 5) logRatioPrior = log(0); // FIXME : log of 0
            //logRatioPrior = logdexp(newValue, 0.001) - logdexp(oldValue, 0.001); 
        	break;

        case 2: // delta - Uniform(-3,3)
            if (newValue > 3 || newValue < -3) logRatioPrior = log(0); // FIXME : log of 0
            break;

        case 3:
        	if (m_lnormPrior[parID]) {  // rUC - LogNormal(0,sd)
                logRatioPrior = logdlnorm(newValue, 0.0,  m_sdLNormPrior[parID]) - logdlnorm(oldValue, 0.0, m_sdLNormPrior[parID]);
            } else {               // rUC - Uniform(0,5)
                if (newValue > 5 || newValue < 0) logRatioPrior = log(0); // FIXME : log of 0
            }
        	break;

        case 4:
        	if (m_lnormPrior[parID]) {  // rVC - LogNormal(0,sd)
                logRatioPrior = logdlnorm(newValue, 0.0,  m_sdLNormPrior[parID]) - logdlnorm(oldValue, 0.0, m_sdLNormPrior[parID]);
            } else {               // rVC - Uniform(0,5)
                if (newValue > 5 || newValue < 0) logRatioPrior = log(0); // FIXME : log of 0
            }
        	break;

        case 5:
        	if (m_lnormPrior[parID]) {  // r1VA - LogNormal(0,sd)
                logRatioPrior = logdlnorm(newValue, 0.0,  m_sdLNormPrior[parID]) - logdlnorm(oldValue, 0.0, m_sdLNormPrior[parID]);
            } else {               // r1VA - Uniform(0,5)
                if (newValue > 0.1 || newValue < 0) logRatioPrior = log(0); // FIXME : log of 0
            }
        	break;
        
        case 6:
        	if (m_lnormPrior[parID]) {  // r2VA - LogNormal(0,sd)
                logRatioPrior = logdlnorm(newValue, 0.0,  m_sdLNormPrior[parID]) - logdlnorm(oldValue, 0.0, m_sdLNormPrior[parID]);
            } else {               // r2VA - Uniform(0,5)
                if (newValue > 0.1 || newValue < 0) logRatioPrior = log(0); 
            }
        	break;
        
        case 7:
        	if (m_lnormPrior[parID]) {  // rUA - LogNormal(0,sd)
                logRatioPrior = logdlnorm(newValue, 0.0,  m_sdLNormPrior[parID]) - logdlnorm(oldValue, 0.0, m_sdLNormPrior[parID]);
            } else {               // rUA - Uniform(0,5)
                if (newValue > 5 || newValue < 0) logRatioPrior = log(0); 
            }
        	break;
        
        case 8:
        	if (m_lnormPrior[parID]) {  // rBA - LogNormal(0,sd)
                logRatioPrior = logdlnorm(newValue, 0.0,  m_sdLNormPrior[parID]) - logdlnorm(oldValue, 0.0, m_sdLNormPrior[parID]);
            } else {               // rBA - Uniform(0,5)
                if (newValue > 5 || newValue < 0) logRatioPrior = log(0); 
            }
        	break;
        
        case 9:
        	if (m_lnormPrior[parID]) {  // rQ - LogNormal(0,sd)
                logRatioPrior = logdlnorm(newValue, 0.0,  m_sdLNormPrior[parID]) - logdlnorm(oldValue, 0.0, m_sdLNormPrior[parID]);
            } else {               // rQ - Uniform(0,5)
                if (newValue > 5 || newValue < 0) logRatioPrior = log(0); 
            }
        	break;

        case 10:
        	if (m_lnormPrior[parID]) {  // rCC - LogNormal(0,sd)
                logRatioPrior = logdlnorm(newValue, 0.0,  m_sdLNormPrior[parID]) - logdlnorm(oldValue, 0.0, m_sdLNormPrior[parID]);
            } else {               // rCC - Uniform(0,5)
                if (newValue > 5 || newValue < 0) logRatioPrior = log(0); 
            }
        	break;
        
        case 11:
        	if (m_lnormPrior[parID]) {  // rMC - LogNormal(0,sd)
                logRatioPrior = logdlnorm(newValue, 0.0,  m_sdLNormPrior[parID]) - logdlnorm(oldValue, 0.0, m_sdLNormPrior[parID]);
            } else {               // rPC - Uniform(0,5)
                if (newValue > 5 || newValue < 0) logRatioPrior = log(0); 
            }
        	break;

        case 12:
        	if (m_lnormPrior[parID]) {  // rFC - LogNormal(0,sd)
                logRatioPrior = logdlnorm(newValue, 0.0,  m_sdLNormPrior[parID]) - logdlnorm(oldValue, 0.0, m_sdLNormPrior[parID]);
            } else {               // rPC - Uniform(0,5)
                if (newValue > 5 || newValue < 0) logRatioPrior = log(0); 
            }
        	break;

         case 13:
        	if (m_lnormPrior[parID]) {  // rFC - LogNormal(0,sd)
                logRatioPrior = logdlnorm(newValue, 0.0,  m_sdLNormPrior[parID]) - logdlnorm(oldValue, 0.0, m_sdLNormPrior[parID]);
            } else {               // rPC - Uniform(0,5)
                if (newValue > 5 || newValue < 0) logRatioPrior = log(0); 
            }
        	break;
    }

	// Compute new log likelihood
    m_parameter[parID] = newValue;
    std::vector<double> newLogLikHH(m_nHH, 0.0);
    double newLogLikGlobal = 0.0;
    for (size_t house = 0; house < m_nHH; house++) {
        newLogLikHH[house] = m_data[house].compute_log_lik(m_parameter, m_selected_param, m_max_PCR_detectability, display);        
        newLogLikGlobal += newLogLikHH[house];
    }

	// Update log likelihood
	double Q = newLogLikGlobal - m_globalLogLik + logRatioProposal + logRatioPrior;

	if (log(runif(m_rnd_generator)) < Q) { // FIXME : you just need to draw from an exponential (param=1)
		m_globalLogLik = newLogLikGlobal;
		m_hhLogLik = newLogLikHH;
        m_num_accepted_move[parID]++;
        m_accept[parID] = (1 + (step - 1) * m_accept[parID]) / step;
	} else {
		m_parameter[parID] = oldValue;
        m_accept[parID] = (0 + (step - 1) * m_accept[parID]) / step;
	}

    // Adaptative mcmc
    if (step <= m_adaptive_steps && m_adaptive_steps > 0) {
        double delta = 1.0 - step / m_adaptive_steps;
        double diff = m_accept[parID] - m_opt_acceptance;
        m_rate_for_random_walk[parID] *= 1.0 + delta * diff;
    }

    m_num_proposed_move[parID]++;
}

// Update augmented infection times
void McmcObject::update_augm_infec_times (bool display) {
    int inf;

    // Loop over houses // FIXME : why ++house?
    for (size_t h = 0; h < m_data.size(); ++h) {
        std::vector<int> infected = m_data[h].get_infec_list(); 
        // FIXME : rationalize this call

        // Loop over infected residents within house
        for (size_t i = 0; i < infected.size(); i++) {
            inf = infected[i];

            // Independent sampler from the incubation period
            double oldValue = m_data[h].get_augm_onset_time(inf);
            double newValue = m_data[h].new_infec_time(inf, m_rnd_generator);

            // Log ratio of the proposal (independent proposal) 
            double oldProposal = m_data[h].log_dIncub_g(inf, display);
            double newProposal = m_data[h].log_dIncub_g(inf, display);
            double logRatioProposal = oldProposal-newProposal;

            // Update infection time and person-to-person transmission rates
            m_data[h].set_infec_time(inf, newValue);
            m_data[h].update_infec_profiles(inf, display);   
            
            // Log likelihood
            double newLogLikOfHH = m_data[h].compute_log_lik(m_parameter, m_selected_param, m_max_PCR_detectability, display);
            double currentLogLik = m_hhLogLik[h];
            double differenceLogLik = newLogLikOfHH - currentLogLik;

            if( log(runif(m_rnd_generator)) < differenceLogLik+logRatioProposal ) { // FIXME : draw from exponential instead (param=1, cutoff?)
                m_globalLogLik += differenceLogLik;
                m_hhLogLik[h] = newLogLikOfHH;
                m_accepted_move_data++;

            } else {                
                m_data[h].set_infec_time(inf, oldValue);
                m_data[h].update_infec_profiles(inf);
            }

            m_proposed_move_data++;
        }
    }
}


// Update augmented symptom onsets
// FIXME : these two loops for update augmented should be integrated...
void McmcObject::update_augm_symptom_onset (bool display) {

    // Loop over houses
    for (size_t h = 0; h<m_data.size(); ++h) {
        std::vector<int> infected = m_data[h].get_infec_list();

        // Loop over infected individuals
        for (size_t i = 0; i < infected.size(); i++) {
            int inf = infected[i];

            // Independent sampler from uniform distribution [0,1]
            double oldValue = m_data[h].get_augm_onset_time(inf); 
            double newValue = m_data[h].new_onset_time(inf, m_rnd_generator);

            // Update infection time and person-to-person transmission rates
            m_data[h].set_onset_time(inf, newValue);
            m_data[h].update_infec_profiles(inf, display);   
            
            // Log likelihood
            double newLogLikOfHH = m_data[h].compute_log_lik(m_parameter, m_selected_param, m_max_PCR_detectability, display);
            double currentLogLik = m_hhLogLik[h];
            double differenceLogLik = newLogLikOfHH - currentLogLik;

            if(log(runif(m_rnd_generator)) < differenceLogLik) { // FIXME : draw from exponential instead (param=1, cutoff?)
                m_globalLogLik += differenceLogLik;
                m_hhLogLik[h] = newLogLikOfHH;
                m_accepted_move_data++;

            } else {                
                m_data[h].set_onset_time(inf, oldValue);
                m_data[h].update_infec_profiles(inf);
            }

            m_proposed_move_data++;
        }
    }
}