#include "mcmc_object.hpp"

using namespace std;

/////////////////////////////////////////////////////
//  McmcObject : Class
/////////////////////////////////////////////////////

// ----------------------------
//  Constructors
// ----------------------------

// default 
McmcObject::McmcObject() {
    m_rand_gen.seed(12345);
    m_global_LogLike = 0.0;
    m_iterations = 0;
    m_iter_time_infec = 0;
    m_hh_LogLike = {};
    m_data = {};
    m_adaptive_steps = 0;

    reset_moves();
}

// user 
McmcObject::McmcObject(
    unsigned seed,
    size_t num_iterations,
    size_t num_adaptive_steps,
    std::vector<Household> data,
    size_t n_iter_time_infec
    ) {

    // Copy input vectors
    m_rand_gen.seed(seed);
    m_hh_LogLike.resize(data.size());
    m_iterations = num_iterations;
    m_adaptive_steps = num_adaptive_steps;
    m_iter_time_infec = n_iter_time_infec;
    m_data = data;

    reset_moves();
}

// ----------------------------
//  Methods : Setter
// ----------------------------

// Initialize param values
void McmcObject::init_param_values () {

    for (std::string parID : m_param_to_infer_names) {
        if (_param_infer(parID) == true){
            _param(parID).set_value(_param(parID).generate_from_prior(m_rand_gen)); 
        }
    }
}

// Modification of attributes
void McmcObject::reset_moves () {

    m_num_propos_data = 0;
    m_num_accept_data = 0;
    for (std::string parID : m_param_to_infer_names) {
        _param(parID).reset_accept_propos();
    }
}


// Initialize infection time
void McmcObject::init_augm_data (bool display) {
    
    // loop over houses
    for (size_t h = 0; h < num_households(); h++) {
        // loop over residents within the house
        for (size_t i = 0; i < m_data[h].get_size(); i++) {
            m_data[h].init_onset_time(i, m_rand_gen, display);
            m_data[h].init_infec_time(i, m_rand_gen, display);
        }
        // Compute person-to-person transmission rate within household
        m_data[h].init_infec_profiles(display);
        // FIXME : update also?
    }
}

// Initialize infection time
void McmcObject::init_LogLike (bool display) {
    int nCases = 0;
    
    for (size_t house = 0; house < num_households(); house++) {
        m_hh_LogLike[house] = m_data[house].compute_loglike(display);
        m_global_LogLike += m_hh_LogLike[house];
        nCases += m_data[house].num_symptomatic();
    }
}

// ----------------------------
//  Methods : Computations
// ----------------------------

// Update param value in the mcmc chain
void McmcObject::update_param (std::string parID, size_t step, bool display) {
    bool accept_move = false;

    // Random walk
	double old_value = _param_value(parID);
	double new_value, Q;
    double log_ratio_proposal, log_ratio_prior, new_global_LogLike;
    std::vector<double> new_LogLike_HH(num_households(), 0.0);

    // propose new value
    new_value = _param(parID).metropolis(m_rand_gen);

    // proceed only if the new value is within range, i. e. no infinities
    if (_param(parID).within_range(new_value) == true) {
        log_ratio_proposal = _param(parID).metropolis_log_ratio_proposal();
        log_ratio_prior = _param(parID).metropolis_log_ratio_prior(new_value);

        // Compute new log likelihood
        _param(parID).set_value(new_value);
        new_global_LogLike = 0.0;
        for (size_t house = 0; house < num_households(); house++) {
            new_LogLike_HH[house] = m_data[house].compute_loglike(display);        
            new_global_LogLike += new_LogLike_HH[house];
        }

        // Metropolis
        Q = new_global_LogLike - m_global_LogLike + log_ratio_proposal + log_ratio_prior;

        if (unif_rand(m_rand_gen) < std::exp(Q)) {
            accept_move = true;
            _param(parID).set_num_accept_pp(); 
            // Update log likelihood
            m_global_LogLike = new_global_LogLike;
            m_hh_LogLike = new_LogLike_HH;

        } else {
            _param(parID).set_value(old_value);
        }
    } 
    _param(parID).set_num_propos_pp();

    //  --- Adaptative MCMC ---  //
    _param(parID).set_accept((accept_move + (step - 1) * _param(parID).get_accept()) / step);
    if (m_adaptive_steps > 0 && step <= m_adaptive_steps) {
        double delta = 1.0 - step / m_adaptive_steps;
        double diff = _param(parID).get_accept() - MCMC_OPTIMAL_ACCEPTANCE;
        _param(parID).set_rate_metropolis(1.0 + delta * diff);
    } 
    // --- --- //

}


// Update augmented infection times
void McmcObject::update_augm_infec_times (bool display) {
    double old_value, new_value, Q;
    double old_proposal, new_proposal;
    double old_LogLike_of_HH, new_LogLike_of_HH;

    // Loop over households
    for (size_t h = 0; h < m_data.size(); h++) {
        // Loop over infected residents within house
        for (size_t i = 0; i < m_data[h].get_size(); i++) {

            // Independent sampler from the incubation period
            old_value = m_data[h].get_augm_onset_time(i);
            new_value = m_data[h].new_infec_time(i, m_rand_gen);

            // Log ratio of the proposal (independent proposal) 
            // FIXME : what is this???
            old_proposal = m_data[h].log_incub_proba(i, display);
            new_proposal = m_data[h].log_incub_proba(i, display);

            // Update infection time and person-to-person transmission rates
            m_data[h].set_augm_infec_time(i, new_value);
            
            // Log likelihood
            old_LogLike_of_HH = m_hh_LogLike[h];
            new_LogLike_of_HH = m_data[h].compute_loglike(display);

            Q = new_LogLike_of_HH - old_LogLike_of_HH + old_proposal - new_proposal;

            if (unif_rand(m_rand_gen) < std::exp(Q)) { 
                m_global_LogLike += new_LogLike_of_HH - old_LogLike_of_HH;
                m_hh_LogLike[h] = new_LogLike_of_HH;
                m_num_accept_data++;

            } else {                
                m_data[h].set_augm_infec_time(i, old_value);
            }

            m_num_propos_data++;
        }
    }
}

// Update augmented symptom onsets
// FIXME : these two loops for update augmented should be integrated...
void McmcObject::update_augm_symptom_onset (bool display) {
    double old_value, new_value, Q;
    //double old_proposal, new_proposal;
    double old_LogLike_of_HH, new_LogLike_of_HH;

    // Loop over households
    for (size_t h = 0; h < m_data.size(); h++) {
        // Loop over infected residents of the househould
        for (size_t i = 0; i < m_data[h].get_size(); i++) {

            // Independent sampler from uniform distribution [0, 1]
            old_value = m_data[h].get_augm_onset_time(i); 
            new_value = m_data[h].new_onset_time(i, m_rand_gen);

            // Log ratio of the proposal (independent proposal) 
            // FIXME :why here there's no proposal?

            // Update infection time and person-to-person transmission rates
            m_data[h].set_augm_onset_time(i, new_value);
            
            // Log likelihood
            new_LogLike_of_HH = m_data[h].compute_loglike(display);
            old_LogLike_of_HH = m_hh_LogLike[h];

            Q = new_LogLike_of_HH - old_LogLike_of_HH;

            if (unif_rand(m_rand_gen) < std::exp(Q)) { 
                m_global_LogLike += new_LogLike_of_HH - old_LogLike_of_HH;
                m_hh_LogLike[h] = new_LogLike_of_HH;
                m_num_accept_data++;

            } else {                
                m_data[h].set_augm_onset_time(i, old_value);
            }

            m_num_propos_data++;
        }
    }
}


//----------------------------------------------------------------------
// Run mcmc
//----------------------------------------------------------------------

void McmcObject::run(
    std::string output_file,
    size_t num_sweeps,
    bool display
    ) {

	ofstream output(output_file.c_str());

	size_t num_iter = int(num_iterations() / num_sweeps);
	size_t n_iter_time_infec = get_num_iter_time_infec();

	// Column names of the output file
	std::string colNames = "iteration LogLik ";
	for (std::string parID : m_param_to_infer_names){
        colNames += parID + " " + parID + "_p " + parID + "_a ";
    }
	colNames += "data_p data_a";
	output << colNames << endl;

	// MCMC chain
	cout << "=============MCMC===============" << endl;

    // Initial state
    init_param_values();
	init_augm_data();
	init_LogLike(); 
	cout << "Initial log likelihood: " << global_LogLike() << endl;

    output << "0 " << global_LogLike() << " "; // Log likelihood
	for (std::string parID : m_param_to_infer_names) {
        output << _param_value(parID) << " 0 0 ";
    }
    output << "0 0" << endl;

    // Chain
	for (size_t i = 0; i < num_iter; i++)	{
        cout << "iteration " << i << endl;
	    reset_moves();

		for (size_t j = 0; j < num_sweeps; j++) {
	        for (std::string parID : m_param_to_infer_names) {
				update_param(parID, i * num_sweeps + j + 1, display);
			}
            //cout << "num_iter " << iter << endl;

			// Augmented data
			for (size_t k = 0; k < n_iter_time_infec; k++) {
                update_augm_infec_times(display); // Update loglik at the household level
                update_augm_symptom_onset(display);
			}
            //cout << iteration * num_sweeps + iter << " " << global_LogLike() << endl;
		}

		// Write log likelihood, parameter values, number of proposed/accepted move per parameter in the output file
		output << (i + 1) * num_sweeps << " " << global_LogLike() << " "; // Log likelihood
        for (std::string parID : m_param_to_infer_names) {
            output << _param_value(parID) << " " << _param(parID).get_num_propos() << " " << _param(parID).get_num_accept() << " ";
        }
        output << m_num_propos_data << " " << m_num_accept_data << endl;
	}

	output.close();
}