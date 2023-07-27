#include "mcmc_class.hpp"

using namespace std;

/////////////////////////////////////////////////////
//  MCMC_Class : Class
/////////////////////////////////////////////////////

// ----------------------------
//  Constructors
// ----------------------------

// default 
MCMC_Class::MCMC_Class() {
    m_rand_gen.seed(12345);
    m_global_LogLike = 0.0;
    m_num_iter = 0;
    m_num_augm_iter = 0;
    m_hh_LogLike = {};
    m_data = {};
    m_num_adapt_iter = 0;
    m_display = false;

    reset_moves();
}

// user 
MCMC_Class::MCMC_Class(
    unsigned seed,
    int num_iter,
    int num_adapt_iter,
    std::vector<Household> data,
    int num_augm_iter,
    bool display
    ) {

    // Copy input vectors
    m_rand_gen.seed(seed);
    m_hh_LogLike.resize(data.size());
    m_num_iter = num_iter;
    m_num_adapt_iter = num_adapt_iter;
    m_num_augm_iter = num_augm_iter;
    m_data = data;
    m_display = display;

    reset_moves();
}

// ----------------------------
//  Methods : Setter
// ----------------------------

// Initialize param values
void MCMC_Class::init_param_values () {

    for (std::string parID : m_param_to_infer_names) {
        if (_param_infer(parID) == true){
            _param(parID).set_value(_param(parID).generate_from_prior(m_rand_gen)); 
        }
    }
}

// Modification of attributes
void MCMC_Class::reset_moves () {

    m_num_propos_data = 0;
    m_num_accept_data = 0;
    for (std::string parID : m_param_to_infer_names) {
        _param(parID).reset_moves();
    }
}


// Initialize infection time
void MCMC_Class::init_augm_data () {
    
    // loop over houses
    for (size_t h = 0; h < num_households(); h++) {
        // loop over residents within the house
        for (size_t i = 0; i < m_data[h].get_size(); i++) {
            m_data[h].init_onset_time(i, m_rand_gen, m_display);
            m_data[h].init_infec_time(i, m_rand_gen, m_display);
        }
        // Compute person-to-person transmission rate within household
        m_data[h].init_infec_profiles(m_display);
        // FIXME : update also?
    }
}

// Initialize infection time
void MCMC_Class::init_LogLike () {
    int nCases = 0;
    
    for (size_t house = 0; house < num_households(); house++) {
        m_hh_LogLike[house] = m_data[house].compute_loglike(m_display);
        m_global_LogLike += m_hh_LogLike[house];
        nCases += m_data[house].num_symptomatic();
    }
}

// ----------------------------
//  Methods : Computations
// ----------------------------

// Update param value in the mcmc chain
bool MCMC_Class::update_param (std::string parID) {
    bool move_accepted = false;

    // Random walk
	double old_value = _param_value(parID); // store old value in case of refusal
	double new_value, Q;
    double log_ratio_proposal, log_ratio_prior, new_global_LogLike;
    std::vector<double> new_LogLike_HH(num_households(), 0.0);

    // propose new value
    new_value = _param(parID).propose_move(m_rand_gen);

    // proceed only if the new value is within range, i. e. no infinities
    if (_param(parID).within_range(new_value) == true) {
        log_ratio_proposal = _param(parID).log_ratio_proposal();
        log_ratio_prior = _param(parID).log_ratio_prior(new_value);

        // Compute new log likelihood
        _param(parID).set_value(new_value);
        new_global_LogLike = 0.0;
        for (size_t house = 0; house < num_households(); house++) {
            new_LogLike_HH[house] = m_data[house].compute_loglike(m_display);        
            new_global_LogLike += new_LogLike_HH[house];
        }

        // Metropolis
        Q = new_global_LogLike - m_global_LogLike + log_ratio_proposal + log_ratio_prior;

        if (unif_rand(m_rand_gen) < std::exp(Q)) {      // move accepted
            move_accepted = true;
            _param(parID).set_num_accept_pp(); 
            // Update log likelihood
            m_global_LogLike = new_global_LogLike;
            m_hh_LogLike = new_LogLike_HH;
        
        } else {                                        // move refused
            _param(parID).set_value(old_value);
        }
    } 
    _param(parID).set_num_propos_pp();

    return move_accepted;
}


// Update augmented infection times
void MCMC_Class::update_augm_infec_times () {
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
            old_proposal = m_data[h].log_incub_proba(i, m_display);
            new_proposal = m_data[h].log_incub_proba(i, m_display);

            // Update infection time and person-to-person transmission rates
            m_data[h].set_augm_infec_time(i, new_value);
            
            // Log likelihood
            old_LogLike_of_HH = m_hh_LogLike[h];
            new_LogLike_of_HH = m_data[h].compute_loglike(m_display);

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
void MCMC_Class::update_augm_symptom_onset () {
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
            new_LogLike_of_HH = m_data[h].compute_loglike(m_display);
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

void MCMC_Class::run (std::string output_file, int num_sweeps) {

	ofstream output(output_file.c_str());

    int step;
    bool move_accepted;
	int num_iter = int(get_num_iter() / num_sweeps);
	int num_augm_iter = get_num_augm_iter();

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
	for (int i = 0; i < num_iter; i++)	{
        cout << "iteration " << i << endl;
	    reset_moves();

		for (int j = 0; j < num_sweeps; j++) {
	        for (std::string parID : m_param_to_infer_names) {
                step = i * num_sweeps + j + 1;
				move_accepted = update_param(parID);
                
                // Adaptative MCMC
                _param(parID).adapt_rate(move_accepted, step, m_num_adapt_iter);
			}

			// Augmented data
			for (int k = 0; k < num_augm_iter; k++) {
                update_augm_infec_times(); // Update loglik at the household level
                update_augm_symptom_onset();
			}
		}

		// Print log likelihood, parameter values, n. of proposed/accepted move per parameter in the output file
		output << (i + 1) * num_sweeps << " " << global_LogLike() << " "; // Log likelihood
        for (std::string parID : m_param_to_infer_names) {
            output << _param_value(parID) << " " << _param(parID).get_num_propos() << " " << _param(parID).get_num_accept() << " ";
        }
        output << m_num_propos_data << " " << m_num_accept_data << endl;
	}

	output.close();
}