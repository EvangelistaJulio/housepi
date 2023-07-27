#include "parameter.hpp"

using namespace std;

/////////////////////////////////////////////////////
//  Metropolis_Strategy : Class
/////////////////////////////////////////////////////

Metropolis_Strategy::Metropolis_Strategy (double new_rate, std::vector<double> new_range){
    // check metropolis_rate values
    if (new_rate <= 0. || new_rate >= 1.) {
        rate = 0.5;
        // NOTE : print warning

    } else {
        rate = new_rate;
    }

    // check range values
    if (range.size() != 2 || range[0] >= range[1]) {
        cout << "Provided an inconsistent range for parameter." << endl;
        exit (EXIT_FAILURE);

    } else {
        range = new_range;
    }

    reset_moves();
}

void Metropolis_Strategy::reset_moves () {
	num_propos = 0;
    num_accept = 0;
    accept = 0.0;
}

bool Metropolis_Strategy::within_range (double new_value) {
    return (new_value > range[0]) && (new_value < range[1]);
}

void Metropolis_Strategy::adapt_rate (bool move_accepted, int step, int num_adapt_iter) {
    double delta, diff;

    accept = ((move_accepted + (step - 1) * accept) / step);
    
    if (num_adapt_iter > 0 && step <= num_adapt_iter) {
        delta = 1.0 - step / num_adapt_iter;
        diff = accept - MCMC_OPTIMAL_ACCEPTANCE;
        set_rate(1.0 + delta * diff);
    }
}

// --- Strategy_exp ---
double Strategy_exp::propose_move (std::mt19937_64& rand_gen, double old_value) {
    double new_value = 0.0;
    new_value = old_value * std::exp(rate * normal_rand(rand_gen));

    return new_value;
}

double Strategy_exp::log_ratio_proposal () {
    return rate;
}

// --- Strategy_normal ---
double Strategy_normal::propose_move (std::mt19937_64& rand_gen, double old_value) {
    double new_value = 0.0;
    std::vector<double> args = {0.0, rate}; 
    new_value = old_value + normal_rand(rand_gen, args);

    return new_value;
}

double Strategy_normal::log_ratio_proposal () {
    return 0.;
}

/////////////////////////////////////////////////////
//  Parameter : Class
/////////////////////////////////////////////////////

// ----------------------------
//  Constructors
// ----------------------------

Parameter::Parameter (
    std::string name,
    double value,
    bool infer,
    std::string prior,
    std::vector<double> args,
    std::string metropolis_strategy,
    double metropolis_rate,
    std::vector<double> range
    ) {
    p_name = name;
    p_value = value;    

    // --- inference mode setup --- //
    p_infer = infer;
    if (infer == true) {

        // --- prior function --- //
        if (prior == "normal") {
            p_prior_rand = normal_rand;
            p_prior_log_pdf = log_normal_pdf;

        } else if (prior == "lognorm") {
            p_prior_rand = lognorm_rand;
            p_prior_log_pdf = log_lognorm_pdf;

        } else if (prior == "gamma") {
            p_prior_rand = gamma_rand;
            p_prior_log_pdf = log_gamma_pdf;

        } else if (prior == "unif") {
            p_prior_rand = unif_rand;
            p_prior_log_pdf = nullptr; // FIXME : look for a smarter way

        } else {
            cout << "Unknown prior " << prior << endl;
            exit (EXIT_FAILURE);    
        }

        p_prior = prior;
        p_prior_args = args;

        // --- metropolis strategy --- //
        if (metropolis_strategy == "exp") { 
            _metropolis = new Strategy_exp(metropolis_rate, range);

        } else if (metropolis_strategy == "normal") { 
            _metropolis = new Strategy_normal(metropolis_rate, range);

        } else {
            cout << "Unknown metropolis strategy " << metropolis_strategy << endl;
            exit (EXIT_FAILURE);  
        }
    }

    reset_moves ();
}

// ----------------------------
//  Methods : Getter
// ----------------------------

double Parameter::get_value () const {
    return p_value;
}

// ----------------------------
//  Methods : Setter
// ----------------------------

void Parameter::set_value (double value) {
    p_value = value;
}

void Parameter::set_rate_metropolis (double new_rate) {
    _metropolis->set_rate(new_rate);
}

// ----------------------------
//  Methods : Computation
// ----------------------------

double Parameter::generate_from_prior (std::mt19937_64& rand_gen) {
    return p_prior_rand(rand_gen, p_prior_args);
}

double Parameter::log_ratio_prior (double new_value) {
    return p_prior_log_pdf(new_value, p_prior_args) - p_prior_log_pdf(p_value, p_prior_args);
}

// -----------------------------------------------
// --- aliases for shared Parameter Dictionary ---
// -----------------------------------------------

Parameter _param (std::string parID) {
    return Parameter_Dictionary.at(parID);
}

double _param_value (std::string parID) {
    return Parameter_Dictionary.at(parID).get_value();
}

bool _param_infer (std::string parID) {
    return Parameter_Dictionary.at(parID).p_infer;
}