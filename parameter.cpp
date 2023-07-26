#include "parameter.hpp"

using namespace std;

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
    init(
        name,
        value,
        infer,
        prior,
        args,
        metropolis_strategy,
        metropolis_rate,
        range
    );
}

// Auxiliary method for Parameter class initialization.
void Parameter::init (
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

    //  inference mode setup  /
    p_infer = infer;
    if (infer == true) {

        //  prior function
        if (prior == "lognorm") {
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
        p_metropolis_strategy = metropolis_strategy;

        // check metropolis_rate values
        if (metropolis_rate <= 0.) {
            // FIXME : print warning
            p_metropolis_rate = 0.5;

        } else {
            p_metropolis_rate = metropolis_rate;
        }

        // check range values
        if (range.size() != 2 || range[0] >= range[1]) {
            cout << "Inconsistent range for parameter " << name << endl;
            exit (EXIT_FAILURE);
        } 
    }

    reset_accept_propos ();
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

void Parameter::set_rate_metropolis (double rate) {
    p_metropolis_rate = rate;
}

// ----------------------------
//  Methods : Setter (mcmc internal)
// ----------------------------

void Parameter::reset_accept_propos () {
    _metropolis_num_propos = 0;
    _metropolis_num_accept = 0;
    _metropolis_accept = 0.0;
}

// ----------------------------
//  Methods : Computation
// ----------------------------

bool Parameter::within_range (double value) {
    return (value >= p_range[0]) && (value <= p_range[1]);
}

double Parameter::generate_from_prior (std::mt19937_64& rand_gen) {
    return p_prior_rand(rand_gen, p_prior_args);
}

double Parameter::metropolis (std::mt19937_64& rand_gen) {
    double new_value = 0.0;

    if (p_metropolis_strategy == "exp") { 
        new_value = p_value * std::exp(p_metropolis_rate * norm_rand(rand_gen));

    } else if (p_metropolis_strategy == "norm") { 
        std::vector<double> args = {0.0, p_metropolis_rate}; 
        new_value = p_value + norm_rand(rand_gen, args);

    } else {
        cout << "Unknown metropolis strategy " << p_metropolis_strategy << endl;
        exit (EXIT_FAILURE);  
    }
    return new_value;
}

double Parameter::metropolis_log_ratio_proposal () {
    if (p_metropolis_strategy == "exp") { 
        return p_metropolis_rate;
        
    } else if (p_metropolis_strategy == "norm") { 
        return 0.;

    } else {
        cout << "Unknown metropolis strategy " << p_metropolis_strategy << endl;
        std::exit (EXIT_FAILURE);  
    }
}

double Parameter::metropolis_log_ratio_prior (double new_value) {
    return p_prior_log_pdf(new_value, p_prior_args) - p_prior_log_pdf(p_value, p_prior_args);
}

// aliases
Parameter _param (std::string parID) {
    return Parameter_Dictionary.at(parID);
}

double _param_value (std::string parID) {
    return Parameter_Dictionary.at(parID).get_value();
}

bool _param_infer (std::string parID) {
    return Parameter_Dictionary.at(parID).p_infer;
}