// setup.hpp //

#define DEFAULT_SEED 20210415
#define TAB "\t"

#define CONTINUOS_TIME 

// model + dataset
#define HEALTHY 0
#define SYMPTOMATIC 1
#define NOT_SYMPTOM 2
#define SYMPTOM_UNK_ONSET 3
#define VACC_ONSET 90.
#define MAX_INCUB_PERIOD 14.

// MCMC
#define DEFAULT_NUM_ITER 10
#define DEFAULT_NUM_PARAM_SWEEEPS 5
#define DEFAULT_NUM_AUGM_SWEEPS 5
#define DEFAULT_NUM_ADAPT_STEPS 50
#define DEFAULT_ACCEPT_RATE 0.25

extern bool Verbose;