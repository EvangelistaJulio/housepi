// setp.hpp //

#define DEFAULT_SEED 20210415

// model + dataset
#define NOT_A_CASE 0
#define SYMPTOMATIC 1
#define NOT_SYMPTOM 2
#define SYMPTOM_UNK_ONSET 3
#define VACC_ONSET 90.

// MCMC
#define DEFAULT_NUM_ITER 2
#define DEFAULT_NUM_MCMC_SWEEEPS 20
#define DEFAULT_NUM_AUGM_ITER 5
#define DEFAULT_NUM_ADAPT 2
#define DEFAULT_ACCEPT_RATE 0.25

extern bool Verbose;