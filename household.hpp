#ifndef DEF_HOUSEHOLD__H
#define DEF_HOUSEHOLD__H
#include <string>
#include <vector>
#include <random>

// Constant parameters
// extern double maxPCRDetectability;
extern double mGamma;
extern double vGamma;
extern double shift;
extern double shapeInf;
extern double scaleInf;

#define INFINITY 1000.0

//----------------------------------------------------------------------
// Individual : Class declaration
//----------------------------------------------------------------------

class Individual {
	public:
		// Attributes
		int indid;
		bool is_case = 0;
		int age;
		int sex;
		double onset_time = 0.;
		double augm_onset_time = 0.; 
		int start_follow_up;
		int study_period;
		int vacc_status;
		int isolation;
		int vacc_delay;
		int recovered;
		double infec_time = INFINITY;

		// Constructors
		Individual(int indid, double onset_time, bool is_case, int age, int sex, int start_follow_up, int study_period, int vacc_status, int isolation, int vacc_delay, int recovered);
		Individual() = default;
		// FIXME : deconstructor?

		// Methods
		double infec_profile(double tinf, int study_period);
		double infec_cumulative(double tinf, int study_period);
};

//----------------------------------------------------------------------
// Household : Class declaration
//----------------------------------------------------------------------

class Household {
	// Attributes
	size_t h_size = 0;
	int h_notis_case = 0;
	int h_start_follow_up = -1;
	bool h_contact_pattern = false;
	double h_infec_time = INFINITY;

	std::vector<int> h_confCase;
	std::vector<Individual> h_residents;
	std::vector<double> h_inst_community_infec;
	std::vector<double> h_cum_community_infec;
	std::vector<std::vector<double>> h_cum_Lambda;
	std::vector<std::vector<double>> h_inst_Lambda;

	public:
		// Constructors
		Household(bool contact_pattern, std::vector<double> instval, std::vector<double>cumval);
		Household() = default ;
		// FIXME : deconstructor?
		
		// Getter
		size_t get_size() const;
		size_t num_infec_cases() const;
		size_t num_symptomatic() const;
		std::string get_resident_index() const;
		std::vector<int> get_infec_list() const;
		double get_augm_onset_time(int index) const;

		// Setter
		void add_individual(Individual this_indiv);
		void set_infec_time(int index, double infec_time);
		void set_onset_time(int index, double onset_time);
		void init_inf_time(int index, std::mt19937_64& rnd_generator, bool display = false);
		void init_onset_time(int index, std::mt19937_64& rnd_generator, bool display = false);

		// Computations
		void compute_infec_profiles(bool display = false);
		void update_infec_profiles(int ind, bool display = false);
		double new_infec_time(int index, std::mt19937_64& rnd_generator);
		double new_onset_time(int index, std::mt19937_64& rnd_generator);
		double log_dIncub(int index, bool display = false);
		double log_dIncub_g(int index, bool display = false);
		double compute_log_lik(std::vector<double> parameters, std::vector<int> selectedParam, double maxPCRDetectability, bool display = false);
		double log_pInf(size_t curr, std::vector<double> parameter, std::vector<int> selectedParam, bool display = false);
		double log_S(size_t curr, double t0, std::vector<double> parameter, std::vector<int> selectedParam, bool display = false);
};

#endif