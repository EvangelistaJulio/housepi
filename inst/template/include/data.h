/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <memory>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

#include "tools/types.h"

#include "household.h"

namespace hsp::rates
{
    constexpr const char *REFERENCE = ".reference";  /// Rate reference choice.
    constexpr const char *MISSING = ".missing_rate"; /// Missing rate.
    constexpr const char *CONTACT = ".contact";      /// Contanct rate.
    constexpr double VALUE_REFERENCE = 1.0;          /// Alias for reference value.
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** Setup : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp::progress
{
    class Bar;
}

struct Setup
{
    Setup();
    ~Setup();

    const std::vector<std::string> file_headers = {"steps", "burnin", "thin",
                                                   "adapt", "chains", "threads",
                                                   "id_draws", "seed", "progress"};

    unsigned long n_steps = 1;  /// n. iterations
    unsigned long n_burnin = 0; /// n. burnin iterations
    unsigned long n_thin = 1;   /// n. thinning iterations
    unsigned long n_adapt = 0;  /// n. (finite) adaptive steps
    unsigned long n_chains = 1; /// n. chains
    unsigned n_threads = 1;     /// n. of OMP threads
    std::int32_t seed;          /// random generator seed
    unsigned long progress = 1; /// thinning seconds for progress bar
    unsigned long n_rungs = 0;  /// n. of rungs

    unsigned long c_sample_size = 0;        /// single chain sample size
    unsigned long c_n_sweeps = 0;           /// single chain n. of sweeps
    double inv_T = 1.;                      /// inverse temperature
    std::vector<std::string> id_draws = {}; /// augmented trace to save

    unsigned long chain = 1; /// current MCMC chain
    unsigned long rung = 0;  /// current MCMC thermodynamic rung

    bool is_power_post() const { return n_rungs > 0; }
    std::string fmt_progtxt();

    void set_n_threads();
    void build(const hsp::str_table &input);
    void assert_mcmc_setup();
    void drop_mcmc_adapt();

    // progress bar
    void prog_declare(const std::string &text = "");
    void prog_set_text(const std::string &text = "");
    void prog_open_for_inference();
    void prog_open_for_power_post();
    void prog_update();
    void prog_close();
    void prog_crash();

    /// Set seed for pseudo random number generator.
    void set_prng_seed();
    /// Return pseudo random number generator.
    hsp::rnd_gen &prng();
    /// Sample a uniform distributed variable between 0 and 1.
    double unif() { return unif_(prng()); }

private:
    std::vector<hsp::rnd_gen> prng_parallel_ = {};    // pseudo-random number generators
    std::uniform_real_distribution<double> unif_;     // uniform sampler
    std::unique_ptr<hsp::progress::Bar> progressbar_; // progress bar
};

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  *** DataWide : struct ***  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp::key
{
    class Registry;
}

struct ReportIndividual;
struct ReportHousehold;
struct ReportHouseSizes;

struct DataWide
{
    std::unique_ptr<hsp::key::Registry> pRegistry;
    std::unique_ptr<hsp::key::Registry> p_size_Registry;
    std::unique_ptr<hsp::key::Registry> p_size_vec_Registry;
    std::unique_ptr<hsp::key::Registry> p_double_vec_Registry;
    std::unique_ptr<hsp::key::Registry> p_double_mtx_Registry;

    // list of all infected statuses in the dataset
    std::unordered_set<int> all_infected_statuses = {};
    std::vector<std::vector<double>> smpl_values = {};
    std::vector<size_t> head_class_house_indexes = {};

    std::unique_ptr<ReportIndividual> report_indiv;
    std::unique_ptr<ReportHousehold> report_house;
    std::unique_ptr<ReportHouseSizes> report_hh_size;

    ///////////////////////////
    //////  CONSTRUCTOR  //////
    ///////////////////////////

    DataWide(const std::vector<std::string> &keys,
             const std::vector<std::string> &size_keys,
             const std::vector<std::string> &size_vec_keys,
             const std::vector<std::string> &double_vec_keys,
             const std::vector<std::string> &double_mtx_keys);
    ~DataWide();

    Household &at(size_t H) { return house_vector_[H]; }
    const Household &at(size_t H) const { return house_vector_[H]; }
    Household &at(std::string_view house_id) { return at(house_id_to_index_(house_id)); }
    const Household &at(std::string_view house_id) const { return at(house_id_to_index_(house_id)); }

    //////////////////////
    //////  SETTER  //////
    //////////////////////

    void build(const hsp::str_table &input);
    void multiplicity_build(const hsp::str_table &input, bool is_any_augm);

    /// Initialize cached structeres according to household size.
    void init_house_cache(Household &house);
    /// Deep copy all household cache attributes in the buffer.
    void buffer_all_cache(Household &house);
    /// Deep copy member cache attributes in the buffer.
    void buffer_member_cache(Household &house, size_t i);
    /// Recover all buffered values of cached depedencies.
    void recover_all_buffer(Household &house);
    /// Recover all buffered values of member cached depedencies.
    void recover_member_buffer(Household &house, size_t i);

    void assert_finite_log_lik() const;
    void init_draws(unsigned long sample_size);
    void store_draws(unsigned long smpl_count);

    /////////////////////
    //////  FLAGS  //////
    /////////////////////

    /// Wheter data have been provided.
    bool is_empty() const noexcept { return house_vector_.empty(); }
    /// Assert data_wide is not empty.
    void assert_non_empty();
    /// Warn if all individuals have missing status.
    void warn_all_missing_status();

    /////////////////////////
    //////  ITERATORS  //////
    /////////////////////////

    auto begin() noexcept { return house_vector_.begin(); }
    auto begin() const noexcept { return house_vector_.begin(); }
    auto end() noexcept { return house_vector_.end(); }
    auto end() const noexcept { return house_vector_.end(); }

    Household &restricted(size_t R) { return house_vector_[head_class_house_indexes[R]]; }

private:
    std::vector<Household> house_vector_ = {};
    hsp::dict<size_t> hhid_to_index_ = {};

    /// Get housepi household index from house_id.
    size_t house_id_to_index_(std::string_view house_id) const;
    /// Clear previous data report.
    void helper_build_clear_data_();
    /// Assert that data headers are correct.
    void helper_build_data_headers_(const hsp::str_table &input);
    /// Use the loaded data to fill-in the respective report.
    void build_helper_fill_report_data_();
};

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  *** DataLong : struct ***  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

struct DataLong
{
public:
    std::unique_ptr<hsp::key::Registry> pRegistry; /// pointer to registry

    std::vector<std::string> headers = {};
    std::vector<std::vector<std::vector<double>>> indiv_datalist = {};

    DataLong(const std::vector<std::string> &keys,
             DataWide *ptr);
    ~DataLong();

    void build(const hsp::str_table &input);

private:
    DataWide *pDataWide_; // pointer to external DataWide struct

    /// Clear previous data report.
    void helper_build_clear_data_();
    /// Assert that data headers are correct.
    void helper_build_data_headers_(const hsp::str_table &input);
};
