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
#include <vector>
#include <unordered_set>
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

struct ReportIndividual;
struct ReportHousehold;
struct ReportHouseSizes;

struct DataWide
{
    std::vector<std::string> headers = {}; /// Headers of input data

    // list of all infected statuses in the dataset
    std::unordered_set<int> all_infected_statuses = {};
    std::vector<std::vector<double>> smpl_values = {};
    std::vector<size_t> head_class_house_indexes = {};

    std::unique_ptr<ReportIndividual> report_indiv;
    std::unique_ptr<ReportHousehold> report_house;
    std::unique_ptr<ReportHouseSizes> report_hh_size;

    /////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  DataWide : constructor  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////

    DataWide();
    ~DataWide();

    Household &at(size_t H) { return house_vector_[H]; }
    const Household &at(size_t H) const { return house_vector_.at(H); }
    Household &at(std::string_view house_id) { return at(house_id_to_index_(house_id)); }
    const Household &at(std::string_view house_id) const { return at(house_id_to_index_(house_id)); }

    /////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  DataWide : getter  //////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////

    bool is_empty() const { return house_vector_.empty(); }

    /////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  DataWide : flags  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////

    /// Assert data_wide is not empty.
    void assert_non_empty();
    /// Warn if all individuals have missing status.
    void warn_all_missing_status();

    /////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  DataWide : setter  //////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////

    void build(const hsp::str_table &input);
    void multiplicity_build(const hsp::str_table &input);

    void init_draws(unsigned long sample_size);
    void store_draws(unsigned long smpl_count);

    /////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  DataWide : iterators  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////

    auto begin() { return house_vector_.begin(); }
    auto begin() const { return house_vector_.begin(); }
    auto end() { return house_vector_.end(); }
    auto end() const { return house_vector_.end(); }

    template <typename T>
    class View
    {
        // conditional to select the vector type if constant or not
        using VectorRef = std::conditional_t<std::is_const<T>::value,
                                             const std::vector<Household> &,
                                             std::vector<Household> &>;

    public:
        VectorRef vec;
        const std::vector<size_t> &indexes; /// NOTE: no need to copy here

        View(VectorRef m, const std::vector<size_t> &i)
            : vec(m), indexes(i) {}

        class iterator
        {
        public:
            iterator(View *v, size_t p) : view_(v), pos_(p) {}

            T &operator*() const { return view_->vec[view_->indexes[pos_]]; }
            T *operator->() const { return &view_->vec[view_->indexes[pos_]]; }

            iterator &operator++()
            {
                ++pos_;
                return *this;
            }
            iterator operator++(int)
            {
                iterator tmp = *this;
                ++(*this);
                return tmp;
            }

            bool operator==(const iterator &other) const { return pos_ == other.pos_; }
            bool operator!=(const iterator &other) const { return !(*this == other); }

        private:
            View *view_;
            size_t pos_;
        };

        iterator begin() { return iterator(this, 0); }
        iterator end() { return iterator(this, indexes.size()); }
    };

    View<Household> restricted() { return View<Household>(house_vector_, head_class_house_indexes); }
    View<const Household> restricted() const { return View<const Household>(house_vector_, head_class_house_indexes); }
    Household &restricted(size_t R) { return house_vector_[head_class_house_indexes[R]]; }

private:
    std::vector<Household> house_vector_ = {};
    hsp::dict<size_t> hhid_to_index_ = {};

    size_t house_id_to_index_(std::string_view house_id) const
    {
        auto it = hhid_to_index_.find(house_id);
        if (it != hhid_to_index_.end())
            return it->second;
        throw std::invalid_argument(
            "Index '" + std::string(house_id) + "' did not match any household in the data.");
    }

    /// Clear data report.
    void build_helper_clear_data_();
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
    std::vector<std::string> headers = {};
    std::vector<hsp::dbl_table> indiv_datalist = {};

    DataLong(std::shared_ptr<DataWide> ptr)
        : pDataWide_(ptr) {}
    ~DataLong() = default;

    void build(const hsp::str_table &input);

private:
    std::shared_ptr<DataWide> pDataWide_; // pointer to external DataWide struct
};
