/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__DATA__H
#define HOUSEPI__DATA__H

#include <memory>

#include "tools/types.h"
#include "household.h"

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** Setup : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp
{
    namespace progress
    {
        class Bar;
    }
}

struct Setup
{
    Setup();
    ~Setup();

    const str_vec file_headers = {"steps", "burnin", "thin",
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

    unsigned long c_sample_size = 0; /// single chain sample size
    unsigned long c_n_sweeps = 0;    /// single chain n. of sweeps
    double inv_T = 1.;               /// inverse temperature
    str_vec id_draws = {};           /// augmented trace to save

    unsigned long chain = 1; /// current MCMC chain
    unsigned long rung = 0;  /// current MCMC thermodynamic rung

    bool is_power_post() const { return n_rungs > 0; }
    std::string fmt_progtxt();

    void set_n_threads();
    void build(const u_dict<str_vec> &deftable);
    void assert_mcmc_setup();
    void drop_mcmc_adapt();

    // progress bar
    void prog_declare(const std::string &text = hsp::na::STRING);
    void prog_set_text(const std::string &text = hsp::na::STRING);
    void prog_open_for_inference();
    void prog_open_for_power_post();
    void prog_update();
    void prog_close();
    void prog_crash();

    /// Set seed for pseudo random number generator.
    void set_prng_seed();
    /// Return pseudo random number generator.
    rnd_gen &prng();
    /// Sample a uniform distributed variable between 0 and 1.
    double unif() { return unif_(prng()); }

private:
    std::vector<rnd_gen> prng_parallel_ = {};         // pseudo-random number generators
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
    str_vec headers = {};

    // list of all infected statuses in the dataset
    int_vec all_infected_status_list = {};
    std::vector<dbl_vec> smpl_values = {};
    std::vector<size_t> head_class_house_indexes = {};

    std::unique_ptr<ReportIndividual> report_indiv;
    std::unique_ptr<ReportHousehold> report_house;
    std::unique_ptr<ReportHouseSizes> report_hh_size;

    /////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  DataWide : constructor  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////

    DataWide();
    ~DataWide();

    Household &operator()(size_t H) { return house_vector_[H]; }
    const Household &operator()(size_t H) const { return house_vector_.at(H); }
    Household &operator()(const std::string &house_id) { return (*this)(house_id_to_index_(house_id)); }
    const Household &operator()(const std::string &house_id) const { return (*this)(house_id_to_index_(house_id)); }

    /////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  DataWide : getter  //////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////

    bool is_empty() const { return house_vector_.empty(); }
    Household copy(size_t H) const { return house_vector_[H]; }

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

    void build(const u_dict<str_vec> &deftable);
    void multiplicity_build(const u_dict<str_vec> &deftable, bool is_any_augm = false);

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
    u_dict<size_t> hhid_to_index_ = {};

    size_t house_id_to_index_(const std::string &house_id) const
    {
        if (hhid_to_index_.count(house_id))
            return hhid_to_index_.at(house_id);
        throw std::invalid_argument(
            "Index '" + house_id + "' did not match any household in the data.");
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
    str_vec headers = {};
    std::vector<u_dict<dbl_vec>> indiv_datalist = {};

    DataLong(std::shared_ptr<DataWide> ptr)
        : pDataWide_(ptr) {}
    ~DataLong() = default;

    void build(const u_dict<str_vec> &deftable);

private:
    std::shared_ptr<DataWide> pDataWide_; // pointer to external DataWide struct
};

#endif // HOUSEPI__DATA__H
