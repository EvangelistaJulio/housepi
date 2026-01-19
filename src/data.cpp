/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <chrono>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/global.h"
#include "tools/messages.h"
#include "tools/openmp.h"
#include "tools/progress.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"

#include "tools/keys.h"
#include "tools/keys_helper.h"

#include "attribute.h"
#include "data.h"
#include "data_reports.h"

namespace
{
    using std::string;
    using std::vector;
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** Setup : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

Setup::Setup() : unif_(0., 1.),
                 progressbar_(std::make_unique<hsp::progress::Bar>()) {}
Setup::~Setup() = default;

void Setup::build(const hsp::str_table &input)
{
    const auto time_raw = std::chrono::steady_clock::now().time_since_epoch().count();
    // ensure time_now is 64bit
    const auto time_now = static_cast<std::uint64_t>(time_raw);
    // xor between upper and lower part of time stamp
    const std::uint64_t mixed = time_now ^ (time_now >> 32);

    // default seed for pseudo random number generator
    // note: constrain to 32bit for compatibility with R (mask sign bit)
    seed = static_cast<std::int32_t>(mixed & 0x7FFFFFFF);

    // current chain
    chain = 1;
    // current thermodynamic rung
    rung = 0;

    if (!input.empty())
    {
        // check no multiple entries if internal error
        if (hsp::read::csv_n_rows(input) != 1)
            throw std::runtime_error(
                "Internal error: multiple entries found in MCMC setup definition.");

        // check entries
        const auto &headers = hsp::extract_keys(input);
        hsp::str::assert_no_unused_entries(headers, file_headers);

        // Random Initialization of Parameters
        if (input.count("progress") && !hsp::na::is(input.at("progress")[0]))
            progress = hsp::str::toul(input.at("progress")[0]);

        // N. of Threads
        if (input.count("threads") && !hsp::na::is(input.at("threads")[0]))
            n_threads = hsp::str::toi(input.at("threads")[0], "threads");

        // PRNG Seed
        if (input.count("seed") && !hsp::na::is(input.at("seed")[0]))
            seed = hsp::str::toi32(input.at("seed")[0], "seed");

        // N. of Chains
        if (input.count("chains") && !hsp::na::is(input.at("chains")[0]))
            n_chains = hsp::str::toul(input.at("chains")[0], "chains");

        // N. of Steps
        if (input.count("steps") && !hsp::na::is(input.at("steps")[0]))
            n_steps = hsp::str::toul(input.at("steps")[0], "steps");

        // N. of Adaptive Steps
        if (input.count("adapt") && !hsp::na::is(input.at("adapt")[0]))
            n_adapt = hsp::str::toul(input.at("adapt")[0], "adpat");

        // N. of Thinning Steps
        if (input.count("thin") && !hsp::na::is(input.at("thin")[0]))
            n_thin = hsp::str::toul(input.at("thin")[0], "thin");

        // N. of Steps to "Burn-In"
        if (input.count("burnin") && !hsp::na::is(input.at("burnin")[0]))
            n_burnin = hsp::str::toul(input.at("burnin")[0], "burnin");

        // Individual ID for Augmented Trace
        if (input.count("id_draws") && !hsp::na::is(input.at("id_draws")[0]))
            id_draws = hsp::read::parse_list(input.at("id_draws")[0]);
    }

    set_n_threads(); // define number of threads
    set_prng_seed(); // pseudo-random number generator seeding
}

/// Set the number of threads for parrell computations (default is 1).
void Setup::set_n_threads()
{
    if (n_threads == 0)
        throw std::invalid_argument("N. of threads was not greater than 0.");

#ifdef _OPENMP
    unsigned int max_threads = omp_get_num_procs();

    if (n_threads > max_threads)
    {
        hsp::msg::warning(
            "N. of threads reduced to " + std::to_string(max_threads) + " available threads.");
        n_threads = max_threads;
    }

    // declare the global number of omp threads
    omp_set_num_threads(n_threads);

#else
    if (n_threads != 1)
        hsp::msg::warning(
            "Multi-threading was disabled. OpenMP is not supported on current machine.");

    n_threads = 1;
#endif
}

/// Assign seed for random generators.
/// Concomitantely create different generators for different parallel threads
void Setup::set_prng_seed()
{
    hsp::rnd_gen aux_prng;
    aux_prng.seed(seed);
    std::uniform_int_distribution<unsigned long> seed_generator(0, 0XFFFFFFFF);

    prng_parallel_.reserve(n_threads);
    for (unsigned n = 0; n < n_threads; ++n)
        prng_parallel_[n].seed(seed_generator(aux_prng));
}

/// Return pseudo random number generator.
hsp::rnd_gen &Setup::prng()
{
#ifdef _OPENMP
    return prng_parallel_[omp_get_thread_num()];
#else
    return prng_parallel_[0];
#endif
}

void Setup::assert_mcmc_setup()
{
    if (inv_T < 0. || inv_T > 1.)
        throw std::range_error(
            "Inverse temperature " + std::to_string(inv_T) + " was out of range [0, 1].");

    if (n_chains < 1)
        throw std::out_of_range("Number of chains was not greater than 0.");

    // Interpret thinning equal to 0 has been interpreted as 1.
    if (n_thin == 0)
        n_thin = 1;

    if (n_steps < 1)
        throw std::out_of_range(
            "N. of steps was not greater than 0.");
    if (n_burnin >= n_steps)
        throw std::out_of_range(
            "Burn-in must was not smaller than n. of steps.");
    if (n_thin > n_steps)
        throw std::out_of_range(
            "Thinning was not be smaller than (or equal to) n. of steps. ");

    c_sample_size = (n_steps - n_burnin) / n_thin;
    // adjust steps and sample size
    c_sample_size += (n_steps - n_burnin) % n_thin > 0;
    n_steps = n_burnin + c_sample_size * n_thin;

    if (c_sample_size == 0)
        throw std::out_of_range(
            "Choice of steps, burn-in and thin generated an empty sample.");

    // assign number of sweeps
    // NOTE: because of rounding I can't put all in the same line
    c_n_sweeps = n_adapt / n_thin;
    c_n_sweeps += n_burnin / n_thin;
    c_n_sweeps += c_sample_size;
}

// Drop adaptation steps.
void Setup::drop_mcmc_adapt()
{
    if (n_adapt > 0)
        hsp::msg::warning("Dropped adaptation steps: no parametric MCMC proposals were found.");
    n_adapt = 0;

    // re-assign number of sweeps
    c_n_sweeps = n_burnin / n_thin;
    c_n_sweeps += c_sample_size;
}

/// The text to print in the Markov chian progress bar.
string Setup::fmt_progtxt()
{
    if (is_power_post())
        return "power post " + std::to_string(rung) + "/" + std::to_string(n_rungs);
    if (n_chains > 1)
        return "sample " + std::to_string(chain) + "/" + std::to_string(n_chains);
    return "sample";
}

void Setup::prog_declare(const string &text)
{
    prog_set_text(text);
    progressbar_->set_thintime(progress);
}
void Setup::prog_set_text(const string &text) { progressbar_->set_text(hsp::na::is(text) ? fmt_progtxt() : text); }
void Setup::prog_open_for_inference() { progressbar_->open(c_n_sweeps * n_chains); }
void Setup::prog_open_for_power_post() { progressbar_->open(c_n_sweeps * n_chains * (n_rungs + 1)); }
void Setup::prog_update() { progressbar_->update(); }
void Setup::prog_close() { progressbar_->close(); }
void Setup::prog_crash() { progressbar_->crash(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** Household : class ***  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Link an individual to the household.
void Household::init_member(size_t index,
                            const string &id,
                            int status,
                            const hsp::dict<double> &data)
{
    // emplace new individual
    members_.emplace_back(index, house_size_, id, status, data);

    // update status lists within the household
    if (members_.back().is_infected_status())
        infected_list_.push_back(house_size_);
    else if (members_.back().is_healthy_status())
        healthy_list_.push_back(house_size_);
    else
        missing_list_.push_back(house_size_);

    ++house_size_;
}

/// Change_status of individual i to new status (without controls).
///	NOTE: this introduces augmented variables -> it requires RJMCMC in case of inference.
void Household::set_status_(Individual &member, int new_status)
{
    if (new_status == member.status())
        return;

    const size_t i = static_cast<size_t>(member); // get within house index

    if (member.is_infected_status())
        infected_list_.erase(std::find(infected_list_.begin(), infected_list_.end(), i));
    else if (member.is_healthy_status())
        healthy_list_.erase(std::find(healthy_list_.begin(), healthy_list_.end(), i));
    else
        missing_list_.erase(std::find(missing_list_.begin(), missing_list_.end(), i));

    member.status_ = new_status; // change member status
    if (member.is_infected_status())
        infected_list_.push_back(i);
    else if (member.is_healthy_status())
        healthy_list_.push_back(i);
    else
        missing_list_.push_back(i);
}

/// Associate a pointer to parameter value for the contact from->to.
void Household::init_contact_rates(size_t from, size_t to, const Attribute &x)
{
    // initialize contact rate matrix if necessary
    if (contact_rates_.empty())
        contact_rates_.resize(
            house_size_,
            vector<Attribute>(house_size_, Attribute(hsp::rates::VALUE_REFERENCE)));

    contact_rates_[from][to] = x;
}

/// Assert that current log-likelihood id not a NaN.
void Household::assert_log_lik_not_NaN() const
{
    if (hsp::na::is(obs_log_lik_))
        throw std::runtime_error("Observation log-likelihood returned NaN.");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Household : cache attributes  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize cached structeres according to household size.
void Household::init_cache()
{
    const size_t hh_size = size();

    // size_t
    cache_size_t_vec_.resize(hsp::key::Registry<HHsize_tKey>::instance().size());
    for (auto &v : cache_size_t_vec_)
        v.assign(hh_size, 0);
    buffer_size_t_vec_ = cache_size_t_vec_;

    // FLAG
    // double
    cache_dbl_vec_.resize(hsp::key::Registry<HHdoubleKey>::instance().size());
    for (auto &v : cache_dbl_vec_)
        v.assign(hh_size, 0.);
    buffer_dbl_vec_ = cache_dbl_vec_;

    // matrix
    cache_dbl_mtx_.resize(hsp::key::Registry<HHmatrixKey>::instance().size());
    for (auto &v : cache_dbl_mtx_)
        v.assign(hh_size, std::vector<double>(hh_size, 0.));
    buffer_dbl_mtx_ = cache_dbl_mtx_;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** DataWide : struct ***  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

DataWide::DataWide() : report_indiv(std::make_unique<ReportIndividual>()),
                       report_house(std::make_unique<ReportHousehold>()),
                       report_hh_size(std::make_unique<ReportHouseSizes>()) {}
DataWide::~DataWide() = default;

/// Assert data_wide is not empty.
void DataWide::assert_non_empty()
{
    if (report_house->total == 0)
        throw std::length_error("Unable to proceed: no data was found.");
}
/// Warn if all individuals have missing status.
void DataWide::warn_all_missing_status()
{
    if (report_indiv->total == report_indiv->missing)
        hsp::msg::warning("All individuals in the data had missing status.");
}

/// Assign data after pre-processing.
void DataWide::build(const hsp::str_table &input)
{
    if (input.empty())
        return;

    // drop all previous data and prepare for new
    build_helper_clear_data_();

    // extract headers
    for (const auto &[col, vec] : input)
        headers.push_back(col);
    // assess the presence of house IDs
    if (!(headers % "house_id"))
        hsp::msg::warning(
            "Ignored household structure. Missing data column 'house_id'.");
    // assess individual IDs are unique
    if (headers % "id" && input.at("id").size() != hsp::str::unique(input.at("id")).size())
        hsp::msg::error(
            "Provided individual IDs are not unique.");

    // safe load of households
    report_indiv->total = input.at(headers[0]).size();
    for (size_t i = 0; i < report_indiv->total; ++i)
    {
        string id = std::to_string(i + 1);                          // default id
        string house_id = std::to_string(house_vector_.size() + 1); // default house_id
        int status = hsp::status::MISSING;
        hsp::dict<double> data = {};

        // read individual data
        for (const auto &key : headers)
            if (key == "id")
            {
                id = input.at(key)[i]; // user id
                if (hsp::na::is(id))
                    throw std::invalid_argument(
                        "Individual in row n." + std::to_string(i + 2) + " had missing " + key + ".");
            }
            else if (key == "house_id")
            {
                house_id = input.at(key)[i]; // user house_id
                if (hsp::na::is(house_id))
                    throw std::invalid_argument(
                        "Individual in row n." + std::to_string(i + 2) + " had missing " + key + ".");
            }
            else if (key == "status")
            {
                // note: NA status is a missing status
                if (!hsp::na::is(input.at("status")[i]))
                    status = hsp::str::toi(input.at("status")[i], "individual status");
            }
            else
            {
                data.emplace(key, hsp::str::tod(input.at(key)[i], "data '" + key + "'"));
            }

        // if household index doesn't exist, add new household
        if (!hhid_to_index_.count(house_id))
        {
            // map data index to household number
            hhid_to_index_.emplace(house_id, house_vector_.size());
            // create a new empty household
            house_vector_.emplace_back(house_id, house_vector_.size());
        }
        // extract the household reference
        auto &house = house_vector_.at(hhid_to_index_.at(house_id));

        // assign the individual the respective household
        house.init_member(i, id, status, data);
    }
    report_house->total = house_vector_.size();

    headers = headers - vector<string>({"id", "house_id", "status"});

    // update report for data
    build_helper_fill_report_data_();
}

/// Use the loaded data to fill-in the respective report.
void DataWide::build_helper_clear_data_()
{
    headers.clear();
    house_vector_.clear();
    hhid_to_index_.clear();
    all_infected_status_list.clear();
    report_indiv = std::make_unique<ReportIndividual>();
    report_house = std::make_unique<ReportHousehold>();
    report_hh_size = std::make_unique<ReportHouseSizes>();
}

/// Use the loaded data to fill-in the respective report.
void DataWide::build_helper_fill_report_data_()
{
    // set report counts
    for (const auto &house : house_vector_)
    {
        // list of unique infected statuses
        for (const auto &member : house.infected())
            if (!(all_infected_status_list % member.status()))
                all_infected_status_list.push_back(member.status());

        // observed quantities in-house
        const unsigned short house_size = house.size();
        const unsigned short n_infec = house.n_infected_status();
        const unsigned short n_miss = house.n_missing_status();
        const unsigned short n_susc = house.n_susceptible_status();

        report_indiv->infectious += n_infec;
        report_indiv->missing += n_miss;
        report_indiv->susceptible += n_susc;

        report_house->no_cases += (n_infec == 0);
        report_house->at_least_1_case += (n_infec > 0);
        report_house->at_least_2_cases += (n_infec > 1);

#define UPDATE_REPORT_SIZE(Report, Value, Symbol) \
    if (Report.first Symbol Value)                \
    {                                             \
        Report.first = Value;                     \
        Report.second = 1;                        \
    }                                             \
    else if (Report.first == Value)               \
        ++Report.second;

        // max number of house cases
        UPDATE_REPORT_SIZE(report_hh_size->max_cases, n_infec, <);
        // max number of missing status
        UPDATE_REPORT_SIZE(report_hh_size->max_miss, n_miss, <);
        // minimum house number of members
        UPDATE_REPORT_SIZE(report_hh_size->min_size, house_size, >);
        // maximum house number of members
        UPDATE_REPORT_SIZE(report_hh_size->max_cases, house_size, <);
    }
}

/// Build multiplicities.
/// UPDATE: I should find a way to deal internally with multiplicity.
void DataWide::multiplicity_build(const hsp::str_table &input, bool is_any_augm)
{
    head_class_house_indexes.clear();
    auto n_houses = report_house->total;
    if (is_any_augm || (input.size() == 0))
    {
        // note: multiplicities are not compatible with data augmentation
        if (!input.empty())
            hsp::msg::warning(
                "Ignored multiplicity optimization: incompatible with data augmentation.");

        // default choice for restricted household indexes
        head_class_house_indexes.resize(n_houses);
        for (size_t H = 0; H < n_houses; ++H)
        {
            house_vector_[H].set_multiplicity(1, H);
            head_class_house_indexes[H] = H;
        }
    }
    else
    {
        // UPDATE: this error is purely internal
        if (input.at("house_id").size() != n_houses)
            throw std::out_of_range(
                "Multiplicity definitions had different dimensions than data.");

        for (size_t r = 0; r < n_houses; ++r)
        {
            auto &house = (*this)(input.at("house_id")[r]);
            auto H0 = (*this)(input.at("class_id")[r]).house_index();
            auto m = hsp::str::toul(input.at("multiplicity")[r]);

            house.set_multiplicity(m, H0);

            if (house.is_multiplicity_head_class())
                head_class_house_indexes.emplace_back(house.house_index());
        }
    }
}

/// Initialize posterior sample of parameters.
void DataWide::init_draws(unsigned long sample_size)
{
    // sample over households
    smpl_values.clear();
    smpl_values.resize(head_class_house_indexes.size());
    for (size_t R = 0; R < head_class_house_indexes.size(); ++R)
    {
        smpl_values[R].clear();
        smpl_values[R].reserve(sample_size);
        smpl_values[R].assign(sample_size, hsp::na::DOUBLE);
    }
}

/// Store current sample of log-likelihood for restricted households.
void DataWide::store_draws(unsigned long smpl_count)
{
    // store sample for restricted households
    for (size_t R = 0; R < head_class_house_indexes.size(); ++R)
    {
        auto &house = house_vector_[head_class_house_indexes[R]];
        double x = house.get_obs_log_lik() / house.multiplicity();
        for (const auto &member : house.infected())
            x += member.get_augm_model_log_lik();
        smpl_values[R][smpl_count] = x;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : DataLong : build  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Assign data after pre-processing.
void DataLong::build(const hsp::str_table &input)
{
    if (input.empty())
        return;

    // data sequence headers
    headers.clear();
    for (const auto &x : input)
        headers.push_back(x.first);

    // check for individual IDs
    if (!(headers % "id"))
        throw std::invalid_argument(
            "Individual identifier column 'id' was missing.");
    // remove "id"
    headers = headers - vector<string>({"id"});

    // individual map for id -> index
    // note: this could be moved elsewhere if required
    hsp::dict<size_t> id_to_index = {};
    size_t c = 0;
    for (const auto &house : *pDataWide_)
        for (const auto &member : house)
        {
            if (id_to_index.count(member.id()))
                throw std::invalid_argument(
                    "There were multiple individuals with id equal to '" + member.id() + "'.");
            id_to_index.emplace(member.id(), c++);
        }

    // data-long allocation
    indiv_datalist.clear();
    indiv_datalist.resize(pDataWide_->report_indiv->total);

    for (auto &xmap : indiv_datalist)
        for (const auto &key : headers)
            xmap[key].clear();

    // load individual data list
    const auto n_rows = input.at(headers[0]).size();
    for (size_t r = 0; r < n_rows; ++r)
    {
        const auto &id = input.at("id")[r];
        if (hsp::na::is(id))
            throw std::invalid_argument(
                "Invalid id found in row n." + std::to_string(r + 2) + ".");
        if (!id_to_index.count(id))
        {
            // UPDATE: these errors should be more self-explicative
            if (pDataWide_->is_empty())
                throw std::invalid_argument(
                    "Tried to define data long format before wide format.");
            throw std::invalid_argument(
                "Id '" + id + "' from data long format had no match in loaded data wide format.");
        }

        // append (since vector is initialized to empty)
        auto &xdict = indiv_datalist.at(id_to_index.at(id));
        for (const auto &key : headers)
            xdict.at(key).emplace_back(
                hsp::str::tod(input.at(key).at(r), "data long format '" + key + "'"));
    }

    // assign pointers to relative data-long to all individuals
    c = 0;
    for (auto &house : *pDataWide_)
        for (auto &member : house)
            member.init_pointer_to_data_long(&indiv_datalist.at(c++));
}
