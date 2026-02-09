/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/keys.h"
#include "tools/loglike.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/openmp.h"
#include "tools/operators.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"

#include "data_wide.h"
#include "attribute.h"
#include "tools/data_reports.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** Household : class ***  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Emplace a new individual in the household.
void Household::init_member(size_t index,
                            const std::string &id,
                            int status,
                            const std::vector<double> &data)
{
    size_t hhindex = size();
    members_.emplace_back(index, hhindex, id, status, data);

    // update status lists within the household
    if (members_.back().is_infected_status())
        infected_list_.emplace(hhindex);
    else if (members_.back().is_healthy_status())
        healthy_list_.emplace(hhindex);
    else
        missing_list_.emplace(hhindex);
}

/// Change_status of individual i to new status (without controls).
///	NOTE: this introduces augmented variables -> it requires RJMCMC in case of inference.
void Household::helper_set_status_(Individual &member, int new_status)
{
    if (new_status == member.status())
        return;

    size_t i = member.hhindex(); // get within house index

    if (member.is_infected_status())
        infected_list_.erase(i);
    else if (member.is_healthy_status())
        healthy_list_.erase(i);
    else
        missing_list_.erase(i);

    member.status_ = new_status; // change member status
    if (member.is_infected_status())
        infected_list_.emplace(i);
    else if (member.is_healthy_status())
        healthy_list_.emplace(i);
    else
        missing_list_.emplace(i);
}

/// Associate a pointer to parameter value for the contact from->to.
void Household::init_contact_rates(size_t from, size_t to, const Attribute &x)
{
    // initialize contact rate matrix if necessary
    if (contact_rates_.empty())
        contact_rates_.assign(size(), std::vector<Attribute>(size(), hsp::rates::VALUE_REFERENCE));
    contact_rates_[from][to] = x;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** DataWide : struct ***  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

DataWide::DataWide(
    const std::vector<std::string> &keys,
    const std::vector<std::string> &size_keys,
    const std::vector<std::string> &size_vec_keys,
    const std::vector<std::string> &double_vec_keys,
    const std::vector<std::string> &double_mtx_keys)
    : pRegistry(std::make_unique<hsp::key::Registry>(keys)),
      p_size_Registry(std::make_unique<hsp::key::Registry>(size_keys)),
      p_size_vec_Registry(std::make_unique<hsp::key::Registry>(size_vec_keys)),
      p_double_vec_Registry(std::make_unique<hsp::key::Registry>(double_vec_keys)),
      p_double_mtx_Registry(std::make_unique<hsp::key::Registry>(double_mtx_keys)) {}
DataWide::~DataWide() = default;

/// Assign data after pre-processing.
void DataWide::build(const hsp::str_table &input)
{
    pRegistry->clear();
    p_size_Registry->clear();
    p_size_vec_Registry->clear();
    p_double_vec_Registry->clear();
    p_double_mtx_Registry->clear();

    if (input.empty())
        return;

    helper_build_clear_data_();        // drop previous data and lists
    helper_build_data_headers_(input); // extract headers
    report_indiv.total = hsp::read::csv_n_rows(input);

    // store households
    for (size_t i = 0; i < report_indiv.total; ++i)
    {
        std::string id = std::to_string(i + 1);                          // default id
        std::string house_id = std::to_string(house_vector_.size() + 1); // default house_id
        int status = hsp::status::MISSING;

        if (input.count("id"))
        {
            id = input.at("id")[i]; // user id
            if (hsp::na::is(id))
                throw std::invalid_argument(
                    "Individual in row n." + std::to_string(i + 2) + " had missing 'id'.");
        }
        if (input.count("house_id"))
        {
            house_id = input.at("house_id")[i]; // user house_id
            if (hsp::na::is(house_id))
                throw std::invalid_argument(
                    "Individual in row n." + std::to_string(i + 2) + " had missing 'house_id'.");
        }
        if (input.count("status"))
        {
            // note: NA status is a missing status
            if (!hsp::na::is(input.at("status")[i]))
                status = hsp::str::toi(input.at("status")[i], "individual status");
        }

        // read individual data
        std::vector<double> data(pRegistry->size());
        for (size_t d = 0; d < data.size(); ++d)
        {
            const std::string datastring = pRegistry->as_string(d);
            data[d] = hsp::str::tod(input.at(datastring)[i], "data '" + datastring + "'");
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
        Household &house = house_vector_.at(hhid_to_index_.at(house_id));

        // assign the individual the respective household
        house.init_member(i, id, status, data);
    }

    // update report for data
    build_helper_fill_report_data_();
}

/// Get housepi household index from house_id.
size_t DataWide::house_id_to_index_(std::string_view house_id) const
{
    auto it = hhid_to_index_.find(house_id);
    if (it != hhid_to_index_.end())
        return it->second;
    throw std::invalid_argument(
        "Index '" + std::string(house_id) + "' did not match any household in the data.");
}

/// Clear previous data report.
void DataWide::helper_build_clear_data_()
{
    house_vector_.clear();
    hhid_to_index_.clear();
    all_infected_statuses.clear();
    report_indiv = ReportIndividual();
    report_house = ReportHousehold();
    report_hh_size = ReportHouseSizes();
}

/// Assert that data headers are correct.
void DataWide::helper_build_data_headers_(const hsp::str_table &input)
{
    std::vector<std::string> headers = hsp::read::extract_keys(input);

    // assess the presence of house IDs
    if (std::find(headers.begin(), headers.end(), "house_id") == headers.end())
        hsp::msg::warning("Ignored household structure. Missing data column 'house_id'.");

    // assess individual IDs are unique
    if (std::find(headers.begin(), headers.end(), "id") != headers.end() &&
        input.at("id").size() != hsp::str::unique(input.at("id")).size())
        hsp::msg::error("Provided individual IDs are not unique.");

    // remove standard headers
    headers = headers - std::vector<std::string>({"id", "house_id", "status"});

    // pedantic check for declared data keys
    hsp::str::assert_no_missing_entries(headers, pRegistry->get_all());

    // add optional keys to registry
    for (const std::string &datastring : headers)
        if (!pRegistry->count(datastring))
            pRegistry->add(datastring);
}

/// Use the loaded data to fill-in the respective report.
void DataWide::build_helper_fill_report_data_()
{
    report_house.total = house_vector_.size();

    // set report counts
    for (const Household &house : house_vector_)
    {
        // list of unique infected statuses
        for (size_t i : house.infected_hhindexes())
            all_infected_statuses.emplace(house.member(i).status());

        // observed quantities in-house
        size_t house_size = house.size();
        size_t n_infec = house.n_infected_status();
        size_t n_miss = house.n_missing_status();
        size_t n_susc = house.n_susceptible_status();

        report_indiv.infectious += n_infec;
        report_indiv.missing += n_miss;
        report_indiv.susceptible += n_susc;

        report_house.no_cases += (n_infec == 0);
        report_house.at_least_1_case += (n_infec > 0);
        report_house.at_least_2_cases += (n_infec > 1);

#define UPDATE_REPORT_SIZE(Report, Value, Symbol) \
    if (Report.first Symbol Value)                \
    {                                             \
        Report.first = Value;                     \
        Report.second = 1;                        \
    }                                             \
    else if (Report.first == Value)               \
        ++Report.second;

        // max number of house cases
        UPDATE_REPORT_SIZE(report_hh_size.max_cases, n_infec, <);
        // max number of missing status
        UPDATE_REPORT_SIZE(report_hh_size.max_miss, n_miss, <);
        // minimum house number of members
        UPDATE_REPORT_SIZE(report_hh_size.min_size, house_size, >);
        // maximum house number of members
        UPDATE_REPORT_SIZE(report_hh_size.max_cases, house_size, <);
    }
}

/// Build multiplicities.
/// UPDATE: I should find a way to deal internally with multiplicity.
void DataWide::multiplicity_build(const hsp::str_table &input, bool is_any_augm)
{
    head_class_house_indexes.clear();
    size_t n_houses = report_house.total;
    if (is_any_augm || input.size() == 0)
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
            Household &house = at(input.at("house_id")[r]);
            size_t H0 = at(input.at("class_id")[r]).house_index();
            unsigned long m = hsp::str::toul(input.at("multiplicity")[r]);

            house.set_multiplicity(m, H0);

            if (house.is_multiplicity_head_class())
                head_class_house_indexes.emplace_back(house.house_index());
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  DataWide : cache attributes  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void DataWide::init_house_cache(Household &house)
{
    size_t hh_size = house.size();

    // size_t
    {
        size_t n = p_size_Registry->size();
        house.cache_size_.assign(n, 0);
        house.buffer_size_.assign(n, 0);
    }

    // size_t vector
    {
        size_t n = p_size_vec_Registry->size();
        house.cache_size_vec_.resize(n);
        house.buffer_size_vec_.resize(n);
        for (size_t e = 0; e < n; ++e)
        {
            house.cache_size_vec_[e].assign(hh_size, 0);
            house.buffer_size_vec_[e].assign(hh_size, 0);
        }
    }

    // double vector
    {
        size_t n = p_double_vec_Registry->size();
        house.cache_dbl_vec_.resize(n);
        house.buffer_dbl_vec_.resize(n);
        for (size_t e = 0; e < n; ++e)
        {

            house.cache_dbl_vec_[e].assign(hh_size, 0);
            house.buffer_dbl_vec_[e].assign(hh_size, 0);
        }
    }

    // double matrix
    {
        size_t n = p_double_mtx_Registry->size();
        house.cache_dbl_mtx_.resize(n);
        house.buffer_dbl_mtx_.resize(n);
        for (size_t e = 0; e < n; ++e)
        {

            house.cache_dbl_mtx_[e].assign(hh_size, std::vector<double>(hh_size, 0.));
            house.buffer_dbl_mtx_[e].assign(hh_size, std::vector<double>(hh_size, 0.));
        }
    }
}

void DataWide::buffer_all_cache(Household &house)
{
    std::copy(house.cache_size_.begin(),
              house.cache_size_.end(),
              house.buffer_size_.begin());

    for (size_t e = 0; e < house.cache_size_vec_.size(); ++e)
        std::copy(house.cache_size_vec_[e].begin(),
                  house.cache_size_vec_[e].end(),
                  house.buffer_size_vec_[e].begin());

    for (size_t e = 0; e < house.cache_dbl_vec_.size(); ++e)
        std::copy(house.cache_dbl_vec_[e].begin(),
                  house.cache_dbl_vec_[e].end(),
                  house.buffer_dbl_vec_[e].begin());

    for (size_t e = 0; e < house.cache_dbl_mtx_.size(); ++e)
        for (size_t i = 0; i < house.cache_dbl_mtx_[e].size(); ++i)
            std::copy(house.cache_dbl_mtx_[e][i].begin(),
                      house.cache_dbl_mtx_[e][i].end(),
                      house.buffer_dbl_mtx_[e][i].begin());
}

void DataWide::buffer_member_cache(Household &house, size_t i)
{
    std::copy(house.cache_size_.begin(),
              house.cache_size_.end(),
              house.buffer_size_.begin());

    for (size_t e = 0; e < house.cache_size_vec_.size(); ++e)
        house.buffer_size_vec_[e][i] = house.cache_size_vec_[e][i];

    for (size_t e = 0; e < house.cache_dbl_vec_.size(); ++e)
        house.buffer_dbl_vec_[e][i] = house.cache_dbl_vec_[e][i];

    for (size_t e = 0; e < house.cache_dbl_mtx_.size(); ++e)
    {
        std::copy(house.cache_dbl_mtx_[e][i].begin(),
                  house.cache_dbl_mtx_[e][i].end(),
                  house.buffer_dbl_mtx_[e][i].begin());

        for (size_t j = 0; j < house.size(); ++j)
            house.buffer_dbl_mtx_[e][j][i] = house.cache_dbl_mtx_[e][j][i];
    }
}

void DataWide::recover_all_buffer(Household &house)
{
    house.cache_size_.swap(house.buffer_size_);
    house.cache_size_vec_.swap(house.buffer_size_vec_);
    house.cache_dbl_vec_.swap(house.buffer_dbl_vec_);
    house.cache_dbl_mtx_.swap(house.buffer_dbl_mtx_);
}

void DataWide::recover_member_buffer(Household &house, size_t i)
{
    house.cache_size_.swap(house.buffer_size_);

    for (size_t e = 0; e < house.cache_size_vec_.size(); ++e)
        house.cache_size_vec_[e][i] = house.buffer_size_vec_[e][i];

    for (size_t e = 0; e < house.cache_dbl_vec_.size(); ++e)
        house.cache_dbl_vec_[e][i] = house.buffer_dbl_vec_[e][i];

    for (size_t e = 0; e < house.cache_dbl_mtx_.size(); ++e)
    {
        house.cache_dbl_mtx_[e][i].swap(house.buffer_dbl_mtx_[e][i]);
        for (size_t j = 0; j < house.size(); ++j)
            if (j != i)
                house.cache_dbl_mtx_[e][j][i] = house.buffer_dbl_mtx_[e][j][i];
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  DataWide : posterior draws  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void DataWide::assert_finite_log_lik() const
{
    for (size_t R = 0; R < head_class_house_indexes.size(); ++R)
    {
        const Household &house = house_vector_[head_class_house_indexes[R]];

        for (size_t i : house.infected_hhindexes())
            if (house.member(i).augm_log_lik_.is_inf())
                throw std::runtime_error(
                    "Initial augmented log-likelihood computation returned infinity for " +
                    house.member(i).err_name() + ".");

        if (house.obs_log_lik_.is_inf())
            throw std::runtime_error(
                "Initial log-likelihood computation returned infinity for " + house.err_name() + ".");
    }
}

/// Initialize posterior sample of parameters.
void DataWide::init_draws(unsigned long sample_size)
{
    // sample over households
    smpl_values.clear();
    smpl_values.reserve(head_class_house_indexes.size());
    for (size_t R = 0; R < head_class_house_indexes.size(); ++R)
        smpl_values.emplace_back(sample_size, hsp::na::DOUBLE);
}

/// Store current sample of log-likelihood for restricted households.
/// UPDATE: multiplicity should not be multiplied within the compute_log_lik!
void DataWide::store_draws(unsigned long smpl_count)
{
    // store sample for restricted households
    for (size_t R = 0; R < head_class_house_indexes.size(); ++R)
    {
        Household &house = house_vector_[head_class_house_indexes[R]];
        double x = house.obs_log_lik_.value();
        for (size_t i : house.infected_hhindexes())
            x += house.member(i).augm_log_lik_.value();
        smpl_values[R][smpl_count] = x;
    }
}
