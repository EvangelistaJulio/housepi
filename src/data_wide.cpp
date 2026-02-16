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
#include "function.h"
#include "tools/data_reports.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** Household : class ***  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Household::init_member(size_t index,
                            const std::string &id,
                            int status,
                            const std::vector<double> &data)
{
    size_t hhindex = size();
    members_.emplace_back();

    Individual &last = members_.back();
    last.index_ = index;
    last.hhindex_ = hhindex;
    last.id_ = id;
    last.status_ = status;
    last.data_vec_ = data;

    // update status lists within the household
    if (last.is_infected_status())
        infected_list_.emplace_back(hhindex);
    else if (last.is_healthy_status())
        healthy_list_.emplace_back(hhindex);
    else
        missing_list_.emplace_back(hhindex);
}

void Household::helper_set_status_(Individual &member, int new_status)
{
    if (new_status == member.status())
        return;

    size_t i = member.hhindex(); // get within house index

    if (member.is_infected_status())
        infected_list_.erase(std::find(infected_list_.begin(), infected_list_.end(), i));
    else if (member.is_healthy_status())
        healthy_list_.erase(std::find(healthy_list_.begin(), healthy_list_.end(), i));
    else
        missing_list_.erase(std::find(missing_list_.begin(), missing_list_.end(), i));

    member.status_ = new_status; // change member status
    if (member.is_infected_status())
        infected_list_.emplace_back(i);
    else if (member.is_healthy_status())
        healthy_list_.emplace_back(i);
    else
        missing_list_.emplace_back(i);
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
            house_vector_.emplace_back();
            Household &last = house_vector_.back();
            last.house_id_ = house_id;
        }
        // extract the household reference
        Household &house = house_vector_.at(hhid_to_index_.at(house_id));

        // assign the individual the respective household
        house.init_member(i, id, status, data);
    }

    // update report for data
    build_helper_fill_report_data_();
}

size_t DataWide::house_id_to_index_(std::string_view house_id) const
{
    auto it = hhid_to_index_.find(house_id);
    if (it != hhid_to_index_.end())
        return it->second;
    throw std::invalid_argument(
        "Index '" + std::string(house_id) + "' did not match any household in the data.");
}

void DataWide::helper_build_clear_data_()
{
    house_vector_.clear();
    hhid_to_index_.clear();
    report_indiv = ReportIndividual();
    report_house = ReportHousehold();
    report_hh_size = ReportHouseSizes();
}

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
        pRegistry->add(datastring);
}

void DataWide::build_helper_fill_report_data_()
{
    report_house.total = house_vector_.size();

    // set report counts
    for (const Household &house : house_vector_)
    {
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

void DataWide::multiplicity_build(const hsp::str_table &input, bool is_any_augm)
{
    size_t n_houses = report_house.total;

    mode_ = Mode::Unique;
    multiplicity_.clear();
    multiplicity_class_.clear();

    // initialize household indexes
    for (Household &house : house_vector_)
        house.restricted_ = hhid_to_index_[house.house_id()];

    // note: multiplicities are not compatible with data augmentation
    if (input.size() == 0 || is_any_augm)
        return;

    // UPDATE: this error is purely internal
    if (input.at("house_id").size() != n_houses)
        throw std::out_of_range(
            "Multiplicity definitions had different dimensions than data.");

    mode_ = Mode::Multiple;
    multiplicity_.reserve(n_houses);
    multiplicity_class_.reserve(n_houses);

    hsp::dict<size_t> hhid_to_restricted = {};
    for (size_t r = 0; r < n_houses; ++r)
    {
        const std::string &house_id = input.at("house_id")[r];
        const std::string &house_id_0 = input.at("class_id")[r];
        unsigned long m = hsp::str::toul(input.at("multiplicity")[r]);

        if (house_id == house_id_0)
        {
            hhid_to_restricted.emplace(house_id, hhid_to_restricted.size());
            multiplicity_.push_back(m);
            multiplicity_class_.push_back(hhid_to_index_[house_id]);
        }
        at(house_id).restricted_ = hhid_to_restricted[house_id_0];
    }
    multiplicity_.shrink_to_fit();
    multiplicity_class_.shrink_to_fit();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  DataWide : cache attributes  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void DataWide::allocate_all_cache()
{
    cached_.resize(restricted_size());
    buffer_.resize(restricted_size());

    for (size_t R = 0; R < restricted_size(); ++R)
    {
        Household &house = restricted(R);
        size_t hh_size = house.size();
        Cached_ &curr = cached_[R];
        Cached_ &prev = buffer_[R];

        // size_t
        curr.size_.assign(p_size_Registry->size(), 0);
        prev.size_.assign(p_size_Registry->size(), 0);

        // size_t vector
        curr.size_vec_.resize(p_size_vec_Registry->size());
        prev.size_vec_.resize(p_size_vec_Registry->size());
        for (size_t e = 0; e < p_size_vec_Registry->size(); ++e)
        {
            curr.size_vec_[e].assign(hh_size, 0);
            prev.size_vec_[e].assign(hh_size, 0);
        }

        // double vector
        curr.dbl_vec_.resize(p_double_vec_Registry->size());
        prev.dbl_vec_.resize(p_double_vec_Registry->size());
        for (size_t e = 0; e < p_double_vec_Registry->size(); ++e)
        {

            curr.dbl_vec_[e].assign(hh_size, 0);
            prev.dbl_vec_[e].assign(hh_size, 0);
        }

        // double matrix
        curr.dbl_mtx_.resize(p_double_mtx_Registry->size());
        prev.dbl_mtx_.resize(p_double_mtx_Registry->size());
        for (size_t e = 0; e < p_double_mtx_Registry->size(); ++e)
        {
            curr.dbl_mtx_[e].assign(hh_size, std::vector<double>(hh_size, 0.));
            prev.dbl_mtx_[e].assign(hh_size, std::vector<double>(hh_size, 0.));
        }
    }
}

void DataWide::buffer_member_cache(Household &house, size_t i)
{
    Cached_ &curr = cached_[house.restricted_];
    Cached_ &prev = buffer_[house.restricted_];
    prev.size_ = curr.size_;

    for (size_t e = 0; e < curr.size_vec_.size(); ++e)
        prev.size_vec_[e][i] = curr.size_vec_[e][i];

    for (size_t e = 0; e < curr.dbl_vec_.size(); ++e)
        prev.dbl_vec_[e][i] = curr.dbl_vec_[e][i];

    for (size_t e = 0; e < curr.dbl_mtx_.size(); ++e)
        for (size_t j = 0; j < house.size(); ++j)
        {
            prev.dbl_mtx_[e][i][j] = curr.dbl_mtx_[e][i][j];
            if (i != j)
                prev.dbl_mtx_[e][j][i] = curr.dbl_mtx_[e][j][i];
        }
}

void DataWide::recover_member_buffer(Household &house, size_t i)
{
    Cached_ &curr = cached_[house.restricted_];
    Cached_ &prev = buffer_[house.restricted_];
    curr.size_.swap(prev.size_);

    for (size_t e = 0; e < prev.size_vec_.size(); ++e)
        curr.size_vec_[e][i] = prev.size_vec_[e][i];

    for (size_t e = 0; e < prev.dbl_vec_.size(); ++e)
        curr.dbl_vec_[e][i] = prev.dbl_vec_[e][i];

    for (size_t e = 0; e < prev.dbl_mtx_.size(); ++e)
    {
        curr.dbl_mtx_[e][i].swap(prev.dbl_mtx_[e][i]);
        for (size_t j = 0; j < house.size(); ++j)
            if (j != i)
                curr.dbl_mtx_[e][j][i] = prev.dbl_mtx_[e][j][i];
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  DataWide : posterior draws  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void DataWide::assert_finite_log_lik() const
{
    for (size_t R = 0; R < restricted_size(); ++R)
    {
        const Household &house = restricted(R);

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

void DataWide::init_draws(unsigned long sample_size)
{
    // sample over households
    smpl_values.clear();
    smpl_values.reserve(restricted_size());
    for (size_t R = 0; R < restricted_size(); ++R)
        smpl_values.emplace_back(sample_size, hsp::na::DOUBLE);
}

void DataWide::store_draws(unsigned long smpl_count)
{
    // store sample for restricted households
    for (size_t R = 0; R < restricted_size(); ++R)
    {
        Household &house = restricted(R);
        double x = house.obs_log_lik_.value();
        for (size_t i : house.infected_hhindexes())
            x += house.member(i).augm_log_lik_.value();
        smpl_values[R][smpl_count] = x;
    }
}

void DataWide::assert_non_empty() const
{
    if (report_house.total == 0)
        throw std::length_error(
            "Unable to proceed: no data was found.");
}

void DataWide::assert_not_all_missing_status() const
{
    if (report_indiv.total == report_indiv.missing)
        throw std::length_error(
            "All individuals in the data had missing status.");
}

void DataWide::assert_no_missing_status() const
{
    if (report_indiv.missing > 0)
        throw std::length_error(
            "Individual with missing status are incompatible with current model.");
}
